// SARS-CoV-2 reference genome
params.reference='https://raw.githubusercontent.com/nextstrain/nextclade_data/master/data/datasets/sars-cov-2/references/MN908947/versions/2022-04-08T12%3A00%3A00Z/files/reference.fasta'
// Nextstrain Clade definition
params.virusprop='https://raw.githubusercontent.com/nextstrain/nextclade_data/master/data/datasets/sars-cov-2/references/MN908947/versions/2022-04-08T12%3A00%3A00Z/files/virus_properties.json'
// Result folder
params.results="results"
// FastQ directory
params.fastqs="data/"
// Number of nt to trim from reads
params.trim=0
// First clade to compare
params.clade1="-1"
// Second clade to compare
params.clade2="-1"
params.help=false

if (params.help) {
    def usage = "geva/covseqclade usage:\nnextflow run geva/covseqclade --fastq <fastq directory> --results <result directory> [--clade1 <first clade to compare> --clade2 <second clade to compare>]"
    log.info usage
    exit 0
}


results=file(params.results)
reference=file(params.reference)
trim=params.trim

fastqs=Channel.fromFilePairs(params.fastqs+'/*_R{1,2}.fastq.gz', flat: true)
dsrc=Channel.fromFilePairs(params.fastqs+'/*_R{1,2}.fastq.dsrc', flat: true)
virusprop=file(params.virusprop)
clade1=params.clade1
clade2=params.clade2

results.mkdir()

// Indexing Reference genome
process bwaIndexGenome {
    label 'bwaindex'

    input:
    file reference

    output:
    tuple val(reference.name), file("${reference}*") into genomeIndex

    script:
    """
    bwa index ${reference}
    """
}

// If input files are dsrc compressed
process dsrc2gz {
	label 'dsrc'
	tag "${name}"

	input:
	tuple val(name), file(r1), file(r2) from dsrc

	output:
	tuple val(name), file("${r1.baseName}.gz"), file("${r2.baseName}.gz") into dsrcgz

	script:
	"""
	dsrc d -t${(int)(task.cpus/2)} -s ${r1} | pigz -p ${(int)(task.cpus/2)} -c -  > ${r1.baseName}.gz
	dsrc d -t${(int)(task.cpus/2)} -s ${r2} | pigz -p ${(int)(task.cpus/2)} -c -  > ${r2.baseName}.gz
	"""
}

samplefiles=fastqs.concat(dsrcgz).filter{ name, f1, f2 -> f1.size()>0 }

// Trim reads
process trimGalore {
    label 'trimgalore'
    tag "${name}"

    input:
    tuple val(name), file(r1),file(r2) from samplefiles
    val trim

    output:
    tuple val(name), file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into galore
    file "*_report.txt" into trimmedlog
    
    script:
    def trimparam= trim>0 ? "--clip_R1 $trim --clip_R2 $trim" : ""
    """
    trim_galore -j ${task.cpus} --phred33 $trimparam \
      --paired \
      $r1 $r2
    """
}

// Mapping reads with bwa mem
process mapping {
    publishDir "${results}/bams", mode: 'link'

    label 'bwa'
    tag "${name}"

    input:
    tuple val(name), file(r1), file(r2) from galore
    tuple val(indexName), file(indexFiles) from genomeIndex
    val reference

    output:
    tuple val(name), file("*.bam"), file("*.bam.bai") into mapped

    script:
    """
     bwa mem -R "@RG\\tID:${name}\\tLB:lib_${name}\\tSM:${name}\\tPL:Illumina" -t ${task.cpus} ${indexName} $r1 $r2 \
       | samtools sort -m ${Math.max(5,task.memory.toGiga()-15)}G -o ${name}.bam
     samtools index ${name}.bam
     samtools flagstat ${name}.bam > ${name}_info.log
    """
}

// Counting nucleotides at each genome position
process countNucleotides {
    label 'bamreadcount'
    
    input:
    tuple val(name), val(bam), val(bai) from mapped
    file reference

    output:
    tuple val(name), file("${name}_readcounts.txt") into readcounts1, readcounts2, readcounts3

    script:
    """
    bam-readcount -q 20 -f ${reference} ${bam} -w 1 > ${name}_readcounts.txt
    """
}

// Extract positions that are specific of clade1 or clade2
// Using Nextstrain clade definition
process extractSpecificMutations {
    label 'python'

    input:
    val clade1
    val clade2
    file virusprop
    file reference

    output:
    file "specific_mutations.tsv" into specificmutations

    when:
    clade1!="-1" && clade2!="-1"

    script:
    """
    specific_mutations.py -v ${virusprop} -c "$clade1" -d "$clade2" -r $reference -o specific_mutations.tsv
    """
}

// Extract positions of all mutations defining all clades
// Using Nextstrain clade definition file
process extractCladePositions {
    label 'python'
    
    input:
    file virusprop
    
    output:
    file "*mutations.tsv" into lineagepositions mode flatten
    
    script:
    """
    lineage_positions.py -v ${virusprop} -o mutations.tsv
    """
}

// 
process computeFrequencies {

    label 'python'

    input:
    tuple val(name), file(readcounts), file(muts) from readcounts3.combine(lineagepositions)

    output:
    tuple val(name), file("${name}_freq.txt") into frequencies1, frequencies2

    script:
    """
    computefreq.py -i ${readcounts} -c ${muts} -o freq_tmp.txt
    awk '{if(NR>1){print "$name\\t" \$0}}' freq_tmp.txt > ${name}_freq.txt
    """
}

process plotMutations {
    publishDir "${results}/plots/", mode: 'link'

    label 'r'

    input:
    file freq from frequencies1.collectFile(){n,f -> [n+"_freq.txt",f]}

    output:
    file "Coinfection_${freq.baseName}*.svg" optional true

    script:
    """
    #!/usr/bin/env Rscript
    library(ggplot2)
    library(reshape2)
    palette=c("#fec00f","#0e72b2", "#1dac4b")

    rcount=read.table("${freq}",header=F)
    colnames(rcount)=c("Sample","Position","Count","Lineage","Total")
    rcount\$Frequency=(rcount[,3])/(rcount\$Total)
    poscount=table(unique(rcount[,c("Position","Lineage")])\$Position)
    rcount\$Specific=FALSE
    rcount[rcount\$Position%in%row.names(poscount[poscount==1]),"Specific"]=TRUE
    rcount\$HighCov=factor(rcount\$Total>100, levels=c(TRUE,FALSE))

    if(length(rcount\$Total)==0){
	quit(save="no",status=0)
    }
    svg("Coinfection_${freq.baseName}_clades.svg",width = 20,height=10)
    ggplot(rcount,aes(x=as.factor(Position), y=Frequency,fill=as.factor(Specific),alpha=HighCov,color=as.factor(Specific)))+
    geom_bar(stat="identity",size = 0)+
    scale_color_manual(values = palette) + 
    scale_fill_manual(values = palette) +
    scale_alpha_manual(values=c(1, 0.2)) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=4), axis.text.y= element_text(size=5))+
    facet_grid(Lineage~Sample)
    dev.off()
    """
}

process computePairFrequencies {

    label 'python'

    input:
    tuple val(name), file(readcounts) from readcounts1
    file sm from specificmutations

    output:
    tuple val(name), file("${name}_freq.txt") into pairfrequencies1, pairfrequencies2

    script:
    """
    computepairfreq.py -i ${readcounts} -c ${sm} -o freq_tmp.txt
    awk '{if(NR==1){print "Sample\\t" \$0}else{print "$name\\t" \$0}}' freq_tmp.txt > ${name}_freq.txt
    """
}

process plotPairMutations {
    publishDir "${results}/plots", mode: 'copy'

    label 'r'

    input:
    tuple val(name), file(freq) from pairfrequencies1

    output:
    file "Cladeseq_pair_${freq.baseName}.svg"

    script:
    """
    #!/usr/bin/env Rscript
    library(ggplot2)
    library(reshape2)
    palette=c("#fec00f","#0e72b2", "#1dac4b")

    rcount=read.table("${freq}",header=T)
    #rcount\$clade1=(rcount[,3])/(rcount[,3]+rcount[,4])
    #rcount\$clade2=(rcount[,4])/(rcount[,3]+rcount[,4])
    rcount\$clade1=(rcount[,3])/(rcount\$Total)
    rcount\$clade2=(rcount[,4])/(rcount\$Total)
    #rcount\$clade1=(rcount[,3])
    #rcount\$clade2=(rcount[,4])
    rcount=rcount[rcount\$Total>100,]

    agg=melt(rcount, id=c("Sample","Position"),measure.vars = c("clade1","clade2"), value.name = "Frequency", variable.name = "VirusStrain")
    levels(agg\$VirusStrain)[levels(agg\$VirusStrain)=="clade1"] <- colnames(rcount)[3]
    levels(agg\$VirusStrain)[levels(agg\$VirusStrain)=="clade2"] <- colnames(rcount)[4]

    svg("Cladeseq_pair_${freq.baseName}.svg",width = 12,height=6)
    ggplot(agg,aes(x=as.factor(Position), y=Frequency, fill=VirusStrain))+
              geom_bar(stat="identity",color="black")+
              scale_color_manual(values = palette) + 
              scale_fill_manual(values = palette) +
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
              facet_grid(Sample~.)
    dev.off()
    """
}

pairfrequencies2.map{n,f->f}.collectFile(name: 'all_pair_freqs.txt',skip:1, keepHeader: true).subscribe{
    	file -> file.copyTo(results)
}
frequencies2.map{n,f ->f}.collectFile(name: 'all_freqs.txt',skip:1, keepHeader: true).subscribe{
    	file -> file.copyTo(results)
}
