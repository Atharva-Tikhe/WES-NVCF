log.info """
    PREPROCESSING PIPELINE
    raw data output - $params.raw_data
    fastqc output - $params.fastqc_out
    repaired files - $params.repaired
"""
process fastqc {
    tag "Running fastqc on sample - $pair_id"
    publishDir params.fastqc_out, failOnError: false
    
    input:
    tuple val(pair_id), path(files)
    
    output:

    path "${pair_id}_{1,2}_fastqc.{html,zip}", emit: files

    script:
    """fastqc -q $files"""

}


process multiqc {
    tag "Running multiqc - $files"
    publishDir params.multiqc_out, overwrite: false, failOnError: false

    input:
    path files 

    output:

    path "multiqc_report.html"

    script:
    """
    multiqc $files
    """
    

}

process trimming {
  tag "Running trim galore - $files"
  publishDir params.trimming_out, failOnError: false 
  
  input:
  tuple val(pair_id), path(files)

  output:
  
  path "${pair_id}_{1,2}_trimmed.fq.gz"
  path "${pair_id}_{1,2}.fq.gz_trimming_report.txt"
  path "${pair_id}_{1,2}_trimmed_fastqc.{html,zip}", emit: trimfqc

  script:
  """
  trim_galore --fastqc --cores 5 -q 25 $files 
  """
}


process trim_multiqc {
    tag "Running multiqc - $files"
    publishDir params.multiqc_out, overwrite: false, failOnError: false

    input:
    
    path files 

    output:
    path "trim_multiqc_report.html"

    script:
    """
    multiqc $files
    mv multiqc_report.html trim_multiqc_report.html
    """

}

process repair_reads {
    tag "Repairing reads for $pair_id"
    publishDir params.repaired, overwrite:true, failOnError: false
    
    input:
    tuple val(pair_id), path(files)

    output:
    
    path("${files[0]}_fixed.fq")
    path("${files[1]}_fixed.fq")

    script:
    """
    cp -r $baseDir/tools/bbmap ./
    ./bbmap/repair.sh in1=${files[0]} in2=${files[1]} out1=${files[0]}_fixed.fq out2=${files[1]}_fixed.fq outsingle=single.fq
    rm -r ./bbmap
    """


}

process alignment {
    tag "Alignment on $files"
    publishDir params.align_out, failOnError: false
    input:

    tuple val(pair_id), path(files)

    output:
    // stdout
    
    path "${pair_id}.sam", emit: samfile
    
    shell:
    """
    bwa mem $params.genome_prefix ${files[0]} ${files[1]} -o ${pair_id}.sam 
    """
    
}


process convertToBam {
    publishDir params.align_out, failOnError: false
    input:

    path sam
    
    output:
    
    path "${sam}.bam", emit: firstBam
    
    script:
    """
    samtools view -Sb $sam > ${sam}.bam
    """
}

process sortBamWithIndex {
    publishDir params.align_out, failOnError: false
    input:

    path bam
    output:
    
    path "${bam}-sorted.bam", emit: sortedBam
    path "${bam}-sorted.bam.csi"

    script:
    """
    samtools sort --write-index -O BAM -o ${bam}-sorted.bam $bam 
    """
}

include { addReadGroups; markDups; removeDups; BQSR; applyBQSR; renameFinal } from './post_alignment.nf'


workflow {
    // do initial raw data QC - generate report using multiqc
    // fastqc.output.signal is a mock dependency (check state dependency in docs)
    raw_ch = channel.fromFilePairs(params.raw_data)
    fastqc(raw_ch)
    // fastqc_out_ch = channel.fromPath(params.fastqc_out, checkIfExists: true)
    // multiqc(fastqc.output.signal, fastqc_out_ch)
    multiqc(fastqc.output.files)
    // trimming
    trim_ch = channel.fromFilePairs(params.raw_data)
    trimming(trim_ch)
    
    trim_multiqc(trimming.output.trimfqc)
    
    
    // repairing reads
    trimmed_ch = channel.fromFilePairs(params.trimmed)
    repair_reads(trimmed_ch)
    repaired_ch = channel.fromFilePairs(params.repaired_files)
    alignment(repaired_ch)
    convertToBam(alignment.output.samfile)
    sortBamWithIndex( convertToBam.output.firstBam)

    
    // preprocessing the reads (bams -> add RGs -> remove Dups -> BQSR -> rename final file)

    addReadGroups(sortBamWithIndex.output.sortedBam)
    markDups(addReadGroups.output.rg_bam)
    removeDups(markDups.output.markedBam)


    genome = channel.fromPath(params.genome_prefix)
    known_sites = channel.fromPath(params.known_sites)
    BQSR(known_sites,removeDups.output.remDups)
    applyBQSR(BQSR.output.recalTable, removeDups.output.remDups)


}
