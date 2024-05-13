process addReadGroups {
    tag "Adding read group to $bam"
    input:
    path bam

    output:
    path "${bam}_rg.bam", emit: rg_bam

    script:
    """
    gatk AddOrReplaceReadGroups I=${bam} O=${bam}_rg.bam RGID=4 RGLB=lib1 RGPL=MGI RGPU=unit1 RGSM=20
    """
}

process markDups {
    tag "Marking duplicates in $bam"
    publishDir params.post_align_out, failOnError: false
    input: 
    path bam

    output:
    path "${bam}_deduped.bam", emit: markedBam
    path "dup_metrics.txt"
    
    
    script:
    """ 
    
    gatk MarkDuplicates -I ${bam} -O ${bam}_deduped.bam -M dup_metrics.txt --TAGGING_POLICY All
    
    """

}

process removeDups {
    tag "Removing duplicates in $bam"
    publishDir params.post_align_out, failOnError: false
    input:
    path bam

    output:
    path "${bam}_unduped.bam", emit: remDups
    path "dup_rem_metrics.txt"

    script:
    """
    
    gatk MarkDuplicates -I ${bam} -O ${bam}_unduped.bam -M dup_rem_metrics.txt --REMOVE_DUPLICATES true 
    
    """


}

process BQSR {
    tag "Running BQSR on $bam"

    input:
    
    // path genome
    path known_sites
    path bam

    output:

    path "recal.table.tsv", emit: recalTable

    script:
    """
    mkdir Genome
    ln -s $baseDir/required_datasets/Genome/* ./Genome/
    ln -s $baseDir/required_datasets/known_sites/common_all_20170710.vcf.gz.tbi ./
    gatk BaseRecalibrator -I ${bam} --known-sites ${known_sites} -O recal.table.tsv -R ./Genome/GRCh37.fa
    
    """
}

process applyBQSR {
    publishDir params.post_align_out, failOnError: false
    input:

    path recal_table
    path bam

    output:

    path "${bam}_bqsr.bam", emit: bqsrBam

    script:
    """
    
    gatk ApplyBQSR -I ${bam} --bqsr-recal-file ${recal_table} -O ${bam}_bqsr.bam
    
    """



}

