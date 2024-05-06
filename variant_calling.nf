process gatk_vc {
  tag "Running haplotypecaller on $bam"

  publishDir params.gatk_htc
  input:
    path bam

  output:
    path "${params.sample}.vcf", emit: vcf
    path "${params.sample}.vcf.idx"
  script:
  """
    mkdir Genome
    ln -s $baseDir/required_datasets/Genome/* ./Genome/
    gatk --java-options "-Xmx10g" HaplotypeCaller -I ${bam} -O ${params.sample}.vcf -R ./Genome/GRCh37.fa

  """
}

process variantRecal {
  tag "running variant recalibrator on $vcf"
  
  publishDir params.post_vc
  input:
  path vcf
  
  output:
  path "output.recal", emit: recal
  path "output.tranches", emit: tranches

  script:
  """ 
  mkdir Genome
  mkdir vqsr-resources
  ln -s $baseDir/required_datasets/vqsr-resources/* ./vqsr-resources/
  ln -s $baseDir/required_datasets/Genome/* ./Genome/

  gatk VariantRecalibrator \
  -R ./Genome/GRCh37.fa \
  -V ${vcf} \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./vqsr-resources/hapmap_3.3.b37.vcf.gz \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 ./vqsr-resources/1000G_omni2.5.b37.vcf.gz \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 ./vqsr-resources/1000G_phase1.snps.high_confidence.b37.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./vqsr-resources/dbsnp_138.b37.vcf.gz \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP \
  -O output.recal \
  --tranches-file output.tranches \
  --dont-run-rscript
  
  """
}

process indexRecal {
  publishDir params.post_vc
  input:
  path recal
  output:
  val ready, emit:signal
  path "${recal}.idx"
  script:
  """ 
  gatk IndexFeatureFile -I $recal 
  
  """
  """ 
  int signal = 0
  Println signal
  """
}

process applyVQSR {
  tag "Applying vqsr to $vcf"
  publishDir params.post_vc, pattern: "*.vcf"

  input:
  val signal
  path vcf
  path recal
  path tranches

  output:
  path "${vcf}_vqsr_out.vcf", emit: readyVCF
  
  when {success indexRecal} 

  script:
  """
  mkdir Genome
  ln -s $baseDir/required_datasets/Genome/* ./Genome/
  ln -s $baseDir/$params.sample_dir/post_vc/${recal}.idx ./

  gatk ApplyVQSR \
   -R ./Genome/GRCh37.fa \
   -V ${vcf} \
   -O ${vcf}_vqsr_out.vcf \
   -ts-filter-level 99.0 \
   --tranches-file $tranches \
   --recal-file $recal \
   -mode SNP
  """
}

process modifyConfig {
  
  input:
  path vcf
  script:
  """ 
  #!/bin/python3

  import configparser
  config = configparser.ConfigParser()
  config_filename = '$baseDir/tools/InterVar/config.ini'
  config.read(config_filename)
  print(config.has_section('InterVar'))
  try:
    config.set('InterVar', 'inputfile', '$params.post_vc/$vcf')
    config.set('InterVar', 'outfile', '${vcf}_annotated.vcf')
    print('updated filename')
  except configparser.NoSectionError:
    pass
  print('done')

  config.write(open(config_filename, 'w'))
  
  """
}


process annotateIntervar {
  tag "Annotating $vcf with InterVar"
  publishDir params.annotations
  input:
  path vcf
  output:
  path "InterVar/${vcf}_annotated.vcf.hg19_multianno.txt"
  path "InterVar/${vcf}_annotated.vcf.hg19_multianno.txt.grl_p"
  path "InterVar/${vcf}_annotated.vcf.hg19_multianno.txt.intervar"

  script:
  """
  ln -s $params.intervar/ ./ 
  cd InterVar
  python3 Intervar.py -c config.ini
  
  """
}

process renameFinal{
    script:
    template "$baseDir/helper_scripts/rename_bam.sh"
}

workflow {
    renameFinal()
    bam_ch = channel.fromPath(params.post_align_bam)
    gatk_vc(bam_ch)
    variantRecal(gatk_vc.output.vcf)
    indexRecal(variantRecal.output.recal)
    applyVQSR(indexRecal.output.signal, gatk_vc.output.vcf, variantRecal.output.recal, variantRecal.output.tranches)
    modifyConfig(applyVQSR.output.readyVCF)
    annotateIntervar(applyVQSR.output.readyVCF)
}