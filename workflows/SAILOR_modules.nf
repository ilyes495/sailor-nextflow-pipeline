process SPLIT_STRANDS {
  tag "$sample"
  cpus = 8
  memory = 32.GB

  input:
  tuple \
    val(sample), \
    path(bam), \
    path(bai)

  output:
  tuple \
    val(sample), \
    path("${bam.baseName}.fwd.bam"), \
    path("${bam.baseName}.rev.bam")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  $projectDir/py/split_strands.py -i ${bam} --reverse-strand --output-forward ${bam.baseName}.fwd.bam --output-reverse ${bam.baseName}.rev.bam
  """
}

process SORT {
  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 64.GB : task.attempt == 2 ? 128.GB : 256.GB }
  tag "$sample"
  cpus = 16
  memory = 64.GB

  input:
  tuple val(sample), path(bam)

  output:
  tuple val(sample), path("${bam.baseName}.sorted.bam")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  samtools sort -@ $task.cpus -n -o namesort.bam ${bam}
  samtools fixmate -@ $task.cpus -m namesort.bam fixmate.bam
  samtools sort -@ $task.cpus -o ${bam.baseName}.sorted.bam fixmate.bam
  """
}

process RMDUP {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 64.GB : task.attempt == 2 ? 128.GB : 256.GB }
  tag "$sample"
  cpus = 16
  memory = 32.GB

  input:
  tuple val(sample), path(bam)

  output:
  tuple val(sample), path("${bam.baseName}.rmdup.bam")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  samtools markdup -@ $task.cpus  -r ${bam} ${bam.baseName}.rmdup.bam
  """
}

process FILTER_READS {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 80.GB : task.attempt == 2 ? 128.GB : 256.GB }
  tag "$sample"
  cpus = 8
  memory = 64.GB

  input:
  tuple val(sample), path(bam)

  output:
  tuple val(sample), path("${bam.baseName}.readfiltered.bam")

  publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  // Extract input params from a module config file
  def junction_overhang = params.junction_overhang ?: 10
  def edge_mutation = params.edge_mutation ?: 5
  def non_ag = params.non_ag ?: 3
  def reverse_stranded_library = params.reverse_stranded_library ?: false
  def ct = params.ct ?: false
  def gt = params.gt ?: false

  if (ct && gt) {
    error "Parameters 'ct' and 'gt' cannot both be true."
  }

  script:
  """
  $projectDir/py/filter_reads.py --input ${bam}\
   --junction_overhang ${junction_overhang} \
   --edge_mutation ${edge_mutation} \
   --non_ag_threshold ${non_ag} \
  ${reverse_stranded_library ? '--reverse-strand' : ''} \
  ${ct ? '--ct' : ''} ${gt ? '--gt' : ''} \
  --output ${bam.baseName}.readfiltered.bam
  """
}

process MPILEUP {
  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2
  memory { task.attempt == 1 ? 160.GB : task.attempt == 2 ? 200.GB : 256.GB }

  tag "$sample"
  cpus = 32
  memory = 32.GB
  time = '12h'
  queue = 'cpu'


  input:
  tuple val(sample), path(bam)

  output:
  tuple val(sample), path("${bam.baseName}.gbcf")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  bcftools mpileup --threads $task.cpus -d 1000 -E -f ${params.genome} -p -O u -I ${bam}  -o ${bam.baseName}.gbcf
  """
}

process CALL_SNVS {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 80.GB : task.attempt == 2 ? 128.GB : 256.GB }

  tag "$sample"
  cpus = 16
  memory = 64.GB
  time = '12h'


  input:
  tuple val(sample), path(mpileup)

  output:
  tuple val(sample), path("${mpileup.baseName}.vcf")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  bcftools view --threads $task.cpus -v snps -O v -o ${mpileup.baseName}.vcf ${mpileup}
  ls -lh
  exit 0
  """
}

process FORMAT_VARIANTS {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 64.GB : task.attempt == 2 ? 128.GB : 256.GB }

  tag "$sample"
  cpus = 8
  memory = 32.GB

  input:
  tuple val(sample), path(gbcf)

  output:
  tuple val(sample), path("${gbcf.baseName}.formatted.vcf")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  bcftools call --threads $task.cpus -c -O v -A -o ${gbcf.baseName}.formatted.vcf ${gbcf}
  """
}

process FILTER_VARIANTS {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 64.GB : task.attempt == 2 ? 128.GB : 256.GB }

  tag "$sample"
  cpus = 8
  memory = 32.GB

  input:
  val(is_reverse)
  tuple val(sample), path(vcf)

  output:
  tuple val(sample), path("${vcf.baseName}.varfiltered.vcf")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  
  def ct = params.ct ?: false
  def gt = params.gt ?: false
  def dp = params.dp ?: "DP4"
  def min_variant_coverage = params.min_variant_coverage ?: 5


  script:
  """
  echo "is_reverse: ${is_reverse}"
  $projectDir/py/filter_variants.py -i ${vcf} \
  -m ${min_variant_coverage} \
  --dp ${dp} \
  ${is_reverse ? '--reverse-split' : ''} \
  ${ct ? '--ct' : ''} ${gt ? '--gt' : ''} \
  -o ${vcf.baseName}.varfiltered.vcf
  """
}

process FILTER_KNOWN_SNP {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 64.GB : task.attempt == 2 ? 128.GB : 256.GB }

  tag "$sample"
  cpus = 8
  memory = 32.GB

  input:
  tuple val(sample), path(vcf)

  output:
  tuple val(sample), path("${vcf.baseName}.snpfiltered.vcf")

  // publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  $projectDir/py/filter_known_snp.py --input ${vcf} --known ${params.known_snp} --output ${vcf.baseName}.snpfiltered.vcf
  """
}

process RANK_EDITS {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 3
  memory { task.attempt == 1 ? 64.GB : task.attempt == 2 ? 128.GB : 256.GB }

  tag "$sample"
  cpus = 8
  memory = 32.GB

  input:
  tuple val(sample), path(vcf)

  output:
  tuple val(sample), path("${vcf.baseName}.ranked.conf"), path("${vcf.baseName}.ranked.bed")

  publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  $projectDir/py/rank_edits.py -i ${vcf} \
  -c ${params.edit_fraction} \
  -a ${params.alpha} \
  -b ${params.beta} \
  ${params.keep_all_edited ? '--keep-100-edited' : ''} \
  ${params.ct ? '--ct' : ''} ${params.gt ? '--gt' : ''} \
  --output ${vcf.baseName}.ranked.conf
  """
}

process MERGE_BEDS {

  tag "$sample"
  cpus = 4
  memory = 16.GB

  input:
    tuple val(sample), path(fwd_bed), path(rev_bed)

  output:
    tuple val(sample), path("${sample}.sorted.bed")
  
  publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  cat ${fwd_bed} ${rev_bed} > ${sample}.bed
  sort -k1,1 -k2,2n ${sample}.bed > ${sample}.sorted.bed
  """
}

process ANNOTATE_BED {

  tag "$sample"
  cpus = 4
  memory = 16.GB

  input:
  tuple val(sample), path(bed)

  output:
  tuple val(sample), path("${sample}_annotated_final.bed")
  
  publishDir "${params.output_dir}/sailor/${sample}", mode: 'copy'

  script:
  """
  
  echo "annotate bed file"
  bedtools intersect -a ${bed} -b ${params.gtf} -wa -wb > ${sample}.annotated.bed
  cut -f1,2,3,4,5,6,9,15 ${sample}.annotated.bed > ${sample}.annotated.filtered.bed
  awk -F'\\t' 'BEGIN {OFS=FS} {match(\$8, /gene_id "([^"]+)";/, Geneid); match(\$8, /gene_name "([^"]+)";/, Genename); print \$1, \$2, \$3, \$4, \$5, \$6, \$7, Geneid[1], Genename[1]}' ${sample}.annotated.filtered.bed > ${sample}_annotated_final.bed

  rm ${sample}.annotated.bed
  rm ${sample}.annotated.filtered.bed
  """
}

workflow PROCESS_STRAND {

  take:
    ch
    is_reverse

  main:
    // Link all subprocesses
    sorted_ch = SORT(ch)
    rmdup_ch = RMDUP(sorted_ch)
    filtered_ch = FILTER_READS(rmdup_ch)
    mpileup_ch = MPILEUP(filtered_ch)
    snvs_ch = CALL_SNVS(mpileup_ch)
    formatted_ch = FORMAT_VARIANTS(snvs_ch)
    varfiltered_ch = FILTER_VARIANTS(is_reverse, formatted_ch)
    snpfiltered_ch = FILTER_KNOWN_SNP(varfiltered_ch)
    ranked_ch = RANK_EDITS(snpfiltered_ch)
    

  emit:
    ranked_ch
}




