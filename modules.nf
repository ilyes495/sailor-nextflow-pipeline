
/*
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */

process PREPARE_GENOME_SAMTOOLS { 
  tag "$genome.baseName"
 
  input: 
    path genome_dir
    path genome
 
  output: 
    path "${genome}.fai"
  
  publishDir "${genome_dir}/", mode: 'copy'

  script:
  """
  samtools faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process PREPARE_GENOME_PICARD {
  tag "$genome.baseName"
  label 'mem_xlarge'
  
  

  input:
    path genome_dir
    path genome
  output:
    path "${genome.baseName}.dict"

  publishDir "${genome_dir}/", mode: 'copy'

  script:
  """
  gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
  """
}


/*
 * Process 1C: Create STAR genome index file.
 */

process PREPARE_STAR_GENOME_INDEX {
  tag "$genome.baseName"

  

  input:
    path genome_dir
    path genome
  output:
    path "${genome.baseName}_index"

  publishDir "${genome_dir}/", mode: 'copy'
  script:
  """  
  mkdir ${genome.baseName}_index

  STAR --runMode genomeGenerate \
       --genomeDir ${genome.baseName}_index \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}

process PREPARE_VCF_FILE {
  tag "$variantsFile.baseName"

  input: 
    path genome_dir
    path variantsFile
    path denylisted

  output:
    tuple \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
  publishDir "${genome_dir}/", mode: 'copy'
  
  script:  
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${denylisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}



process CLEAN_FASTQ{
    tag "$sample"
    label "mem_xlarge"
    
    cpus = 8
    memory = 20.GB
    time = '48h'
    
    input:
        tuple val(sample), val(lane), val(file_name), path(read1), path(read2)
    
    output:
        tuple \
            val(sample), \
            val(lane), \
            val(file_name), \
            path("${lane}_${file_name}_R1_cleaned.fastq.gz"), \
            path("${lane}_${file_name}_R2_cleaned.fastq.gz"), \
            path("${lane}_${file_name}.json"), \
            path("${lane}_${file_name}.html")

    publishDir "${params.output_dir}/fastp", pattern: "*.{json,html}", mode: 'copy'
    
    script:
    """
    fastp \
        --in1 $read1 \
        --in2 $read2 \
        --out1 ${lane}_${file_name}_R1_cleaned.fastq.gz \
        --out2 ${lane}_${file_name}_R2_cleaned.fastq.gz \
        --json ${lane}_${file_name}.json \
        --html ${lane}_${file_name}.html \
        --report_title ${lane}_${file_name} \
        --thread $task.cpus
    
    """
}

 /*
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process RNASEQ_MAPPING_STAR {
  
  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  
  cpus = 8
  memory = 65.GB
  time = '48h'
  scratch = false

  tag "${lane}_${file_name}"

  input:
    path genome
    path genomeDir
    tuple val(sample), val(lane), val(file_name), path(read1), path(read2)

  output:
    tuple \
      val(sample), \
      val(lane), \
      val(file_name), \
      path("${lane}_${file_name}.sortedByCoord.out.bam"), \
      path("${lane}_${file_name}.sortedByCoord.out.bam.bai")
      
      

  script:
  """
   STAR \
    --runThreadN 16 \
    --genomeDir $genomeDir \
    --readFilesIn $read1 $read2 \
    --readFilesCommand zcat \
    --outSAMattributes NH HI AS nM NM MD \
    --outSAMattrRGline ID:${lane}_${file_name} LB:library PL:illumina PU:${lane}_${file_name} SM:${lane}_${file_name} \
    --outFileNamePrefix ${lane}_${file_name}. \
    --genomeLoad NoSharedMemory \
    --outSAMmode Full \
    --outStd Log \
    --outBAMsortingBinsN 60 \
    --limitBAMsortRAM 31000000000


    sleep 5m

    # Sort BAM file
    samtools sort -@ 8 -o ${lane}_${file_name}.sortedByCoord.out.bam ${lane}_${file_name}.Aligned.out.sam
    samtools index ${lane}_${file_name}.sortedByCoord.out.bam

    rm ${lane}_${file_name}.Aligned.out.sam

  """  
  // # Run Qualimap
  //   #qualimap bamqc -bam ${lane}_${file_name}.sortedByCoord.out.bam -gtf $genome -outdir qualimap_${lane}_${file_name}
}

process RUN_QUALIMAP {
  tag 'qualimap'

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

  cpus = 8
  memory = 20.GB
  time = '2h'
  scratch true

  publishDir "${params.output_dir}/qualimap", mode: 'copy'

  input:
  tuple \
    val(sample), \
    val(lane), \
    val(file_name), \
    path(bam), \
    path(bai)

  output:
  tuple \
    val(sample), \
    val(lane), \
    val(file_name), \
    path("qualimap_${lane}_${file_name}")

  script:
  """

  qualimap bamqc --java-mem-size=20G -bam $bam -outdir qualimap_${lane}_${file_name}
  """
}

process RUN_MULTIQC {
  tag 'multiqc'

  publishDir "${params.output_dir}/multiqc", mode: 'copy'

  input:
  path qualimap_results

  output:
  path '*'

  script:
  """
  multiqc ${qualimap_results.join(' ')} -o .
  """
}

process MERGE_BAMS {
    tag "$sample"

    errorStrategy 'ignore' 

    scratch '/scratch/kharasm/baalii'

    cpus = 16
    memory = 80.GB
    time = '48h'

    input:
    tuple val(sample), path(bams)

    output:
    tuple val(sample), path("${sample}_merged.bam"), path("${sample}_merged.bam.bai")

    publishDir "${params.output_dir}/bam_files/${sample}", mode: 'copy'

    script:
    """
    # Merge BAM files from different lanes
    samtools merge -@ 16 -f ${sample}_merged.bam ${bams.join(' ')}

    # Index the merged BAM file
    samtools index ${sample}_merged.bam
    """
}


 /*
 * Process 6: Process to format the results
 *            the process runs an R script that gathers the files from 3_variant step
 */


process FORMAT_OUTPUT {

  cpus = 8
  memory = 20.GB
  time = '20h'


  publishDir "${params.output_dir}/hypertribe/", mode: 'copy'

  input:
    val gtf_path
    tuple val(samples), path(vcf_path_list, stageAs: "?/*"), path(bam_path_list, stageAs: "?/*"), path(bai_path_list, stageAs: "?/*")
  


  output:
    path '4_formatted_output.csv'
  //path '*.rds'

  script:
  """
  echo \$(pwd)
  echo ${samples}
  Rscript $baseDir/scripts/base_editing/4_format_output.R \
  --working_dir \$(pwd) \
  --gtf_path $gtf_path \
  --samples $samples \
  --bam_path_list ${bam_path_list.join(',')} \
  --vcf_path_list ${vcf_path_list.join(',')} \
  --genome_build ${params.genome_build} \
  --editor ${params.editor} ${workflow.resume ? '--resume' : ''}
  """
}

/*
 * Process 7: Base editing stats and comparisons 
 *            the process runs a R script that does differential analysis based on sample sheet
 */

process BE_STATS {

  cpus = 8
  memory = 20.GB
  time = '20h'

  publishDir "${params.output_dir}/hypertribe/DES/", mode: 'copy'
  tag "BE_STATS"

  input:
  tuple val(meta), path (formatted_file)
  

  output:
  tuple val(meta),path ('*csv')

  script:
  def ctrl      = "${meta.ctrl}"
  def condition = "${meta.cond}"
  """
  Rscript $baseDir/scripts/base_editing/5_test_differential.R \
     --formatted_file ${formatted_file} \
     --cond ${condition} \
     --ctrl ${ctrl}
  """

}

/* 
* PART 2: RNA-seq expression MODULES 
*/

/*
 * Process 1 : RNA seq counts from merged bams 
 *                 using htseq counts         
 */

process RNASEQ_COUNT {
  tag "$sample"
  label "mem_xlarge"

  cpus = 8
  memory = 40.GB
  time = '48h'

  publishDir "${params.output_dir}/deseq2/count/", mode: 'copy'

  input:
    path gtf_path
    tuple val(sample), path(bam), path(bai)

  output:
    tuple val(sample), path('*counts')

  script:
  """
  htseq-count \
  -f bam \
  -r pos \
  -s reverse \
  -a 10 \
  -i gene_id \
  -m union \
  -n ${task.cpus} \
  ${bam} ${gtf_path} > "${sample}.counts"
  """
}
process RNASEQ_COUNT_FEATURECOUNTS {
  tag "${feature}_${sample}"
  label "mem_xlarge"

  cpus = 8
  memory = 80.GB
  time = '48h'

  publishDir "${params.output_dir}/deseq2/featurecounts/", mode: 'copy'

  input:
    path gtf_path
    val(feature)
    path(bams)


  output:
    tuple path('*counts')

  script:
  """
  featureCounts -p --countReadPairs \
  -t ${feature} \
  -a ${gtf_path} \
  -o ${feature}.counts \
  ${bams.join(' ')}
  """
}

/*
 * Process 2 : Gathering all the sample counts from htseq  
 *             using R script into a single file      
 */

process FORMAT_OUTPUT_COUNTS {

  cpus = 8
  memory = 20.GB
  time = '24h'

  

  publishDir "${params.output_dir}/deseq2/", mode: 'copy'
  tag "FORMAT_OUTPUT_COUNTS"

  input:
    val gtf_path
    tuple val(samples), path(counts_path_list, stageAs: "?/*")
  

  output:
    path '4_formatted_output_counts.csv'
  //path '*.rds'

  script:
  """

  echo \$(pwd)

  Rscript $baseDir/nextflow_scripts/deseq2/4_format_output.R \
  --gtf_path $gtf_path \
  --samples ${samples.join(',')} \
  --counts_path_list ${counts_path_list.join(',')} \
  --genome_build ${params.genome_build} 
  """

}

/*
 * Process 3 : Differential gene expression using DESeq2     
 */

process DESEQ2_RNA_SEQ {

  cpus = 8
  memory = 20.GB
  time = '12h'

  publishDir "${params.output_dir}/deseq2/DEG/", mode: 'copy'
  tag "DESEQ2"

  input:
  tuple val(meta), path (counts_file)

  output:
  tuple val(meta), path ('*csv')

  script:
  def ctrl      = "${meta.ctrl}"
  def condition = "${meta.cond}"
  """
  Rscript $baseDir/nextflow_scripts/deseq2/5_test_differential.R \
     --counts_file ${counts_file} \
     --cond ${condition} \
     --ctrl ${ctrl}
  """
}




process PLOTS_RNA_SEQ {

  publishDir "${params.output_dir}/deseq2/plots/", mode: 'copy'
  tag "PLOTS"

  cpus = 8
  memory = 20.GB
  time = '12h'

  input:
    val (threshold)
    tuple val(meta), path (deseq_file)
    val (formatted_file)
    val (analysis)

  output:
    tuple val(meta), path ('*jpeg') 

  script:
  def ctrl      = "${meta.ctrl}"
  def condition = "${meta.cond}"
  """
  Rscript $baseDir/nextflow_scripts/plots_pipeline.R \
     --deseq_file ${deseq_file} \
     --formatted_file ${formatted_file} \
     --cond ${condition} \
     --ctrl ${ctrl} \
     --analysis ${analysis} \
     --padj_threshold ${params.padj_threshold} \
     --diff_threshold ${threshold} \
     --selection_name "NULL" \
     --project_name ${params.project_name}
  """
}
