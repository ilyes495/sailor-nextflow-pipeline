/*
 * Copyright (c) 2020, Seqera Labs.
 * Copyright (c) 2017-2019, Centre for Genomic Regulation (CRG).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

/*
 * SAILOR Nextflow Pipeline
 *
 * A Nextflow pipeline for RNA editing analysis using SAILOR methodology
 * Integrates RNA-seq mapping, editing detection, and differential expression analysis
 *
 * Based on SAILOR: https://github.com/YeoLab/SAILOR
 * Adapted for Nextflow workflow management
 */

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */

// Default parameters - modify these or provide via command line
params.genome_dir = null
params.genome = null
params.genome_index = null
params.genomeDir = null
params.gtf = null
params.known_snp = null
params.variants = null
params.variants_index = null

// Analysis parameters
params.genome_build = 'hg38'
params.editor = 'APOBEC'
params.input_dir = null
params.output_dir = "./results"

// Sample sheets
params.deseq2_test_sheet = "$projectDir/assets/sample_sheet_deseq2.csv"
params.sample_sheet = "$projectDir/assets/sample_sheet.csv"

// Workflow flags
params.workflow_DE = "DESEQ2"
params.workflow_SAILOR = "SAILOR"

// Thresholds
params.diff_threshold_deseq2 = 1
params.padj_threshold = 0.05
params.project_name = "SAILOR_analysis"



log.info """\
S A I L O R   N E X T F L O W   P I P E L I N E
=============================================
  projectDir   : ${projectDir}
  genome       : ${params.genome}
  genome_index : ${params.genome_index}
  genomeDir    : ${params.genomeDir}
  variants     : ${params.variants}
  gtf          : ${params.gtf}
  known_snp    : ${params.known_snp}
  input_dir    : ${params.input_dir}
  results      : ${params.output_dir}
  genome_build : ${params.genome_build}
  editor       : ${params.editor}
  project_name : ${params.project_name}

"""

/*
 * Import modules
 */
include {
  CLEAN_FASTQ;
  RNASEQ_MAPPING_STAR;
  RUN_QUALIMAP;
  RUN_MULTIQC;
  MERGE_BAMS;
  } from './modules.nf'

include { DESEQ2 }     from './workflows/DESEQ2.nf'
include { SAILOR }     from './workflows/SAILOR.nf'


/*
 * main pipeline logic
 */


def group_per_sample = { channel ->
  channel.groupTuple(by: [0])
}

workflow {
      
      // ** -- Read the sample sheet
      deseq_ch = join_csv(file(params.deseq2_test_sheet))
      // deseq_ch.view()
      
      // ** -- Read the sample names sheet 
      samplesheet_ch = sample_names(file(params.sample_sheet, checkIfExists: true))
      // samplesheet_ch.view()
      
      // PART 1: STAR RNA-Seq Mapping

      CLEAN_FASTQ(samplesheet_ch)
      // CLEAN_FASTQ.out.view()

      fastp_map_ch = CLEAN_FASTQ.out
            .map { id, sample_lane, file_name, files1, files2,files3, files4  ->
            return [id, sample_lane,file_name, files1, files2]}

      rnaseq_map_ch = RNASEQ_MAPPING_STAR(
            params.genome,
            params.genomeDir,
            fastp_map_ch)

      qualimap_ch = RUN_QUALIMAP(
        rnaseq_map_ch)
      

      // // // PART 2: Run MultiQC on all qualimap results
      // RUN_MULTIQC(qualimap_ch.collect { it[3] })
  
      // // Merge BAM files from different lanes
      merge_ch_input = rnaseq_map_ch.groupTuple(by: [0]).map { s ->
        [s[0],s[3]]}
      

      // merge_ch_input.view { "merge_ch_input: ${it}"}

      MERGE_BAMS(merge_ch_input)
      merged_bams_ch = MERGE_BAMS.out

      // ** -- Run Differential Expression Analysis
      // If DESEQ2 is selected, run the DESEQ2 workflow
      if (params.workflow_DE == "DESEQ2") {
        DESEQ2(merged_bams_ch,deseq_ch)
      }

      
      // // *-------- BASE EDITING SCRIPTS WORKFLOW -------* //
      // if (params.workflow_HY == "HYPERTRIBE") {
      //   HYPERTRIBE(merged_bams_ch,deseq_ch)
      // }

      // // // // *-------- SAILOR WORKFLOW -------* //
      if (params.workflow_SAILOR == "SAILOR") {
        SAILOR(merged_bams_ch)
      }
      
  }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def join_csv(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
           def condition  = row['Condition']
           def control = row['Control']
           return [condition,control]
        }.map{condition, control ->
          def meta = [
                   "cond" : condition,
                   "ctrl" : control 
          ]
          return[meta]}
         

}


def sample_names(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
           def group_name  = row['Group']
           def sample_ID = row['Sample']
           def file_name = row['File_name']
           def fastq_1 = "${params.input_dir}${row['Folder_name']}/${row['fastq_1']}"  // Combine input directory 
           def fastq_2 = "${params.input_dir}${row['Folder_name']}/${row['fastq_2']}" 
           return [group_name,sample_ID,file_name,fastq_1,fastq_2]
        }
         
}