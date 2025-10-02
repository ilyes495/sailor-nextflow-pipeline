include { SPLIT_STRANDS; 
PROCESS_STRAND as fwd_PROCESS_STRAND;
PROCESS_STRAND as rev_PROCESS_STRAND;
MERGE_BEDS;
ANNOTATE_BED;
} from './SAILOR_modules.nf'

workflow SAILOR {

  take:
    merge_bams_ch

  main:
    
    // Split strands
    split_ch = SPLIT_STRANDS(merge_bams_ch)

    fwd_ch = split_ch.map { sample, fwd_bam, rev_bam -> tuple(sample, fwd_bam) }
    rev_ch = split_ch.map { sample, fwd_bam, rev_bam -> tuple(sample, rev_bam) }

    // Process forward and reverse strands
    fwd_results = fwd_PROCESS_STRAND(fwd_ch, false)
    rev_results = rev_PROCESS_STRAND(rev_ch, true)

    // Merge results
    merged_results = fwd_results.combine(rev_results, by: 0).map { it -> 
      tuple(it[0], it[2], it[4])
    }
    // merged_results.view()
    
    merge_beds_ch = MERGE_BEDS(merged_results)

    // Annotate bed files
    ANNOTATE_BED(merge_beds_ch)



}

