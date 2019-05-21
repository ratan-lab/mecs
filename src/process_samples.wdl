import "process_sample.wdl" as sl

workflow process_samples {
  File inputs
  String scriptdir

  String bwa
  String reference   
  String sambamba
  String samtools
  String tmp_dir
  File target
  String vep
  String vepdata
  String assembly
 
  String? bwa_options
  String bwaopts = select_first([bwa_options, "-Y"])

  Int? mem_sambamba
  Int memsambamba = select_first([mem_sambamba, 6])
  
  Int? num_threads
  Int threads = select_first([num_threads, 1])

  Int? first_keep
  Int firstkeep = select_first([first_keep, 31])

  Int? min_reads_per_umi
  Int minreadsperumi = select_first([min_reads_per_umi, 10])

  Int? min_coverage_per_umi
  Int mincoverageperumi = select_first([min_coverage_per_umi, 5])

  Float? min_agreement
  Float minagreement  = select_first([min_agreement, 0.9])

  Float? max_Ns_in_consensus
  Float maxNsinconsensus = select_first([max_Ns_in_consensus, 0.1])

  Int? max_template_length
  Int maxtemplatelength = select_first([max_template_length, 300])

  Int? min_coverage_per_loci
  Int mincoverageperloci = select_first([min_coverage_per_loci, 500])

  Float? min_vaf_as_germline
  Float minvafasgermline = select_first([min_vaf_as_germline, 0.35])

  Array[Array[File]] samples = read_tsv(inputs)
  
  scatter(sample in samples) {
    call sl.process_sample {
      input:
        sample = sample[0],
        lib = sample[1],
        fq = sample[2],
        mq = sample[3],
        bwa = bwa,
        reference = reference,
        sambamba = sambamba,
        samtools = samtools,
        target = target,
        tmp_dir = tmp_dir,
        scriptdir = scriptdir,
        num_threads = threads,
        first_keep = firstkeep,
        min_reads_per_umi = minreadsperumi,
        min_coverage_per_umi = mincoverageperumi,
        min_agreement = minagreement,
        max_Ns_in_consensus = maxNsinconsensus,
        max_template_length = maxtemplatelength
    }
  }

  call identify_muts {
    input:
      samtools = samtools,
      reference = reference,
      target = target,
      scriptdir = scriptdir,
      vep = vep,
      vepdata = vepdata,
      assembly = assembly,
      bams = process_sample.output_bam,
      bais = process_sample.output_bai,
      min_coverage_per_loci = mincoverageperloci,
      min_vaf_as_germline = minvafasgermline
  }

  output {
    File variants = identify_muts.variants
  }
}

task identify_muts {
  Array[File] bams
  Array[File] bais
  String scriptdir
  String samtools
  String reference
  File target
  String vep
  String vepdata
  String assembly
    
  Int min_coverage_per_loci
  Float min_vaf_as_germline

  String ftb = "$"

  command <<<
    bamnames=$(for f in ${sep=' ' bams}; do name=$(basename ${ftb}{f}); echo ${ftb}{name%".consensus.bam"}; done) 

    ${samtools} mpileup -BQ0 -d 10000000000000 -f ${reference} -l ${target} \
      ${sep=' ' bams} \
    | ${scriptdir}/call_variants -c ${min_coverage_per_loci} -g ${min_vaf_as_germline} ${ftb}{bamnames} \
    > variants.txt

    cat variants.txt \
    | ${scriptdir}/prepare_for_vep > variants.vep.input.txt

    perl ${vep} -i variants.vep.input.txt --dir ${vepdata} \
       --assembly ${assembly} --offline --symbol --canonical 

    cat variant_effect_output.txt| ${scriptdir}/filter_vep_results > variants.ann.txt
  >>>

  output {
    File variants = "variants.ann.txt"
  }
}
