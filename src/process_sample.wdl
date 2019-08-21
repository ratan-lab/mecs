workflow process_sample {
  # information from the samplesheet
  String sample
  String lib
  String? lane_num
  String lane = select_first([lane_num, "1"])
  String? platform_name
  String platform = select_first([platform_name, "illumina"])

  String bwa
  String fq 
  String mq 
  String reference   
  String sambamba
  String samtools
  String pear
  String tmp_dir
  File target
  String scriptdir

  String? bwa_options
  String bwaopts = select_first([bwa_options, "-Y"])

  Int? mem_sambamba
  Int memsambamba = select_first([mem_sambamba, 6])
  
  String? out_prefix
  String outprefix = select_first([out_prefix, sample+"_"+lib])

  String? read_group
  String rg = "@RG\\tID:" + sample + "_" + lib + "_" + lane + "\\tSM:" + sample + "\\tPU:" + lane + "\\tPL:" + platform
  String readgroup = select_first([read_group, rg])

  Int? num_threads
  Int threads = select_first([num_threads, 1])
  
  Int? first_keep
  Int firstkeep = select_first([first_keep, 1])

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

  call merge_pairs {
    input:
      pear = pear,
      fq = fq,
      mq = mq,
      threads = threads,
      first_keep = firstkeep
  }

  call align_sort_index {
    input:
      bwa = bwa,
      bwaopts = bwaopts,
      merged = merge_pairs.merged,
      fq = merge_pairs.unmerged_ff,
      mq = merge_pairs.unmerged_rf,
      memsambamba = memsambamba,
      outprefix = outprefix,
      readgroup = readgroup,
      reference = reference,
      sambamba = sambamba,
      threads = threads,
      tmp_dir = tmp_dir,
      samtools = samtools,
      scriptdir = scriptdir
  }

  call get_consensus_reads {
    input:
      scriptdir = scriptdir,
      outprefix = outprefix,
      target = target,
      threads = threads,
      bamfile = align_sort_index.alignments,
      bamfile_index = align_sort_index.alignments_index,
      min_reads_per_umi = minreadsperumi,
      min_coverage_per_umi = mincoverageperumi,
      min_agreement = minagreement,
      max_Ns_in_consensus = maxNsinconsensus,
      max_template_length = maxtemplatelength
  }

  call align_sort_index_fa {
    input:
      bwa = bwa,
      bwaopts = bwaopts,
      fa = get_consensus_reads.consensus,
      memsambamba = memsambamba,
      outprefix = outprefix,
      readgroup = readgroup,
      reference = reference,
      sambamba = sambamba,
      threads = threads,
      tmp_dir = tmp_dir
  }

  output {
    File output_bam = align_sort_index_fa.alignments
    File output_bai = align_sort_index_fa.alignments_index
  }   
}

task merge_pairs {
  String pear
  String fq
  String mq
  Int threads
  Int first_keep

  command <<<
    gzip -dc ${fq} \
    | awk -v n=${first_keep} '{if((NR-1)%4==0) {split($2,a,":"); printf "%s:%s\n", $1, a[4]}  else if(((NR-2)%4==0) || ((NR-4)%4==0)) {printf "%s\n", substr($1,n)} else {print $1}}'\
    > read_1.fq

     gzip -dc ${mq} \
    | awk -v n=${first_keep} '{if((NR-1)%4==0) {split($2,a,":"); printf "%s:%s\n", $1, a[4]}  else if(((NR-2)%4==0) || ((NR-4)%4==0)) {printf "%s\n", substr($1,n)} else {print $1}}'\
    > read_2.fq

    ${pear} -f read_1.fq -r read_2.fq -j ${threads} -o merge -p 0.05

    rm read_1.fq read_2.fq
  >>>

  output {
    File merged = "merge.assembled.fastq"
    File unmerged_ff = "merge.unassembled.forward.fastq"
    File unmerged_rf = "merge.unassembled.reverse.fastq"
  }
}

task align_sort_index {
  String bwa
  String bwaopts
  String merged
  String fq
  String mq
  Int memsambamba
  String outprefix
  String readgroup
  String reference
  String sambamba
  String samtools
  Int threads
  String tmp_dir
  String scriptdir

  command <<<
    set -e 

    mkdir -p ${tmp_dir}

    ${bwa} mem -t ${threads} ${bwaopts} -R "${readgroup}" \
      -K 100000000 ${reference} ${merged} \
    | ${sambamba} view -f bam -S -t ${threads} /dev/stdin \
    | ${samtools} view -F 2820 -b \
      -o ${tmp_dir}/${outprefix}.ns.alignment.bam
    echo "alignment finished"

    # sort the alignments 
    ${sambamba} sort -m ${memsambamba}G --tmpdir=${tmp_dir} \
      -t ${threads} \
      -o ${outprefix}.alignment.bam \
      ${tmp_dir}/${outprefix}.ns.alignment.bam
    echo "mapped: coord sort finished"

    # index the alignments
    ${sambamba} index -t ${threads} \
      ${outprefix}.alignment.bam
    echo "mapped: index finished"

    # delete the temp directory
    rm ${tmp_dir}/${outprefix}.ns.alignment.bam
  >>>

  output {
    File alignments = "${outprefix}.alignment.bam"
    File alignments_index = "${outprefix}.alignment.bam.bai"
  }
}

task get_consensus_reads {
  String scriptdir
  String outprefix
  File target 
  File bamfile
  File bamfile_index
  Int threads
  Int min_reads_per_umi
  Int min_coverage_per_umi
  Float min_agreement
  Float max_Ns_in_consensus
  Int max_template_length

  String tmpfile = outprefix + ".tmp"
  String outfile = outprefix + ".fa"

  command <<<
    ${scriptdir}/get_consensus_reads -r ${min_reads_per_umi} -c ${min_coverage_per_umi} -a ${min_agreement} -n ${max_Ns_in_consensus} -t ${max_template_length} -p ${threads} ${bamfile} ${target} ${tmpfile} > ${outfile}
  >>>

  output {
    File consensus = "${outfile}"
  }
}

task align_sort_index_fa {
  String bwa
  String bwaopts
  String fa
  Int memsambamba
  String outprefix
  String readgroup
  String reference
  String sambamba
  Int threads
  String tmp_dir

  command <<<
    set -e 

    mkdir -p ${tmp_dir}

    ${bwa} mem -t ${threads} ${bwaopts} -R "${readgroup}" \
      -K 100000000 ${reference} ${fa} \
    | ${sambamba} view -f bam -S \
      -o ${tmp_dir}/${outprefix}.ns.consensus.bam -t ${threads} \
      /dev/stdin 
    echo "alignment finished"

    # sort the alignments 
    ${sambamba} sort -m ${memsambamba}G --tmpdir=${tmp_dir} \
      -t ${threads} \
      -o ${outprefix}.consensus.bam \
      ${tmp_dir}/${outprefix}.ns.consensus.bam
    echo "mapped: coord sort finished"

    # index the alignments
    ${sambamba} index -t ${threads} \
      ${outprefix}.consensus.bam
    echo "mapped: index finished"

    # delete the temp directory
    rm ${tmp_dir}/${outprefix}.ns.consensus.bam
  >>>

  output {
    File alignments = "${outprefix}.consensus.bam"
    File alignments_index = "${outprefix}.consensus.bam.bai"
  }
}

