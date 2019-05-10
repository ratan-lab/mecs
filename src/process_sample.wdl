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

task align_sort_index {
  String bwa
  String bwaopts
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

  Int first_keep=31

  command <<<
    set -e 

    mkdir -p ${tmp_dir}

    ${bwa} mem -t ${threads} ${bwaopts} -R "${readgroup}" \
      -K 100000000 ${reference} \
      <(gzip -dc ${fq} |  awk -v n=${first_keep} '{if((NR-1)%4==0) {split($2,a,":"); printf "%s:%s\n", $1, a[4]}  else if(((NR-2)%4==0) || ((NR-4)%4==0)) {printf "%s\n", substr($1,n)} else {print $1}}') \
      <(gzip -dc ${mq} | awk -v n=${first_keep} '{if((NR-1)%4==0) {split($2,a,":"); printf "%s:%s\n", $1, a[4]}  else if(((NR-2)%4==0) || ((NR-4)%4==0)) {printf "%s\n", substr($1,n)} else {print $1}}') \
    | ${sambamba} view -f bam -S -t ${threads} /dev/stdin \
    | ${samtools} fixmate -O BAM - ${tmp_dir}/${outprefix}.ns.alignment.bam 
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

  String tmpfile = outprefix + ".tmp"
  String outfile = outprefix + ".fa"

  command <<<
    ${scriptdir}/get_consensus_reads ${bamfile} ${target} ${tmpfile} ${threads} > ${outfile}
  >>>

  output {
    File consensus = "${outfile}"
  }
}

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

  call align_sort_index {
    input:
      bwa = bwa,
      bwaopts = bwaopts,
      fq = fq,
      mq = mq,
      memsambamba = memsambamba,
      outprefix = outprefix,
      readgroup = readgroup,
      reference = reference,
      sambamba = sambamba,
      threads = threads,
      tmp_dir = tmp_dir,
      samtools = samtools
  }

  call get_consensus_reads {
    input:
      scriptdir = scriptdir,
      outprefix = outprefix,
      target = target,
      threads = threads,
      bamfile = align_sort_index.alignments,
      bamfile_index = align_sort_index.alignments_index
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
