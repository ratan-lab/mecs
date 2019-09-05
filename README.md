# Mutation calling for Error-Corrected Sequencing

An implementation of the bioinformatics methods detailed at https://www.jove.com/video/57509/rare-event-detection-using-error-corrected-dna-and-rna-sequencing to call variants from error corrected DNA sequences.

Multiple samples are sequenced to high coverage using targeted capture, and paired-end fastq files are generated for each of them. The following analysis is done on each one of them to generate a BAM file for each of them. The final SNP calls are generated based on the error-profile that is based on all the alignment files.

For each sample:
* Trim off the first 30 nucleotides of each demultiplexed read to remove oligo sequences from the gene panel.
* Align the reads to the human reference genome using BWA
* Process the alignment to identify and align reads that share the same UMIs to one another to form read families.
* Perform de-duplication and error-correction using the following recommended parameters.
  * Use ≥5 read pairs in the same read family.
  * Compare nucleotide at every position across all reads in the same read family, and generate a consensus nucleotide if there is at least 90% concordance among the reads for the particular nucleotide. Call an N if there is less than 90% agreement for the nucleotide position.
  * Discard consensus reads that have >10% of the total number of consensus nucleotides being called as N.
* Align all retained consensus reads locally to either hg19 or hg38 human reference genome using researcher’s preferred aligner(s) such as Bowtie2 and BWA.
* Process aligned reads with Mpileup using parameters –BQ0 –d 10,000,000,000,000 to remove coverage thresholds to ensure a proper pileup output regardless of VAF.
* Filter out positions with less than 500x (user-specified) consensus read coverage.
* Use binomial distribution to call single nucleotide variants (SNPs) in retained data from Step 2.5.7 with the following parameters. The binomial statistic will be based on a genomic position-specific error model. Each genomic position is modeled independently after summing out the error rates of all samples for that particular position. Following the example:

  Probability of nucleotide profile at a given genomic position, p
  ∑ Variant RF2 ∑ Total RFs
  = 26/255505
  = 0.000101759
  Binomial probability of 24 variant RFs out of 35911 total RFs, P(X ≥ x) in sample K
  = 1 - binomial(24, 35911, 0.000101759)
  = 2.26485E-13

NOTE: For each genomic position queried, there would be three possible mutational changes (i.e.,A>T, A>C, A>G), and each of which would be represented as background artifact. Somatic events that are significantly different from the background after Bonferroni correction are retained. In the example shown in Table 1, the number of tests performed was 11, hence a Bonferroni corrected p-value ≤0.00454545 (0.05/11) was required to call an event as statistically significant.
Somatic events are required to be present in both replicates from the same specimen; otherwise, regard them as false positives.

## Requirements

### Tools & Frameworks
1. BWA
2. Sambamba
3. SAMtools
4. PEAR (Paired-End reAd mergeR)
5. snpEff
6. Cromwell (https://github.com/broadinstitute/cromwell)
7. Java
8. Python

### Additional python libraries
1. scipy
3. pysam
4. numpy

### R libraries
1. tidyverse

## Running a simulation test

Create a test dataset. ${reference_fasta} refers to a fasta sequence of hg19. 

```
cd tests
./simulate_fragments ${reference_fasta} > fragments.fa
./simulate_pe_reads
gzip read_1.fq
gzip read_2.fq
```

Now run the pipeline on the simulated dataset after setting the values in the input json file. The value corresponding to process_samples.inputs in the JSON fiel should be a file with the columns that refer to the sample name, library name, absolute path to the first read file, and absolute path to the second read file. The variable ${cromwell} should point to the JAR for cromwell (https://github.com/broadinstitute/cromwell). An example configuration file is included, but should be modified so that it point to the correct files and directories for the user.

```
java -Xmx4g -Dconfig.file=../src/local.conf -jar ${cromwell} run ../src/process_samples.wdl --inputs process_samples.json
```

The output file of the workflow is named 'variants.ann.txt' and can be found in the cromwell-executions folder. Lets make a soft link to it.

```
ln -s `find . -name "variants.ann.txt" ` .
```

Lets create a simple plot which show the sequencing depth and VAF of the mutations in this dataset, and whether we found them or not.

```
./plot_mutations.R
```

Take a look at the file assessment.pdf. Most of the variants we miss should be the ones with extremely low VAF.

