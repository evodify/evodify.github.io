---
layout: post
title: Genomic variant calling pipeline
date: 2018-08-30 15:36:50 +02:00
categories: Pipelines
tags: dog, genomics, python, variant
image: /assets/posts/2018-08-30-genomic-variant-calling-pipeline/Genomic-variant-calling-pipeline.jpeg
alt: genomic variant calling pipeline
description: An automated genomic variant calling pipeline becomes essential when a project scales to hundreds of genomes. Here is my genomic variant calling pipeline.
---

<p>As probably any beginner, I used to process my genomic data with manual interference at every step. So, I would submit mapping jobs for all samples on a computing cluster, when they all done I would submit mark duplicates jobs, etc. Moreover, I would also manually write sbatch scripts (my cluster <a href="http://www.uppmax.uu.se/" target="_blank">UPPMAX</a> uses the Slurm Workload Manager). It was not efficient.</p>

<p>Well, I used replacements (with <code>sed</code>) and loops (with <code>for i in x; do ...</code>) to reduce the amount of work, but there were many manual steps. I managed to process 24-31 small <em>Capsella</em> genomes (~200Mb) this way during my PhD projects. Now, I work with the dog genome which is much bigger (~2.5Gb) and I also need to analyze many more samples (82 genomes at the moment). So, I had to write this genomic variant calling pipeline to make my workflow as automatic as possible.</p>

<!--more-->

<p>This genomic variant calling pipeline was written for the dog data, but you can modify it for any other organism with a reference genome (see below). If you work with non-model organisms, I recommend you also check my <a href="{{ site.baseurl }}/gatk-in-non-model-organism/">GATK tutorial for non-model organisms</a>.</p>

<h2>Requirements</h2>

<p>For this genomic variant calling pipeline you would need the following software installed:</p>
<ul>
<li><a href="http://bio-bwa.sourceforge.net/" target="_blank">BWA</a> or <a href="http://www.well.ox.ac.uk/project-stampy" target="_blank">Stampy</a> - mapping.</li>
<li><a href="http://samtools.sourceforge.net/" target="_blank">Samtools</a> - SAM/BAM files manipulation.</li>
<li><a href="https://broadinstitute.github.io/picard/">Picard</a> - mark duplicates.</li>
<li><a href="http://qualimap.bioinfo.cipf.es/" target="_blank">QualiMap</a> (optional) - check the mapping quality.</li>
<li><a href="https://software.broadinstitute.org/gatk/" target="_blank">GATK</a> - genotype calling.</li>
</ul>

<h2>Genomic variant calling pipeline</h2>

<h3>Main steps</h3>

<p>This genomic variant calling pipeline includes the following steps:</p>
<ol>
<li>Mapping to the reference.</li>
<li>Merging BAM files of different lanes.</li>
<li>Mark duplicates.</li>
<li>Base Quality Score Recalibration (BQSR).</li>
<li>Check mapping quality (optional).</li>
<li>Genotype each sample in the GVCF mode.</li>
<li>Merge gVCF files from different samples.</li>
<li>Joint genotyping of all samples.</li>
</ol>

<h3>How to use</h3>

<p>This genomic variant calling pipeline consists of two files <em>reads-to-VCF_dog.py</em>  (main script) and <em>sbatch.py</em> (all the functions) that can be downloaded from <a href="https://github.com/evodify/to-generate-sbatch" target="_blank">my Github repository</a>.</p>

<p>To run it, you need a list of FASTQ files with R1 only names. You can obtain such a list by listing the content of the folder with your FASTQ files:</p>

```bash
ls -l fastq/*.fastq.gz | grep R1 | sed 's/ \+/ /g' | cut -d " " -f 9 > R1.txt
```

<p>Next, you run the <em>reads-to-VCF.py</em> script using <em>R1.txt</em> as an input.</p>

```bash
python reads-to-VCF.py -i R1.txt -r canFam3.fa -p snic2017 -n 20 -t 5-00:00:00 -l \$SNIC_TMP -a BWA -m 4 -k "dogs.557.publicSamples.ann.chrAll.PASS.vcf.gz,00-All_chrAll.vcf.gz"

```

<p>You can see the description of each option by running:</p>

```bash
python reads-to-VCF.py -h
```

<p>I only would like to mention that you can use either <code>BWA</code> or <code>stampy</code> aligners for the mapping step (option <code>-a</code>). In case of <code>stampy</code>, you can also specify the per base pair divergence from the reference (<code>-d</code>).</p>

<p>Among other options, <code>snic2017</code> is the project name that is needed for the Slurm on Uppmax, I am not sure if you need it in other systems.</p>

<p>I also recommend providing the maximum available number of cores (<code>-n</code>) for parallel jobs, and maximum available RAM per core (<code>-m</code>).</p>

<p>I also request 5 days (<code>-t</code>) for each main step to run on a cluster.</p>

<p>For many intermediate steps, I make use of the temporary storage <code>$SNIC_TMP</code> (you need to add <code>\</code> before the <code>$</code> sign because it is a special character). This helps to keep my account disk use low. Moreover, it makes the pipeline more efficient because <code>$SNIC_TMP</code> is located on the same node where the job is running. You need to provide your temporary directory with the option <code>-l</code>.</p>

<p>Assuming the content of <code>R1.txt</code> is the following</p>

```bash
WCRO84_S23_L005_R1_test.fastq.gz
WCRO84_S23_L006_R1_test.fastq.gz
WCRO86_S21_L005_R1_test.fastq.gz
WCRO86_S21_L007_R1_test.fastq.gz
```

<p>The script <em>reads-to-VCF.py</em> will generate these files:</p>

```bash
GVCF_all.sbatch
WCRO84_L005_map.sh
WCRO84_L006_map.sh
WCRO84_gVCF.sh
WCRO84_mergeMarkDuplBQSR.sh
WCRO84_qualimap.sh
WCRO86_L005_map.sh
WCRO86_L007_map.sh
WCRO86_gVCF.sh
WCRO86_mergeMarkDuplBQSR.sh
WCRO86_qualimap.sh
reads-to-VCF_wolf.sbatch
```

<p>All <code>*.sh</code> files contain commands to execute jobs. The file is <em>reads-to-VCF.sbatch</em> contains the BASH script to submit these jobs to the cluster with the dependencies between jobs.</p>

<p>You place all these files in the same directory and submit all the jobs with:</p>

```bash
sh reads-to-VCF.sbatch
```

<p>If everything worked fine, you will see this output:</p>
<div class="image">
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2018-08-30-genomic-variant-calling-pipeline/genomic-variant-calling-pipeline-sh-reads-to-VCF.jpeg" alt="genomic variant calling pipeline execution" /><figcaption class="aligncenter"> Execution of reads-to-VCF.sbatch</figcaption>
</figure></div>
<p>If you check the submitted jobs (with <code>jobinfo</code>), you will see the job dependencies (in red):</p>
<div class="image">
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2018-08-30-genomic-variant-calling-pipeline/genomic-variant-calling-pipeline-jobinfo-1024x146.jpeg" alt="submitted jobs of the genomic variant calling pipeline" /><figcaption class="aligncenter"> Lits of submitted jobs (click to enlarge)</figcaption>
</figure></div>
<p>Following the job IDs, you can see that, for instance, the job <em>WCRO84_mergeMarkDuplBQSR</em> won't start until the two mapping jobs <em>WCRO84_L005_map</em> and <em>WCRO84_L006_map</em> are not finished successfully. When all these jobs finish, you will get a set of <a href="https://gatkforums.broadinstitute.org/gatk/discussion/11004/gvcf-genomic-variant-call-format" target="_blank">gVCF</a> files, which you can merge and jointly genotype with <em>GVCF_all.sbatch</em>:</p>

```bash
sbatch GVCF_all.sbatch
```

<p><strong>Thus, you run only 4 commands and you get from a set if FASTQ files to one VCF file.</strong></p>

<p>I do not make the last command to execute automatically because it depends on all samples in this genomic variant calling pipeline. If your pipeline runs on many samples it is very likely that some samples will fail simply due to chance. For example, when I run this pipeline on 82 samples, 3 of them failed due to node failure.</p>

<p>If some of the jobs fail, the jobs dependent on it will be canceled.</p>

<p>To check if any job failed, run (assuming you use Slurm):</p>

```bash
finishedjobinfo | grep FAILED
```

<p>Then investigate why it failed. If it is not your mistake and it failed due to node failure, restart failed jobs by extracting and executing the code for the failed jobs from <em>reads-to-VCF.sbatch</em>.</p>

<p>I also suggest to list all the final BAM and gVCF files and check their size before executing <em>GVCF_all.sbatch</em>:</p>

```bash
ls -l *.gz *.bam
```

<p>If some samples have extremely large or extremely small files, check their <code>*.fastq, *.sh *.err</code> and <code>*.out</code> files to figure out what the problem is.</p>

<p>Also, you may want to modify <code>GVCF_all.sbatch</code> for your needs. For example, you may want to genotype mitochondrial and sex chromosomes with the control of ploidy (e.g. <code>-ploidy 1 -L chrX -L chrM</code>).</p>

<p>Generally, there should not be any problem and you will get your VCF file with as little effort as 4 commands.</p>

<p>Note, this genomic variant calling pipeline was designed and tested on <a href="http://www.uppmax.uu.se/" target="_blank">UPPMAX</a> with the Slurm Workload Manager, if you run it on another cluster, you may need to modify it. Below, I provide more details for each step and describe the hardcoded parts.</p>

<h3>Description of each step</h3>

<h4>1. Mapping to the reference</h4>

<p>The very first step of this genomic variant calling pipeline is to map all reads to the reference genome. Many of my samples were sequenced on different lanes. So, to control for this factor and to make the mapping process faster, I perform the mapping of each lane separately. An example command of the BWA mapping with subsequent sorting and converting to BAM by Samtools looks like this:</p>

```bash
bwa mem -t 20 -M -R '@RG\tID:WCRO84_L005\tPL:illumina\tLB:WCRO84\tSM:WCRO84' canFam3.fa WCRO84_S23_L005_R1_test.fastq.gz WCRO84_S23_L005_R2_test.fastq.gz | samtools sort -O bam -@ 20 -m 4G > WCRO84_L005.bam
```

<p>where <code>WCRO84</code> is the sample name, and <code>L005</code> is the lane number. These names are changed by the <em>reads-to-VCF.py</em> script. You can read about the options of BWA and Samtools in the documentation of these programs.</p>

<p>This step can be parallelized and I use 20 threads here.</p>

<p>If you choose to map with Stampy, the mapping command will look like this:</p>

```bash
stampy.py -g canFam3 -h canFam3 --substitutionrate=0.002 -t 20 -o $SNIC_TMP/WCRO84_L005_stampy.bam --readgroup=ID:WCRO84_L005  -M WCRO84_S23_L005_R1_test.fastq.gz WCRO84_S23_L005_R2_test.fastq.gz
samtools sort -O bam -@ 20 -m 4G $SNIC_TMP/WCRO84_L005_stampy.bam > $SNIC_TMP/WCRO84_L005.bam

java -Xmx4G -jar picard.jar AddOrReplaceReadGroups \
I=$SNIC_TMP/WCRO84_S23_L005.bam \
O=WCRO84_S23_L005.bam \
RGID=WCRO84_S23_L005 \
RGLB=WCRO84 \
RGPL=illumina \
RGPU=WCRO84_S23_L005 \
RGSM=WCRO84
```

<p>The last command adds read-group tags. I was not able to past these tags correctly with <code>--read-group</code> option in stampy. So, I added this extra step.</p>

<p>Mapping with Stampy is preferable if you <a href="{{ site.baseurl }}/gatk-in-non-model-organism/">map reads that are divergent from the reference</a>. For example, I use Stampy to map wolf reads to the dog reference genome with the account for divergence from the reference (<code>-d 0.002</code>).</p>

<p>There are separate files generated for each mapping job. The example above would be named as <em>WCRO84_L005_map.sh</em>.</p>

<h4>2. Merge lanes, mark duplicates, BQSR</h4>

<p>All these three steps are performed as a sequence in one cluster job on a single core. So, they all are listed in one file, e.g. <em>WCRO84_mergeMarkDuplBQSR.sh</em>.</p>

<h5>Merge lanes</h5>

<p>The BAM files of different lanes are then merged for each samples:</p>

```bash
ls WCRO84_L00?.bam | sed 's/  / /g;s/ /\n/g' > WCRO84_bam.list
samtools merge -b WCRO84_bam.list $SNIC_TMP/WCRO84_merged.bam
xargs -a WCRO84_bam.list rm
```

<p>The first line creates a list of BAM files to merge. The second line performs merging. The last line removed the BAM files that have been merged to free some disk space. The file <em>WCRO84_bam.list</em> can be checked to make sure all the BAM files were included during the merging step.</p>

<h5>Mark duplicates</h5>

<p>This step identifies reads that are likely artificial duplicates originated during sequencing workflow. These reads are not removed but they will be ignored by default during the variant calling process.</p>

```bash
java -Xmx4G -Djava.io.tmpdir=$SNIC_TMP -jar picard.jar MarkDuplicates \
VALIDATION_STRINGENCY=LENIENT \
METRICS_FILE=WCRO84_merged_markDupl_metrix.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=15000 \
INPUT=$SNIC_TMP/WCRO84_merged.bam \
OUTPUT=$SNIC_TMP/WCRO84_merged_markDupl.bam

rm $SNIC_TMP/WCRO84_merged.bam
```

<p>The <code>rm</code> cleans unneeded files.</p>

<h5>Base Quality Score Recalibration (BQSR)</h5>

<p>This part is the slowest in the script. It applies some machine learning to detect and correct systematic errors in the base quality scores.</p>

<p>You need a good reference database of SNPs and indels for this step. I have used some of the files which I cannot share, but if you also run this genomic variant calling pipeline for dog data, you can download these two publicly available variant files:</p>

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/eva/PRJEB24066/dogs.557.publicSamples.ann.vcf.gz*
wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/dog_9615/VCF/00-All.vcf.gz*

```

<p>I also removed all filtered positions, because BQSR needs only the high confidence sites:</p>

```bash
zcat 00-All.vcf | awk '{if (/^#/) {print $0} else {print "chr"$0}}' | bgzip > 00-All_chrAll.vcf
zcat dogs.557.publicSamples.ann.vcf.gz | awk '{if (/^#/) {print $0} else if ($7=="PASS") {print "chr"$0}}' | bgzip > dogs.557.publicSamples.ann.chrAll.PASS.vcf
tabix 00-All_chrAll.vcf & tabix dogs.557.publicSamples.ann.chrAll.PASS.vcf

```

<p>To perform BQSR, you run these commands:</p>

```bash
# Generate the first pass BQSR table file
gatk --java-options "-Xmx4G"  BaseRecalibrator \
-R canFam3.fa \
-I $SNIC_TMP/WCRO84_merged_markDupl.bam \
--known-sites 00-All_chrAll.vcf.gz \
--known-sites dogs.557.publicSamples.ann.chrAll.PASS.vcf.gz \
-O WCRO84_merged_markDupl_BQSR.table

# Apply BQSR
gatk --java-options "-Xmx4G"  ApplyBQSR \
-R canFam3.fa \
-I $SNIC_TMP/WCRO84_merged_markDupl.bam \
-bqsr WCRO84_merged_markDupl_BQSR.table \
-O WCRO84_merged_markDupl_BQSR.bam

rm $SNIC_TMP/WCRO84_merged_markDupl.bam

# Generate the second pass BQSR table file
gatk --java-options "-Xmx4G"  BaseRecalibrator \
-R canFam3.fa \
-I WCRO84_merged_markDupl_BQSR.bam \
--known-sites 00-All_chrAll.vcf.gz \
--known-sites dogs.557.publicSamples.ann.chrAll.PASS.vcf.gz \
-O WCRO84_merged_markDupl_BQSR2.table

# Plot the recalibration results
gatk --java-options "-Xmx4G"  AnalyzeCovariates \
-before WCRO84_merged_markDupl_BQSR.table \
-after WCRO84_merged_markDupl_BQSR2.table \
-plots WCRO84_merged_markDupl_BQSR.pdf
```

<p>As I have mentioned it in my <a href="{{ site.baseurl }}/gatk-the-best-practice-for-genotype-calling-in-a-non-model-organism/"> GATK: the best practice for non-model organisms</a>, the BQSR can actually hurt if your variant reference database is not good enough. So, check the plots in <em>WCRO84_merged_markDupl_BQSR.pdf</em> to make sure your recalibration worked correctly. Mine worked fine and it looked like this:</p>
<div class="image">
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2018-08-30-genomic-variant-calling-pipeline/genomic-variant-calling-pipeline-BQSR.jpeg" alt="Some of the BQSR plots" /><figcaption class="aligncenter"> Some of the BQSR plots</figcaption>
</figure></div>
<h4>3. Check mapping quality (Optional)</h4>

<p>This step is not essential for the genomic variant calling pipeline, but it useful to check the quality of your BAM files. With Qualimap you get a lot of useful information such as coverage, number of mapped/unmapped reads, mapping quality, etc. An example of the file name for this step is <em>WCRO84_qualimap.sh</em>.</p>

<p>The command is the following:</p>

```bash
qualimap bamqc -nt 20 --java-mem-size=4G -bam WCRO84_merged_markDupl_BQSR.bam -outdir WCRO84_qualimap
```
<div class="image">
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2018-08-30-genomic-variant-calling-pipeline/genomic-variant-calling-pipeline-qualimap.jpeg" alt="Qualimap plot showing coverage" /><figcaption class="aligncenter"> One of the Qualimap plots</figcaption>
</figure></div>
<p>The report will be in an HTML file which you can view in your browser. For really many samples, you can also batch-process the TXT files with results. For example, to extract the mean coverage information, run</p>

```bash
grep "mean coverageData" */genome_results.txt
```

<p>At this stage, the work with BAM files is finished and the pipeline proceeds to the genotyping steps.</p>

<h4>4. Call variants per-sample (gVCF)</h4>

<p>This step runs the GATK HaplotypeCaller in GVCF mode:</p>

```bash
gatk --java-options "-Xmx4G" HaplotypeCaller \
-R canFam3.fa \
-ERC GVCF \
-I WCRO84_merged_markDupl_BQSR.bam \
-O WCRO84_merged_markDupl_BQSR.g.vcf.gz

tabix WCRO84_merged_markDupl_BQSR.g.vcf.gz

mkdir WCRO84
mv WCRO84_* WCRO84
```

<p>This enables to generate intermediate files, which are then used for joint genotyping in an efficient way. <code>tabix</code> is necessary to index the compressed file for the GATK compatibility in the next steps.</p>

<p>At this step, I also create sample specific directories and move all the files into the sample's own directory. This helps to keep the order when you work with many files.</p>

<p>An example of the file name for this step is <em>WCRO84_gVCF.sh</em>.</p>

<h4>5. Joint genotyping of all samples.</h4>

<p>This step merges all gVCF files from multiple samples:</p>

```bash
gatk --java-options "-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CombineGVCFs \
-R canFam3.fa \
-V WCRO84/WCRO84_merged_markDupl_BQSR.g.vcf.gz \
-V WCRO86/WCRO86_merged_markDupl_BQSR.g.vcf.gz \
-O GVCF_merged_markDupl_BQSR.g.vcf.gz

tabix GVCF_merged_markDupl_BQSR.g.vcf.gz
```

<p><em>In the example command above, I provide only 2 samples to keep it simple. Obviously, there will be more samples in a real run.</em></p>

<p>The final command performs joint-call of SNPs and indels:</p>

```bash
gatk --java-options "-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
-R reference/canFam3.fa \
-V GVCF_merged_markDupl_BQSR.g.vcf.gz \
-O GVCF_merged_markDupl_BQSR.vcf.gz

tabix GVCF_merged_markDupl_BQSR.vcf.gz
```

<p>The file name of the script that includes both these commands is <em>GVCF_all.sbatch</em>.</p>

<h3>Script to submit all these jobs</h3>

<p>All the job described above are submitted by executing this script (2 samples example):</p>

```bash
#!/bin/sh

WCRO84_L005=$(sbatch -A snic2017 -p core -n 20 -t 5-00:00:00 -J WCRO84_L005_map -e WCRO84_L005_map.err -o WCRO84_L005_map.out WCRO84_L005_map.sh | cut -d " " -f 4)
echo "WCRO84_L005_map.sh under the ID $WCRO84_L005 has been submitted"

WCRO84_L006=$(sbatch -A snic2017 -p core -n 20 -t 5-00:00:00 -J WCRO84_L006_map -e WCRO84_L006_map.err -o WCRO84_L006_map.out WCRO84_L006_map.sh | cut -d " " -f 4)
echo "WCRO84_L006_map.sh under the ID $WCRO84_L006 has been submitted"

WCRO84_mergeMarkDuplBQSR=$(sbatch --dependency=afterok:$WCRO84_L005:$WCRO84_L006 -A snic2017 -p core -n 1 -t 5-00:00:00 -J WCRO84_mergeMarkDuplBQSR -e WCRO84_mergeMarkDuplBQSR.err -o WCRO84_mergeMarkDuplBQSR.out WCRO84_mergeMarkDuplBQSR.sh | cut -d " " -f 4)
echo "WCRO84_mergeMarkDuplBQSR.sh under the ID $WCRO84_mergeMarkDuplBQSR has been submitted"

WCRO84_gVCF=$(sbatch --dependency=afterok:$WCRO84_mergeMarkDuplBQSR -A snic2017 -p core -n 1 -t 5-00:00:00 -J WCRO84_gVCF -e WCRO84_gVCF.err -o WCRO84_gVCF.out WCRO84_gVCF.sh | cut -d " " -f 4)
echo "WCRO84_gVCF.sh under the ID $WCRO84_gVCF has been submitted"

WCRO84_qualimap=$(sbatch --dependency=afterok:$WCRO84_mergeMarkDuplBQSR -A snic2017 -p core -n 1 -t 5-00:00:00 -J WCRO84_qualimap -e WCRO84_qualimap.err -o WCRO84_qualimap.out WCRO84_qualimap.sh | cut -d " " -f 4)
echo "WCRO84_qualimap.sh under the ID $WCRO84_qualimap has been submitted"

WCRO86_L005=$(sbatch -A snic2017 -p core -n 20 -t 5-00:00:00 -J WCRO86_L005_map -e WCRO86_L005_map.err -o WCRO86_L005_map.out WCRO86_L005_map.sh | cut -d " " -f 4)
echo "WCRO86_L005_map.sh under the ID $WCRO86_L005 has been submitted"

WCRO86_L007=$(sbatch -A snic2017 -p core -n 20 -t 5-00:00:00 -J WCRO86_L007_map -e WCRO86_L007_map.err -o WCRO86_L007_map.out WCRO86_L007_map.sh | cut -d " " -f 4)
echo "WCRO86_L007_map.sh under the ID $WCRO86_L007 has been submitted"

WCRO86_mergeMarkDuplBQSR=$(sbatch --dependency=afterok:$WCRO86_L005:$WCRO86_L007 -A snic2017 -p core -n 1 -t 5-00:00:00 -J WCRO86_mergeMarkDuplBQSR -e WCRO86_mergeMarkDuplBQSR.err -o WCRO86_mergeMarkDuplBQSR.out WCRO86_mergeMarkDuplBQSR.sh | cut -d " " -f 4)
echo "WCRO86_mergeMarkDuplBQSR.sh under the ID $WCRO86_mergeMarkDuplBQSR has been submitted"

WCRO86_gVCF=$(sbatch --dependency=afterok:$WCRO86_mergeMarkDuplBQSR -A snic2017 -p core -n 1 -t 5-00:00:00 -J WCRO86_gVCF -e WCRO86_gVCF.err -o WCRO86_gVCF.out WCRO86_gVCF.sh | cut -d " " -f 4)
echo "WCRO86_gVCF.sh under the ID $WCRO86_gVCF has been submitted"

WCRO86_qualimap=$(sbatch --dependency=afterok:$WCRO86_mergeMarkDuplBQSR -A snic2017 -p core -n 1 -t 5-00:00:00 -J WCRO86_qualimap -e WCRO86_qualimap.err -o WCRO86_qualimap.out WCRO86_qualimap.sh | cut -d " " -f 4)
echo "WCRO86_qualimap.sh under the ID $WCRO86_qualimap has been submitted"
```

<p><em>Again, I assume and describe a cluster system with the Slurm Workload Manager, if your cluster is different, you may need to modify this script.</em></p>

<p>It works pretty simply. The script submits all mapping jobs (e.g. <code>WCRO84_L00?_map.sh</code>) as independent jobs, store their job ID in variables (e.g. <code>$WCRO84_L005, $WCRO84_L006</code>) and uses these IDs as dependencies for the subsequent <code>WCRO84_mergeMarkDuplBQSR.sh</code> job (<code>--dependency=afterok:</code>). In its turn, <code>WCRO84_gVCF.sh</code> uses the ID of <code>WCRO84_mergeMarkDuplBQSR.sh</code> as a dependency, etc.</p>

<h2>Hardcoded parts</h2>

<p>Most of the variables of this genomic variant calling pipeline  can be changed with the command line options of <em>reads-to-VCF.py.</em> Let me remind you that you can see the available options with:</p>

```bash
python reads-to-VCF.py -h
```

<p>However, a few parts are hardcoded and you may need to modify them.</p>

<p>I also need to explain here that the script <em>reads-to-VCF.py</em> depends on <a href="https://github.com/evodify/to-generate-sbatch/blob/master/sbatch.py" target="_blank">sbatch.py</a> module that is also located in <a href="https://github.com/evodify/to-generate-sbatch/blob/master/reads-to-VCF_dog.py" target="_blank">my Github repository</a>. I store all the functions in that module, so if you need to edit some hardcoded parts, you are most likely need to edit <em>sbatch.py</em>.</p>

<p>1. You may need to modify the functions <code>writeMapBWAJob</code> and <code>writeMapStampyJob</code> if your FASTQ file names don't follow this structure:</p>

```bash
WCRO84_S1_L005_R1_001.fastq.gz
```

<p>Alternatively, you can rename your files to fit this structure.</p>

<p>2. Uppmax Slurm outputs the job ID in a sentence: <code>Submitted batch job 4724521</code> So, I had to use the command <code>cut -d " " -f 4</code> to extract the ID. If your cluster simply outputs the ID, you can remove the <code>cut</code> command from the function <code>writeSbatchScript</code>.</p>

<p>3. Possibly, you also need to adjust <code>MAX_FILE_HANDLES_FOR_READ_ENDS_MAP</code> in the function <code>writeMarkDuplJob</code> according to your system. As a rule, it should be little smaller than the output of the <code>ulimit -n</code> command.</p>

<p>4. <code>picard.jar</code> is assumed to be located in the working directory. You can either edit the code (functions <code>writeMapStampyJob</code> , <code>writeMarkDuplJob</code>) to specify the full path to <code>picard.jar</code>, or copy <code>picard.jar</code> to your working directory before submitting the jobs.</p>

<h2>Final thoughts</h2>

<p>I hope you will use this genomic variant calling pipeline in your workflow. In my case of 82 dog samples with the mean coverage of x30, running the whole pipeline took little more than a week. Given that it is automatic and I was able to do some other work during this time, it is not a lot.</p>

<p>I know that such a pipeline can be made even more efficient and reliable with <a href="https://www.nextflow.io/" target="_blank">Nextflow</a> and <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">SnakeMake</a>. But I found out about these workflow management systems after I wrote this pipeline. Moreover, I first thought to write this genomic variant calling pipeline in <a href="{{ site.baseurl }}/processing-genomic-data-apache-spark-big-data-tutorial/">Apache Spark</a>, but Spark is not mainstream yet and Spark version of the GATK is still in beta. So, I hope I will update this pipeline with a more advanced framework in the future.</p>
