# LycAdmixMosaic
Notes, scripts and analyses for my study of admixture and ancestry in *Lycaeides* from PoolSeq data

This is a fork of the contemporary evolution data set and analyses that excludes the time dimension. See [Pool_DNA_plates](https://drive.google.com/drive/folders/1U4AsshyMvlySNtODuSLWo0dYDH_rgho4) for details on samples and sample concentrations.

# Data
The raw pool-seq data are currently in /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/. The set of samples I am using here, which are from a single (initial) round of sequencing are as follows:

| Population | Samples (N)| Nominal taxon |
|------------|------------|---------------|
| ABM | ABM20 (48) | L. melissa |
| BCR | BCR17 (48), BCR17rep (48) | JH (admixed) |
| BHP | BHP19 (48) | L. melissa |
| BKM | BKM19 (33) | Warners (admixed) |
| BTB | BTB17 (48), BTB17rep (48) | JH (admixed) | 
| CLH | CLH19 (36) | White Mt (admixed) |
| CP | CP19 (48), CP19rep (48) | Sierra (admixed) |
| EP | EP19 (48), EP19rep (48) | Warners (admixed) |
| GNP | GNP17 (56) | L. iads |
| HJ | HJ20 (48) | L. melissa |
| HNV | HNV17 (48) | JH (admixed) |
| LS | LS19 (48) | L. anna |
| MEN | MEN12 (10) | L. argus (France) |
| MR |  MR20 (48) | Sierra (admixed) |
| MTU | MTU20 (48) | L. melissa |
| SBW | SBW18 (20) | L. idas (Alaska) | 
| SHC | SHC11 (46) | L. anna ricei |
| SIN | SIN10 (48) | L. melissa |
| SUV | SUV20 (51) | L. melissa |
| TBY | TBY11 (24) | L. idas sublivens | 
| TIC | TIC19 (48) | Sierra (admixed) |
| VE | VE20 (48) | L. melissa |
| YG | YG20 (48) | L. anna |

Replicates (rep) are the same DNA extractions re-pooled and sequenced as a distinct sample. MEN = outgroup, *Plebejus argus*. The data were generated and cleaned up by BGI (with `soapnuke`). Here is the report from BGI:  [BGI_F22FTSUSAT0310-01_LYCgpswR_report_en.pdf](https://github.com/zgompert/LycSpaceTimePoolSeq/files/9940314/BGI_F22FTSUSAT0310-01_LYCgpswR_report_en.pdf).

# DNA Sequence Alignment

I am aligning the DNA sequence data to the updated (based on PacBio) *L. melissa* genome. I am using `bwa-mem2` for this, which is basically just a sped up version of `bwa mem` that also works directly with gzipped files [https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2). I am using `bwa-mem2` version 2.0pre2. 

First, I (re)indexed the reference genome.

```bash
## index genome with bwa-mem2
/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta
```
Then, I set up the alignment. The submission scrip is (from /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts):

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bwa-mem2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignments

perl BwaMemFork.pl ../F22FTSUSAT0310-01_LYCgpswR/soapnuke/clean/*/*1.fq.gz 
```

Which runs the following:

```perl
#!/usr/bin/perl
#
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";

FILES:
foreach $fq1 (@ARGV){
	$pm->start and next FILES; ## fork
	$fq2 = $fq1;
	$fq2 =~ s/_1\.fq\.gz/_2.fq.gz/ or die "failed substitution for $fq1\n";
        $fq1 =~ m/clean\/([A-Za-z0-9]+)/ or die "failed to match id $fq1\n";
	$ind = $1;
	$fq1 =~ m/([A-Za-z_\-0-9]+)_1\.fq\.gz$/ or die "failed match for file $fq1\n";
	$file = $1;
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:Lyc-"."$ind\\tLB:Lyc-"."$ind\\tSM:Lyc-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$file.bam - && samtools index -@ 2 $ind"."_$file.bam\n";

	$pm->finish;
}
```
I am pipping the results on to `samtools` (version 1.16) to compress, sort and index the alignments.

I then merged the bam files for each population using `samtools` version 1.16. This was donw with the following shell script:

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=merge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignment

perl ../Scripts/MergeFork.pl 
```

Which runs


```perl
#!/usr/bin/perl
#
# merge alignments for each population sample with samtools version XX 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

open(IDS,"pids.txt");
while(<IDS>){
	chomp;
	push(@IDs,$_);
}
close(IDS);

FILES:
foreach $id (@IDs){
	$pm->start and next FILES; ## fork
        system "samtools merge -c -p -o Merged/$id.bam $id"."_*.bam\n";
	system "samtools index -@ 2 Merged/$id.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```

pids.txt lists all of the population IDs.

# Removing PCR duplicates

I am using `samtools` (version 1.16) to remove PCR duplicates. This is more efficient than `PicardTools` and performs similarly, see [Ebbert 2016](https://link.springer.com/article/10.1186/s12859-016-1097-3). I am following the standard prtocol from [`samtools`](https://www.htslib.org/doc/samtools-markdup.html). I am using the default option (same as `-m t`) to measure positions based on template start/end. And I am using `-r` to not just mark but remove duplicates. The submission script is:

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=dedup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)


cd /scratch/general/nfs1/dedup

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts/RemoveDupsFork.pl *bam
```

Which runs

```perl
#!/usr/bin/perl
#
# PCR duplicate removal with samtools
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $bam (@ARGV){
	$pm->start and next FILES; ## fork
	$bam =~ m/^([A-Za-z0-9]+)/ or die "failed to match $bam\n";
	$base = $1;
	system "samtools collate -o co_$base.bam $bam /scratch/general/nfs1/dedup/t$bam\n";
	system "samtools fixmate -m co_$base.bam fix_$base.bam\n";
	system "samtools sort -o sort_$base.bam fix_$base.bam\n";
	## using default definition of dups
	## measure positions based on template start/end (default). = -m t
	system "markdup -T /scratch/general/nfs1/dedup -r sort_$base.bam dedup_$base.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```

# Variant Calling

I called variants with `bcftools` (version 1.16). I did not perform INDEL realignment (it doesn't seem to be important) and this approach is not aware I have pooled data. Nonetheless, this should be a solid approach for estimating allele frequencies and generating alignment data for the tree-based analyses. 

For this, I ran the following submission script,

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bcf_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /scratch/general/nfs1/dedup

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts/BcfForkLg.pl chrom*list 
```

which runs

```perl
#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#


use Parallel::ForkManager;
my $max = 26;
my $pm = Parallel::ForkManager->new($max);

my $genome ="/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";




foreach $chrom (@ARGV){
	$pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
	$out = "o_lycpool_chrom$1";
	system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o $out"."vcf\n";
	$pm->finish;

}

$pm->wait_all_children;
```
Note that each chromosome (big scaffold) is being processed separately (chrom*list). 

The variant data are in `/uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/SpecGenomVars`.

I filtered the vcf file with `GATK` version (4.1.4.1), keeping only those with mapping quality > 30, depth > 1350 and bias scores less than +- 3. This uses [VarFiltFork2.pl](VarFiltFork2.pl), which is also shown below.

```perl
#!/usr/bin/perl
#
# filter vcf files 
#

use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

foreach $vcf (@ARGV){ ## takes vcf.gz
        $pm->start and next; ## fork
        $o = $vcf;
        $o =~ s/filt_// or die "failed sub $o\n";
        $in = "b_$vcf";
        system "gunzip -c $vcf | bgzip > $in\n";
        system "tabix $in\n";
        system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar IndexFeatureFile -I $in\n";
        system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration -R /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta -V $in -O fff_$o --filter-name \"bqbz\" --filter-expression \"BQBZ > 3.0 || BQBZ < -3.0\" --filter-name \"mqbz\" --filter-expression \"MQBZ > 3.0 || MQBZ < -3.0\" --filter-name \"rpbz\" --filter-expression \"RPBZ > 3.0 || RPBZ < -3.0\" --filter-name \"depth\" --filter-expression \"DP < 1350\" --filter-name \"mapping\" --filter-expression \"MQ < 30\" --verbosity ERROR\n";
        system "bgzip -d fff_$o\n";

        $pm->finish;

}

$pm->wait_all_children;
```

Next, I extracted the allele depths (count of each allele for each SNP) from the filtered vcf files. This script also drops INDELS and multiallelic data.

```bash
#!/usr/bin/bash
#
# extract allele depth AD from biallelic SNPs that passed filtering 
#

for f in fff*vcf
do
	echo "Processing $f"
	out="$(echo $f | sed -e 's/vcf/txt/')"
	echo "Output is ad1_$out"
	grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\1/g' > ad1_$out   
	grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\2/g' > ad2_$out
done

```

This creates allele depth files for each allele (ad1* and ad2*) and chromosome, which I can use for downstream analyses. 

I also grapped the SNP information (alleles):

```bash
#!/usr/bin/bash
#
# extract alleles from biallelic SNPs that passed filtering 
#

for f in fff*vcf
do
	echo "Processing $f"
	out="$(echo $f | sed -e 's/vcf/txt/')"
	echo "Output is snps_$out"
	grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | cut -f 4,5 > snps_$out &   
done
```

All downstream analyses will be in `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/`.

# Population genetic structure

As a first pass, I summarized patterns of population structure based on allele frequencies for chromosome and sets of chromosomes (autosomes vs Z). My initial analyses, which include some estimates of Fst as well, are in [pcaFst.R](pcaFst.R), wheras code for a summary figure with PCAs for all autosomes vs Z and a map (with a blank spot for a tree) are in [mkPCAfig.R](mkPCAfig.R). Neither is final. I want to pay attention to how exactly I get allele frequeny estiamtes and to any additional filtering for the PCA (e.g., focus on SNPs with some minimal coverage in all populations). This is all in the subdirectory `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/LycAdmix/GenData`.

# Time-calibrated phylogenetic tree

I am working on a time-calibrated tree of *Lycaeides* with `Beast` (version 2.7.7). I first tried this by sampling a multiplocus haplotype from the allele frequencies (as I am still currently doing for `Caster`, see below), but I think this was inflating the terminal branch lengths. I thus instead decided to generate a consensus sequence for each population, whereby I chose the most common allele for each SNPs. This is all done from the allele depth information and I insert N when fewer than 5 reads were observed. This was done with the [mkBeastDat.R](mkBeastDat.R) script (in `GenData`), which outputs the max_chrom*.fasta alignment files (the script for `Caster` makes similar files without the max). I moved these to the `Beast` subdirectory. I then used [SubAlign.pl](SubAlign.pl) to drop the replicate population samples, which creates the sub_max_chom*.fasta files.

`Beast` cannot handle the full set of SNPs. I thus subset the SNP data from the alignments. This involves first generating a filtered and reduced set of SNPs with [GetSNPSubstMax.R](GetSNPSubstMax.R), which uses coverted (to numeric) alignments generated by [mkNumericFasta.pl](mkNumericFasta.pl). I specifically retain 0.035 percent of the SNPs from each chromosome and only those with no missing data. Also at this point, I dropped chromsome 21.2 (the small bit that might be part of 21). This outputs the keepSNPs_max_chrom files. I then create the combined (across chromosomes) fasta alignment file with just the subset of SNPs, lyc_genomemax.fasta. This was done with [SubSetFasta.pl](SubSetFasta.pl). I then converted this to a nexus format alignment for `Beast`.

```bash
seqmagick convert --output-format nexus --alphabet dna lyc_genomemax.fasta lyc_genomemax.nex
```

The alignment file includes 5408 characters (SNPs), 23 taxa and no missing data. I experimented with `Beast` to assess the effects of different priors and setting on perfomance, especially in terms of how to appropriately use the time calibrarions. This is all entered via `beauti`.

```bash
ml beast
beauti
```

Part of my approach involves entering the (approximate) number of invariant sites. For this, I am computing the total number of each base (A, C, G and T) in the genome, reducing this to 0.035 of the total (to acount for my subsampling) and subtracting of the SNPs. The current version doesn't use the final SNP set, so I will need to update this slightly. I get the counts (for the genome and SNP set) with [countBases.pl](countBases.pl).

I then add the conts by appending "Orig" to the name, "lyc_genomemax' in the xml from `beauti`,  and adding the following just after the data (number of constant A, C, G and T):

```xml
    <data id='lyc_genomemax' spec='FilteredAlignment' filter='-' data='@lyc_genomemaxOrig' constantSiteWeights='50190 28322 28284 50114'>
    </data>
```

I am using two sources for the calibration dates.

1. [Vila et al. 2011](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2010.2213) focuses on the colonization of the Americas by *Polyommatus* blues. This is older and has more limited data (a few mitochondrial and nuclear genes) but includes a good sample of old and new world *Lycaeides*. The TMRCA for *Lycaeides* (old and new world) was estimated at 2.4 million years.

2. [Kawahara et al. 2023](https://www.nature.com/articles/s41559-023-02041-9) uses many more genes and species for a global phylogeny of butterflies. With that siad, this only has *L. melissa*, *L. anna* and *L. idas* from the contiguous USA, and thus the TMRCA for these likely underestimates that of our set of *Lycaeides* (which includes Alaska, and this appears to be the deepest split in the ingroup for our tree). Figure S1 in the paper gives the TMRCA for *Lycaeides* as 1.29 million years and for *Lycaeides* + *Plebejus argus* as 5.22 million years (this is our outgroup, their sample is from Japan, ours from France). My most recent runs set normal priors on both of these nodes based on this data. For *Lycaeides*, I used mean = 1.845 and SD = 0.283, which gives 95 density intervals of 1.29 to 2.4 (my dates). This is worth thinking about a bit more, but is likely about right. I used 5.22 for the mean for *Lycaeides* + *P. argus* and the same SD, 0.283, to reflect a similar level of uncertainty.

Other than that, here are my main thoughts on choices of priors:

- GTR substitution model with relative rates estimated and gamma rate heterogeneity approximated by four categoreis (everything else default). This is flexible and we seem to have sufficient data to estimate the parameters. So, I am happy with this choice.

- I considered two clock models: the Optimised Relaxed Clock (ORC) and the Random Local Clock (RLC). The  ORC allows for rate variation and works well whether the rates and branch lengths are correlated or not . I want to read and think more about this. As noted in the help, the RLC allows the clock model to behave like a strict clock (if there are no rate changes) and a relaxed clock (if there are many rate changes). The number of rate changes is sampled during MCMC, and is drive by the data and prior. From this perspective, it is superior to the strict clock and relaxed clock models. However, it may suffer from convergence issues. The latter was an issue but I solved it (see below).

- I am **estimating the clock rate**, which is important given my calibration information. This requires messing with settings in the file menu to be allowed. When I don't do this, I get very wonky results. I could instead maybe set it to 0.0029 based on the *Heliconius* mutation rate[ Keightley et al. 2015](https://academic.oup.com/mbe/article/32/1/239/2925597), but estimating it and using the calibration nodes is more sensible (I am using 0.0029 as the starting value).

- I considered the Coalescent Bayesian Skyline prior (BSP) and the Extended Coalescent Bayesian Skyline Prior (EBSP). Either is appropriate in cases where you have multiple populations or individuals from the same species. The two approaches differ in terms of where (at nodes or anywhere) effective population size changes. I left the rest of the priors at their defaults, but this is worth revisiting.

I ran full analyses with three different combinations: ORC and BSP [lyc_wgs_max.xml](lyc_wgs_max.xml), ORC and EBSP [lyc_wgs_max_ebsp.xml](lyc_wgs_max_ebsp.xml), and RLC and BSP [lyc_wgs_max_ranlc.xml](lyc_wgs_max_ranlc.xml). All three gave the same topology and mostly similar branch lengths. I decided to focus on the RLC with BSP as I think the RLC is better dealing with the SNP-based nature of the data (all was derived from allele frequencies) and yielding more sensible branch lengths for the shallower divergences. This initially gave me some trouble with mixing (unlike the others). I solved this by using [coupled MCMC](https://github.com/nicfel/CoupledMCMC?tab=readme-ov-file) with four hot chains and one cold chain (delta temp = 0.025). I ran six coupled MCMC runs (each with the five aforementioned chains), each comprising 500,000,000 iterations. Log samples (parameter values) were stored every 1000 iterations, whereas tree samples were stored every 5000. 


I ran these with `Beast`:

```bash
ml beast
beast -prefix ch1 lyc_wgs_max_ranlc.xml
beast -prefix ch2 lyc_wgs_max_ranlc.xml 
beast -prefix ch3 lyc_wgs_max_ranlc.xml 
beast -prefix ch4 lyc_wgs_max_ranlc.xml 
beast -prefix ch5 lyc_wgs_max_ranlc.xml  
beast -prefix ch6 lyc_wgs_max_ranlc.xml 

```

I then used `logcombiner` (version 2.7.5) to combine the trees and logs. I applied the 20% burnin to each chain at this stage and furtehr thinned the samples (mostly to reduce file sizes), such that log files containe every 5000th sample and treee files every 20,000th. The outfiles for the main analysis (RLC with BSP) are: combined_ranlc_bsp.log  (480,006 samples) and combined_ranlc_bsp.trees (120,006 trees). I used `tracer` (version 1.7) to then check the effective sample sizes, all were higher than 200, most were much higher.

I use `treeannotator` to make the consensus tree (median heights), mctree.combined-wgs_max_ranlc_bsp.tre. I visualized the tree with:

```bash
java -jar FigTree_v1.4.4/lib/figtree.jar mctree.combined-wgs_max_ranlc_bsp.tre. 
```
I saved tree  plots with time estimates and HPD estiamtes on these.

As noted, my end focus was on the RLC BSP model, but I also have results from the other runs. These did not use coupled MCMC (they mixed well without it). The ORC with BSP results are in  lyc_wgs*, lyc_wgs_max_c2lyc_wgs_max*, lyc_wgs_max_c3lyc_wgs_max*,  and lyc_wgs_max_c5lyc_wgs_max* (four chains) and combined in combined_orclock_bsp* (10% burnin, no thinning beyond the initial), and the single chaing for RLC with EBSP is lyc_wgs_max_ebsp*. These give the same topology.

# Quantifying treeness and identifiying putative cases of admixture with Treemix

Next, I used `treemix` (version 1.13) to quantify how deviations from strict bifurcating trees varied across the genome, how this varied (for autosomes) with chromosome size (and thus with average rates of recombination), and to identify the most noteworthy instances of putative admixture for further evaluation with `caster`.

I first created infiles for `treemix` for each chromosome using the read count data. This was done in the `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/LycAdmix/GenData` directory with the [mkTreeMixin.R](mkTreeMixin.R). The rest of these analyses are in the `TreeMx` subdirectory. 

I ran `treemix` for 0 to xx migration (admixture) events on the full data set from each of the 23 chromosomes.

```bash
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=tmix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


cd /uufs/chpc.utah.edu/common/home/gompert-group5/projects/LycAdmix/TreeMx

perl run_max_treemix.pl treemix_in_ch*gz
```
Which runs [run_max_treemix.pl](run_max_treemix.pl)

```perl
#!/usr/bin/perl
#
## version of treemix fork script that obtains the ML result across $N runs


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

$N = 20;

foreach $fi (@ARGV){
	$fi =~ m/(\d+)/ or die "failed here $fi\n";
	$ch = $1;
	foreach $m (0..8){
		$pm->start and next;
		$out = "tro_ch$ch"."_m$m";
		system "treemix -i $fi -o $out -k 100 -m $m -root MEN12\n";
		open(IN, "$out\.llik");
		$a = <IN>;
		$a = <IN>;
		chomp($a);
		$a =~ m/:\s+(\-[0-9\.]+)/ or die "can't find the ll: $a\n";
		$ll = $1;
		close(IN);
		print "starting ll = $ll\n";
		foreach $i (1..$N){
			$tout = "temp_tro_ch$ch"."_m$m";
			system "treemix -i $fi -o $tout -k 100 -m $m -root MEN12\n";
			open(IN, "$tout\.llik");
			$a = <IN>;
			$a = <IN>;
			chomp($a);
			$a =~ m/:\s+(\-[0-9\.]+)/ or die "can't find the ll: $a\n";
			$llt = $1;
			if($llt > $ll){
				$ll = $llt;
				system "mv $tout\.treeout.gz $out\.treeout.gz\n";
				system "mv $tout\.modelcov.gz $out\.modelcov.gz\n";
				system "mv $tout\.vertices.gz $out\.vertices.gz\n";
				system "mv $tout\.cov.gz $out\.cov.gz\n";
				system "mv $tout\.covse.gz $out\.covse.gz\n";
				system "mv $tout\.edges.gz $out\.edges.gz\n";
				system "mv $tout\.llik $out\.llik\n";
			}
		}
		$pm->finish;
	}
}

$pm->wait_all_children;
```
Importantly, this runs the full ML analysis 20 times for each set of conditions are retains the results with the highest likelihood out of the set of runs. In all cases, I set *P. argus* as the root (MEN12) and used 100 SNPs per block for estimation of the covariance matrix.

Next, I summarized the results in R. 


# Examining genome-wide heterogeneity in relationships with Caster

I am using `caster` (version v1.20.2.5)--Coalescence-aware Alignment-based Species Tree EstimatoR--to quantify variation in the best supported tree topologies across the genome, first among chromosomes and then in a focused manner across chromosomes for subsets of focal taxa (putative cases of admixture or tree discordance). See [Zhang et al. 2025](https://www.science.org/doi/abs/10.1126/science.adk9688).

I am working with the same most likely sequence alignments that I began with for Beast (prior to subsetting) and then once again getting rid of the replicat samples for tree construction. I am not further subsetting the data but am seperately inferring the tree for each chromosome. This is all being done in `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/LycAdmix/Caster`.

```perl
#!/usr/bin/perl
#

foreach $i (1..23){

	print "ASTER-Linux/bin/caster-site -i sub_max_chrom$i.fasta -o cout_max_$i --root MEN --thread 24\n";
	system "ASTER-Linux/bin/caster-site -i sub_max_chrom$i.fasta -o cout_max_$i --root MEN --thread 24\n";
}
```
I am currently using the cater-site model, but I want to look back over the paper and consider the pair of sites model (these essentialy correspond with different models of DNA sequence evolution).

I am then using `ape` (version 5.8) to plot the 23 trees (while rotate around nodes to maximize visula similarity for comparison). See [plotTrees.R](plotTrees.R).

# Window-based analyses with Caster

I have now run a number of sliding window analyses with Caster. For each set of 4 taxa (A, B, C and outgroup) I compute scores in 10 kb windows, then average over sets of 5 windows to plot normalized (sum to 1) scores across the genome. Here is what I have so far:

| A | B | C | P(A+B) | P(A+C) | P(B+C) | Graph | Direcotry |
|---|---|---|--------|--------|--------|-------|-----------|
| BHP | SIN | YG | 0.834 | 0.010 | 0.064 | ![winBHPxSINxYG](https://github.com/user-attachments/assets/2a993fdb-4494-4e98-aeec-610e5f57d686) | winBHPxSINxYG |
| BTB | GNP | SIN | 0.146 | 0.736 | 0.118 | ![winBTBxGNPxSIN](https://github.com/user-attachments/assets/7f1d597c-e83c-430e-b2b8-d85234a2566b) | winBTBxGNPxSIN |
| CLH | TIC | BHP | 0.382 | 0.442 | 0.277 | ![winCHLxTICxBHP](https://github.com/user-attachments/assets/1a56de44-5681-4d95-96d0-6d081f12cf91) | winCLHxTICxBHP |
| CLH | YG | BHP | 0.168 | 0.750 | 0.082 | ![winCHLxYGxBHP](https://github.com/user-attachments/assets/fc27c381-4f2f-4896-9d7b-0cce354658b5) | winCLHxYGxBHP |
| CP | YG | BHP | 0.414 | 0.491 | 0.095 | ![winCPxYGxBHP](https://github.com/user-attachments/assets/f38d5fb5-aca0-4b5a-adb1-24b45d01b4ce) | winCPxYGxBHP |
| CP | YG | TIC | 0.254 | 0.674 | 0.071 | ![winCPxYGxTIC](https://github.com/user-attachments/assets/28e6732d-0014-4efd-90ff-7116dac8215d) | winCPxYGxTIC |
| EP | SHC | TIC | 0.372 | 0.518 | 0.110 | ![winEPxSHCxBHP](https://github.com/user-attachments/assets/ff6b141d-b7f9-40af-a45a-f8e3c0c94d20) | winEPxSHCxTIC |
| HJ | VE | SIN | 0.635 | 0.247 | 0.117 | ![winHJxVExSIN](https://github.com/user-attachments/assets/7fc5995d-1773-4a3d-bb94-c4e99358fb97) | winHJxVExSIN |
| HNV | GNP | SIN | 0.277 | 0.488 | 0.235 | ![winHNVxGNPxSIN](https://github.com/user-attachments/assets/fc9f3599-ab57-4905-a537-c177b54260b8)  | winHNVxGNPxSIN |
| LS | YG | SIN | 0.967 | 0.018 | 0.015 | ![winLSxYGxSIN](https://github.com/user-attachments/assets/54b91e31-3c4a-4e60-8a60-b0298bf59ae5) | winLSxYGxSIN |
| MR | CP | YG | 0.337 | 0.537 | 0.126 | ![winMRxCPxYG](https://github.com/user-attachments/assets/9ac3e06a-eaad-41db-a5dc-a84b3e3f2a55) | winMRxCPxYG |
