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

`Beast` cannot handle the full set of SNPs. I thus subset the SNP data from the alignments. This involves first generating a filtered and reduced set of SNPs with [GetSNPSubstMax.R](GetSNPSubstMax.R), which uses coverted (to numeric) alignments generated by [mkNumericFasta.pl](mkNumericFasta.pl). I specifically retain 0.025 percent of the SNPs from each chromosome and only those with no missing data. Also at this point, I dropped chromsome 21.2 (the small bit that might be part of 21). This outputs the keepSNPs_max_chrom files. I then create the combined (across chromosomes) fasta alignment file with just the subset of SNPs, lyc_genomemax.fasta. This was done with [SubSetFasta.pl](SubSetFasta.pl). I then converted this to a nexus format alignment for `Beast`.

```bash
seqmagick convert --output-format nexus --alphabet dna lyc_genomemax.fasta lyc_genomemax.nex
```

I have been playing with `Beast` to assess the effects of different priors and setting on perfomance, especially in terms of how to appropriately use the time calibrarions. This is all entered via `beauti`.

```bash
ml beast
beauti
```

Part of my approach involves entering the (approximate) number of invariant sites. For this, I am computing the total number of each base (A, C, G and T) in the genome, reducing this to 0.025 of the total (to acount for my subsampling) and subtracting of the SNPs. The current version doesn't use the final SNP set, so I will need to update this slightly. I get the counts (for the genome and SNP set) with [countBases.pl](countBases.pl).

I then add the conts by appending "Orig" to the name, "lyc_genomemax' in the xml from `beauti`,  and adding the following just after the data (number of constant A, C, G and T):

```xml
    <data id='lyc_genomemax' spec='FilteredAlignment' filter='-' data='@lyc_genomemaxOrig' constantSiteWeights='35769 20210 20276 35745'>
    </data>
```

I am using two sources for the calibration dates.

1. [Vila et al. 2011](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2010.2213) focuses on the colonization of the Americas by *Polyommatus* blues. This is older and has more limited data (a few mitochondrial and nuclear genes) but includes a good sample of old and new world *Lycaeides*. The TMRCA for *Lycaeides* (old and new world) was estimated at 2.4 million years.

2. [Kawahara et al. 2023](https://www.nature.com/articles/s41559-023-02041-9) uses many more genes and species for a global phylogeny of butterflies. With that siad, this only has *L. melissa*, *L. anna* and *L. idas* from the contiguous USA, and thus the TMRCA for these likely underestimates that of our set of *Lycaeides* (which includes Alaska, and this appears to be the deepest split in the ingroup for our tree). Figure S1 in the paper gives the TMRCA for *Lycaeides* as 1.29 million years and for *Lycaeides* + *Plebejus argus* as 5.22 million years (this is our outgroup, their sample is from Japan, ours from France). My most recent runs set normal priors on both of these nodes based on this data. For *Lycaeides*, I used mean = 1.845 and SD = 0.283, which gives 95 density intervals of 1.29 to 2.4 (my dates). This is worth thinking about a bit more, but is likely about right. I used 5.22 for the mean for *Lycaeides* + *P. argus* and the same SD, 0.283, to reflect a similar level of uncertainty (again worth thinking about).

Other than that, here are my main thoughts and current choices for priors (which seem to be working okay, though the terminal branches are still probably too long, but their is also uncertainty):

- HKY substitution model with kappa estimated and gamma rate heterogeneity approximated by four categoreis (everything else default). I probably want this or GTR; it is likely shallow enough for HKY. Read a bit more on this.

- I am currently using the Optimised Relaxed Clock, which is more flexible than the Strict Clock but a bit less flexible than the Random Local Clock. This basically allows rates to change along the tree but in a correlated way. I want to read and think more about this. I am **estimating the clock rate**, which is important given my calibration information. This requires messing with settings in the file menu to be allowed. When I don't do this, I get very wonky results. I could instead maybe set it to 0.0029 based on the *Heliconius* mutation rate[ Keightley et al. 2015](https://academic.oup.com/mbe/article/32/1/239/2925597), but estimating it and using the calibration nodes is probably the way to go.

- I am using the Coalescent Bayesian Skyline prior, which makes sense for cases where you have multiple populations or individuals from the same species. I might also consider the Extended version of this. They differ in terms of where (at nodes or wherever) Ne changes. I left the rest of the priors at their defaults, but this is worth revisiting.

- I am using the defaults from MCMC so far.

I then run `Beast` with:

```bash
ml beast
beast lyc_genomemax.nex
```

I use `treeannotator` to make the consensus tree (median heights) but should also run some MCMC diagostics. I have been visualizing the tree with:

```bash
java -jar FigTree_v1.4.4/lib/figtree.jar mctree.combined-wgs_max.tre 
```

The latter is the output tree from `treeannotator`.

# Quantifying treeness and identifiying putative cases of admixture with Treemix

# Examining genome-wide heterogeneity in relationships with Caster
