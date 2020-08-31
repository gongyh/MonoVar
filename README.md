## Overview ##

**Monovar_NBfork** is a single nucleotide variant (SNV) detection and genotyping algorithm for single-cell DNA sequencing data. It takes a list of bam files as input and outputs a vcf file containing the detected SNVs.

## Dependencies ##
* The fork works with python 2.7 or 3.x
* Python: NumPy v1.8.1 ([http://www.numpy.org/]()), SciPy v0.14.0 ([http://www.scipy.org/]()), Pysam v0.8.1 ([https://code.google.com/p/pysam/]())

## Installation ##

Clone the Monovar repository: 

```
#!python
git clone git@bitbucket.org:hamimzafar/monovar.git
cd monovar
```
Install the Monovar python package:

```
#!python

sudo python setup.py install
```


## Usage ##
The program requires multiple bam files. The bam files should be sorted by coordinates. The raw sequence reads in .fastq format should be aligned to a reference genome with the help of an aligner program (e.g., BWA ([http://bio-bwa.sourceforge.net/]())). Aligner like BWA generates sam files containing aligned reads. The sam files can be converted to compressed bam files using ```samtools view``` command (see Samtools manual for details [http://www.htslib.org/doc/samtools.html]()). 

We have included three sample bam files in the folder examples. To run Monovar, a reference genome file is also needed. Assuming indexed reference genome file to be ref.fa and present in the examples directory, go to the examples directory and run Monovar on the provided bam files as follows:

```
#!python

samtools mpileup -B -d10000 -f ref.fa -q 40 -b filenames.txt | monovar.py -p 0.002 -a 0.2 -t 0.05 -m 2 -f ref.fa -b filenames.txt -o output.vcf
```

> ## NOTE! Do not use the '-Q0' argument, it will inflate the False Positive rates massively!

The arguments of Monovar are as follows:

```
#!python

-i: Pileup file (optional, if no Bam file is provided)
-b: Text file containing the full path for each Bam file. One file per line.
-f: Reference genome file.
-o: Output file.
-t: Threshold to be used for variant calling (Recommended value: 0.05)
-p: Offset for prior probability for false-positive error (Recommended value: 0.002)
-a: Offset for prior probability for allelic drop out (Default value: 0.2)
-m: Number of threads to use in multiprocessing (Default value: 1)
-c: Flag indicating whether to use Consensus Filter (CF) or not (Possible values: 0, 1; Default Value: 1; if 1 then CF is used, otherwise not used)  
-d: Flag indicating debugging mode/no threading (1: enabled)
# Newly added arguments:
-i: pileup file (instead of stdin)
-mpd: Maximum pileup depth to take into account (Default value: 10000)
-mrd: Minimum read depth required for SNV calling (per cell-locus) (Default value = 1)
-th: Heterozygosity rate theta (Default value: 0.001)
-d: Debbuging mode: Turn of threading

```
We recommend using cutoff 40 for mapping quality when using ```samtools mpileup```. To use the probabilistic realignment for the computation of Base Alignment Quality, drop the ```-B``` while running ```samtools mpileup```.
