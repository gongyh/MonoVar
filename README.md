## Overview ##

**MonoVar** is a single nucleotide variant (SNV) detection and genotyping algorithm for single cell DNA sequencing data. It takes a list of bam files as input and outputs a vcf file containing the detected SNVs.

## Dependencies ##

* Python: NumPy, SciPy
* Samtools (included in the external folder)

## Installation ##

Clone the MonoVar repository: 

```
#!python
git clone git@bitbucket.org:hamimzafar/monovar.git
cd monovar

```
Install the MonoVar python package:

```
#!python

sudo python setup.py install
```

Add the samtools binary to the PATH:

```
#!python

export PATH=$PATH:$CURR_DIR/external/samtools/samtools
```

## Usage ##
The program requires multiple bam files. We have included three sample bam files in the folder examples. To run MonoVar, a reference genome file is also needed. Assuming indexed reference genome file to be ref.fa and present in the examples directory, run MonoVar on the provided bam files as follows:

```
#!python

samtools mpileup -BQ0 -d10000 -f ref.fa -q 40 -b filenames.txt | python ../src/monovar.py -n 3 -p 0.002 -a 0.2 -t 0.05 -m 2 -f ref.fa -b filenames.txt -o output.vcf
```
 



### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact