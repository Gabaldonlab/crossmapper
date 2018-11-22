![alt text](https://image.ibb.co/bs7fAV/logos.png)

## General descripiton of Crossmapper

Crossmapper is an automated bioinformatics pipeline for asessing the rate of read crossmapping when two or more organisms are sequenced as one sample. The software can be used for planning such kind of experimental setups as dual- or multiple RNA-seq (mainly for host-pathogen, symbiont and cohabitant interaction studies), metagenomics studies, sequencing and analysis of hybrid species, allele-specific expression studies, and can be extended for the use in large sequencing facilities for resource optimization.

Based on in-silico read simulation and back-mapping to the original genomes of sequenced organisms, Crossmapper allows the users to assess the rate of incorrect unique and multimapped reads to non-corresponding genomes and thus helps to optimize the sequencing parameters such as the read lenght, paired/single-end, mapping parameters, etc., prior to performing the actual sequencing experiment.

### Developers

Crossmapper is developed by 

Hrant Hovhannisyan (Toni Gabaldon's group, Centre for Genomic Regulation, Barcelona, Spain, contact email: grant.hovhannisyan@gmail.com) and

Ahmed Hafez (Biotechvana, Valencia, Spain, contact email: ). 

Both are a part of [EU funded Innovative Training Network OPATHY - From Omics to Patient: Improving Diagnostics of Pathogenic Yeasts](https://www.opathy.eu/).

### Licensing
Crossmapper is distribued under the [GNU GENERAL PUBLIC LICENSE v3](https://github.com/GrantHov/crossmap/blob/master/LICENSE).


## Installation and setup
We have implemented the Crossmapper in Python 3.6 as an Anaconda package. Thus, the Crossmapper installation is as easy as any other Anaconda package, without a need to solve any dependency issues. 

If anaconda or miniconda package managers are installed, simply run:

`conda install -c bioconda crossmapper`

to install crossmapper and all its dependencies.

The examplified procedure of Crossmapper installation and optional specific enviroment creation can be found in our step-by-step tutorial. 

## Usage and options

The basic usage arguments of Crossmapper can be found by typing `crossmap.py -h` in the command line:
```
usage: crossmap.py [-h] [-v] {DNA,RNA} ...

-- crossmap.py Software

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

SimulationType:
  {DNA,RNA}      Simulation type. Choose to simulate either DNA or RNA data
    DNA          Simulate DNA data
    RNA          Simulate RNA data
```
Thus, Crossmapper has two running options: DNA and RNA. By running, for example, `crossmap.py RNA -h`, the user can see the parameters for RNA mode. Optional arguments are the same for DNA mode, but on the bottom of the help page the user can find the arguments specific only for RNA mode.

```
usage: crossmap.py RNA [-h] -g GENOMES GENOMES [-t THREADS] [-e ERROR]
                       [-d OUTER_DIST] [-s S_DEV]
                       (-N N_READ N_READ | -C COVERAGE COVERAGE)
                       [-rlay {SE,PE,both}] [-rlen READ_LENGTH] [-r MUT_RATE]
                       [-R INDEL_FRACTION] [-X INDEL_EXTEND] [-S RANDOM_SEED]
                       [-A DISCARD_AMBIG] [-hapl] [-o OUT_DIR]
                       [-max_mismatch Int] -a ANNOTATIONS ANNOTATIONS

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of cores to be used for all multicore-
                        supporting steps (default: 1)
  -e ERROR, --error ERROR
                        Base error rate (default: 0.02)
  -d OUTER_DIST, --outer_dist OUTER_DIST
                        Outer distance between the two reads. For example, in
                        case of 2x50 reads, d=300 and s=0 the mates will be
                        200 bp apart from each other. (default: 500)
  -s S_DEV, --s_dev S_DEV
                        Standard deviation of outer distance. (default: 30)
  -N N_READ N_READ, --N_read N_READ N_READ
                        The number of reads/read pairs to generate. This
                        paremeter can not be used alongside with -C (default:
                        None)
  -C COVERAGE COVERAGE, --coverage COVERAGE COVERAGE
                        Generate the number of reads that reaches the
                        specified coverage. Coverage is calculated as:C =
                        N*rlen/L, where L is the length of the
                        genome/transcriptome (default: None)
  -rlay {SE,PE,both}, --read_layout {SE,PE,both}
                        Specify the read configuration - single-end (SE),
                        paired-end (PE), or both (both). If chosen 'both', the
                        software will make separate analysis with each
                        configuration (default: SE)
  -rlen READ_LENGTH, --read_length READ_LENGTH
                        Specify the read length. Choose from the possible read
                        lengths available for Illumina machines:
                        25,50,75,100,125,150,300. The user can either enter a
                        specific length, or specify a (!) COMMA-SEPARATED (no
                        spaces are allowed between commas) list of desired
                        read lengths. In the latter case, the software will
                        perform the analysis for all specified values
                        separatelly and will report mapping statistics in a
                        form of a graph (default: 50)
  -r MUT_RATE, --mut_rate MUT_RATE
                        Mutation rate. (default: 0.001)
  -R INDEL_FRACTION, --indel_fraction INDEL_FRACTION
                        Fraction of indels. (default: 0.015)
  -X INDEL_EXTEND, --indel_extend INDEL_EXTEND
                        Probability of an indel to be extended. (default: 0.3)
  -S RANDOM_SEED, --random_seed RANDOM_SEED
                        Seed for random generator. (default: -1)
  -A DISCARD_AMBIG, --discard_ambig DISCARD_AMBIG
                        Disgard if the fraction of ambiguous bases is higher
                        than this number. (default: 0.05)
  -hapl, --haplotype_mode
                        Haplotype mode. If specified, the haploid mutations
                        will be simulated instead of diploid. (default: False)
  -o OUT_DIR, --out_dir OUT_DIR
                        Specify the output directory for crossmap output
                        files. (default: crossmap_out)

Required Arguments:
  -g GENOMES GENOMES, --genomes GENOMES GENOMES
                        Specify the genome files in fasta format. Enter genome
                        names separated by whitespace. NOTE: Keep the same
                        order of listing for gtf/gff files (default: None)

Mapper and annotation Arguments:
  Arguments specific to STAR Mapper

  -max_mismatch Int, --outFilterMismatchNmax Int
                        From STAR manual: alignment will be output only if it
                        has no more mismatches than this value (default: 10)
  -a ANNOTATIONS ANNOTATIONS, --annotations ANNOTATIONS ANNOTATIONS
                        Specify the gtf/gff files. Enter the file names
                        separated by whitespace. NOTE: Keep the same order of
                        listing as for genome files (default: None)
```


## Step-by-Step tutorial

We have created a [Step-by-Step tutorial](https://github.com/GrantHov/crossmap/wiki/Step-by-step-Crossmapper-usage-tutorial)  with detailed description of options and workflow steps performed by Crossmapper. We guide the user from the initial software installation and setup untill the interpretations of the obtained results. 


## Contact and reporting

The users of Crossmapper can report bugs/issues in our [Github Issues section](https://github.com/GrantHov/crossmap/issues), or alternativelly can ask questions, report bugs and suggest new features for Crossmapper in our [google group](https://groups.google.com/forum/#!forum/crossmapper-users).
