[![Anaconda-Server Badge](https://anaconda.org/gabaldonlab/crossmapper/badges/platforms.svg)](https://anaconda.org/gabaldonlab/crossmapper)

![alt text](https://image.ibb.co/bs7fAV/logos.png)



## General descripiton of Crossmapper

Crossmapper is an automated bioinformatics pipeline for asessing the rate of read crossmapping when two or more organisms are sequenced as one sample. The software can be used for planning such kind of experimental setups as dual- or multiple RNA-seq (mainly for host-pathogen, symbiont and cohabitant interaction studies), metagenomics studies, sequencing and analysis of hybrid species, allele-specific expression studies, and can be extended for the use in large sequencing facilities for resource optimization.

Based on in-silico read simulation and back-mapping to the original genomes of sequenced organisms, Crossmapper allows the users to assess the rate of incorrect unique and multimapped reads to non-corresponding genomes and thus helps to optimize the sequencing parameters such as the read lenght, paired/single-end, mapping parameters, etc., prior of performing the actual sequencing experiment.

### Developers

Crossmapper is developed by:

Hrant Hovhannisyan (Comparative Genomics group, Centre for Genomic Regulation, Barcelona, Spain, contact email: grant.hovhannisyan@gmail.com) and

Ahmed Hafez (Biotechvana, Valencia, Spain, contact email: ah.hafez@gmail.com ). 

Both are part of [EU funded Innovative Training Network OPATHY - From Omics to Patient: Improving Diagnostics of Pathogenic Yeasts](https://www.opathy.eu/).

### Licensing
Crossmapper is distribued under the [GNU GENERAL PUBLIC LICENSE v3](https://github.com/GrantHov/crossmap/blob/master/LICENSE).


## Installation and setup
We have implemented the Crossmapper in Python 3.6 as an Anaconda package. Thus, the Crossmapper installation is as easy as any other Anaconda package, without a need to solve any dependency issues. 

If anaconda or miniconda package managers are installed, simply run:

```bash
conda install -c gabaldonlab -c bioconda crossmapper
```

to install crossmapper and all its dependencies.

The examplified procedure of Crossmapper installation and optional specific enviroment creation can be found in our step-by-step tutorial. 

## Usage and options

The basic usage arguments of Crossmapper can be found by typing `crossmapper -h` in the command line:
```
usage: crossmapper [-h] [-v] {DNA,RNA} ...

-- crossmapper Software

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

SimulationType:
  {DNA,RNA}      Simulation type. Choose to simulate either DNA or RNA data
    DNA          Simulate DNA data
    RNA          Simulate RNA data
```
Crossmapper has two running options: DNA and RNA. By running, for example, `crossmapper RNA -h`, the user can see the parameters for RNA mode. Optional arguments are the same for DNA mode, but on the bottom of the help page the user can find the arguments specific only for RNA mode.

```
usage: crossmapper RNA [-h] -g fasta [fasta ...] [-t int] [-e float] [-d int]
                       [-s int]
                       (-N int [int ...] | -C float/int [float/int ...])
                       [-rlay {SE,PE,both}] [-rlen int] [-r float] [-R float]
                       [-X float] [-S int] [-AMB float] [-hapl] [-o PATH]
                       [-gb] [-rc] [-gn name [name ...]]
                       [--mapper-template PATH] [-max_mismatch_per_len float]
                       [-bact_mode] [-max_mismatch int] -a gtf [gtf ...]
                       [-star_tmp PATH]

optional arguments:
  -h, --help            show this help message and exit
  -t int, --threads int
                        Number of cores to be used for all multicore-
                        supporting steps (default: 1)
  -e float, --error float
                        Base error rate (default: 0.02)
  -d int, --outer_dist int
                        Outer distance between the two reads. For example, in
                        case of 2x50 reads, d=300 and s=0 the mates will be
                        200 bp apart from each other. (default: 500)
  -s int, --s_dev int   Standard deviation of outer distance. (default: 30)
  -N int [int ...], --N_read int [int ...]
                        The number of reads/read pairs to generate. This
                        parameter can not be used alongside with -C (default:
                        None)
  -C float/int [float/int ...], --coverage float/int [float/int ...]
                        Generate the number of reads that reaches the
                        specified coverage. Coverage is calculated as:C =
                        N*rlen/L, where L is the length of the
                        genome/transcriptome (default: None)
  -rlay {SE,PE,both}, --read_layout {SE,PE,both}
                        Specify the read configuration - single-end (SE),
                        paired-end (PE), or both (both). If chosen 'both', the
                        software will make separate analysis with each
                        configuration (default: SE)
  -rlen int, --read_length int
                        Specify the read length. Choose from the possible read
                        lengths available for Illumina machines:
                        50,75,100,125,150,300. The user can either enter a
                        specific length, or specify a (!) COMMA-SEPARATED (no
                        spaces are allowed between commas) list of desired
                        read lengths. In the latter case, the software will
                        perform the analysis for all specified values
                        separatelly and will report mapping statistics in a
                        form of a graph (default: 50)
  -r float, --mut_rate float
                        Mutation rate. (default: 0.001)
  -R float, --indel_fraction float
                        Fraction of indels. (default: 0.015)
  -X float, --indel_extend float
                        Probability of an indel to be extended. (default: 0.3)
  -S int, --random_seed int
                        Seed for random generator. (default: -1)
  -AMB float, --discard_ambig float
                        Disgard if the fraction of ambiguous bases is higher
                        than this number. (default: 0.05)
  -hapl, --haplotype_mode
                        Haplotype mode. If specified, the haploid mutations
                        will be simulated instead of diploid. (default: False)
  -o PATH, --out_dir PATH
                        Specify the output directory for crossmap output
                        files. (default: crossmap_out)
  -gb, --groupBarChart  Use a grouped bar chart in the output report instead
                        of individual bar chart. (default: False)
  -rc, --reportCrossmapped
                        Report all cross mapped reads into csv file. (default:
                        False)
  -gn name [name ...], --genome_names name [name ...]
                        Specify names of the genomes. The names will appear in
                        the report file. (default: None)
  --mapper-template PATH, --mapper-template PATH
                        --mapper-template (default: None)

Required Arguments:
  -g fasta [fasta ...], --genomes fasta [fasta ...]
                        Specify the genome files in fasta format. Enter genome
                        names separated by whitespace. NOTE: Keep the same
                        order of listing for gtf/gff files (default: None)

Mapper and annotation Arguments:
  Arguments specific to STAR Mapper

  -max_mismatch_per_len float, --outFilterMismatchNoverReadLmax float
                        From STAR manual: alignment will be output only if its
                        ratio of mismatches to *read* length is less than or
                        equal to this value: for 2x100b, max number of
                        mismatches is 0.04*200=8 for the paired read.
                        (default: 0.04)
  -bact_mode, --bacterial_mode
                        This option prohibits spliced alignments for STAR and
                        it can be used for mapping bacterial data. (default:
                        False)
  -max_mismatch int, --outFilterMismatchNmax int
                        From STAR manual: alignment will be output only if it
                        has no more mismatches than this value (default: 10)
  -a gtf [gtf ...], --annotations gtf [gtf ...]
                        Specify the gtf/gff files. Enter the file names
                        separated by whitespace. NOTE: Keep the same order of
                        listing as for genome files (default: None)
  -star_tmp PATH, --star_temp_dir PATH
                        Specify a full path to a local temprorary directory,
                        where all intermediate files of STAR will be written.
                        This option can be used when running Crossmaper from a
                        local machine in a server or cluster with SAMBA
                        connection. (default: ./TMPs)
                        
```


## Step-by-Step tutorial

We have created a [Step-by-Step tutorial](https://github.com/Gabaldonlab/crossmapper/wiki/Crossmapper-Tutorial)  with detailed description of options and workflow steps performed by Crossmapper. We guide the user from the initial software installation and setup untill the interpretations of the obtained results. 

## Precomputed values

To save time for users, we have precomputed the cross-mapping values for some widely sequences species, namelly human combined with each of the following species: *[M. musculus](https://www.dropbox.com/s/18ne29pr3l7xcjt/report_human_mouse.html?dl=0), [A. thaliana](https://www.dropbox.com/s/4yy02qhnlv9z5xz/report_human_arabidopsis.html?dl=0), [D. melanogaster](https://www.dropbox.com/s/b0h4agubz4r3b48/report_human_fly.html?dl=0), [C. elegans](https://www.dropbox.com/s/d6bzsepljo0cxlr/report_nematode.html?dl=0), [S. cerevisiae](https://www.dropbox.com/s/bnsryeg7rgbv4ll/report_human_yeast.html?dl=0)* and *[S. typhimurium](https://www.dropbox.com/s/30kmknnhkf4ullk/report_human_salmonella.html?dl=0)* (please download the files to view correctly). We will expand these list of precomputed values in future.


## Contact and reporting

The users of Crossmapper can report bugs/issues in our [Github Issues section](https://github.com/Gabaldonlab/crossmapper/issues), or alternativelly can ask questions, report bugs and suggest new features for Crossmapper in our [google group](https://groups.google.com/forum/#!forum/crossmapper-users).
