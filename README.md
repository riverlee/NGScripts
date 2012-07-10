NGScripts--A collection of scripts for NGS data analysis.
=========

## Table of contents
1. get_reads_lane_info.pl
2. vcfStat.pl


##1. get_reads_lane_info.pl
To fetch the assay information (run, lane etc.) from a fastq format file.
The information is based on the head line of each read. Currently it only
supports the format comes from Illumina sequencing.

```
Usage: perl get_reads_lane_info.pl -i/--input <inputfile> -c/--casava18 -h/--help
          -i/--input    input fastq file
          -c/--casava18 indicate the header is in Casava 1.8 format 
          -h/--help     help page
--------------------------------------------------------------------------------------------
There are two kinds of header format
1) Default header format(see at http://en.wikipedia.org/wiki/FASTQ_format)
e.g., @HWUSI-EAS100R:6:73:941:1973#0/1
========================================
HWUSI-EAS100R	the unique instrument name
6               flowcell lane
73              tile number within the flowcell lane
941             'x'-coordinate of the cluster within the tile
1973            'y'-coordinate of the cluster within the tile
#0              index number for a multiplexed sample (0 for no indexing)
/1              the member of a pair, /1 or /2 (paired-end or mate-pair reads only)

2) Casava1.8 format
e.g., @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
==============================================================
EAS139    the unique instrument name
136       the run id
FC706VJ   the flowcell id
2         flowcell lane
2104      tile number within the flowcell lane
15343     'x'-coordinate of the cluster within the tile
197393    'y'-coordinate of the cluster within the tile
1         the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
Y         Y if the read fails filter (read is bad), N otherwise
18        0 when none of the control bits are on, otherwise it is an even number
ATCACG    index sequence

Output format
1) output of a default header
=============================
#Instrument LaneID   TotalReads
@HWI-ST829  156      93702907

2) output of a casava1.8 header
===============================
#Instrument RunID LaneID   TotalReads  FailedReads
@HWI-ST829  156     3       93702907    11125894
```
##2. vcfStat.pl
To get a summary of SNPs from a VCF file. Information includes Ti/Tv, novel SNPs etc.
[Example](https://github.com/riverlee/NGScripts/raw/master/vcf_report.html)

```
Usage: perl vcfStat.pl -i/--input <inputfile> -d/--dbsnp <dbsnp file> -n/--nonsynonymous <non-synonymousfile> -h/--help -nf --qual --filter --gq --dp
          -i/--input          input vcffile (*required)
          -d/--dbsnp          input dbsnp file either in vcf format or a 'list'(chr pos [rsid] [ref] [alt]), default is vcf format (*required)
          -df/--dbsnpformat   if this option is provided, it means the dbsnp file in a 'list' format (optional)
          -n/--nonsynoymous   non-synonymous snps file annotated by ANNOVAR. The default format assumes the 2nd, 4th and 5th columns
                              are 'snptype(synonymous, non-synonymous etc.)','chr' and 'pos'. If option provided with '-nf', it means 
                              the input is in a 'list' format (chr pos [rsid] [ref] [alt]) (*required)
          -nf                 if this option is provided, it means the non-synonymous file in a 'list' format (optional)
          -nsoft              if this option provided, 'non-synonymous' SNPs not include those 'stopgain' and 'stoploss' ones. Default inclues those two types. (optional)
          --qual              'qual' filter for the input vcf file (qual < given cutoff will be skipped). The qual filter is the 6th column in a vcf.(optional)
          --filter            'filter' filter for the input vcf file, only FILTER name eq to 'filter' in the 7th columns will be kept. Default is 'PASS'. 
                              Give '*' to make no filter. Multile filter names could be give by 'filter1|filter2|filter3'.(optional)
          --gq                Genotype quality filter for each sample. If a given samples' GQ less than cutoff, the genotype will be masked as './.'. (optional)
          --dp                Depth filter for each sample. If a given samples' DP less than cutoff, the genotype will be masked as './.'. (optional)
          -o/--out            output filename for the report, default is vcf_report.html (optional)
          -h/--help           help page
```











