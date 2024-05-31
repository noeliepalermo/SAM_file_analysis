# Bioinformatics Master 1 project - UE HAI724I Système

# SAM_file_analysis

Project done for the Master in Bioinformatics of Montpellier, France (2022-2024) - Analysis of the different components of a SAM file.

The script sam_analysis.py can also be directly downloaded from the GitHub repositorie of Fatima-Zahra ABANI:https://github.com/abanifatimazahra/M1_Bioinformatic_mapping_sam and executed under Linux(Ubuntu 22.04.1 LTS) with Python 3.10.6.

## README 

Our bioinformatic project consisted in the processing and analysis of the data resulting from the alignment of the sequencing data.
Variation in sequencing data types (paired reads, single reads...) leads to variation in alignment data. These data, which are stored in a SAM format file, can be summarized as follows:

* Read's name : NAME of the request sequence.
* FLAG : Combination of FLAGs at bit level.
* Reference : reference sequence name. Ref is set to '*' when we have no reference.
* Reference's position : The first base of a reference sequence has coordinate 1. POS is set to 0 for an unmapped read.
* Mapping quality : Basic ASCII plus 33 quality.
* CIGAR : Concise Idiosyncratic Gapped Alignment Report. This parameter is formed by an alternative succession of numbers and 	letters that describe the insertions, deletions, matches and mismatches that result from the alignment of each read

-To help you use the Script we have developed, we invite you to read this description and respect the conditions of use.

## SCRIPT EXECUTION

-You must give the name of the SAM file you want to analyze as a first parameter after the execute command of this script.

-This script can only process files :
	With SAM format (Sequence Alignment/Map): the file in question must have the ".sam" extension.
	which are given by the user and exist in his directory.
	Not empty.

-The analysis begins by differentiating between mapped reads (which have a reference and a position values) and unmapped reads (which have a reference equal to '*' and a position equal to 0).

-If you want to perform quality control of data in the sam file based on a specific quality limit value, you should define it as a second parameter.
Otherwise, this quality check will be performed with a default quality value of 20.

-The analysis of the flag will be based on a list that we created on the site https://broadinstitute.github.io/picard/explain-flags.html from the data available on the platform https://www.samformat.info/sam-format-flag.

-The CIGAR check is based on the evaluation of the number of matches with the reference genome, the more matches we have, we say that the read has been well aligned.

## RESULTS
-The results of the analysis will be in the form of several files that contain :
	The referenced reads
	The non referenced reads
	The reads with a quality superior or equal to the input value (default 20)
	The reads with a quality strictly inferior to the input value (default 20)
	The reads paired with the correct FLAG references
	The reads no paired with the incorrect FLAG references
	The reads with a CIGAR superior or equal to 1M
	The reads with a CIGAR different to 1M
	
-The results of the analysis are printed on the terminal as:
```
Total number of reads
Number of referenced
Number of reads with a quality superior to  20 and reads with a quality 
Number of paired reads
Number of aligned reads
```
## TEST
								
-The application of our script on the file mapping.sam allowed us to obtain the following results for 264 reads :
```
Total number of reads:  264
Number of referenced reads:  264
Number of reads with a quality superior to  20 :  238
Number of paired reads:  264
Number of aligned reads:  264
```
