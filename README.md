# Irinotecan_resistance
Evolution of drug resistance in cancer cells involves generation of numerous mutations in non-coding genome that reduces the chances of DNA breaks
Authors
Santosh Kumar1†, Valid Gahramanov1†, Julia Yaglom1, Shivani Patel1†, Lukasz Kaczmarczyk1, Ivan Alexandrov2, Gabi Gerlitz1, Mali Salmon-Divon1, Michael Y. Sherman1*

Affiliations 
1.	Department of Molecular Biology, Ariel University, Israel-40700.
2.	Research Center of Biotechnology of the Russian Academy of Sciences, Moscow, Russia-119071.
†Equal contribution


Dear reader,
Welcome!!!

The text is intended to give information about the source codes that has been used to conclude the manuscript quoted.
This section is dedicated to provide information about the files that are present in the folders.

List of files and usage:
1.Abstract.md: It contains the abstract for the manuscript.
2.Irino_allgenes.csv: File is the output for the analysis of shRNA screen of Irinotecan. This file need to be used as source file for analysis
3.Screening_data.R: R codes to find out sensitizers and protectors in shRNA screeen.
4.Irino_sensitive_protective_gene_list.csv: Output from the R codes run on Irino_allgenes.csv.
5.Linux_commands_Transcriptome_analysis: Commands that were used in linux server to pre-process the RNA reads and obtain featurecount file.
6.featureCounts.txt: File that is output and modified to be used in R environment.
7.Transcriptome_pipeline.R: Codes required to compute differential expressed genes using featureCounts.txt file as input.
8.Linux_commands_DNA_analysis: It contains all the representative codes that are necessary to compute mutation using GATK. Perform post analysis such as RepeatMasker, demultiplexing of barcodes.
9.triple_mut.txt: Data for common mutations obtained form VCF file. It shows common mutations among the three resistant mutant clones.
10.Karyoplot_FigS3.R: Code to obtain graph for supplment FigS3.

Feel free to contact "kumars@ariel.ac.il" or "shashisantosh2007@gmail.com" or +918600668019/+972559944894 for any clarification on codes troubleshooting.

Thank you!!!

Santosh Kumar
PhD final year, Ariel University, Ariel, Israel
#*******End******
