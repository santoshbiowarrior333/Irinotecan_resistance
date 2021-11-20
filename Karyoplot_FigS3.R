#Evolution of drug resistance in cancer cells involves generation of numerous mutations in non-coding genome that reduces the chances of DNA breaks
#Authors:
#Santosh Kumar1†, Valid Gahramanov1†, Julia Yaglom1, Shivani Patel1†, Lukasz Kaczmarczyk1, Ivan Alexandrov2, Gabi Gerlitz1, Mali Salmon-Divon1, Michael Y. Sherman1*
#Affiliations 
#1.	Department of Molecular Biology, Ariel University, Israel-40700.
#2.	Research Center of Biotechnology of the Russian Academy of Sciences, Moscow, Russia-119071.

#Set directory

#BiocManager::install("karyoploteR")
library("karyoploteR")

#Set_directory
setwd("~/Documents/FOCI_Study_Publication")
snps<-read.table("triple_mut.txt",sep="\t",header=T,blank.lines.skip=TRUE,
                 comment.char = "#")
View(snps)
colnames(snps)<-c("chr","start","end","refallele","altallele","qual", "filter","info","format")

#make a GRanges with your data (we need to repeat column 2 as start and end for this to work)
snps.gr <- toGRanges(snps[,c(1,3,3)],comment.char = "#")
available.genomes()
#Begin plotting
kp <- plotKaryotype("hg38",cex = 1.1, chromosomes = "chr1") #<- use the genome you need, if not human


#kpPlotDensity(kp, data=snps.gr, col="blue")
kp <- kpPlotDensity(kp, data=snps.gr, col="blue", window.size = 2e6)
#To add the axis, first get the max density value
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2, cex=1,r0=0, r1=0.4)
kpAddBaseNumbers(kp, add.units = T, cex = 0.8)

#top1: Correlated top1_cutting site reported elsewere (Ref GEO accession number in manuscript)
setwd("~/Documents/PhD studies 2018_Ariel University/Proposal Graphs/Paper Irinotecan")

top1<-read.table("top1.txt",sep="\t",header=F,blank.lines.skip=TRUE,
                 comment.char = "#")
View(top1)
colnames(top1)<-c("chr","start","end")

#make a GRanges with your data (we need to repeat column 2 as start and end for this to work)
top1.gr <- toGRanges(top1[,c(1,2,3)],comment.char = "#")
available.genomes()
#Begin plotting
kp1 <- plotKaryotype("hg38",cex = 1.1) #<- use the genome you need, if not human
#kpPlotDensity(kp, data=snps.gr, col="blue")
kp1 <- kpPlotDensity(kp1, data=top1.gr, col="dark green", window.size = 2e6)

#To add the axis, first get the max density value
max.density <- kp1$latest.plot$computed.values$max.density
kpAxis(kp1, ymin=0, ymax=max.density, numticks = 2, cex=1,r0=0, r1=0.4)
kpAddBaseNumbers(kp1, add.units = T, cex = 0.8)

library(ggplot2)
library(cowplot)

dev.off()
