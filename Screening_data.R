#Evolution of drug resistance in cancer cells involves generation of numerous mutations in non-coding genome that reduces the chances of DNA breaks
#Authors:
#Santosh Kumar1†, Valid Gahramanov1†, Julia Yaglom1, Shivani Patel1†, Lukasz Kaczmarczyk1, Ivan Alexandrov2, Gabi Gerlitz1, Mali Salmon-Divon1, Michael Y. Sherman1*
#Affiliations 
#1.	Department of Molecular Biology, Ariel University, Israel-40700.
#2.	Research Center of Biotechnology of the Russian Academy of Sciences, Moscow, Russia-119071.

#Set Directory
setwd("~/Documents/PhD studies 2018_Ariel University/Genetic_Screen_Nav_FX")
#Load library
library(dplyr)
library(plyr)
#Files were prepared using custom software that calculates the ratio of barcodes between control and treated samples.
#File is available in the repository
Irino_allgenes <- read.table("~/Documents/PhD studies 2018_Ariel University/Genetic_Screen_Nav_FX/Irino_allgenes.csv", dec=".", quote=",", sep = ",", stringsAsFactors = T)
View(Irino_allgenes)
Irino_allgenes_Subset<- Irino_allgenes[,2:5]

View(Irino_allgenes_Subset)
#Renaming of columns
colnames(Irino_allgenes_Subset)<-c("Genes","Before", "After","Rep_Treat_Control")
Irino_allgenes_matrix<-as.data.frame(Irino_allgenes_Subset)

Control<-Irino_allgenes_matrix[!duplicated(Irino_allgenes_matrix$Genes), ]
fre<-count(Control, vars = "Genes")
bind_control<-merge(Control,fre, by="Genes")
Frequency4<-bind_control %>% select(1,4,5)
colnames(Frequency4)<-c("Genes","value","freq")

View(Irino_allgenes_matrix)

#plot simple representation

library(ggplot2)
bp <- ggplot(Irino_allgenes_matrix, aes(x=Genes, y=Irino_allgenes_matrix$Rep_Treat_Control)) +
  geom_point(aes(colour = cut(Rep_Treat_Control, c(-Inf, 0.3, 3, Inf))),
             size = 1, alpha=0.3) +
  scale_color_manual(name = "Genetic_Screen_Effect",
                     values = c("(-Inf,0.3]" = "red",
                                "(0.3,3]" = "black",
                                "(3, Inf]" = "blue"),
                     labels = c("Sensitizers", "No-Effect", "Protector"))+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank())+
  theme(legend.position="bottom")+
  scale_y_continuous(trans = 'log10')+ 
  xlab("shRNA Barcodes")+
  ylab("Log10(Enrichment of shRNA barcodes)")


library(dplyr)
library(plyr)
#Filter all rows which is less or equal to 0.3 (Sensitizer)
result<-Irino_allgenes_matrix[Irino_allgenes_matrix[, "Rep_Treat_Control"] <= 0.3,]
View(result)
#result0.3_wilcox<-wilcox.test(result$Before, result$After, paired = TRUE, alternative = "two.sided")
#result0.3_T.test<-t.test(result$Before, result$After, paired = TRUE, alternative = "two.sided")
frequency<-count(result, vars = "Genes")

#To avarage all similar rows after cutoff
library(data.table)
keys_sensitizer <- colnames(result)[!grepl('Rep_Treat_Control',colnames(result))]
Irino_Average <- as.data.table(result)
Irino_Average_Filter<-Irino_Average[,list(Average= mean(Rep_Treat_Control)),keys_sensitizer]
View(Irino_Average_Filter)
bind_sensitizer<-merge(Irino_Average_Filter,frequency, by="Genes")
View(bind_sensitizer)
frequency2<-bind_sensitizer[bind_sensitizer$freq>=3,]

library(tidyverse)
Frequency3<-frequency2 %>% select(1,4,5)
Sensitizers <- ddply(frequency2, 'Genes', summarize, value =mean(Average), freq=head(freq,1))

Sensitizers$Genetic_Screen_Effect<-rep("Sensitizers")
#results <- lapply(Sensitizers$value, t.test)
Sensitizers$log2.change<-log2(Sensitizers$value)
nrow(Sensitizers)
Sensitizers<-Sensitizers[!duplicated(Sensitizers$Genes), ]
Protector$Genetic_Screen_Effect<-rep("Protectors")
Sensitizers$Genetic_Screen_Effect<-rep("Sensitizers")

allinone<- as.data.frame(bind_rows(Protector,Sensitizers))
nrow(allinone)
nrow(unique(allinone))
allinone<- allinone %>% select(1,4)
write.csv(allinone, "allinone.csv", quote =F , row.names = T)

####********For protector*********#########


#Filter all rows which is greater or equal to 3(Protector)
result2<-Irino_allgenes_matrix[Irino_allgenes_matrix[, "Rep_Treat_Control"] >= 3,]
View(result2)
#result2_3_wilcox<-wilcox.test(result2$Before, result2$After, paired = TRUE, alternative = "two.sided")
#result2_3_T.test<-t.test(result2$Before, result2$After, paired = TRUE, alternative = "two.sided")

freq<-count(result2, vars = "Genes")

#To avarage all similar rows after cutoff
keys_protector <- colnames(result2)[!grepl('Rep_Treat_Control',colnames(result2))]
Irino_Average_protector <- as.data.table(result2)
Irino_Average_Filter_protector<-Irino_Average_protector[,list(Average= mean(Rep_Treat_Control)),keys_protector]
View(Irino_Average_Filter_protector)
bind_protector<-merge(Irino_Average_Filter_protector,freq, by="Genes")
View(bind_protector)

freq_Protector<-bind_protector[bind_protector$freq>=3,]
freq_Protector$ttest<- t.test(freq_Protector$Before, freq_Protector$After)$p.value
freq_Protector<-freq_Protector[!duplicated(freq_Protector$Genes), ]
freq_Protector$log2.Fold<-log2(freq_Protector$Average)

library(tidyverse)
Frequency4<-freq_Protector %>% select(1,4,5)
Protector <- ddply(Frequency4, 'Genes', summarize, value =mean(Average), freq=head(freq,1))

Protector$log2.change<-log2(Protector$value)
Protector<-Protector[!duplicated(Protector$Genes), ]

nrow(Protector)
Protector$Genetic_Screen_Effect<-rep("Protectors")
bind_rows(Sensitizers,Protector)

Irinotecan_combined<-bind_rows(Sensitizers,Protector)

as.data.frame(Irinotecan_combined)
View(Irinotecan_combined)
nrow(Irinotecan_combined)
Irinotecan<-Irinotecan_combined[!duplicated(Irinotecan_combined$Genes), ]
nrow(Irinotecan)
row.names(Irinotecan)<- Irinotecan$Genes
write.csv(Irinotecan, "Irino_sensi_prot.csv",quote = F)
View(Irinotecan_combined)
library(ggplot2)
library(ggrepel)

# A basic scatterplot with color depending on Species
p1<-ggplot(Irinotecan_combined, aes(x=Genes, y=value, color=Genetic_Screen_Effect)) + 
  geom_point(size=2, alpha=1/3)
p1+ scale_y_continuous(trans = 'log10')+ 
  xlab("Genes")+
  ylab("Log10(Enrichment of shRNA barcodes)")+
  theme_classic()+
  theme(legend.position="bottom")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank())+
  annotate(geom = "text", x = 250, y = 0.02400178, label = "H2AFX", size=2.5)+
  annotate(geom = "text", x = 460, y = 0.19625866, label = "POLE3", size=2.5)+
  annotate(geom = "text", x = 461, y = 0.03635632, label = "POLE4", size=2.5)+
  annotate(geom = "text", x = 335, y = 0.1173669, label = "LIG4", size=2.5)+
  annotate(geom = "text", x = 500, y = 0.13453690, label = "RAD1", size=2.5)+
  annotate(geom = "text", x = 503, y = 0.05944874, label = "RAD9A", size=2.5)+
  annotate(geom = "text", x = 308, y = 0.05838512, label = "KAT5", size=2.5)


