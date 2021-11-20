#Evolution of drug resistance in cancer cells involves generation of numerous mutations in non-coding genome that reduces the chances of DNA breaks
#Authors:
#Santosh Kumar1†, Valid Gahramanov1†, Julia Yaglom1, Shivani Patel1†, Lukasz Kaczmarczyk1, Ivan Alexandrov2, Gabi Gerlitz1, Mali Salmon-Divon1, Michael Y. Sherman1*
  #Affiliations 
#1.	Department of Molecular Biology, Ariel University, Israel-40700.
#2.	Research Center of Biotechnology of the Russian Academy of Sciences, Moscow, Russia-119071.

#load the library
library(edgeR)
library(limma)

#Set Directory
setwd("~/Documents/JULIA")
#Load data
Santosh_data <- read.delim("featureCounts.txt",header = T)
#View data
View(Santosh_data)
#Rename rows and data mining
row.names(Santosh_data)<-Santosh_data[,1]
Santosh_data_1<-Santosh_data[,-1]
#View polished data
View(Santosh_data_1)
#Print class
print(lapply(Santosh_data_1,class))

#Set conditions
cond <- c("SCC1_Control", "SCC1_Control","MSC1", "MSC1","MSC2","MSC2","MSC3","MSC3")
conditions <- factor(cond, levels=c("SCC1_Control","MSC1","MSC2","MSC3"))
#Set design
design <- model.matrix(~0+conditions)
#Naming of columns
colnames(design) <- levels(conditions)

#Make DGE object
dge <- DGEList(counts=Santosh_data_1, group=factor(cond))
#Filter to remove rows where sum of count for all samples is less than 300
filter<-rowSums(dge$counts) > 300
#Keep rest data
keep <- filterByExpr(dge,design)
dge <- dge[filter, , keep.lib.sizes=FALSE]

View(dge$counts)
counts<-dge$counts
#Normalize the output
normdge <- calcNormFactors(dge)
#Plot MDS to observe the distance between the samples and control (correlation)
plotMDS(dge,main="After filtration of Reads <400",plot = T)

View(normdge$counts)

#Apply voom function
v <- voom(dge,design)
dataE <- v$E
dataF<-v$targets
#Fit the data on linear model
fit2 <- lmFit(v,design)
View(fit2)

#Follow Protocol for analyzing as for a single factor
cont.matrix <- makeContrasts(MSC1-SCC1_Control,MSC2-SCC1_Control,MSC3-SCC1_Control,levels=design)

fit2 <- contrasts.fit(fit2, cont.matrix)
fit2 <- eBayes(fit2)
#MA plots
limma::plotMA(fit2, coef = 1)
limma::plotMA(fit2, coef = 2)
limma::plotMA(fit2, coef = 3)

#Volcano plots (example plots)
volcanoplot(fit2, coef = 1,main="MSC1vsControl")
volcanoplot(fit2, coef = 2,main="MSC2vsControl")
volcanoplot(fit2, coef = 3,main="MSC3vsControl",style = "p-value", highlight = 0, hl.col="blue",
            xlab = "Log2 Fold Change")

#topTable(fit2, number=50 or inf for all)

MSC1<-topTable(fit2, number = Inf, coef=1, adjust="BH",p.value = 0.05)
MSC2<-topTable(fit2, number = Inf, coef=2, adjust="BH",p.value = 0.05)
MSC3<-topTable(fit2, number = Inf, coef=3, adjust="BH",p.value = 0.05, lfc = 1)

#Call for differential expression of genes
genes_differentially_expressed_full_table <- decideTests(fit2, p.value=0.005)

View(genesAll1)

#Annotation to gene name and symbol
library(org.Hs.eg.db)
Average_exp<-dataE
symbols <- mapIds(org.Hs.eg.db, keys=row.names(Average_exp), keytype="ENSEMBL", column="SYMBOL")
Avg_exp <- cbind(Average_exp, symbols)
View(MSC1_name)
symbol1<-mapIds(org.Hs.eg.db, keys=row.names(Avg_exp), keytype="ENSEMBL", column = "GENENAME")
Avg_exp <- cbind(Avg_exp, symbol1)
View(Avg_exp)
write.csv(Avg_exp, file="Avg_exp.csv")

#******Follow same codes for rest two samples******

###########End of codes##########