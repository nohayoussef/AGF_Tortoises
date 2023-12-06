library(readxl)
library(phyloseq)
library(ape)
library(plyr)
library(vegan)
library(rbiom)
library(ggplot2)
library(ggpubr)
library(forcats)
library(metagMisc)

#Fig1B,D,E
# create a table with the headers Genera, Age, and Abundance
# This table can be created in excel
data=read.table("/Users/nohayoussef/Desktop/Abundances.txt", header=TRUE)
dist <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/TortReads/Dist/data_sorted.txt", header=TRUE)
Occurrence=read.table("/Users/nohayoussef/Desktop/Occurrences.txt", header=TRUE)
AnimalOrd <-c("Cow", "Deer", "Goat", "Horse", "Sheep", "OtherMammals", "Tortoise")
GeneraOrd <-c("NY36", "NY54", "NY56")
ggplot(dist, aes(x=Group, y=Dist))+
  geom_boxplot(outlier.shape = NA)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ForRare <-read.table("/Users/nohayoussef/Desktop/Rarefaction.txt", header=TRUE, row.names=1)
rarecurve(ForRare, step=1, xlab = "Sample Size", ylab = "Species", label=FALSE)

Coverage <-phyloseq_coverage(PhyseqT)
Coverage
ggplot(Occurrence, aes(x=Animal, y=NY36))+
  geom_bar(stat="identity")+
  scale_x_discrete(limits = AnimalOrd)
ggplot(Occurrence, aes(x=Animal, y=NY56))+
  geom_bar(stat="identity")+
  scale_x_discrete(limits = AnimalOrd)
ggplot(Occurrence, aes(x=Animal, y=NY54))+
  geom_bar(stat="identity")+
  scale_x_discrete(limits = AnimalOrd)



# grouped boxplot Alpha
sig54 = compare_means(formula=NY54~ Animal, data=data, method="wilcox.test")
sig54
sig56 = compare_means(formula=NY56~ Animal, data=data, method="wilcox.test")
sig56
sig36 = compare_means(formula=NY36~ Animal, data=data, method="wilcox.test")
sig36


ggplot(data, aes(x=Animal, y=NY54)) + 
  geom_boxplot() + scale_y_log10()
stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=45)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggplot(data, aes(x=Animal, y=NY56)) + 
  geom_boxplot() + scale_y_log10()
stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=45)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggplot(data, aes(x=Animal, y=NY36)) + 
  geom_boxplot() + scale_y_log10()
stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=45)+
  theme(axis.text.x = element_text(angle=90, hjust=1))


Fungi <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/FungiCopies.txt", header=TRUE)
n_fun <- function(x){
  return(data.frame(y = 0.95*log10(50),
                    label = length(x)))
}

prettyLogs <- function(x){
  pretty_range <- range(x[x > 0])
  pretty_logs <- 10^(-10:10)
  log_index <- which(pretty_logs < pretty_range[2] &
                       pretty_logs > pretty_range[1])
  log_index <- c(log_index[1]-1,log_index, log_index[length(log_index)]+1)
  pretty_logs_new <-  pretty_logs[log_index]
  return(pretty_logs_new)
}

fancyNumbers <- function(n){
  nNoNA <- n[!is.na(n)]
  x <-gsub(pattern = "1e",replacement = "10^",
           x = format(nNoNA, scientific = TRUE))
  exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))
  
  base <- ifelse(exponents == 0, "1", ifelse(exponents == 1, "10","10^"))
  exponents[base == "1" | base == "10"] <- ""
  textNums <- rep(NA, length(n))
  textNums[!is.na(n)] <- paste0(base,exponents)
  
  textReturn <- parse(text=textNums)
  return(textReturn)
}

p= ggplot(data=Fungi, aes(x=Animal, y=Fungi))+
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6, fill = "lightgrey") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
  scale_y_log10(limits = c(0.01, 10000),
                expand = expansion(mult = c(0, 0)),
                labels=fancyNumbers,
                breaks=prettyLogs) +
  annotation_logticks(sides = c("rl")) +
  xlab("Animal") +
  ylab("Copies/ng DNA")
p                

Bacteria <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/Copies_Bacteria.txt", header=TRUE)  

p2= ggplot(data=Bacteria, aes(x=Animal, y=Bacteria))+
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6, fill = "lightgrey") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
  scale_y_log10(limits = c(100000, 100000000000),
                expand = expansion(mult = c(0, 0)),
                labels=fancyNumbers,
                breaks=prettyLogs) +
  annotation_logticks(sides = c("rl")) +
  xlab("Animal") +
  ylab("Copies/ng DNA")
p2

qPCR <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/qPCR.txt", header=TRUE)  
Fungi <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/Fungi.txt", header=TRUE)
Bact <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/Bacteria.txt", header=TRUE)

p3= ggplot(data=qPCR, aes(x=Animal, y=qPCR, fill=Type)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  scale_y_log10(limits = c(10, 100000000000)) +
  xlab("Animal") +
  ylab("Copies/g feces")
p3

p4= ggplot(data=Fungi, aes(x=Animal, y=qPCR)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  scale_y_log10(limits = c(10, 1000000)) +
  xlab("Animal") +
  ylab("Copies/g feces")+
  stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=1000000)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
p4
AOV=aov(qPCR~Animal, data=Fungi)
summary(AOV)
TukeyHSD(AOV)

p5= ggplot(data=Bact, aes(x=Animal, y=qPCR)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  scale_y_log10(limits = c(1000000, 10000000000)) +
  xlab("Animal") +
  ylab("Copies/g feces")+
  stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=10000000000)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
p5
pwc <- Fungi %>% 
  emmeans_test(qPCR ~ Animal, p.adjust.method = "bonferroni") 
pwc



## create phyloseq object with the genera tree
otu_matT <-read_excel("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Physeq_tort2.xlsx", sheet="OTU")
tax_matT <- read_excel("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Physeq_tort2.xlsx", sheet="taxon")
MetaT <-read_excel("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Physeq_tort2.xlsx", sheet="Samples")
otu_matT <- otu_matT %>%
  tibble::column_to_rownames("#OTU ID")
tax_matT <- tax_matT %>%
  tibble::column_to_rownames("#OTU ID")
MetaT <- MetaT  %>%
  tibble::column_to_rownames("Sample")
otu_matT <- as.matrix(otu_matT)
tax_matT <- as.matrix(tax_matT)
OTUT = otu_table(otu_matT, taxa_are_rows = TRUE)
TAXT = tax_table(tax_matT)
samplesT = sample_data(MetaT)
TreeT <-ape::read.tree(file="~/Desktop/GradStudents/Carrie/FecesReptilesPaper/TortReads/Genera_mafft.Aligned2.tree")
PhyseqT <-phyloseq(OTUT, TAXT, samplesT, TreeT)
PhyseqT


##Fig2A
Alpha=read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/TortReads/AlphaT.txt", header=TRUE)
sigObs = compare_means(formula=Observed~ Animal, data=Alpha, method="wilcox.test")
sigObs
sigShannon = compare_means(formula=Shannon~ Animal, data=Alpha, method="wilcox.test")
sigShannon
sigSimpson = compare_means(formula=Simpson~ Animal, data=Alpha, method="wilcox.test")
sigSimpson
sigInvSimpson = compare_means(formula=InvSimpson~ Animal, data=Alpha, method="wilcox.test")
sigInvSimpson


ggplot(Alpha, aes(x=Animal, y=Observed)) + 
  geom_boxplot() +
  stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=45)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggplot(Alpha, aes(x=Animal, y=Shannon)) + 
  geom_boxplot() +
  stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=3.5)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggplot(Alpha, aes(x=Animal, y=Simpson)) + 
  geom_boxplot() +
  stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggplot(Alpha, aes(x=Animal, y=InvSimpson)) + 
  geom_boxplot() +
  stat_compare_means(aes(group=Animal), label="p.signif", hide.ns = TRUE, label.y=12)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
#Fig2B
ordiT = ordinate(PhyseqT, method="PCoA", distance="Wunifrac")
P1 = plot_ordination(PhyseqT, ordiT, "samples", color="HostClass")
P1=P1+geom_point(size=2)+stat_ellipse(lwd=0.5, level=0.95)
P1
write.table(P1$data, file="/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/TortReads/DPCoA.txt")
#Fig2C
P2 = plot_ordination(PhyseqT, ordiT, "samples", color="Animal")
P2=P2+geom_point(size=2)+stat_ellipse(lwd=0.5, level=0.95)
P2
#Fig2D
ordiT2 = ordinate(PhyseqT, method="DPCoA", distance="Wunifrac")

P3 = plot_ordination(PhyseqT, ordiT2, "biplot", color="Animal")
P3=P3+geom_point(size=2)+stat_ellipse(lwd=0.5, level=0.95)
P3

## Calculate distance matrices based on the community composition
TUnifrac_w <-phyloseq::distance(PhyseqT, method="wunifrac")
Factors <-read.table("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Factors.txt", header=TRUE)
Adonis4 <-adonis(TUnifrac_w ~HostClass, Factors)
Adonis4$aov

##Fig4
Cazy <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/CAZY_forBoxPlots.txt", header=TRUE)
CEs <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/CEs.txt", header=TRUE)

p = ggplot(CEs, aes(CE, Value, fill=Group)) +
  geom_boxplot()
p
pwc <- CEs %>% 
  group_by(CE) %>%
  emmeans_test(Value ~ Group, p.adjust.method = "bonferroni") 
pwc
write.table(pwc, file="/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/CEs_Wilcoxon.txt")
GHs <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/GHs.txt", header=TRUE)
p2 = ggplot(GHs, aes(GH, Value, fill=Group)) +
  geom_boxplot()
p2
pwc2 <- GHs %>% 
  group_by(GH) %>%
  emmeans_test(Value ~ Group, p.adjust.method = "bonferroni") 
pwc2
write.table(pwc2, file="/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/GHs_Wilcoxon.txt")


PLs <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/PLs.txt", header=TRUE)
p3 = ggplot(PLs, aes(PL, Value, fill=Group)) +
  geom_boxplot()
p3
pwc3 <- PLs %>% 
  group_by(PL) %>%
  emmeans_test(Value ~ Group, p.adjust.method = "bonferroni") 
pwc3
write.table(pwc3, file="/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/PLs_Wilcoxon.txt")

#FigS4A
Transcripts <-read.table("/Users/nohayoussef/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/TranscriptsStats.txt", header=TRUE)

AOV5=aov(Average_length~Class, data=Transcripts)
summary(AOV5)

AOV7=aov(Average_GC~Class, data=Transcripts)
summary(AOV7)


p1 = ggplot(Transcripts, aes(Class, Average_length, fill=Class)) +
  geom_boxplot()

p2 = ggplot(Transcripts, aes(Class, Average_GC, fill=Class)) +
  geom_boxplot()

#FigS2
## create phyloseq object with the genera tree
otu_matZ <-read_excel("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Paper/Figures/FinalForPNAS/ForSendingOut/Rebuttal_Nat_Commun/ZooCompRed.xlsx", sheet="OTU")
tax_matZ <- read_excel("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Paper/Figures/FinalForPNAS/ForSendingOut/Rebuttal_Nat_Commun/ZooCompRed.xlsx", sheet="taxon")
MetaZ <-read_excel("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Paper/Figures/FinalForPNAS/ForSendingOut/Rebuttal_Nat_Commun/ZooCompRed.xlsx", sheet="Sample")
otu_matZ <- otu_matZ %>%
  tibble::column_to_rownames("#OTU ID")
tax_matZ <- tax_matZ %>%
  tibble::column_to_rownames("#OTU ID")
MetaZ <- MetaZ  %>%
  tibble::column_to_rownames("Sample")
otu_matZ <- as.matrix(otu_matZ)
tax_matZ <- as.matrix(tax_matZ)
OTUZ = otu_table(otu_matZ, taxa_are_rows = TRUE)
TAXZ = tax_table(tax_matZ)
samplesZ = sample_data(MetaZ)

TreeZ <-ape::read.tree(file="~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Paper/Figures/FinalForPNAS/ForSendingOut/Rebuttal_Nat_Commun/Genera_mafft.Aligned.tree")
Tree
PhyseqZ <-phyloseq(OTUZ, TAXZ, samplesZ, TreeZ)
PhyseqZ
ordiZ = ordinate(PhyseqZ, method="PCoA", distance="Wunifrac")
P1 = plot_ordination(PhyseqZ, ordiZ, "samples", color="Class")
P1=P1+geom_point(size=2)+stat_ellipse(lwd=0.5, level=0.95)
P1
P2 = plot_ordination(PhyseqZ, ordiZ, "samples", color="Family")
P2=P2+geom_point(size=2)+stat_ellipse(lwd=0.5, level=0.95)
P2

AlphaZ <-estimate_richness(PhyseqZ, measures=c("Observed", "Shannon", "Simpson", "InvSimpson"))
AlphaZ1=plot_richness(PhyseqZ, x="Class", measures="Shannon", color = "Class")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1,vjust=1,size=12))
AlphaZ2=plot_richness(PhyseqZ, x="Class", measures="Simpson", color = "Class")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1,vjust=1,size=12))
AlphaZ3=plot_richness(PhyseqZ, x="Class", measures="InvSimpson", color = "Class")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1,vjust=1,size=12))
AlphaZ4=plot_richness(PhyseqZ, x="Class", measures="Observed", color = "Class")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1,vjust=1,size=12))

Class_anova.Shannon = aov(AlphaZ$Shannon ~ sample_data(PhyseqZ)$Class)
summary(Class_anova.Shannon)
Class_anova.Shannon
Class_anova.Simpson = aov(AlphaZ$Simpson ~ sample_data(PhyseqZ)$Class)
summary(Class_anova.Simpson)
Class_anova.Simpson
Class_anova.InvSimpson = aov(AlphaZ$InvSimpson ~ sample_data(PhyseqZ)$Class)
summary(Class_anova.InvSimpson)
Class_anova.InvSimpson
Class_anova.Obs = aov(AlphaZ$Observed ~ sample_data(PhyseqZ)$Class)
summary(Class_anova.Obs)
Class_anova.Obs

ZUnifrac_w <-phyloseq::distance(PhyseqZ, method="wunifrac")
Factors <-read.table("~/Desktop/GradStudents/Carrie/FecesReptilesPaper/NonMammalian_TranscriptomeWork/Paper/Figures/FinalForPNAS/ForSendingOut/Rebuttal_Nat_Commun/Factors.txt", header=TRUE)
## adonis (PERMANOVA)
Adonis1 <-adonis(ZUnifrac_w ~Class, Factors)
Adonis1$aov

