library(tidyverse)
require(reshape2)
library(ggpubr)
library(seqinr)

dir.create("04_VariantsEffects")
path="04_VariantsEffects"

#required data:
load("nscore_df_ADan.RData")


#
peptide_seq<-"EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY"
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], '\n', n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "*")

min<-min(singles_stops$nscore_c)
max<-max(singles_stops$nscore_c)
cols <- c(colorRampPalette(c( "brown3", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(-min+max)*100)-0.5))

### SNVs

nt_seq<-"GAAGCCAGCAATTGTTTCGCAATTCGGCATTTTGAAAACAAATTTGCCGTGGAAACTTTAATTTGTTTTAATTTGTTCTTGAACAGTCAAGAAAAACATTAT"

SNVs_maker<-function(nt_seq){
  nt_seq<-c(strsplit(nt_seq, '')[[1]])
  nucleotides<-c("G", "A", "C", "T")
  HGVSc<-c()
  nt_seq_SNVs<-c()
  aa_seq_SNVs<-c()
  for (i in 1:length(nt_seq)){
    new_seq<-nt_seq
    for (nt in nucleotides){
      if (nt_seq[i] != nt){
        HGVSc<-c(HGVSc, paste0("c.", i+729, nt_seq[i], ">", nt))
        new_seq[i]<-nt
        nt_seq_SNVs<-c(nt_seq_SNVs, paste(new_seq, collapse=""))
        aa_seq<-paste(translate(new_seq), collapse="")
        aa_seq_SNVs<-c(aa_seq_SNVs, paste(aa_seq, collapse=""))
      }
    }
    df<-data.frame("HGVSc"=HGVSc, "nt_seq"=nt_seq_SNVs, "aa_seq"=aa_seq_SNVs)
  }
  return(df)
}

SNVs<-SNVs_maker(nt_seq)

##################################################################################

#####  heatmap (contains stops)

#add syn
positions<-(1:34)
syn.df<-data.frame(
  "aa_seq"="EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY",
  "WT_AA"= peptide_seq,
  "Mut"= peptide_seq,
  "Pos"= positions, 
  "sigma"=0, 
  "nscore_c"=0,
  "ID"="syn",
  "mean_count"=NA
  
)

heatmap_df<-rbind(singles_stops[,c("aa_seq","WT_AA", "Mut", "Pos", "sigma", "nscore_c", "ID", "mean_count")], syn.df)

heatmap_df$WT<-""
heatmap_df[heatmap_df$ID=="syn",]$WT<-"ADan"

heatmap_df$SNVs<-""
heatmap_df[(heatmap_df$aa_seq %in% SNVs$aa_seq) & (heatmap_df$ID != "syn"), "SNVs"]<-"|"

heatmap_df$Comb_mut<-""
heatmap_df[heatmap_df$ID=="E-1-K",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="S-3-N",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="C-5-Y",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="F-6-L",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="A-7-S",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="A-7-T",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="I-8-V",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="R-9-Q",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="K-14-E",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="A-16-V",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="V-17-M",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="L-20-S",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="L-20-V",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="I-21-V",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="S-29-G",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="K-32-T",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="H-33-Y",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="Y-34-H",]$Comb_mut<-"|"


p_heatmap<-ggplot(heatmap_df)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)), fill=nscore_c), color='white', size=1)+
  geom_text(aes(Pos, Mut, label=WT), size=3)+
  geom_text(aes(Pos, Mut, label=SNVs), color="black", size=3, nudge_x=0.32, nudge_y=0.37, angle=45, fontface="bold")+
  geom_text(aes(Pos, Mut, label=Comb_mut), color="black", size=3, nudge_x=0.26, nudge_y=0.32, angle=45, fontface="bold")+
  scale_x_continuous(breaks=seq(1,34), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  labs(x="ADan peptide", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=16))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60",breaks=c(-7, -4, -1, 1))
p_heatmap

ggsave(p_heatmap, path=path, file="p_heatmap_nscore.jpg", width=15, height=8)

##### Violin plot per position
subsp_median_df<-as.data.frame(heatmap_df %>% group_by(Pos) %>% dplyr::summarise(median_p=median(nscore_c)))
heatmap_df<-left_join(heatmap_df, subsp_median_df)

p_violin_p<-ggplot(heatmap_df, aes(x=as.factor(Pos), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", size=0.2, aes(fill=median_p))+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  scale_x_discrete(breaks=seq(1,34), labels = peptide_seq_pos, expand=c(0,0))+
  theme_bw()+
  labs(x="ADan peptide", y="Nucleation score", fill="Median NS")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-7, -4, -1, 1))

p_violin_p

ggsave(p_violin_p, path=path, file="p_violinplotp_nscore.jpg", width=15, height=8)

##### violin per change

subsm_median_df<-as.data.frame(heatmap_df %>% group_by(Mut) %>% dplyr::summarise(median_m=median(nscore_c)))
heatmap_df<-left_join(heatmap_df, subsm_median_df)

p_violin_m<-ggplot(heatmap_df, aes(x=factor(Mut, levels=rev(vectorAA)), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", size=0.2, aes(fill=median_m))+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  theme_bw()+
  labs(x="ADan peptide", y="Nucleation score", fill="Median NS")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-7, -4, -1, 1))+ 
  coord_flip() 

p_violin_m

ggsave(p_violin_m, file="p_violinplot_mut.jpg",width=4, height=8)

#### stacked bars
# Considering stops
fdr_categories<-singles_stops[,c("aa_seq", "Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10")]

fdr_categories$category<-"WT-like"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS-"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+"

fdr_categories$SNVs<-F
fdr_categories[fdr_categories$aa_seq %in% SNVs$aa_seq, "SNVs"]<-T

categories_p <- fdr_categories %>% group_by(Pos,category) %>% dplyr::summarise(Freq_p=n()) 

levels_cat = c("NS-", "WT-like", "NS+")
colors_cat<-c("#DF9292", "#F2F2F2", "#7979BE")

p_categories_p<-ggplot(categories_p, aes(fill=factor(category, levels=levels_cat), x=factor(Pos), y=Freq_p))+
  scale_x_discrete(breaks=seq(1,34), labels = peptide_seq_pos, expand=c(0,0))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "ADan peptide", y="Counts")

p_categories_p

ggsave(p_categories_p, path=path, file="p_categoriesp_nscore.jpg", width=15, height=8)

# changes from

categories_from <- fdr_categories %>% group_by(WT_AA, category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) )


levels_cat = c("NS-", "WT-like", "NS+")
colors_cat<-c("#DF9292", "#F2F2F2", "#7979BE")

p_categories_from<-ggplot(categories_from, aes(fill=factor(category, levels=levels_cat), x=factor(WT_AA), y=freq))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "WT AA", y="Frequency")

p_categories_from

ggsave(p_categories_from, path=path, file="p_categoriesfrom_nscore.jpg", width=5, height=4)

######

categories_m <- fdr_categories %>% group_by(Mut, category) %>% dplyr::summarise(Freq_m=n()) 

p_categories_m<-ggplot(categories_m, aes(fill=factor(category, levels=levels_cat), y=factor(Mut, levels=rev(vectorAA)), x=Freq_m))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.title.y = element_blank(),
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "Counts")

p_categories_m

ggsave(p_categories_m, path=path, file="p_categoriesm_nscore.jpg", width=15, height=8)
##################

p_hv <- ggarrange(p_heatmap + rremove("x.text"), p_violin_m + rremove("y.text"), p_categories_m + rremove("y.text"), 
                  ncol=3, align = "h",  legend="none", widths = c(1, 0.2, 0.2), common.legend=TRUE) 

p_hv

ggsave(p_hv, path=path, file="p_heatmap_violinm_categoreism.jpg", width=15, height=8)


### Arrange plots
p_heatmap <- p_heatmap + theme(axis.title.x = element_blank())
p_hvc<-ggarrange(p_heatmap + rremove("x.text"), p_violin_p + rremove("x.text"), p_categories_p, 
                 nrow=3, align = "v",  legend="none", heights = c(1, 0.3, 0.3))

p_hvc

ggsave(p_hvc, path=path, file="p_heatmap_violin_categories_nscore.jpg", width=11, height=12)

###

# Effect of SNVs on nucleation:
SNV_df<-singles_stops

SNV_df$category<-"WT-like"
SNV_df[SNV_df$low_sigma == F,]$category<- "unknown"
SNV_df[(SNV_df$p.adjust<0.1 & SNV_df$nscore_c<0),]$category<- "NS-"
SNV_df[(SNV_df$p.adjust<0.1 & SNV_df$nscore_c>0),]$category<- "NS+"

SNV_df$SNVs<-""
SNV_df[(SNV_df$aa_seq %in% SNVs$aa_seq) & (SNV_df$ID != "syn"), "SNVs"]<-"|"

SNV_df$Comb_mut<-""
SNV_df[SNV_df$ID=="E-1-K",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="S-3-N",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="C-5-Y",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="F-6-L",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="A-7-S",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="A-7-T",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="I-8-V",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="R-9-Q",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="K-14-E",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="A-16-V",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="V-17-M",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="L-20-S",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="L-20-V",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="I-21-V",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="S-29-G",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="K-32-T",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="H-33-Y",]$Comb_mut<-"|"
SNV_df[SNV_df$ID=="Y-34-H",]$Comb_mut<-"|"

# reported SNPs
SNV_df %>% group_by(Comb_mut, category) %>% dplyr::summarise(Freq_p=n())

# all SNVs
SNV_df %>% group_by(SNVs, category) %>% dplyr::summarise(Freq_p=n())

##################################################################################
####  heatmap FDR

singles_stops<-singles_stops[singles_stops$low_sigma == T | singles_stops$sig_10 == T, ]

heatmap_fdr<-singles_stops[,c("aa_seq", "Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10", "low_sigma")]

#add syn
syn.df<-data.frame(
  "aa_seq"="EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY",
  "WT_AA"= peptide_seq,
  "Mut"= peptide_seq,
  "Pos"= c(1:34), 
  "nscore_c"=0,
  "ID"="syn",
  "p.adjust"=NA,
  "category_10"="WT-like",
  "low_sigma"=FALSE
)

heatmap_fdr<-rbind(heatmap_fdr, syn.df)

heatmap_fdr$category<-"WT-like"

heatmap_fdr[(heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c<0),]$category<- "NS- 25%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c<0),]$category<- "NS- 10%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c<0),]$category<- "NS- 5%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c<0),]$category<- "NS- 1%"

heatmap_fdr[(heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 25%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 10%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 5%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 1%"

heatmap_fdr$WT<-""
heatmap_fdr[heatmap_fdr$ID=="syn",]$WT<-"ADan"

heatmap_fdr$SNVs<-""
heatmap_fdr[(heatmap_fdr$aa_seq %in% SNVs$aa_seq) & (heatmap_fdr$ID != "syn"), "SNVs"]<-"|"

heatmap_fdr$Comb_mut<-""
heatmap_fdr[heatmap_fdr$ID=="E-1-K",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="S-2-N",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="C-5-Y",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="F-6-L",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="A-7-S",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="A-7-T",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="I-8-V",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="R-9-Q",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="K-14-E",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="A-16-V",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="V-17-M",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="L-20-S",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="L-20-V",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="I-21-V",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="S-29-G",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="K-32-T",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="H-33-Y",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="Y-34-H",]$Comb_mut<-"|"


levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")

colors<-c( "NS- 1%"="#CD3333", 
           "NS- 5%"="#D66262", 
           "NS- 10%"="#DF9292", 
           "NS- 25%"="#E8C2C2",
           "WT-like"="#F2F2F2", 
           "NS+ 25%"="#B5B5D8",
           "NS+ 10%"="#7979BE",
           "NS+ 5%"="#3C3CA4",
           "NS+ 1%"="#00008B")

p_heatmap_fdr<-ggplot(heatmap_fdr)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)),fill=factor(category, levels=levels)), size=1, color='white')+
  geom_text(aes(Pos, Mut, label=WT), size=3)+
  geom_text(aes(Pos, Mut, label=SNVs), color="black", size=3, nudge_x=0.32, nudge_y=0.37, angle=45, fontface="bold")+
  geom_text(aes(Pos, Mut, label=Comb_mut), color="black", size=3, nudge_x=0.26, nudge_y=0.32, angle=45, fontface="bold")+
  scale_fill_manual("Category (FDR)", values= colors, labels=levels)+
  scale_x_continuous(breaks=seq(1,34), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  labs(x="", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=16))
p_heatmap_fdr

ggsave(p_heatmap_fdr,path=path, file="p_heatmap_FDR.jpg", width=15, height=8)

#####
fdr_categories<-singles_stops[,c("Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10")]

fdr_categories$category<-"WT-like"
fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c<0),]$category<- "NS- 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS- 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c<0),]$category<- "NS- 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c<0),]$category<- "NS- 1%"

fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c>0),]$category<- "NS+ 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+ 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c>0),]$category<- "NS+ 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c>0),]$category<- "NS+ 1%"


categories_fdr <- fdr_categories %>% group_by(Pos, category) %>% dplyr::summarise(Freq=n()) 

levels_cat_fdr = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")
colors_cat_fdr<-c( "NS- 1%"="#CD3333", 
                   "NS- 5%"="#D66262", 
                   "NS- 10%"="#DF9292", 
                   "NS- 25%"="#E8C2C2",
                   "WT-like"="#F2F2F2", 
                   "NS+ 25%"="#B5B5D8",
                   "NS+ 10%"="#7979BE",
                   "NS+ 5%"="#3C3CA4",
                   "NS+ 1%"="#00008B")

p_categories_fdr<-ggplot(categories_fdr, aes(fill=factor(category, levels=levels_cat_fdr), x=factor(Pos), y=Freq))+
  scale_x_discrete(breaks=seq(1,34), labels = peptide_seq_pos, expand=c(0,0))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color="black", size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation FDR(%)", values=colors_cat_fdr)+
  labs(x= "ADan peptide", y="Counts")

p_categories_fdr

ggsave(p_categories_fdr,path=path, file="p_categories_FDR.jpg", width=15, height=8)

####
p_hfdrc<-ggarrange(p_heatmap_fdr + rremove("x.text"), p_categories_fdr, 
                   nrow=2, align = "v",  legend="right", heights = c(1, 0.4),
                   common.legend = TRUE)

p_hfdrc

ggsave(p_hfdrc, path=path, file="p_heatmap_fdr_categories.jpg", width=14, height=10)

#####


#### stacked bars - all variants
# Considering stops
fdr_categories_all<-singles_stops[,c("Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10", "low_sigma")]

fdr_categories_all$category<-"WT-like"
fdr_categories_all[fdr_categories_all$low_sigma == F,]$category<- "unknown"
fdr_categories_all[(fdr_categories_all$p.adjust<0.1 & fdr_categories_all$nscore_c<0),]$category<- "NS-"
fdr_categories_all[(fdr_categories_all$p.adjust<0.1 & fdr_categories_all$nscore_c>0),]$category<- "NS+"

categories_p_all <- fdr_categories_all %>% group_by(category) %>% dplyr::summarise(Freq_p=n()) 
categories_p_all$library<-"ADan"
categories_p_all$percentage<-round((categories_p_all$Freq_p/sum(categories_p_all$Freq_p))*100, 2)

levels_cat = c("NS-", "WT-like", "NS+", "unknown")
colors_cat<-c("#DF9292", "#F2F2F2", "#7979BE", "white")

p_categories_p_all<-ggplot(categories_p_all, aes(fill=factor(category, levels=levels_cat), x=library, y=percentage))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8, color="darkgrey")+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "", y="% of variants")+
  annotate("text", x="ADan", y=12.5, label="NS+", size=5)+
  annotate("text", x="ADan", y=40, label="WT-like", size=5)+
  annotate("text", x="ADan", y=78, label="NS-", size=5)+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        plot.margin = unit(c(0,0,0,0), 'cm'))

p_categories_p_all

ggsave(p_categories_p_all, path=path, file="p_categories_p_all.jpg", width=4, height=4)

######
fdr_categories_all$region<-"flanking"
fdr_categories_all[fdr_categories_all$Pos>19 & fdr_categories_all$Pos<26, "region"]<-"core"

categories_p_grouped<-fdr_categories_all %>% group_by(region, category) %>% dplyr::summarise(Freq_p=n()) 
categories_p_grouped$freq_norm<-0
categories_p_grouped[categories_p_grouped$region == "core", "freq_norm"]<-categories_p_grouped[categories_p_grouped$region == "core", "Freq_p"]/6
categories_p_grouped[categories_p_grouped$region == "flanking", "freq_norm"]<-categories_p_grouped[categories_p_grouped$region == "flanking", "Freq_p"]/28
categories_p_grouped
