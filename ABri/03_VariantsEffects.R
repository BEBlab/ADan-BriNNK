library(tidyverse)
require(reshape2)
library(ggpubr)
library(seqinr)

dir.create("03_VariantsEffects")
path="03_VariantsEffects"

#required data:
load("nscore_df_ABri.RData")

#
peptide_seq<-"EASNCFAIRHFENKFAVETLICSRTVKKNIIEEN"
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

nt_seq<-"GAAGCCAGCAATTGTTTCGCAATTCGGCATTTTGAAAACAAATTTGCCGTGGAAACTTTAATTTGTTCTCGAACAGTCAAGAAAAACATTATTGAGGAAAAT"

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
singles_stops<-singles_stops[singles_stops$low_sigma == T ,]

#####  heatmap (contains stops)


#add syn
positions<-(1:34)
syn.df<-data.frame(
  "aa_seq" = "EASNCFAIRHFENKFAVETLICSRTVKKNIIEEN",
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
heatmap_df[heatmap_df$ID=="syn",]$WT<-"ABri"

heatmap_df$FD<-""
heatmap_df[heatmap_df$ID=="R-24-L",]$FD<-"FCD"
heatmap_df[heatmap_df$ID=="R-24-S",]$FD<-"FKD"


heatmap_df$SNVs<-""
heatmap_df[(heatmap_df$aa_seq %in% SNVs$aa_seq) & (heatmap_df$ID != "syn"), "SNVs"]<-"|"
heatmap_df[heatmap_df$Pos==24,]$SNVs<-""
heatmap_df[heatmap_df$ID=="R-24-L",]$SNVs<-"|"
heatmap_df[heatmap_df$ID=="R-24-S",]$SNVs<-"|"
heatmap_df[heatmap_df$ID=="R-24-G",]$SNVs<-"|"
heatmap_df[heatmap_df$ID=="R-24-*",]$SNVs<-"|"
heatmap_df[heatmap_df$ID=="R-24-C",]$SNVs<-"|"
heatmap_df[heatmap_df$ID=="R-24-W",]$SNVs<-"|"

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
heatmap_df[heatmap_df$ID=="K-27-R",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="K-28-R",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="N-29-H",]$Comb_mut<-"|"
heatmap_df[heatmap_df$ID=="I-31-T",]$Comb_mut<-"|"

p_heatmap<-ggplot(heatmap_df)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)), fill=nscore_c), color='white', size=1)+
  geom_text(aes(Pos, Mut, label=FD, color=FD), size=3, fontface="bold")+
  geom_text(aes(Pos, Mut, label=WT), size=3)+
  geom_text(aes(Pos, Mut, label=SNVs), color="black", size=3, nudge_x=0.32, nudge_y=0.37, angle=45, fontface="bold")+
  geom_text(aes(Pos, Mut, label=Comb_mut), color="black", size=3, nudge_x=0.26, nudge_y=0.32, angle=45, fontface="bold")+
  scale_color_manual(values=c("white", "darkgreen", "brown", "black"))+
  scale_x_continuous(breaks=seq(1,34), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  labs(x="ABri peptide", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=16))+ 
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-2, 0, 2, 4, 6, 8))
p_heatmap

ggsave(p_heatmap, path=path, file="p_heatmap_nscore.jpg", width=15, height=8)


##################################################################################
####  heatmap FDR

heatmap_fdr<-singles_stops[,c("aa_seq","Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10", "low_sigma")]

#add syn
syn.df<-data.frame(
  "aa_seq" = "EASNCFAIRHFENKFAVETLICSRTVKKNIIEEN",
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
heatmap_fdr[heatmap_fdr$ID=="syn",]$WT<-"ABri"

heatmap_fdr$SNVs<-""
heatmap_fdr[(heatmap_fdr$aa_seq %in% SNVs$aa_seq) & (heatmap_fdr$ID != "syn"), "SNVs"]<-"|"
heatmap_fdr[heatmap_fdr$ID=="R-24-L",]$SNVs<-"|"
heatmap_fdr[heatmap_fdr$ID=="R-24-S",]$SNVs<-"|"
heatmap_fdr[heatmap_fdr$ID=="R-24-G",]$SNVs<-"|"
heatmap_fdr[heatmap_fdr$ID=="R-24-*",]$SNVs<-"|"
heatmap_fdr[heatmap_fdr$ID=="R-24-C",]$SNVs<-"|"
heatmap_fdr[heatmap_fdr$ID=="R-24-W",]$SNVs<-"|"

heatmap_fdr$FD<-""
heatmap_fdr[heatmap_fdr$ID=="R-24-L",]$FD<-"FCD"
heatmap_fdr[heatmap_fdr$ID=="R-24-S",]$FD<-"FKD"

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
heatmap_fdr[heatmap_fdr$ID=="K-27-R",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="K-28-R",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="N-29-H",]$Comb_mut<-"|"
heatmap_fdr[heatmap_fdr$ID=="I-31-T",]$Comb_mut<-"|"


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
  geom_text(aes(Pos, Mut, label=FD, color=FD), size=3, fontface="bold")+
  geom_text(aes(Pos, Mut, label=WT), size=3)+
  geom_text(aes(Pos, Mut, label=SNVs), color="black", size=3, nudge_x=0.32, nudge_y=0.37, angle=45, fontface="bold")+
  geom_text(aes(Pos, Mut, label=Comb_mut), color="black", size=3, nudge_x=0.26, nudge_y=0.32, angle=45, fontface="bold")+
  scale_fill_manual("Category (FDR)", values= colors, labels=levels)+
  scale_color_manual(values=c("white", "darkgreen", "brown", "black"))+
  scale_x_continuous(breaks=seq(1,34), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  labs(x="ABri peptide", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=16))
p_heatmap_fdr

ggsave(p_heatmap_fdr,path=path, file="p_heatmap_FDR.jpg",width=15, height=8)

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
  labs(x= "ABri peptide", y="Counts")

p_categories_fdr

ggsave(p_categories_fdr,path=path, file="p_categories_FDR.jpg",width=15, height=8)

####
p_hfdrc<-ggarrange(p_heatmap_fdr + rremove("x.text"), p_categories_fdr, 
                   nrow=2, align = "v",  legend="right", heights = c(1, 0.4),
                   common.legend = TRUE)

p_hfdrc

ggsave(p_hfdrc, path=path, file="p_heatmap_fdr_categories.jpg", width=14, height=10)

