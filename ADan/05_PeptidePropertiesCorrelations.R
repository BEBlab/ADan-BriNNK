library(Peptides)
library(tidyverse)
require(ggplot2)
require(ggpubr)
require(ggrepel)
library(gridExtra)


dir.create("05_ScoresCorrelations")
path="05_ScoresCorrelations"

##  import required data 

name<-"ADan"

load(paste0("nscore_df_", name, ".RData"))

# Let's work only with those trustable sequences
singles_stops<-singles_stops[singles_stops$low_sigma == TRUE | singles_stops$sig_10 == TRUE,]
# Drop the nonsense
singles<-singles_stops[singles_stops$STOP == F,]
singles<-singles[,c("ID", "aa_seq", "nscore_c")]

wt<-data.frame("ID"="ADan", "aa_seq"="EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY", "nscore_c"=0.03251531)

singles<-rbind(singles, wt)

###
# Merge with scores calculated with localCider and some chunk of IDRLab/calvadoss (python code - jupyter notebook)
# and ALBATROSS (colab)

localcider_calvadoss_scores<-read_csv("ADan_df_localCider_calvadoss.csv", show_col_types = FALSE)
singles<-left_join(singles, localcider_calvadoss_scores, by="ID")

albatross_scores<-read_csv("ADan_df_ALBATROSS.csv", show_col_types = FALSE)
albatross_scores<-rename(albatross_scores, ID = fasta_header, length = lengths) 
albatross_scores<-albatross_scores[,c("ID", "length", "radius_of_gyration", "end_to_end_distance", 
                                      "scaling_exponent", "asphericity", "prefactor")]
singles<-left_join(singles, albatross_scores, by="ID")

### 
# Add correlation with predictors: Tango, CamSol, AMYPRED, Polyphen
# Tango
tango_scores<-read_csv("ADan_df_tango_score.csv", show_col_types = FALSE)
tango_scores<-tango_scores[,c("ID", "Tango_score")]
singles<-left_join(singles, tango_scores, by="ID")
# Camsol
camsol_scores<-read_csv("ADan_df_camsol_score.csv", show_col_types = FALSE)
camsol_scores<-camsol_scores[,c("ID", "CamSol_score")]
singles<-left_join(singles, camsol_scores, by="ID")
# AMYPred-FRL
amypred_scores<-read_csv("ADan_df_amypred_score.csv", show_col_types = FALSE)
amypred_scores<-amypred_scores[,c("ID", "Amypred_score")]
singles<-left_join(singles, amypred_scores, by="ID")


###
# Scores calculation for the whole sequence - singles, insertions and deletions
# Measurments of Hydrophobicity (KyteDooltitle), Boman score, aIndex
# charge, aaDescriptors(Cruciani properties, Fasgai vectors, Kidera Factors,
# pI, ProtFP, stScales, vheScales, zScales..)

for (i in 1:nrow(singles)){
  singles[i, "Hydrophobicity(KD)"]<-hydrophobicity(singles[i, "aa_seq"], scale = "KyteDoolittle")
  singles[i, "Boman"]<-boman(singles[i, "aa_seq"])
  singles[i, "aIndex"]<-aIndex(singles[i, "aa_seq"])
  singles[i, "instaIndex"]<-instaIndex(singles[i, "aa_seq"])
  singles[i, "Charge"]<-charge(singles[i, "aa_seq"], pH = 7, pKscale = "Lehninger")
  singles[i, "pI"]<-pI(singles[i, "aa_seq"], pKscale = "EMBOSS")
  
  cruciani<-crucianiProperties(singles[i, "aa_seq"])
  singles[i, colnames(cruciani)]<-cruciani
  fasgai<-fasgaiVectors(singles[i, "aa_seq"])
  singles[i, colnames(fasgai)]<-fasgai
  kidera<-kideraFactors(singles[i, "aa_seq"])
  singles[i, colnames(kidera)]<-kidera
  msw<-mswhimScores(singles[i, "aa_seq"])
  singles[i, colnames(msw)]<-msw
  protfp<-protFP(singles[i, "aa_seq"])
  singles[i, colnames(protfp)]<-protfp
  tscales<-tScales(singles[i, "aa_seq"])
  singles[i, colnames(tscales)]<-tscales
  vhs<-vhseScales(singles[i, "aa_seq"]) 
  singles[i, colnames(vhs)]<-vhs
  zscales<-zScales(singles[i, "aa_seq"]) 
  singles[i, colnames(zscales)]<-zscales
  blosum<-blosumIndices(singles[i, "aa_seq"])
  singles[i, colnames(blosum)]<-blosum
}

singles<-as.data.frame(singles)
# Scores per aa change 

load("AA_descriptors.RData")

# Keep only those scores that are "approved"
AA_descriptors<-AA_descriptors[AA_descriptors$Status %in% c("Accepted", "Added"),]

# Iterate over the different type of scores
# Sum the value of each aa in a sequence divided by the length of the sequence
for (i in AA_descriptors$Identifier){
  score<-AA_descriptors[AA_descriptors$Identifier == `i`,]
  # Iterate over each aa_seq
  for (j in 1:nrow(singles)){
    peptide_seq<-unlist(strsplit(singles[j, "aa_seq"], ""))
    values<-c()
    for (aa in peptide_seq){
      values<-c(values, as.numeric(score[aa]))
    }
    singles[j, i]<-sum(values)/length(values)
  }
}

###########################################################################################################
###########################################################################################################
### All the variants
# Define a correlation vector that keeps all the correlation values
# Vector for whole peptide
corr_vector=c()


to_correlate<-c("Hydrophobicity(KD)", "Boman", "aIndex", "Charge", "instaIndex", "pI", colnames(vhs), colnames(blosum),
                colnames(cruciani), colnames(fasgai), colnames(kidera), colnames(msw), colnames(protFP), colnames(tscales),
                "radius_of_gyration", "end_to_end_distance", "scaling_exponent", "asphericity", "prefactor",
                "fK", "fR", "fE", "fD", "fARO", "Mean_lambda", "SHD", "Omega_ARO", "SCD", "kappa", "FCR", "NCPR",
                "Tango_score", "CamSol_score", "Amypred_score", AA_descriptors$Identifier)

for (i in to_correlate){
  # Whole sequence
  correlation<-cor.test(singles$nscore_c,
                        singles[[i]], use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

### Whole sequence  
corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Identifier", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}

# Keep the significant ones
corr_text_sig<-corr_text[(corr_text$pvalue<0.05),]
corr_text_sig<-arrange(corr_text_sig, corr)
corr_text_sig$R2<-corr_text_sig$corr^2
corr_text_sig<-arrange(corr_text_sig, desc(R2))

# Plot the top 20 correlations
corr_bar<-ggplot(head(corr_text_sig, 20), aes(x=Identifier, y=corr))+
  geom_bar(stat="identity")+
  ylim(-1,1)+theme_bw()+
  annotate("text", label="NS vs others - Whole Peptide", x=-Inf, y=Inf, vjust=2, hjust=-.1)+
  labs(x="", y="Pearson Correlation")+
  theme(axis.text=element_text(angle = 40, vjust = 0.5, hjust=1),
        title=element_text(size=14))
corr_bar

ggsave(corr_bar, file="ScoresCorr_bars.jpg", width = 6, height = 4, path=path)


###
# Plot Scatter plot of predictors + hydrophobicity

#to_plot<-c("Tango_score", "CamSol_score", "Amypred_score", "Hydrophobicity(KD)")

singles <- rename(singles, "Tango Score" = "Tango_score", "Amypred Score" = "Amypred_score",
                           "CamSol Score \n (Solubility)" = "CamSol_score")

corr_text[corr_text$Identifier == "Tango_score", "Identifier"]<-"Tango Score"
corr_text[corr_text$Identifier == "Amypred_score", "Identifier"]<-"Amypred Score"
corr_text[corr_text$Identifier == "CamSol_score", "Identifier"]<-"CamSol Score \n (Solubility)"

to_plot<-c("Tango Score", "CamSol Score \n (Solubility)", "Amypred Score", "Hydrophobicity(KD)")

scatter_plots<-list()
n<-1
for (i in to_plot){
  scatter_p<-ggplot(singles)+
    geom_point(aes(x=nscore_c, y=.data[[i]]), color="grey80", size=2)+
    geom_point(data=singles[singles$ID == "ADan",], aes(x=nscore_c, y=.data[[i]]), color="red", size=3, alpha=0.5)+
    labs(color="", shape="", x="Nucleation Score")+
    annotate("text", label=paste0("R=", round(corr_text[corr_text$Identifier == i,]$corr, 2)), x=-5.5, y=Inf, vjust=2, size=5)+
    annotate("text", label=paste0("p=", round(corr_text[corr_text$Identifier == i,]$pvalue, 2)), x=-5.5, y=Inf, vjust=3.5, size=5)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 16),
          axis.text = element_text(size=16),
          legend.text= element_text(size=14),
          plot.title = element_text(hjust = 0.5))
  scatter_plots[[n]]<-scatter_p
  n=sum(n, 1)
}

corr_scatter_merge<-ggarrange(plotlist=scatter_plots,
                              ncol=2, nrow=2, align = "hv", 
                              common.legend = TRUE, legend = "top")
corr_scatter_merge

ggsave(corr_scatter_merge, file="ScoresCorr_scatter_region.jpg", width = 10, height = 8, path=path)

#

save(singles, corr_text, corr_text_sig,  file = paste0(name, "_scores_df.RData"))

#

#load(paste0(name, "_scores_df.RData"))
