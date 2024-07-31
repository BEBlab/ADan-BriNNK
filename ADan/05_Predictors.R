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
# Measurments of Hydrophobicity (KyteDooltitle)


for (i in 1:nrow(singles)){
  singles[i, "Hydrophobicity(KD)"]<-hydrophobicity(singles[i, "aa_seq"], scale = "KyteDoolittle")
 }

singles<-as.data.frame(singles)
# Scores per aa change 

load("AA_descriptors.RData")


###########################################################################################################
### All the variants
# Define a correlation vector that keeps all the correlation values
# Vector for whole peptide
corr_vector=c()

singles <- rename(singles, "Tango Score" = "Tango_score", "Amypred Score" = "Amypred_score",
                  "CamSol Score \n (Solubility)" = "CamSol_score")

to_correlate<-c("Hydrophobicity(KD)", "Tango Score", "Amypred Score",
                "CamSol Score \n (Solubility)")

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

###
# Plot Scatter plot of predictors + hydrophobicity


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

