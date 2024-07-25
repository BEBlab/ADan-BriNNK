library(pROC)
library(ggpubr)
library(tidyverse)

dir.create("04_ROCs")
path="04_ROCs"

load("NNK_all_df.RData")

hydrophobicity<-read_tsv("Bri2NNK_all_df_hydrophobicity.tsv")
hydrophobicity<-rename(hydrophobicity, aa_seq = sequence, Hydrophobicity = `Hydrophobicity (Kite-Doolittle)`)
hydrophobicity<-hydrophobicity[,c("aa_seq", "Hydrophobicity")]

# Tango run with default parameters (T = 298, pH=7.2, ionic = 0.1, not protected terminus)
tangoscore<-read_tsv("Bri2NNK_all_df_tango.txt")
tangoscore<-rename(tangoscore, aa_seq = Sequence, Tango_score = Aggregation)
tangoscore<-tangoscore[,c("aa_seq", "Tango_score")]

# CamSol run with default parameters
camsolscore<-read_tsv("Bri2NNK_all_df_camsol.txt")
camsolscore<-rename(camsolscore, aa_seq = Name, CamSol_score = `protein variant score`)
camsolscore<-camsolscore[,c("aa_seq", "CamSol_score")]

# Amypred run with default parameters
amypred<-read_tsv("Bri2NNK_all_df_amypred.txt")
amypred<-rename(amypred, aa_seq = Name, Amypred_score = `Probability`)
amypred<-amypred[,c("aa_seq", "Amypred_score")]

# Aggrescan 
aggrescan<-read_tsv("Bri2NNK_all_df_1_aggrescan.txt")
names<-aggrescan$`Sequence Name`
aggrescan <- as.data.frame(t(aggrescan[,-1]))
colnames(aggrescan) <- names
aggrescan$aa_seq <- factor(row.names(aggrescan))

for (i in (2:3)){
  aggrescan_x<-read_tsv(paste0("Bri2NNK_all_df_", i,"_aggrescan.txt"))
  names<-aggrescan_x$`Sequence Name`
  aggrescan_x<-as.data.frame(t(aggrescan_x[,-1]))
  colnames(aggrescan_x)<-names
  aggrescan_x$aa_seq<-factor(row.names(aggrescan_x))
  aggrescan<-rbind(aggrescan, aggrescan_x)
}

aggrescan<-rename(aggrescan, Aggrescan_score = `Area of the profile Above Threshold (AAT):`)
aggrescan<-aggrescan[,c("aa_seq", "Aggrescan_score")]


#

NNK_all_df<-rename(NNK_all_df, aa_seq = full_seq)

df_list<-list(NNK_all_df, hydrophobicity, 
              tangoscore, camsolscore, amypred, 
              aggrescan)

NNK_all_df<-df_list %>% reduce(full_join, by="aa_seq")
NNK_all_df<-NNK_all_df[!is.na(NNK_all_df$mode_seed_bh),]

#

roc.list <- roc(NNK_all_df$mode_seed_bh ~ Tango_score + CamSol_score + 
                  Amypred_score + Aggrescan_score + Hydrophobicity, 
                data = NNK_all_df,
                ci=T, boot.n=1000)


ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

aucs<-data.frame("label"=c("Tango_score", "CamSol_score", "Amypred_score",
                           "Aggrescan_score","Hydrophobicity"),
                 "AUClabel"=c(paste0("Tango ", round(roc.list$Tango_score$auc[1], 4), " ± ", round(roc.list$Tango_score$auc-roc.list$Tango_score$ci[1], 4)), 
                              paste0("CamSol ", round(roc.list$CamSol_score$auc[1], 4), " ± ", round(roc.list$CamSol_score$auc-roc.list$CamSol_score$ci[1], 4)),
                              paste0("Amypred ", round(roc.list$Amypred_score$auc[1], 4), " ± ", round(roc.list$Amypred_score$auc-roc.list$Amypred_score$ci[1], 4)),
                              paste0("Aggrescan ", round(roc.list$Aggrescan_score$auc[1], 4), " ± ", round(roc.list$Aggrescan_score$auc-roc.list$Aggrescan_score$ci[1], 4)),
                              paste0("Hydrophobicity ", round(roc.list$Hydrophobicity$auc[1], 4), " ± ", round(roc.list$Hydrophobicity$auc-roc.list$Hydrophobicity$ci[1], 4))))

p <- ggroc(roc.list, size=.2) + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(labels=aucs$AUClabel, values=c("#F0E442", "#CC79A7", "#56B4E9", "#E69F00", "#999999"))+
  labs(x="False Positive Rate", y="True Positive Rate", color="")+
  coord_cartesian(clip = "off")+
  theme_classic() +
  theme(
    plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
    legend.position = c(0.7, 0.25),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fill_ribbon <- c("#F0E442", "#CC79A7", "#56B4E9", "#E69F00", "#999999")

for(i in 1:length(roc.list)) {
  p <- p + geom_ribbon(
    data = dat.ci.list[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = fill_ribbon[i],
    alpha = 0.3,
    inherit.aes = F) 
} 


p

ggsave(p, file="roc_plot_predictors.jpg", width = 4, height = 4, path=path)


# Extremes

NNK_all_df_extremes<-NNK_all_df[NNK_all_df$category_fdr == "non-nucleators" | NNK_all_df$category_fdr == "Top10 nucleators",]

roc.list <- roc(NNK_all_df_extremes$mode_seed_bh ~ Tango_score + CamSol_score + 
                  Amypred_score + Aggrescan_score + Hydrophobicity, 
                data = NNK_all_df_extremes,
                ci=T, boot.n=1000)


ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

aucs<-data.frame("label"=c("Tango_score", "CamSol_score", "Amypred_score",
                           "Aggrescan_score","Hydrophobicity"),
                 "AUClabel"=c(paste0("Tango ", round(roc.list$Tango_score$auc[1], 4), " ± ", round(roc.list$Tango_score$auc-roc.list$Tango_score$ci[1], 4)), 
                              paste0("CamSol ", round(roc.list$CamSol_score$auc[1], 4), " ± ", round(roc.list$CamSol_score$auc-roc.list$CamSol_score$ci[1], 4)),
                              paste0("Amypred ", round(roc.list$Amypred_score$auc[1], 4), " ± ", round(roc.list$Amypred_score$auc-roc.list$Amypred_score$ci[1], 4)),
                              paste0("Aggrescan ", round(roc.list$Aggrescan_score$auc[1], 4), " ± ", round(roc.list$Aggrescan_score$auc-roc.list$Aggrescan_score$ci[1], 4)),
                              paste0("Hydrophobicity ", round(roc.list$Hydrophobicity$auc[1], 4), " ± ", round(roc.list$Hydrophobicity$auc-roc.list$Hydrophobicity$ci[1], 4))))

p <- ggroc(roc.list, size=.2) + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(labels=aucs$AUClabel, values=c("#F0E442", "#CC79A7", "#56B4E9", "#E69F00", "#999999"))+
  labs(x="False Positive Rate", y="True Positive Rate", color="")+
  coord_cartesian(clip = "off")+
  theme_classic() +
  theme(
    plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
    legend.position = c(0.7, 0.25),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fill_ribbon <- c("#F0E442", "#CC79A7", "#56B4E9", "#E69F00", "#999999")

for(i in 1:length(roc.list)) {
  p <- p + geom_ribbon(
    data = dat.ci.list[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = fill_ribbon[i],
    alpha = 0.3,
    inherit.aes = F) 
} 


p

ggsave(p, file="roc_plot_predictors_extremes.jpg", width = 4, height = 4, path=path)
