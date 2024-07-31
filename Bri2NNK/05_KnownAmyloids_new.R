library(tidyverse)
library(ggpubr)
require(stringr)
require("ggrepel")
library(ggExtra)


dir.create("05_KnownAmyloids")
path="05_KnownAmyloids"


###

load("NNK_all_df.RData")

vectorAA <- c("A","V","L","M","I","F","Y","W","K","R","D","E","S","T","N","Q","H","G", "P","C")

###

NNK_all_df_index<-read_tsv("Bri2NNK_all_df_aaindex.tsv")
NNK_all_df_index<-rename(NNK_all_df_index, full_seq = index)
NNK_all_df_index$name<-""

human_proteome<-read_csv("human_proteome_aaindex.csv")
human_proteome$group<-"Human Proteome"

#human_amyloids<-read_csv("human_amyloids_aaindex.csv")
#human_amyloids$group<-"Human Amyloids"


NNK_all_df<-left_join(NNK_all_df, NNK_all_df_index, by="full_seq")
NNK_all_df$group<-NNK_all_df$category_fdr


##
# Hydropathy index (Kyte-Doolittle, 1982): KYTJ820101
# Mean polarity (Radzicka-Wolfenden, 1988): RADA880108
# Polarity (Grantham, 1974): GRAR740102
# Average relative probability of beta-sheet (Kanehisa-Tsong, 1980): KANM800102
# Beta-sheet propensity derived from designed sequences (Koehl-Levitt, 1999): KOEP990102
# Average flexibility indices (Bhaskaran-Ponnuswamy, 1988): BHAR880101
# Net charge (Klein et al., 1984): KLEP840101


merge_dfs<-function(scores_interest){
  columns_interest<-c("aa_seq", "group", "name",
                      scores_interest)
  
  NNK_all_df<-NNK_all_df[, columns_interest]
  
  human_proteome<-human_proteome[, columns_interest]
  
  #human_amyloids<-human_amyloids[, columns_interest]
  
  all_sequences<-rbind(NNK_all_df, human_proteome)
  
  all_sequences$label<-NA
  all_sequences[all_sequences$name == "AB42", "label"]<-"AB42"
  all_sequences[all_sequences$name == "RIPK1", "label"]<-"RIPK1"
  all_sequences[all_sequences$name == "Tau", "label"]<-"Tau"
  all_sequences[all_sequences$name == "HNRNPA2", "label"]<-"HNRNPA2"
  
  return(all_sequences)
}  



known_amyloids_nnk_plot<-function(scores_interest, x, y, x_name, y_name){
    all_sequences<-merge_dfs(scores_interest)
    
    knownamyloids_scatter<-ggplot(all_sequences, 
                                  aes_string(x=x, y=y, color="group"))+
    geom_point(aes(color=factor(group, levels=c("Human Proteome", "nucleators", "non-nucleators", "Top10 nucleators"), labels=c("Human Proteome", "nucleators", "non-nucleators", "Top10 nucleators"))),
                 alpha=0.2, size=1)+
    geom_density_2d(aes(color=factor(group, levels=c("Human Proteome", "nucleators", "non-nucleators", "Top10 nucleators"), labels=c("Human Proteome", "nucleators", "non-nucleators", "Top10 nucleators"))), 
                    bins = 100, linewidth=1, alpha=0.5)+
    scale_color_manual(values=c("grey80", "#DF9292", "#7979BE", "darkblue"))+
    labs(x=x_name, y=y_name, color="")+
    coord_cartesian(clip = "off")+
    theme_classic() +
    theme(plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
          legend.position = c(0.25, 0.8),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank(),
          legend.spacing.x = unit(0, "cm")
    )
  
  print(knownamyloids_scatter)
  ggsave(knownamyloids_scatter, file=paste0("knownamyloids_scatter_", x, "_", y, "_NNK.jpg"), width = 4, height = 4, path=path)

  knownamyloids_scatter_marginal<-ggMarginal(knownamyloids_scatter, type="density",  groupColour = TRUE, groupFill = TRUE)
  print(knownamyloids_scatter_marginal)
  
  ggsave(knownamyloids_scatter_marginal, file=paste0("knownamyloids_scatter_", x, "_", y, "_marginal_NNK.jpg"), width = 4.2, height = 4.2, path=path)

}

# Hydrophobicity vs Beta-sheet
known_amyloids_nnk_plot(c("KYTJ820101", "KANM800102"), "KYTJ820101", "KANM800102", "Hydrophobicity \n(Kyte-Doolitle)", "Beta-sheet proability \n(Kanehisa-Tsong)")


