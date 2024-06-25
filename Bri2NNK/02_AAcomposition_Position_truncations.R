library(tidyverse)
require(reshape2)
library(ggpubr)
require(stringr)
library(ggsignif)
require("ggrepel")


setwd("C:/Users/mmartin/OneDrive - IBEC/Documentos/Postdoc-BEBlab/Lab/LabProjects/BRI2_ADan_ABri/Sequencing_Analysis/Bri2_NNK_all_240417")


#required data:
dir.create("02_AAcomposition_Position_Truncations")


##
# Length of extension:
#seq_len<-11
for (seq_len in 2:11){
  name<-"NNK_all_df.RData"
  load(name)
  
  NNK_all_df<-NNK_all_df[NNK_all_df$aa_len == seq_len,]
  
  #
  dir.create(paste0("02_AAcomposition_Position_Truncations/02_AAcomposition_Position_", seq_len))
  path=paste0("02_AAcomposition_Position_Truncations/02_AAcomposition_Position_", seq_len)
  
  
  # Define important vectors and data frames with AA info
  vectorAA <- c("A","V","L","M","I","F","Y","W","K","R","D","E","S","T","N","Q","H","G", "P","C")
  
  AA_type<-data.frame("AA"= vectorAA,
                      "name_AA"=c("Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                                  "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                                  "Asparagine","Glutamine","Histidine","Glycine", "Proline", "Cysteine"),
                      "type"=c(rep("aliphatic",5), rep("aromatic",3), rep("positive",2), rep("negative",2), rep("polar",5), "glycine", "proline", "cysteine"))
  
  color_AAtype=c("aliphatic"="darkgrey",
                 "aromatic"="#9A703EFF",
                 "negative"="#EE0011FF",
                 "positive"="#0C5BB0FF", 
                 "polar"="#15983DFF", 
                 "glycine"="grey30", 
                 "cysteine"="#FEC10BFF",
                 "proline"="black")
  
  
  
  ################################################################################################################
  # All sequences in the input (including those that doesn't have output counts)
  
  NNK_all_df$pos_1<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 1)
  NNK_all_df$pos_2<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 2)
  NNK_all_df$pos_3<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 3)
  NNK_all_df$pos_4<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 4)
  NNK_all_df$pos_5<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 5)
  NNK_all_df$pos_6<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 6)
  NNK_all_df$pos_7<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 7)
  NNK_all_df$pos_8<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 8)
  NNK_all_df$pos_9<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 9)
  NNK_all_df$pos_10<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 10)
  NNK_all_df$pos_11<-sapply(strsplit(NNK_all_df$aa_seq_cor, ""), '[', 11)
  
  
  nucl_df<-NNK_all_df
  
  pos_1<-table(nucl_df$pos_1)/length(nucl_df$pos_1)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_1))){
      pos_1[aa]<-0}}
  pos_1<-pos_1[order(names(pos_1))]
  
  pos_2<-table(nucl_df$pos_2)/length(nucl_df$pos_2)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_2))){
      pos_2[aa]<-0}}
  pos_2<-pos_2[order(names(pos_2))]
  
  pos_3<-table(nucl_df$pos_3)/length(nucl_df$pos_3)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_3))){
      pos_3[aa]<-0}}
  pos_3<-pos_3[order(names(pos_3))]
  
  pos_4<-table(nucl_df$pos_4)/length(nucl_df$pos_4)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_4))){
      pos_4[aa]<-0}}
  pos_4<-pos_4[order(names(pos_4))]
  
  pos_5<-table(nucl_df$pos_5)/length(nucl_df$pos_5)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_5))){
      pos_5[aa]<-0}}
  pos_5<-pos_5[order(names(pos_5))]
  
  pos_6<-table(nucl_df$pos_6)/length(nucl_df$pos_6)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_6))){
      pos_6[aa]<-0}}
  pos_6<-pos_6[order(names(pos_6))]
  
  pos_7<-table(nucl_df$pos_7)/length(nucl_df$pos_7)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_7))){
      pos_7[aa]<-0}}
  pos_7<-pos_7[order(names(pos_7))]
  
  pos_8<-table(nucl_df$pos_8)/length(nucl_df$pos_8)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_8))){
      pos_8[aa]<-0}}
  pos_8<-pos_8[order(names(pos_8))]
  
  pos_9<-table(nucl_df$pos_9)/length(nucl_df$pos_9)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_9))){
      pos_9[aa]<-0}}
  pos_9<-pos_9[order(names(pos_9))]
  
  pos_10<-table(nucl_df$pos_10)/length(nucl_df$pos_10)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_10))){
      pos_10[aa]<-0}}
  pos_10<-pos_10[order(names(pos_10))]
  
  pos_11<-table(nucl_df$pos_11)/length(nucl_df$pos_11)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_11))){
      pos_11[aa]<-0}}
  pos_11<-pos_11[order(names(pos_11))]
  
  
  aa_frequencies_pos_all<-rbind(pos_1, pos_2, pos_3, pos_4, pos_5,
                                pos_6, pos_7, pos_8, pos_9, pos_10, pos_11)
  
  aa_frequencies_pos_all<-as.data.frame(aa_frequencies_pos_all)
  
  aa_frequencies_pos_all<-t(aa_frequencies_pos_all)
  
  co_all=melt(aa_frequencies_pos_all)
  co_all<-co_all[co_all$Var1 != "*",]
  
  # Get rid of positions with null values
  for (i in 1:11){
    if (length(unique(co_all[co_all$Var2 == paste0("pos_", i), "value"])) == 1){
      co_all<-co_all[co_all$Var2 != paste0("pos_", i),]
    }
  }
  
  
  freq_heatmap_all<-ggplot(co_all, aes(factor(Var2, labels=(23:(22+seq_len))), factor(Var1, levels=rev(vectorAA))))+
    geom_tile(aes(fill=value), color="white", size=1)+
    scale_fill_continuous(high="black", low="grey95", limit=c(0, 0.25))+
    labs(x="Position", y="Amino acid", fill="Frequency", title="AA frequencies per position in all sequences")+
    coord_cartesian(clip = "off")+
    theme_minimal()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"))
  
  print(freq_heatmap_all)
  
  ggsave(freq_heatmap_all, file="AA_frequency_all.jpg", width = 14, height = 8, path=path)
  
  #
  p_scatter_all<-ggplot(co_all, aes(x=factor(Var2, labels=23:(22+seq_len)), y=value, group=1))+
    geom_point()+
    geom_line()+
    ylim(c(-0.15, 0.15))+
    facet_wrap(~factor(Var1, levels=vectorAA))+
    labs(x="Position", y="Frequency")+
    theme_bw()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
      legend.position = c(0.85, 0.2),
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      axis.line=element_line()
    )
  p_scatter_all
  
  ggsave(p_scatter_all, file="AA_frequency_all_pos.jpg", width = 14, height = 8, path=path)
  
  ####
  #aa_frequencies_pos_all<-aa_frequencies_pos_all[,1:seq_len]
  
  ################################################################################################################
  # Work with nucleators
  
  nucleators_nucleators<-NNK_all_df[NNK_all_df$category_fdr == "nucleators",]
  
  nucleators_nucleators$pos_1<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 1)
  nucleators_nucleators$pos_2<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 2)
  nucleators_nucleators$pos_3<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 3)
  nucleators_nucleators$pos_4<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 4)
  nucleators_nucleators$pos_5<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 5)
  nucleators_nucleators$pos_6<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 6)
  nucleators_nucleators$pos_7<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 7)
  nucleators_nucleators$pos_8<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 8)
  nucleators_nucleators$pos_9<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 9)
  nucleators_nucleators$pos_10<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 10)
  nucleators_nucleators$pos_11<-sapply(strsplit(nucleators_nucleators$aa_seq_cor, ""), '[', 11)
  
  nucl_df<-nucleators_nucleators
  
  pos_1<-table(nucl_df$pos_1)/length(nucl_df$pos_1)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_1))){
      pos_1[aa]<-0}}
  pos_1<-pos_1[order(names(pos_1))]
  
  pos_2<-table(nucl_df$pos_2)/length(nucl_df$pos_2)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_2))){
      pos_2[aa]<-0}}
  pos_2<-pos_2[order(names(pos_2))]
  
  pos_3<-table(nucl_df$pos_3)/length(nucl_df$pos_3)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_3))){
      pos_3[aa]<-0}}
  pos_3<-pos_3[order(names(pos_3))]
  
  pos_4<-table(nucl_df$pos_4)/length(nucl_df$pos_4)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_4))){
      pos_4[aa]<-0}}
  pos_4<-pos_4[order(names(pos_4))]
  
  pos_5<-table(nucl_df$pos_5)/length(nucl_df$pos_5)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_5))){
      pos_5[aa]<-0}}
  pos_5<-pos_5[order(names(pos_5))]
  
  pos_6<-table(nucl_df$pos_6)/length(nucl_df$pos_6)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_6))){
      pos_6[aa]<-0}}
  pos_6<-pos_6[order(names(pos_6))]
  
  pos_7<-table(nucl_df$pos_7)/length(nucl_df$pos_7)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_7))){
      pos_7[aa]<-0}}
  pos_7<-pos_7[order(names(pos_7))]
  
  pos_8<-table(nucl_df$pos_8)/length(nucl_df$pos_8)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_8))){
      pos_8[aa]<-0}}
  pos_8<-pos_8[order(names(pos_8))]
  
  pos_9<-table(nucl_df$pos_9)/length(nucl_df$pos_9)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_9))){
      pos_9[aa]<-0}}
  pos_9<-pos_9[order(names(pos_9))]
  
  pos_10<-table(nucl_df$pos_10)/length(nucl_df$pos_10)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_10))){
      pos_10[aa]<-0}}
  pos_10<-pos_10[order(names(pos_10))]
  
  pos_11<-table(nucl_df$pos_11)/length(nucl_df$pos_11)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_11))){
      pos_11[aa]<-0}}
  pos_11<-pos_11[order(names(pos_11))]
  
  aa_frequencies_pos_nucleators_<-rbind(pos_1, pos_2, pos_3, pos_4, pos_5,
                                        pos_6, pos_7, pos_8, pos_9, pos_10, pos_11)
  
  aa_frequencies_pos_nucleators_<-as.data.frame(aa_frequencies_pos_nucleators_)
  
  # not normalized
  aa_frequencies_pos_nucleators<-t(aa_frequencies_pos_nucleators_)
  
  co_nucleators=melt(aa_frequencies_pos_nucleators)
  co_nucleators<-co_nucleators[co_nucleators$Var1 != "*",]
  
  
  for (i in 1:11){
    if (length(unique(co_nucleators[co_nucleators$Var2 == paste0("pos_", i), "value"])) == 1){
      co_nucleators<-co_nucleators[co_nucleators$Var2 != paste0("pos_", i),]
    }
  }
  
  
  freq_heatmap_nucleators<-ggplot(co_nucleators, aes(factor(Var2, labels=(23:(22+seq_len))), factor(Var1, levels=rev(vectorAA))))+
    geom_tile(aes(fill=value), color="white", size=1)+
    scale_fill_gradient2(high="darkblue", mid= "white", low="black", limit=c(0, 0.25))+
    labs(x="Position", y="Amino acid", fill="Frequency", title="AA frequencies per position in nucleators")+
    theme_minimal()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"))
  print(freq_heatmap_nucleators)
  
  ggsave(freq_heatmap_nucleators, file="AA_frequency_nucleators.jpg", width = 14, height = 8, path=path)
  
  #
  p_scatter_nucleators<-ggplot(co_nucleators, aes(x=factor(Var2, labels=23:(22+seq_len)), y=value, group=1))+
    geom_point(color="darkblue")+
    geom_line(color="darkblue")+
    ylim(c(-0.15, 0.15))+
    facet_wrap(~factor(Var1, levels=vectorAA))+
    labs(x="Position", y="Frequency", title="nucleators")+
    theme_bw()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
      legend.position = c(0.85, 0.2),
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      axis.line=element_line()
    )
  p_scatter_nucleators
  
  ggsave(p_scatter_nucleators, file="AA_frequency_nucleators_pos.jpg", width = 14, height = 8, path=path)
  
  
  # Normalized
  
  aa_frequencies_pos_nucleators<-t(aa_frequencies_pos_nucleators_)-aa_frequencies_pos_all
  
  co_nucleators=melt(aa_frequencies_pos_nucleators)
  co_nucleators<-co_nucleators[co_nucleators$Var1 != "*",]
  
  for (i in 1:11){
    if (length(unique(co_nucleators[co_nucleators$Var2 == paste0("pos_", i), "value"])) == 1){
        co_nucleators<-co_nucleators[co_nucleators$Var2 != paste0("pos_", i),]
    }
  }
  
  freq_heatmap_nucleators<-ggplot(co_nucleators, aes(factor(Var2, labels=(23:(22+seq_len))), factor(Var1, levels=rev(vectorAA))))+
    geom_tile(aes(fill=value), color="white", size=1)+
    scale_fill_gradient2(high="darkblue", mid= "white", low="black", limit=c(-0.05, 0.05))+
    labs(x="Position", y="Amino acid", fill="Frequency", title="AA frequencies per position in nucleators")+
    theme_minimal()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"))
  print(freq_heatmap_nucleators)
  
  ggsave(freq_heatmap_nucleators, file="AA_frequency_nucleators_normall.jpg", width = 14, height = 8, path=path)
  
  #
  p_scatter_nucleators<-ggplot(co_nucleators, aes(x=factor(Var2, labels=23:(22+seq_len)), y=value, group=1))+
    geom_point(color="darkblue")+
    geom_line(color="darkblue")+
    ylim(c(-0.1, 0.1))+
    facet_wrap(~factor(Var1, levels=vectorAA))+
    labs(x="Position", y="Frequency", title="nucleators (Normalized)")+
    theme_bw()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
      legend.position = c(0.85, 0.2),
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      axis.line=element_line()
    )
  
  p_scatter_nucleators
  
  ggsave(p_scatter_nucleators, file="AA_frequency_nucleators_posnormall.jpg", width = 14, height = 8, path=path)
  
  ##################################################################################
  # Work with non-nucleators only
  nonnucleators_nucleators<-NNK_all_df[NNK_all_df$category_fdr == "non-nucleators",]
  
  nonnucleators_nucleators$pos_1<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 1)
  nonnucleators_nucleators$pos_2<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 2)
  nonnucleators_nucleators$pos_3<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 3)
  nonnucleators_nucleators$pos_4<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 4)
  nonnucleators_nucleators$pos_5<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 5)
  nonnucleators_nucleators$pos_6<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 6)
  nonnucleators_nucleators$pos_7<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 7)
  nonnucleators_nucleators$pos_8<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 8)
  nonnucleators_nucleators$pos_9<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 9)
  nonnucleators_nucleators$pos_10<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 10)
  nonnucleators_nucleators$pos_11<-sapply(strsplit(nonnucleators_nucleators$aa_seq_cor, ""), '[', 11)
  
  nucl_df<-nonnucleators_nucleators
  
  pos_1<-table(nucl_df$pos_1)/length(nucl_df$pos_1)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_1))){
      pos_1[aa]<-0}}
  pos_1<-pos_1[order(names(pos_1))]
  
  pos_2<-table(nucl_df$pos_2)/length(nucl_df$pos_2)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_2))){
      pos_2[aa]<-0}}
  pos_2<-pos_2[order(names(pos_2))]
  
  pos_3<-table(nucl_df$pos_3)/length(nucl_df$pos_3)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_3))){
      pos_3[aa]<-0}}
  pos_3<-pos_3[order(names(pos_3))]
  
  pos_4<-table(nucl_df$pos_4)/length(nucl_df$pos_4)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_4))){
      pos_4[aa]<-0}}
  pos_4<-pos_4[order(names(pos_4))]
  
  pos_5<-table(nucl_df$pos_5)/length(nucl_df$pos_5)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_5))){
      pos_5[aa]<-0}}
  pos_5<-pos_5[order(names(pos_5))]
  
  pos_6<-table(nucl_df$pos_6)/length(nucl_df$pos_6)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_6))){
      pos_6[aa]<-0}}
  pos_6<-pos_6[order(names(pos_6))]
  
  pos_7<-table(nucl_df$pos_7)/length(nucl_df$pos_7)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_7))){
      pos_7[aa]<-0}}
  pos_7<-pos_7[order(names(pos_7))]
  
  pos_8<-table(nucl_df$pos_8)/length(nucl_df$pos_8)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_8))){
      pos_8[aa]<-0}}
  pos_8<-pos_8[order(names(pos_8))]
  
  pos_9<-table(nucl_df$pos_9)/length(nucl_df$pos_9)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_9))){
      pos_9[aa]<-0}}
  pos_9<-pos_9[order(names(pos_9))]
  
  pos_10<-table(nucl_df$pos_10)/length(nucl_df$pos_10)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_10))){
      pos_10[aa]<-0}}
  pos_10<-pos_10[order(names(pos_10))]
  
  pos_11<-table(nucl_df$pos_11)/length(nucl_df$pos_11)
  for (aa in vectorAA){
    if (!(aa %in% names(pos_11))){
      pos_11[aa]<-0}}
  pos_11<-pos_11[order(names(pos_11))]
  
  aa_frequencies_pos_nonnucleators_<-rbind(pos_1, pos_2, pos_3, pos_4, pos_5,
                                           pos_6, pos_7, pos_8, pos_9, pos_10, pos_11)
  
  aa_frequencies_pos_nonnucleators_<-as.data.frame(aa_frequencies_pos_nonnucleators_)
  
  # not normalized
  aa_frequencies_pos_nonnucleators<-t(aa_frequencies_pos_nonnucleators_)
  
  co_nonnucleators=melt(aa_frequencies_pos_nonnucleators)
  co_nonnucleators<-co_nonnucleators[co_nonnucleators$Var1 != "*",]
  
  for (i in 1:11){
    if (length(unique(co_nonnucleators[co_nonnucleators$Var2 == paste0("pos_", i), "value"])) == 1){
      co_nonnucleators<-co_nonnucleators[co_nonnucleators$Var2 != paste0("pos_", i),]
    }
  }
  
  freq_heatmap_nonnucleators<-ggplot(co_nonnucleators, aes(factor(Var2, labels=(23:(22+seq_len))), factor(Var1, levels=rev(vectorAA))))+
    geom_tile(aes(fill=value), color="white", size=1)+
    scale_fill_gradient2(high="darkred", mid= "white", low="black", limit=c(0, 0.25))+
    labs(x="Position", y="Amino acid", fill="Frequency", title="AA frequencies per position in nonnucleators Nucleators")+
    theme_minimal()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"))
  print(freq_heatmap_nonnucleators)
  
  ggsave(freq_heatmap_nonnucleators, file="AA_frequency_nonseeder.jpg", width = 14, height = 8, path=path)
  
  #
  p_scatter_nonnucleators<-ggplot(co_nonnucleators, aes(x=factor(Var2, labels=23:(22+seq_len)), y=value, group=1))+
    geom_point(color="darkred")+
    geom_line(color="darkred")+
    ylim(c(-0.15, 0.15))+
    facet_wrap(~factor(Var1, levels=vectorAA))+
    labs(x="Position", y="Frequency", title="Non-nucleators")+
    theme_bw()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
      legend.position = c(0.85, 0.2),
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      axis.line=element_line()
    )
  p_scatter_nonnucleators
  
  ggsave(p_scatter_nonnucleators, file="AA_frequency_nonnucleators_pos.jpg", width = 14, height = 8, path=path)
  
  # Normalized
  
  aa_frequencies_pos_nonnucleators<-t(aa_frequencies_pos_nonnucleators_)-aa_frequencies_pos_all
  
  co_nonnucleators=melt(aa_frequencies_pos_nonnucleators)
  co_nonnucleators<-co_nonnucleators[co_nonnucleators$Var1 != "*",]
  
  for (i in 1:11){
    if (length(unique(co_nonnucleators[co_nonnucleators$Var2 == paste0("pos_", i), "value"])) == 1){
      co_nonnucleators<-co_nonnucleators[co_nonnucleators$Var2 != paste0("pos_", i),]
    }
  }
  
  freq_heatmap_nonnucleators<-ggplot(co_nonnucleators, aes(factor(Var2, labels=(23:(22+seq_len))), factor(Var1, levels=rev(vectorAA))))+
    geom_tile(aes(fill=value), color="white", size=1)+
    scale_fill_gradient2(high="darkred", mid= "white", low="black", limit=c(-0.05, 0.05))+
    labs(x="Position", y="Amino acid", fill="Frequency", title="AA frequencies per position in nonnucleators Nucleators")+
    theme_minimal()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"))
  print(freq_heatmap_nonnucleators)
  
  ggsave(freq_heatmap_nonnucleators, file="AA_frequency_nonseedernormall.jpg", width = 14, height = 8, path=path)
  
  #
  p_scatter_nonnucleators<-ggplot(co_nonnucleators, aes(x=factor(Var2, labels=23:(22+seq_len)), y=value, group=1))+
    geom_point(color="darkred")+
    geom_line(color="darkred")+
    ylim(c(-0.1, 0.1))+
    facet_wrap(~factor(Var1, levels=vectorAA))+
    labs(x="Position", y="Frequency", title="Non-nucleators (Normalized)")+
    theme_bw()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
      legend.position = c(0.85, 0.2),
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      axis.line=element_line()
    )
  p_scatter_nonnucleators
  
  ggsave(p_scatter_nonnucleators, file="AA_frequency_nonnucleators_pos_normall.jpg", width = 14, height = 8, path=path)
  
  
  ## Differential heatmaps
  co_all$ID<-paste0(co_all$Var1, "_", co_all$Var2)
  co_all<-rename(co_all, freq_all = value)
  
  co_nucleators$ID<-paste0(co_nucleators$Var1, "_", co_nucleators$Var2)
  co_nucleators<-rename(co_nucleators, freq_nucleators = value)
  
  co_nonnucleators$ID<-paste0(co_nonnucleators$Var1, "_", co_nonnucleators$Var2)
  co_nonnucleators<-rename(co_nonnucleators, freq_nonnucleators = value)
  
  # Using Normalized frequencies
  
  co_nucleators_nonnucleators<-full_join(co_nucleators[,c("ID", "freq_nucleators", "Var1", "Var2")], 
                                         co_nonnucleators[,c("ID", "freq_nonnucleators")], by="ID")
  
  
  co_nucleators_nonnucleators$diff<-co_nucleators_nonnucleators$freq_nucleators-co_nucleators_nonnucleators$freq_nonnucleators
  
  diff_heatmap_nucleators_nonnucleators<-ggplot(co_nucleators_nonnucleators, aes(factor(Var2, labels=(23:(22+seq_len))), factor(Var1, levels=rev(vectorAA))))+
    geom_tile(aes(fill=diff), color="white", size=1)+
    scale_fill_gradient2(high="darkblue", mid= "white", low="darkred", limit=c(-0.15, 0.15))+
    labs(x="Position", y="Amino acid", fill="Frequency \nDifference")+
    coord_cartesian(clip = "off")+
    theme_minimal()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points")
    )
  
  print(diff_heatmap_nucleators_nonnucleators)
  
  ggsave(diff_heatmap_nucleators_nonnucleators, file="diff_nucleators_nonnucleators.jpg", width = 10, height = 8, path=path)
  
  #
  subsp_median_df<-as.data.frame(co_nucleators_nonnucleators %>% group_by(Var2) %>% dplyr::summarise(median_p=median(diff)))
  co_nucleators_nonnucleators<-left_join(co_nucleators_nonnucleators, subsp_median_df, by="Var2")
  
  cols <- c(colorRampPalette(c("darkred", "grey95"))(1000), colorRampPalette("grey95")(1),
            colorRampPalette(c("grey95",  "darkblue"))(1000))
  
  p_violin_p<-ggplot(co_nucleators_nonnucleators, aes(x=factor(Var2, labels=(23:(22+seq_len))), y=diff))+
    geom_hline(yintercept = 0, size=0.1)+
    geom_violin(scale = "width", size=0.2, aes(fill=median_p))+
    geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
    scale_fill_gradientn(colours=cols, limits=c(-0.08, 0.08), na.value = "grey60")+
    #labs(x="Position", y="Enirchment (Difference in Frequency)")+
    labs(x="Position", y="")+
    theme_classic()+
    theme(plot.margin = unit(c(0,0,0,0), 'cm'),
          legend.position = "none"
    )
  
  p_violin_p
  
  ggsave(p_violin_p, file="p_violinplot_pos.jpg", width=9, height=4)
  
  
  #
  subsm_median_df<-as.data.frame(co_nucleators_nonnucleators %>% group_by(Var1) %>% dplyr::summarise(median_m=median(diff)))
  co_nucleators_nonnucleators<-left_join(co_nucleators_nonnucleators, subsm_median_df, by="Var1")
  
  p_violin_m<-ggplot(co_nucleators_nonnucleators, aes(x=factor(Var1, levels=rev(vectorAA)), y=diff))+
    geom_hline(yintercept = 0, size=0.1)+
    geom_violin(scale = "width", size=0.2, aes(fill=median_m))+
    geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
    scale_fill_gradientn(colours=cols, limits=c(-0.08, 0.08), na.value = "grey60")+
    labs(x="", y="")+
    #labs(x="Amino acid", y="Enirchment (Difference in Frequency)")+
    theme_classic()+
    theme(plot.margin = unit(c(0,0,0,0), 'cm'),
          legend.position = "none"
    )+
    coord_flip() 
  
  p_violin_m
  
  ggsave(p_violin_m, file="p_violinplot_AA.jpg", width=4, height=8)
  
  ##################
  
  p_hv <- ggarrange(diff_heatmap_nucleators_nonnucleators + rremove("x.text"), p_violin_m + rremove("y.text"),
                    ncol=2, align = "h",  legend="none", widths = c(1, 0.2), common.legend=TRUE) 
  
  p_hv
  
  ggsave(p_hv, path=path, file="diff_nucleators_nonnucleators_violinAA.jpg", width=12, height=9)
  
  
  #
  diff_heatmap_nucleators_nonnucleators <- diff_heatmap_nucleators_nonnucleators + theme(axis.title.x = element_blank())
  p_hvc<-ggarrange(diff_heatmap_nucleators_nonnucleators + rremove("x.text"), p_violin_p, 
                   nrow=2, align = "v",  legend="none", heights = c(1, 0.2))
  
  p_hvc
  
  ggsave(p_hvc, path=path, file="p_heatmap_violin_categories_nscore.jpg", width=10, height=10)
  
  ###
  
  #
  co_nucleators$group<-"nucleators"
  co_nonnucleators$group<-"Non-nucleators"
  co_nucleators_nonnucleators<-rbind(rename(co_nucleators, freq = freq_nucleators), rename(co_nonnucleators, freq = freq_nonnucleators))
  
  
  test_position_aa <- function(pos_, aa){
    pos_seeder<-table(nucleators_nucleators[[paste0("pos_",pos_)]])
    
    yes_seeder<-pos_seeder
    no_seeder<-sum(pos_seeder)-pos_seeder
    seeder<-rbind(yes_seeder, no_seeder)
    colnames(seeder)<-paste0(colnames(seeder),"_seeder")
    
    pos_nonseeder<-table(nonnucleators_nucleators[[paste0("pos_",pos_)]])
    
    yes_nonseeder<-pos_nonseeder
    no_nonseeder<-sum(pos_nonseeder)-pos_nonseeder
    nonseeder<-rbind(yes_nonseeder, no_nonseeder)
    colnames(nonseeder)<-paste0(colnames(nonseeder),"_nonseeder")
    
    merged_<-cbind(seeder, nonseeder)
    
    if (length(colnames(merged_[,startsWith(colnames(merged_), aa)])) == 2){
      chisq<-chisq.test(merged_[,startsWith(colnames(merged_), aa)])
      return(chisq$p.value)}
    else{return(1)}
  }
  
  
  test_stars <- function(p.value){
    if (p.value < 0.001){
      return("***")
    }
    else if (p.value < 0.01){
      return("**")
    }
    else if (p.value < 0.05){
      return("*")
    }
    else{
      return("")
    }
  }
  ###
  
  co_nucleators_nonnucleators$significance<-""
  for (aa in vectorAA){
    for (pos_ in 1:seq_len){
      co_nucleators_nonnucleators[co_nucleators_nonnucleators$Var1 == aa & 
                                    co_nucleators_nonnucleators$Var2 == paste0("pos_", pos_) &
                                    co_nucleators_nonnucleators$group == "nucleators", 
                                  "significance"] <- test_stars(test_position_aa(pos_, aa))
    }
  }
  
  
  
  p_scatter_nucleators_nonnucleators<-ggplot(co_nucleators_nonnucleators, aes(x=factor(Var2, labels=23:(22+seq_len)), y=freq ,group=group, color=group))+
    geom_point(size=2)+
    geom_text(aes(label=significance), size=5, color="black", vjust=0.65, angle=90, y=0.08)+
    geom_line(size=1)+
    #scale_x_discrete(breaks=c(5,10,15,20), labels=c(5,10,15,20))+
    scale_color_manual(values=c("darkred", "darkblue"), label=c("non-nucleator", "nucleator"))+
    ylim(c(-0.2, 0.2))+
    facet_wrap(~factor(Var1, levels=vectorAA))+
    labs(x="Position", y="Normalized Frequency")+
    theme_bw()+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
      legend.position = "none",
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      axis.line=element_line()
    )
  p_scatter_nucleators_nonnucleators
  
  ggsave(p_scatter_nucleators_nonnucleators, file="AA_frequency_nucleators_nonnucleators_pos_normall_stars.jpg", width = 12, height = 10, path=path)
  
  
  
  # Join to get aa types
  co_nucleators_nonnucleators<-left_join(co_nucleators_nonnucleators, AA_type[,c("AA", "type")], by=c("Var1"="AA"))
  
  # Group by aa type, group and position and sum the frequencies per group
  co_nucleators_nonnucleators_sum<-co_nucleators_nonnucleators %>% group_by(type, group, Var2) %>% mutate(sum_freq=sum(freq))
  
  
  ####################
  
  test_position_types <- function(pos_, type_){
    aa_types_<-AA_type[AA_type$type == type_, "AA"]
    pos_seeder<-table(nucleators_nucleators[[paste0("pos_",pos_)]])
  
    for (i in aa_types_){
      if (is.na(pos_seeder[i])){
        aa_types_<-aa_types_[!aa_types_ %in% c(i)]
      }
    }
    yes_seeder<-sum(pos_seeder[aa_types_])
    no_seeder<-sum(pos_seeder)-yes_seeder
    seeder<-rbind(yes_seeder, no_seeder)
    colnames(seeder)<-paste0(type_, "_seeder")
    
    aa_types_<-AA_type[AA_type$type == type_, "AA"]
    pos_nonseeder<-table(nonnucleators_nucleators[[paste0("pos_",pos_)]])
  
    for (i in aa_types_){
      if (is.na(pos_nonseeder[i])){
        aa_types_<-aa_types_[!aa_types_ %in% c(i)]
      }
    }
    yes_nonseeder<-sum(pos_nonseeder[aa_types_])
    no_nonseeder<-sum(pos_nonseeder)-yes_nonseeder
    nonseeder<-rbind(yes_nonseeder, no_nonseeder)
    colnames(nonseeder)<-paste0(type_, "_nonseeder")
  
    merged_<-cbind(seeder, nonseeder)
  
    chisq<-chisq.test(merged_)
    if (is.na(chisq$p.value)){
      return(1)}
    else{return(chisq$p.value)}
  
  }
  
  
  
  co_nucleators_nonnucleators_sum$significance<-""
  for (type_ in unique(AA_type$type)){
    for (pos_ in 1:seq_len){
      co_nucleators_nonnucleators_sum[co_nucleators_nonnucleators_sum$type == type_ & 
                                        co_nucleators_nonnucleators_sum$Var2 == paste0("pos_", pos_) &
                                        co_nucleators_nonnucleators_sum$group == "nucleators", 
                                      "significance"] <- test_stars(test_position_types(pos_, type_))
    }
  }
  
  
  
  p_scatter_nucleators_nonnucleators_types<-ggplot(co_nucleators_nonnucleators_sum, aes(x=factor(Var2, labels=23:(22+seq_len)), y=sum_freq ,group=group, color=group))+
    geom_point(size=2)+
    geom_line(size=1)+
    geom_text(aes(label=significance), size=5, color="black", vjust=0.65, angle=90, y=0.19)+
    scale_color_manual(values=c("darkred", "darkblue"), label=c("non-nucleator", "nucleator"))+
    ylim(c(-0.21, 0.21))+
    facet_wrap(~factor(type, levels=c("glycine", "proline", "cysteine", "positive", "negative", "polar", "aliphatic", "aromatic")))+  labs(x="Position", y="Normalized Frequency")+
    labs(color="")+
    coord_cartesian(clip = "off")+
    theme_bw()+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
    theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
      legend.position = c(0.85, 0.2),
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      axis.line=element_line()
    )
  p_scatter_nucleators_nonnucleators_types
  
  ggsave(p_scatter_nucleators_nonnucleators_types, file="AA_frequency_nucleators_nonnucleators_types_pos_normall_stars.jpg", width = 8, height = 8, path=path)
  
  #
  
  save(co_nucleators_nonnucleators, co_nucleators_nonnucleators_sum,  file = paste0("freq_", seq_len, ".RData"))
}

################################################################################
# Let's merge the datasets

setwd("C:/Users/mmartin/OneDrive - IBEC/Documentos/Postdoc-BEBlab/Lab/LabProjects/BRI2_ADan_ABri/Sequencing_Analysis/Bri2_NNK_all_240417")

path="02_AAcomposition_Position_Truncations"

### AA
#required data:
load(paste0("freq_11.RData"))
co_nucleators_nonnucleators$aa_len<-11
co_nucleators_nonnucleators_all<-co_nucleators_nonnucleators

seq_len<-2:10

for (i in seq_len){
  load(paste0("freq_", i, ".RData"))
  co_nucleators_nonnucleators$aa_len<-i
  co_nucleators_nonnucleators_all<-rbind(co_nucleators_nonnucleators_all, co_nucleators_nonnucleators)
}

p_scatter_nucleators_nonnucleators<-ggplot(co_nucleators_nonnucleators_all, aes(x=factor(Var2, labels=(23:33)), y=freq ,group=group, color=group))+
  geom_point(size=2)+
  geom_text(aes(label=significance), size=5, color="black", vjust=0.65, angle=90, y=0.12)+
  geom_line(size=1)+
  scale_color_manual(values=c("darkred", "darkblue"), label=c("non-nucleator", "nucleator"))+
  ylim(c(-0.2, 0.2))+
  facet_grid(aa_len~factor(Var1, levels=vectorAA))+
  labs(x="Position", y="Normalized Frequency")+
  theme_bw()+
  theme(
    #plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
    #legend.position = c(1.2, .5),
    strip.text = element_text(size = 12),
    axis.line=element_line()
  )
p_scatter_nucleators_nonnucleators

ggsave(p_scatter_nucleators_nonnucleators, file="AA_frequency_nucleators_nonnucleators_pos_normall_stars_truncations.jpg", width = 40, height = 15, path=path)

### AA type

#required data:
load(paste0("freq_11.RData"))
co_nucleators_nonnucleators_sum$aa_len<-11
co_nucleators_nonnucleators_all<-co_nucleators_nonnucleators_sum

seq_len<-2:10

for (i in seq_len){
  load(paste0("freq_", i, ".RData"))
  co_nucleators_nonnucleators_sum$aa_len<-i
  co_nucleators_nonnucleators_all<-rbind(co_nucleators_nonnucleators_all, co_nucleators_nonnucleators_sum)
}

p_scatter_nucleators_nonnucleators_types<-ggplot(co_nucleators_nonnucleators_all, aes(x=factor(Var2, labels=(23:33)), y=sum_freq ,group=group, color=group))+
  geom_point(size=2)+
  geom_line(size=1)+
  geom_text(aes(label=significance), size=5, color="black", vjust=0.65, angle=90, y=0.19)+
  scale_color_manual(values=c("darkred", "darkblue"), label=c("non-nucleator", "nucleator"))+
  ylim(c(-0.21, 0.21))+
  facet_grid(aa_len~factor(type, levels=c("glycine", "proline", "cysteine", "positive", "negative", "polar", "aliphatic", "aromatic")))+  
  labs(x="Position", y="Normalized Frequency")+
  labs(color="")+
  coord_cartesian(clip = "off")+
  theme_bw()+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(
    #plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
    #legend.position = c(1.2, .5),
    strip.text = element_text(size = 12),
    axis.line=element_line()
  )
p_scatter_nucleators_nonnucleators_types

ggsave(p_scatter_nucleators_nonnucleators_types, file="AA_frequency_nucleators_nonnucleators_types_pos_normall_stars_truncations.jpg", width = 20, height = 15, path=path)

###

p_scatter_nucleators_nonnucleators_types<-ggplot(co_nucleators_nonnucleators_all[co_nucleators_nonnucleators_all$type == "aromatic",], 
                                                 aes(x=factor(Var2, labels=(23:33)), y=sum_freq ,group=group, color=group))+
  geom_point(size=2)+
  geom_line(size=1)+
  geom_text(aes(label=significance), size=5, color="black", vjust=0.65, angle=90, y=0.21)+
  scale_color_manual(values=c("darkred", "darkblue"), label=c("non-nucleator", "nucleator"))+
  ylim(c(-0.23, 0.23))+
  facet_wrap(~aa_len, ncol=5)+  
  labs(x="Position", y="Normalized Frequency", title="Aromatics")+
  labs(color="")+
  coord_cartesian(clip = "off")+
  theme_bw()+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(
    #plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
    #legend.position = c(1.2, .5),
    strip.text = element_text(size = 12),
    axis.line=element_line(),
    plot.title = element_text(hjust = 0.5)
  )
p_scatter_nucleators_nonnucleators_types

ggsave(p_scatter_nucleators_nonnucleators_types, file="AA_frequency_nucleators_nonnucleators_aromatics_pos_normall_stars_truncations.jpg", width = 12, height = 5, path=path)
