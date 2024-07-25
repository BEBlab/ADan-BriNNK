library(tidyverse)
require(reshape2)
library(ggpubr)
require(stringr)
library(ggsignif)
require("ggrepel")
library(lsr)

dir.create("01_AAcomposition")
path="01_AAcomposition"


# Mode value
mode <- function(x, na.rm = FALSE) {
  
  if(na.rm){ #if na.rm is TRUE, remove NA values from input x
    x = x[!is.na(x)]
  }
  
  val <- unique(x)
  return(val[which.max(tabulate(match(x, val)))])
}

#required data:

name<-"NNK_all_df.RData"

load(name)

# Let's see if there's an enrichment in particular amino-acids per position for the different groups

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


# Let's work only with sequences that are 12 aa long
NNK_all_df<-NNK_all_df[NNK_all_df$aa_len == 12,]

plots_aa<-list()
plots_aa_corrected<-list()
for (aa in vectorAA){
  NNK_all_df[aa]<-str_count(NNK_all_df$aa_seq_cor, aa)
  NNK_all_df[paste0(aa, "_perc")]<-NNK_all_df[aa]/NNK_all_df$aa_len*100
  
  y_name<-paste0(aa, "_perc")
  aa_plot<-ggplot(NNK_all_df, aes(x=factor(.data[["mode_seed_bh"]], levels=c(0, 1), labels=c("non-nucleators", "nucleators")), y=.data[[y_name]]))+
    geom_jitter(aes(color=factor(.data[["mode_seed_bh"]], levels=c(0, 1), labels=c("non-nucleators", "nucleators"))), 
                alpha=0.1)+
    geom_boxplot(aes(fill=factor(.data[["mode_seed_bh"]], levels=c(0, 1), labels=c("non-nucleators", "nucleators"))), 
                 alpha=0.8, width=0.5, lwd=0.5, outlier.shape = NA)+
    scale_fill_manual(values=c("darkred", "darkblue"), labels=c("non-nucleators", "nucleators"))+
    scale_color_manual(values=c("darkred", "darkblue"), labels=c("non-nucleators", "nucleators"))+
    geom_signif(comparisons = list(c("non-nucleators", "nucleators")),
                map_signif_level = TRUE, tip_length = 0, vjust = 2)+
    labs(x="", y="%", fill="", color="", title=aa)+
    theme_classic()+
    theme(plot.title = element_text(size=16, hjust = 0.5),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14), 
          axis.title = element_text(size = 14),
          axis.text.x = element_text(color="black", size=14),
          axis.text.y = element_text(color="black", size=14))
  plots_aa[[aa]]<-aa_plot
}


plots_aa_<-ggarrange(plotlist=plots_aa, ncol=5, nrow=4, common.legend = T)
plots_aa_
ggsave(plots_aa_, file="plots_aa.jpg", width = 16, height = 16, path = path)

# 

aa_types<-c("glycine", "proline", "cysteine", "positive", "negative", "polar", "aliphatic", "aromatic")

plots_aa<-list()
plots_aa_corrected<-list()
for (type_ in aa_types){
  NNK_all_df[type_]<-rowSums(NNK_all_df[c(AA_type[AA_type$type == type_, "AA"])])
  NNK_all_df[paste0(type_, "_perc")]<-NNK_all_df[type_]/NNK_all_df$aa_len*100
  
  y_name<-paste0(type_, "_perc")
  aa_plot<-ggplot(NNK_all_df, aes(x=factor(.data[["mode_seed_bh"]], levels=c(0, 1), labels=c("non-nucleators", "nucleators")), y=.data[[y_name]]))+
    geom_jitter(aes(color=factor(.data[["mode_seed_bh"]], levels=c(0, 1), labels=c("non-nucleators", "nucleators"))), 
                alpha=0.1)+
    geom_boxplot(aes(fill=factor(.data[["mode_seed_bh"]], levels=c(0, 1), labels=c("non-nucleators", "nucleators"))), 
                 alpha=0.8, width=0.5, lwd=0.5, outlier.shape = NA)+
    scale_fill_manual(values=c("darkred", "darkblue"), labels=c("non-nucleators", "nucleators"))+
    scale_color_manual(values=c("darkred", "darkblue"), labels=c("non-nucleators", "nucleators"))+
    geom_signif(comparisons = list(c("non-nucleators", "nucleators")),
                map_signif_level = TRUE, tip_length = 0, vjust = 2)+
    labs(x="", y="%", fill="", color="", title=type_)+
    theme_classic()+
    theme(plot.title = element_text(size=16, hjust = 0.5),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14), 
          axis.title = element_text(size = 14),
          axis.text.x = element_text(color="black", size=14),
          axis.text.y = element_text(color="black", size=14))
  plots_aa[[type_]]<-aa_plot
}


plots_aa_<-ggarrange(plotlist=plots_aa, common.legend = T)
plots_aa_
ggsave(plots_aa_, file="plots_aatype.jpg", width = 16, height = 16, path=path)


####
# Scatter plots
columns_interest<-c("aa_seq_cor", "mode_seed_bh", "G_perc", "A_perc", "V_perc", "L_perc",
                    "M_perc", "I_perc", "F_perc", "Y_perc", "W_perc", "K_perc",
                    "R_perc", "D_perc", "E_perc", "S_perc", "T_perc", "C_perc",
                    "N_perc", "Q_perc", "H_perc", "P_perc")
new_df<-NNK_all_df[,columns_interest]
melted_df<-melt(new_df, id=c("aa_seq_cor", "mode_seed_bh"))
non_nucleators_perc<-melted_df[melted_df$mode_seed_bh == 0,]
non_nucleators_perc<-rename(non_nucleators_perc, AA = variable, perc = value)
non_nucleators_perc<-non_nucleators_perc[,c("aa_seq_cor", "AA", "perc")]
non_nucleators_test<-non_nucleators_perc
non_nucleators_perc<-non_nucleators_perc %>% group_by(AA) %>% mutate(mean_perc_ns = mean(perc), std_perc_ns = sd(perc))
non_nucleators_perc<-non_nucleators_perc[, c("AA", "mean_perc_ns", "std_perc_ns")]
non_nucleators_perc<-distinct(non_nucleators_perc)
nucleators_perc<-melted_df[melted_df$mode_seed_bh == 1,]
nucleators_perc<-rename(nucleators_perc, AA = variable, perc = value)
nucleators_perc<-nucleators_perc[,c("aa_seq_cor", "AA", "perc")]
nucleators_test<-nucleators_perc
nucleators_perc<-nucleators_perc %>% group_by(AA) %>% mutate(mean_perc_s = mean(perc), std_perc_s = sd(perc))
nucleators_perc<-nucleators_perc[, c("AA", "mean_perc_s", "std_perc_s")]
nucleators_perc<-distinct(nucleators_perc)

# t-test and cohen's d effect size
non_nucleators_test$perc_ns<-non_nucleators_test$perc
non_nucleators_test<-non_nucleators_test[c("AA", "perc_ns")]
nucleators_test$perc_s<-nucleators_test$perc
nucleators_test<-nucleators_test[c("AA", "perc_s")]

aa_<-c()
pvalues<-c()
freq_nucleators<-c()
freq_nonnucleators<-c()
differences<-c()
effectsizes<-c()

print("t-test and Cohen's:")
for(aa in c("G_perc", "A_perc", "V_perc", "L_perc",
            "M_perc", "I_perc", "F_perc", "Y_perc", "W_perc", "K_perc",
            "R_perc", "D_perc", "E_perc", "S_perc", "T_perc", "C_perc",
            "N_perc", "Q_perc", "H_perc", "P_perc")){
  print("***********************************************************")
  print(aa)
  aa_<-c(aa_, aa)
  test<-t.test(
    x           = nucleators_test[nucleators_test$AA == aa, "perc_s"],
    y           = non_nucleators_test[non_nucleators_test$AA == aa, "perc_ns"],
    alternative = "two.sided",
    mu          = 0,
    var.equal   = TRUE,
    conf.level  = 0.95)
  
  print(paste("t-test p-value:", test$p.value))
  pvalues<-c(pvalues,test$p.value)
  
  freq_nucleators<-c(freq_nucleators, test$estimate[1])
  freq_nonnucleators<-c(freq_nonnucleators, test$estimate[2])
  
  print(paste("means difference:", test$estimate[1]-test$estimate[2]))
  differences<-c(differences,test$estimate[1]-test$estimate[2])
  
  effectsize<-cohensD(x=nucleators_test[nucleators_test$AA == aa, "perc_s"],
                      y=non_nucleators_test[non_nucleators_test$AA == aa, "perc_ns"])
  print(paste("Cohen's d effect-size:", effectsize))
  effectsizes<-c(effectsizes,effectsize)
}

table_test<-data.frame("AA"=aa_,
                       "freq_nucleators"=freq_nucleators,
                       "freq_nonnucleators"=freq_nonnucleators,
                       "t-test pvalue"=pvalues,
                       "Mean difference"=differences,
                       "Cohen's effect size"=effectsizes)
write_tsv(table_test, "Supplementary Table 2.tsv")

# Prepare for plotting

perc_new_df<-merge(nucleators_perc, non_nucleators_perc, by="AA")

perc_new_df<-separate_wider_delim(perc_new_df, AA, delim="_", names=c("AA", "perc"))

perc_new_df<-merge(perc_new_df, AA_type, by="AA")

# Significance according to t-test
no_sign<-c("W", "Y", "E", "T", "Q", "P")


perc_new_df$significance<-"p-value<0.05"
perc_new_df[perc_new_df$AA %in% no_sign, "significance"]<-"ns"

scatter_perc <- ggplot(perc_new_df, aes(x = mean_perc_ns, y = mean_perc_s, color = type, label = AA)) +
  geom_abline(slope=1, linewidth=.5, linetype = "dashed", color = "grey80")+
  geom_point(aes(shape = factor(significance, levels = c("ns", "p-value<0.05"))), size = 2, stroke = 2) +
  geom_text_repel(aes(label = AA), box.padding = 0.5, max.overlaps = Inf, 
                  min.segment.length=0.1, size=4) +
  scale_color_manual(values = color_AAtype, guide=guide_legend()) +
  scale_shape_manual(values = c("ns" = 1, "p-value<0.05" = 19)) +
  labs(color = "", x = "AA Percentage \n(non-nucleator sequences)", y = "AA Percentage \n(nucleator sequences)", shape = "") +
  xlim(1, 11) +
  ylim(1, 11) +
  coord_cartesian(clip = "off")+
  theme_classic() +
  theme(
    plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points")
  )

scatter_perc
ggsave(scatter_perc, file="scatter_AA_perc.jpg", width = 5.5, height = 4, path=path)


#
scatter_perc <- ggplot(perc_new_df, aes(x = mean_perc_ns, y = mean_perc_s, color = type, label = AA)) +
  geom_abline(slope=1, linewidth=.5, linetype = "dashed", color = "grey80")+
  geom_point(aes(shape = factor(significance, levels = c("ns", "p-value<0.05"))), size = 2, stroke = 2) +
  geom_text_repel(aes(label = AA), box.padding = 0.5, max.overlaps = Inf, 
                  min.segment.length=0.1, size=4) +
  scale_color_manual(values = color_AAtype, guide=guide_legend()) +
  scale_shape_manual(values = c("ns" = 1, "p-value<0.05" = 19)) +
  labs(color = "", x = "AA Percentage \n(non-nucleator sequences)", y = "AA Percentage \n(nucleator sequences)", shape = "") +
  xlim(1, 11) +
  ylim(1, 11) +
  coord_cartesian(clip = "off")+
  theme_classic() +
  theme(
    plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points"),
    legend.position = "none"
  )

scatter_perc
ggsave(scatter_perc, file="scatter_AA_perc_nolegend.jpg", width = 4, height = 4, path=path)


####

columns_interest<-c("aa_seq_cor", "mode_seed_bh", "aliphatic_perc", "aromatic_perc",
                    "positive_perc",  "negative_perc",  "polar_perc", "glycine_perc",
                    "proline_perc", "cysteine_perc" )

new_df<-NNK_all_df[,columns_interest]

melted_df<-melt(new_df, id=c("aa_seq_cor", "mode_seed_bh"))

non_nucleators_perc<-melted_df[melted_df$mode_seed_bh == 0,]
non_nucleators_perc<-rename(non_nucleators_perc, AA = variable, perc = value)
non_nucleators_perc<-non_nucleators_perc[,c("aa_seq_cor", "AA", "perc")]
non_nucleators_perc<-non_nucleators_perc %>% group_by(AA) %>% mutate(mean_perc_ns = mean(perc), std_perc_ns = sd(perc))
non_nucleators_perc<-non_nucleators_perc[, c("AA", "mean_perc_ns", "std_perc_ns")]
non_nucleators_perc<-distinct(non_nucleators_perc)


nucleators_perc<-melted_df[melted_df$mode_seed_bh == 1,]
nucleators_perc<-rename(nucleators_perc, AA = variable, perc = value)
nucleators_perc<-nucleators_perc[,c("aa_seq_cor", "AA", "perc")]
nucleators_perc<-nucleators_perc %>% group_by(AA) %>% mutate(mean_perc_s = mean(perc), std_perc_s = sd(perc))
nucleators_perc<-nucleators_perc[, c("AA", "mean_perc_s", "std_perc_s")]
nucleators_perc<-distinct(nucleators_perc)

perc_new_df<-merge(nucleators_perc, non_nucleators_perc, by="AA")

perc_new_df<-separate_wider_delim(perc_new_df, AA, delim="_", names=c("AA", "perc"))


scatter_perc<-ggplot(perc_new_df, aes(x=mean_perc_ns, y=mean_perc_s, color=AA, label=AA))+
  geom_abline(slope=1, linewidth=.5, linetype = "dashed", color = "grey80")+
  geom_point(size=2, shape=19, stroke=2)+
  geom_text_repel(aes(label = AA), box.padding = 0.3, max.overlaps = Inf, 
                  min.segment.length=0.1, size=4) +  scale_color_manual(values=color_AAtype)+
  labs(color="AA type", x="AA Percentage \n(non-nucleator sequences)", y="AA Percentage \n(nucleator sequences)")+
  theme_classic() +
  theme(
    plot.margin = ggplot2::unit(c(1, 1, 1, 1), "points")
  )

scatter_perc
ggsave(scatter_perc, file="scatter_types_perc.jpg", width = 5.5, height = 4, path=path)

##############
#required data:

name<-"NNK_all_df.RData"

load(name)

#NNK_all_df<-NNK_all_df[NNK_all_df$aa_len<12,]

NNK_all_df_grouped_category <- NNK_all_df %>% group_by(aa_len, category_fdr) %>% count(category_fdr)
NNK_all_df_grouped_category <- NNK_all_df_grouped_category %>% group_by(aa_len) %>% mutate(percent = n/sum(n)*100)

#NNK_all_df_grouped_category <- NNK_all_df_grouped_category[NNK_all_df_grouped_category$category_fdr != "non-nucleators",]
NNK_all_df_grouped_category <- rbind(NNK_all_df_grouped_category, data.frame("aa_len"=1, "category_fdr"="Top10 nucleators", n=0, "percent"=0))

bars_plot<-ggplot(NNK_all_df_grouped_category, aes(x=factor(aa_len, levels=(1:12), labels=(23:34)),
                                                   y=percent,
                                                   fill=factor(category_fdr, levels=c("non-nucleators", "nucleators", "Top10 nucleators"))))+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c("#DF9292", "#7979BE", "darkblue"))+
  labs(x="Extension Length", y="Percentage (%)", fill="")+
  theme_classic()+
  theme(legend.position = "none")
bars_plot

ggsave(bars_plot, file="truncations_nucleators_perc.jpg", width = 4, height = 4, path = path)

###

save(NNK_all_df, file = name)

write.csv(NNK_all_df, "NNK_all_df.csv")