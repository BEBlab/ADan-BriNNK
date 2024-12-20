# ThT - Kinetics

library(tidyverse)
library(ggpubr)
library(reshape2)
library(DescTools)
library(ggbreak) 
library(writexl)  


dir.create("Plots")
path="Plots"

########
# Bri2, ABri & ADan

file_name<-"Bri2_ABri_ADan_ThT.txt"

tecan_results<-read_delim(file_name, locale=locale(encoding="latin1"), delim="\t")

################################################################################
tecan_melted<-melt(select(tecan_results, !c("Time [s]", "Temp. [ C]", "Time_hs")),
                   id="Time")
tecan_melted$value<-as.numeric(tecan_melted$value)

tecan_melted$plate_column<-substring(tecan_melted$variable, 2, 3)
tecan_melted$plate_row<-substring(tecan_melted$variable, 1, 1)

tecan_melted[tecan_melted$plate_column == 1, "conc"]<-"blank"
tecan_melted[tecan_melted$plate_column == 2, "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 3, "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 4, "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 5, "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 6, "conc"]<-"0.75uM"
tecan_melted[tecan_melted$plate_column == 7, "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 8, "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 9, "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 10, "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 11, "conc"]<-"0.75uM"

tecan_melted$peptide<-"blank"
tecan_melted[tecan_melted$plate_row %in% c("B", "C", "D") & tecan_melted$plate_column %in% c(2, 3, 4, 5, 6), "peptide"]<-"Bri2"
tecan_melted[tecan_melted$plate_row %in% c("B", "C", "D") & tecan_melted$plate_column %in% c(7, 8, 9, 10, 11), "peptide"]<-"ABri"
tecan_melted[tecan_melted$plate_row %in% c("E", "F", "G") & tecan_melted$plate_column %in% c(7, 8, 9, 10, 11), "peptide"]<-"ADan"

tecan_melted$sample<-paste0(tecan_melted$peptide, "-", tecan_melted$conc)

# Substract the blank and normalize
tecan_melted$value_correct <- tecan_melted$value-min(tecan_melted[tecan_melted$peptide == "blank", "value"])

tecan_melted_grouped <- tecan_melted %>% group_by(Time, sample) %>% mutate(mean_value=mean(value_correct), std_value=sd(value_correct))

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(conc) %>% mutate(value_norm=(value_correct/max(mean_value))*100)

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))

# Raw
df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("6uM", "3uM", "1.5uM", "0.75uM"),]

scatter_tht<-ggplot(df_to_plot, 
                    aes(y=mean_value, x=Time, 
                        color=factor(peptide, levels=c("ADan", "ABri", "Bri2"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("6uM", "3uM", "1.5uM", "0.75uM"), labels=c("6 μM", "3 μM", "1.5 μM", "0.75 μM")))+
  scale_colour_manual(values=c("grey20", "grey50", "grey80"))+
  theme_bw()+
  labs(x="Time (hs)", y="ThT Fluorescence (rfu)", color="")
scatter_tht

ggsave(scatter_tht, file="20241127_Bri2_ABri_ADan_ThT_std.jpg", width=5, height=4, path=path)


# Normalized
scatter_tht_normalized<-ggplot(df_to_plot, 
                             aes(y=mean_value_norm, x=Time, 
                                 color=factor(peptide, levels=c("ADan", "ABri", "Bri2"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_norm-std_value_norm, ymax=mean_value_norm+std_value_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("6uM", "3uM", "1.5uM", "0.75uM"), labels=c("6 μM", "3 μM", "1.5 μM", "0.75 μM")))+
  scale_colour_manual(values=c("grey20", "grey50", "grey80"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized

ggsave(scatter_tht_normalized, file="20241127_Bri2_ABri_ADan_ThT_std_norm.jpg", width=5, height=4, path=path)

#

scatter_both <- ggarrange(scatter_tht, scatter_tht_normalized, common.legend = T)
ggsave(scatter_both, file="20241127_Bri2_ABri_ADan_ThT_both.jpg", width=8, height=4, path=path)

# Save as xlsx
tecan_melted_grouped<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("blank", "6uM", "3uM", "1.5uM", "0.75uM"),]
write_xlsx(tecan_melted_grouped,"Bri2_ABri_ADan_ThT.xlsx")



########
# Bri2 & ABri - 12 uM

# ThT - Kinetics

file_name<-"Bri2_ABri_ThT.txt"

tecan_results<-read_delim(file_name, locale=locale(encoding="latin1"), delim="\t")

################################################################################
tecan_melted<-melt(select(tecan_results, !c("Time [s]", "Temp. [ C]", "Time_hs")),
                   id="Time")
tecan_melted$value<-as.numeric(tecan_melted$value)

tecan_melted$plate_column<-substring(tecan_melted$variable, 2, 3)
tecan_melted$plate_row<-substring(tecan_melted$variable, 1, 1)

tecan_melted[tecan_melted$plate_column == 1, "conc"]<-"blank"
tecan_melted[tecan_melted$plate_column == 3, "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 4, "conc"]<-"12uM"


tecan_melted$peptide<-"blank"
tecan_melted[tecan_melted$plate_row %in% c("E", "F", "G") & tecan_melted$plate_column %in% c(3), "peptide"]<-"Bri2"
tecan_melted[tecan_melted$plate_row %in% c("E", "F", "G") & tecan_melted$plate_column %in% c(4), "peptide"]<-"ABri"

tecan_melted$sample<-paste0(tecan_melted$peptide, "-", tecan_melted$conc)

# Substract the blank and normalize
tecan_melted$value_correct <- tecan_melted$value-min(tecan_melted[tecan_melted$peptide == "blank", "value"])

tecan_melted_grouped <- tecan_melted %>% group_by(Time, sample) %>% mutate(mean_value=mean(value_correct), std_value=sd(value_correct))

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(conc) %>% mutate(value_norm=(value_correct/max(mean_value))*100)

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))

df_to_plot <- tecan_melted_grouped[tecan_melted_grouped$peptide != "blank",]

# Raw
scatter_tht<-ggplot(df_to_plot, 
                    aes(y=mean_value, x=Time, 
                    color=factor(peptide, levels=c("Bri2", "ABri"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc))+
  scale_y_continuous(breaks=c(0, 20000, 40000, 60000))+
  ylim(0, 60000)+
  scale_colour_manual(values=c("grey80", "grey50"))+
  theme_bw()+
  labs(x="Time (hs)", y="ThT Fluorescence (rfu)", color="")
scatter_tht

ggsave(scatter_tht, file="20241129_Bri2_ABri_ThT_std.jpg", width=3, height=2, path=path)

# Normalized to max value per concentration
scatter_tht_normalized<-ggplot(df_to_plot, 
                             aes(y=mean_value_norm, x=Time, 
                                 color=factor(peptide, levels=c("Bri2", "ABri"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_norm-std_value_norm, ymax=mean_value_norm+std_value_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("12uM"), labels=c("12 μM")))+
  scale_y_continuous(breaks=c(20, 40, 60, 80, 100))+
  scale_colour_manual(values=c("grey80", "grey50"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized

ggsave(scatter_tht_normalized, file="20241129_Bri2_ABri_ThT_std_norm.jpg", width=3, height=2, path=path)

# Normalized to each trace

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_correct_trace=value_correct-min(mean_value))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace=mean(value_correct_trace), std_value_trace=sd(value_correct_trace))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_trace_norm=(value_correct_trace/max(mean_value_trace))*100)
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

df_to_plot <- tecan_melted_grouped[tecan_melted_grouped$peptide != "blank",]

# 
scatter_tht_normalized_trace<-ggplot(df_to_plot, 
                               aes(y=mean_value_trace_norm, x=Time, 
                                   color=factor(peptide, levels=c("Bri2", "ABri"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("12uM"), labels=c("12 μM")))+
  #scale_y_continuous(breaks=c(20, 40, 60, 80, 100))+
  scale_colour_manual(values=c("grey80", "grey50"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

ggsave(scatter_tht_normalized_trace, file="20241129_Bri2_ABri_ThT_std_norm_trace.jpg", width=3, height=2, path=path)

#

scatter_all <- ggarrange(scatter_tht, scatter_tht_normalized_trace, common.legend = T, ncol=2)
scatter_all
ggsave(scatter_all, file="20241129_Bri2_ABri_ThT.jpg", width=6, height=3, path=path)

# Save as xlsx
write_xlsx(tecan_melted_grouped,"Bri2_ABri_ThT.xlsx")



########
# ADan vs mutants

file_name<-"ADan_mutants_ThT.txt"

tecan_results<-read_delim(file_name, locale=locale(encoding="latin1"), delim="\t")

################################################################################
tecan_melted<-melt(select(tecan_results, !c("Time [s]", "Temp. [ C]", "Time_hs")),
                   id="Time")
tecan_melted$value<-as.numeric(tecan_melted$value)

tecan_melted$plate_column<-substring(tecan_melted$variable, 2, 3)
tecan_melted$plate_row<-substring(tecan_melted$variable, 1, 1)

tecan_melted[tecan_melted$plate_column == 1, "conc"]<-"blank"
tecan_melted[tecan_melted$plate_column == 2, "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 3, "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 4, "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 5, "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 6, "conc"]<-"0.75uM"
tecan_melted[tecan_melted$plate_column == 7, "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 8, "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 9, "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 10, "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 11, "conc"]<-"0.75uM"

tecan_melted$peptide<-"blank"
tecan_melted[tecan_melted$plate_row %in% c("B", "C", "D") & tecan_melted$plate_column %in% c(2, 3, 4, 5, 6), "peptide"]<-"ADan"
tecan_melted[tecan_melted$plate_row %in% c("E", "F", "G") & tecan_melted$plate_column %in% c(2, 3, 4, 5, 6), "peptide"]<-"L27*"
tecan_melted[tecan_melted$plate_row %in% c("E", "F", "G") & tecan_melted$plate_column %in% c(7, 8, 9, 10, 11), "peptide"]<-"N13G"

tecan_melted$sample<-paste0(tecan_melted$peptide, "-", tecan_melted$conc)


# Substract the blank and normalize
tecan_melted$value_correct <- tecan_melted$value-min(tecan_melted[tecan_melted$peptide == "blank", "value"])

tecan_melted_grouped <- tecan_melted %>% group_by(Time, sample) %>% mutate(mean_value=mean(value_correct), std_value=sd(value_correct))

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(conc) %>% mutate(value_norm=(value_correct/max(mean_value))*100)

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))

#
df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("3uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]

# Raw
scatter_tht<-ggplot(df_to_plot, 
                             aes(y=mean_value, x=Time, 
                                 color=factor(peptide, levels=c("ADan", "N13G", "L27*"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM"), labels=c("3 μM")))+
  scale_colour_manual(values=c("grey20", "lightblue", "blue3"))+
  theme_bw()+
  labs(x="Time (hs)", y="ThT Fluorescence (rfu)", color="")
scatter_tht

ggsave(scatter_tht, file="20241209_ADan_mutants_ThT_std.jpg", width=3, height=2, path=path)


# Normalized to max of each conc
scatter_tht_normalized<-ggplot(df_to_plot, 
                             aes(y=mean_value_norm, x=Time, 
                                 color=factor(peptide, levels=c("ADan", "N13G", "L27*"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_norm-std_value_norm, ymax=mean_value_norm+std_value_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM"), labels=c("3 μM")))+
  scale_colour_manual(values=c("grey20", "lightblue", "blue3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized

ggsave(scatter_tht_normalized, file="20241209_ADan_mutants_ThT_std_norm.jpg", width=3, height=2, path=path)

# Normalized to each trace

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_correct_trace=value_correct-min(mean_value))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace=mean(value_correct_trace), std_value_trace=sd(value_correct_trace))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_trace_norm=(value_correct_trace/max(mean_value_trace))*100)
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("3uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]
# 
scatter_tht_normalized_trace<-ggplot(df_to_plot, 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("ADan", "N13G", "L27*"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM"), labels=c("3 μM")))+
  scale_colour_manual(values=c("grey20", "lightblue", "blue3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

ggsave(scatter_tht_normalized_trace, file="20241209_ADan_mutants_std_norm_trace.jpg", width=3, height=2, path=path)

#

scatter_all <- ggarrange(scatter_tht, scatter_tht_normalized_trace, common.legend = T, ncol=2)
scatter_all
ggsave(scatter_all, file="20241209_ADan_mutants_ThT_all.jpg", width=6, height=3, path=path)

# Save as xlsx
tecan_melted_grouped<-tecan_melted_grouped[tecan_melted_grouped$conc != "12uM",]
write_xlsx(tecan_melted_grouped,"ADan_mutants_ThT.xlsx")

###
# Use loess fit to predict t1/2
# Check loess fitting
scatter_tht_normalized_trace<-ggplot(df_to_plot[df_to_plot$Time<5,], 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("ADan", "N13G", "L27*"))))+
  geom_point()+
  geom_smooth(method="loess")+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM"), labels=c("3 μM")))+
  scale_colour_manual(values=c("grey20", "lightblue", "blue3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace


loess_model_adan <- loess(Time ~ mean_value_trace_norm, data = tecan_melted_grouped[tecan_melted_grouped$Time<4.3 & tecan_melted_grouped$sample == "ADan-3uM",])
t1_2_adan<-predict(loess_model_adan, 50)

loess_model_L27 <- loess(Time ~ mean_value_trace_norm, data = tecan_melted_grouped[tecan_melted_grouped$Time<4.3 & tecan_melted_grouped$sample == "L27*-3uM",])
t1_2_L27<-predict(loess_model_L27, 50)

loess_model_N13G <- loess(Time ~ mean_value_trace_norm, data = tecan_melted_grouped[tecan_melted_grouped$Time<4.3 & tecan_melted_grouped$sample == "N13G-3uM",])
t1_2_N13G<-predict(loess_model_N13G, 50)

# plot t1/2
scatter_tht_normalized_trace<-ggplot(df_to_plot[df_to_plot$Time<4.3,], 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("ADan", "N13G", "L27*"))))+
  geom_vline(xintercept=t1_2_adan, color="grey20")+
  geom_hline(yintercept=50, color="grey20")+
  geom_vline(xintercept=t1_2_L27, color="blue3")+
  geom_hline(yintercept=50, color="blue3")+
  geom_vline(xintercept=t1_2_N13G, color="lightblue")+
  geom_hline(yintercept=50, color="lightblue")+
  geom_point()+
  geom_smooth(method="loess")+
  
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM"), labels=c("3 μM")))+
  scale_colour_manual(values=c("grey20", "lightblue", "blue3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace






########
# Bri2 vs Bri2NNK extensions

file_name<-"Bri2NNK_ThT.txt"

tecan_results<-read_delim(file_name, locale=locale(encoding="latin1"), delim="\t")

################################################################################
tecan_melted<-melt(select(tecan_results, !c("Time [s]", "Temp. [ C]", "Time_hs")),
                   id="Time")
tecan_melted$value<-as.numeric(tecan_melted$value)

tecan_melted$plate_column<-substring(tecan_melted$variable, 2, 3)
tecan_melted$plate_row<-substring(tecan_melted$variable, 1, 1)

tecan_melted[tecan_melted$plate_column == 1, "conc"]<-"blank"

tecan_melted[tecan_melted$plate_column == 2 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 3 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 2 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 3 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 4 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 5 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 4 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 5 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 6 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 7 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 6 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 7 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 8 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 9 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 8 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 9 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"1.5uM"
tecan_melted[tecan_melted$plate_column == 10 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"12uM"
tecan_melted[tecan_melted$plate_column == 11 & tecan_melted$plate_row %in% c("B", "C", "D"), "conc"]<-"6uM"
tecan_melted[tecan_melted$plate_column == 10 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"3uM"
tecan_melted[tecan_melted$plate_column == 11 & tecan_melted$plate_row %in% c("E", "F", "G"), "conc"]<-"1.5uM"

tecan_melted$peptide<-"blank"
tecan_melted[tecan_melted$plate_column == 2 & tecan_melted$plate_row %in% c("B", "C", "D"), "peptide"]<-"EASNCFAIRHFENKFAVETLICS (Bri2)"
tecan_melted[tecan_melted$plate_column == 3 & tecan_melted$plate_row %in% c("B", "C", "D"), "peptide"]<-"EASNCFAIRHFENKFAVETLICS (Bri2)"
tecan_melted[tecan_melted$plate_column == 2 & tecan_melted$plate_row %in% c("E", "F", "G"), "peptide"]<-"EASNCFAIRHFENKFAVETLICS (Bri2)"
tecan_melted[tecan_melted$plate_column == 3 & tecan_melted$plate_row %in% c("E", "F", "G"), "peptide"]<-"EASNCFAIRHFENKFAVETLICS (Bri2)"
tecan_melted[tecan_melted$plate_column == 4 & tecan_melted$plate_row %in% c("B", "C", "D"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG"
tecan_melted[tecan_melted$plate_column == 5 & tecan_melted$plate_row %in% c("B", "C", "D"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG"
tecan_melted[tecan_melted$plate_column == 4 & tecan_melted$plate_row %in% c("E", "F", "G"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG"
tecan_melted[tecan_melted$plate_column == 5 & tecan_melted$plate_row %in% c("E", "F", "G"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG"
tecan_melted[tecan_melted$plate_column == 8 & tecan_melted$plate_row %in% c("B", "C", "D"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSIV"
tecan_melted[tecan_melted$plate_column == 9 & tecan_melted$plate_row %in% c("B", "C", "D"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSIV"
tecan_melted[tecan_melted$plate_column == 8 & tecan_melted$plate_row %in% c("E", "F", "G"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSIV"
tecan_melted[tecan_melted$plate_column == 9 & tecan_melted$plate_row %in% c("E", "F", "G"), "peptide"]<-"EASNCFAIRHFENKFAVETLICSIV"

tecan_melted$sample<-paste0(tecan_melted$peptide, "-", tecan_melted$conc)
tecan_melted$value <- replace_na(tecan_melted$value, max(tecan_melted$value, na.rm = T))


# Substract the blank and normalize
tecan_melted$value_correct <- tecan_melted$value-min(tecan_melted[tecan_melted$peptide == "blank", "value"])

tecan_melted_grouped <- tecan_melted %>% group_by(Time, sample) %>% mutate(mean_value=mean(value_correct), std_value=sd(value_correct))

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(conc) %>% mutate(value_norm=(value_correct/max(mean_value))*100)

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))

#
df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("6uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]

# Raw
scatter_tht<-ggplot(df_to_plot, 
                             aes(y=mean_value, x=Time, 
                                 color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("6uM"), labels=c("6uM")))+
  #facet_wrap(~factor(conc, levels=c("6uM","3uM"), labels=c("6uM","3uM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="ThT Fluorescence", color="")
scatter_tht

ggsave(scatter_tht, file="20241209_Bri2NNK_ThT_std.jpg", width=6, height=2, path=path)

# Normalized to max of each conc
scatter_tht_normalized<-ggplot(df_to_plot, 
                             aes(y=mean_value_norm, x=Time, 
                                 color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_norm-std_value_norm, ymax=mean_value_norm+std_value_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("6uM"), labels=c("6uM")))+
  #facet_wrap(~factor(conc, levels=c("6uM","3uM"), labels=c("6uM","3uM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized

ggsave(scatter_tht_normalized, file="20241209_Bri2NNK_ThT_std_norm.jpg", width=6, height=2, path=path)

# Normalized to each trace

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_correct_trace=value_correct-min(mean_value))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace=mean(value_correct_trace), std_value_trace=sd(value_correct_trace))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_trace_norm=(value_correct_trace/max(mean_value_trace))*100)
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("6uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]

# 
scatter_tht_normalized_trace<-ggplot(tecan_melted_grouped[tecan_melted_grouped$conc %in% c("6uM"),], 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("6uM"), labels=c("6uM")))+
  #facet_wrap(~factor(conc, levels=c("6uM","3uM"), labels=c("6uM","3uM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

ggsave(scatter_tht_normalized_trace, file="20241209_Bri2NNK_std_norm_trace.jpg", width=6, height=2, path=path)

#

scatter_all <- ggarrange(scatter_tht, scatter_tht_normalized_trace, common.legend = T, 
                         ncol=2)
scatter_all
ggsave(scatter_all, file="20241209_Bri2NNK_ThT_all.jpg", width=6, height=3, path=path)

###
# Use loess fit to predict t1/2
# Check loess fitting
scatter_tht_normalized_trace<-ggplot(df_to_plot[df_to_plot$Time<5,], 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_smooth(method="loess")+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("6uM"), labels=c("6 μM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace


loess_model_Bri2NNK1 <- loess(Time ~ mean_value_trace_norm, data = tecan_melted_grouped[tecan_melted_grouped$Time<5 & tecan_melted_grouped$sample == "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG-6uM",])
t1_2_Bri2NNK1<-predict(loess_model_Bri2NNK1, 50)

loess_model_Bri2NNK2 <- loess(Time ~ mean_value_trace_norm, data = tecan_melted_grouped[tecan_melted_grouped$Time<5 & tecan_melted_grouped$sample == "EASNCFAIRHFENKFAVETLICSIV-6uM",])
t1_2_Bri2NNK2<-predict(loess_model_Bri2NNK2, 50)

# plot t1/2
scatter_tht_normalized_trace<-ggplot(df_to_plot[df_to_plot$Time<5,], 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_vline(xintercept=t1_2_Bri2NNK1, color="darkred")+
  geom_hline(yintercept=50, color="darkred")+
  geom_vline(xintercept=t1_2_Bri2NNK2, color="red3")+
  geom_hline(yintercept=50, color="red3")+
  geom_point()+
  geom_smooth(method="loess")+
  
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("6uM"), labels=c("6 μM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

#

# Fig. 4d.
df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("6uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]

# Raw
scatter_tht<-ggplot(df_to_plot, 
                    aes(y=mean_value, x=Time, 
                        color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2,
                alpha=0.05)+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_classic()+
  labs(x="Time (hs)", y="ThT Fluorescence", color="")+
  theme(legend.position = "none")
scatter_tht

ggsave(scatter_tht, file="20241209_Bri2NNK_ThT_std_fig4d.jpg", width=2, height=2, path=path)

# Normalized to max of each conc
scatter_tht_normalized<-ggplot(df_to_plot, 
                               aes(y=mean_value_norm, x=Time, 
                                   color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_norm-std_value_norm, ymax=mean_value_norm+std_value_norm), width=.2,
                alpha=0.05)+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_classic()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")+
  theme(legend.position = "none")
scatter_tht_normalized

ggsave(scatter_tht_normalized, file="20241209_Bri2NNK_ThT_std_norm_fig4d.jpg", width=2, height=2, path=path)

#
both_scatters<-ggarrange(scatter_tht, scatter_tht_normalized, nrow=2)
both_scatters
ggsave(both_scatters, file="20241209_Bri2NNK_ThT_fig4d.jpg", width=2.5, height=4, path=path)


#### 3 uM

#
df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("3uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]

# Raw
scatter_tht<-ggplot(df_to_plot, 
                    aes(y=mean_value, x=Time, 
                        color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM"), labels=c("3uM")))+
  #facet_wrap(~factor(conc, levels=c("6uM","3uM"), labels=c("6uM","3uM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="ThT Fluorescence", color="")
scatter_tht

ggsave(scatter_tht, file="20241209_Bri2NNK_ThT_std_3uM.jpg", width=6, height=2, path=path)

# Normalized to max of each conc
scatter_tht_normalized<-ggplot(df_to_plot, 
                               aes(y=mean_value_norm, x=Time, 
                                   color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_norm-std_value_norm, ymax=mean_value_norm+std_value_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM"), labels=c("3uM")))+
  #facet_wrap(~factor(conc, levels=c("6uM","3uM"), labels=c("6uM","3uM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized

ggsave(scatter_tht_normalized, file="20241209_Bri2NNK_ThT_std_norm_3uM.jpg", width=6, height=2, path=path)

# Normalized to each trace

tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_correct_trace=value_correct-min(mean_value))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace=mean(value_correct_trace), std_value_trace=sd(value_correct_trace))
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(sample) %>% mutate(value_trace_norm=(value_correct_trace/max(mean_value_trace))*100)
tecan_melted_grouped <- tecan_melted_grouped %>% group_by(Time, sample) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("3uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]

# 
scatter_tht_normalized_trace<-ggplot(tecan_melted_grouped[tecan_melted_grouped$conc %in% c("3uM"),],
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

ggsave(scatter_tht_normalized_trace, file="20241209_Bri2NNK_std_norm_trace_3uM.jpg", width=6, height=2, path=path)

#

scatter_all <- ggarrange(scatter_tht, scatter_tht_normalized_trace, common.legend = T, 
                         ncol=2)
scatter_all
ggsave(scatter_all, file="20241209_Bri2NNK_ThT_all_3uM.jpg", width=9, height=3, path=path)

###
# Use loess fit to predict t1/2
# Check loess fitting
scatter_tht_normalized_trace<-ggplot(df_to_plot[df_to_plot$Time<5,], 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_smooth(method="loess")+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace


loess_model_Bri2NNK1 <- loess(Time ~ mean_value_trace_norm, data = tecan_melted_grouped[tecan_melted_grouped$Time<5 & tecan_melted_grouped$sample == "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG-3uM",])
t1_2_Bri2NNK1<-predict(loess_model_Bri2NNK1, 50)

loess_model_Bri2NNK2 <- loess(Time ~ mean_value_trace_norm, data = tecan_melted_grouped[tecan_melted_grouped$Time<5 & tecan_melted_grouped$sample == "EASNCFAIRHFENKFAVETLICSIV-3uM",])
t1_2_Bri2NNK2<-predict(loess_model_Bri2NNK2, 50)

# plot t1/2
scatter_tht_normalized_trace<-ggplot(df_to_plot[df_to_plot$Time<5,], 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_vline(xintercept=t1_2_Bri2NNK1, color="darkred")+
  geom_hline(yintercept=50, color="darkred")+
  geom_vline(xintercept=t1_2_Bri2NNK2, color="red3")+
  geom_hline(yintercept=50, color="red3")+
  geom_point()+
  geom_smooth(method="loess")+
  
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  facet_wrap(~factor(conc, levels=c("3uM")))+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_bw()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

#

# Fig. 4d.
df_to_plot<-tecan_melted_grouped[tecan_melted_grouped$conc %in% c("3uM"),]
df_to_plot<-df_to_plot[df_to_plot$Time < 8,]

# Raw
scatter_tht<-ggplot(df_to_plot, 
                    aes(y=mean_value, x=Time, 
                        color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2,
                alpha=0.05)+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  scale_y_continuous(breaks=c(0, 20000, 40000, 60000))+
  ylim(0, 60000)+
  theme_classic()+
  labs(x="Time (hs)", y="ThT Fluorescence", color="")+
  theme(legend.position = "none")
scatter_tht

ggsave(scatter_tht, file="20241209_Bri2NNK_ThT_std_fig4d_3uM.jpg", width=2, height=2, path=path)

# Normalized to max of each conc
scatter_tht_normalized<-ggplot(df_to_plot, 
                               aes(y=mean_value_norm, x=Time, 
                                   color=factor(peptide, levels=c("EASNCFAIRHFENKFAVETLICS (Bri2)", "EASNCFAIRHFENKFAVETLICSQLIMIYEDRKG", "EASNCFAIRHFENKFAVETLICSIV"))))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_norm-std_value_norm, ymax=mean_value_norm+std_value_norm), width=.2,
                alpha=0.05)+
  scale_colour_manual(values=c("grey80", "darkred", "red3"))+
  theme_classic()+
  labs(x="Time (hs)", y="Normalized\n ThT Fluorescence", color="")+
  theme(legend.position = "none")
scatter_tht_normalized

ggsave(scatter_tht_normalized, file="20241209_Bri2NNK_ThT_std_norm_fig4d_3uM.jpg", width=2, height=2, path=path)

#
both_scatters<-ggarrange(scatter_tht, scatter_tht_normalized, nrow=2)
both_scatters
ggsave(both_scatters, file="20241209_Bri2NNK_ThT_fig4d_3uM.jpg", width=2.5, height=4, path=path)


# Save as xlsx
tecan_melted_grouped<-tecan_melted_grouped[tecan_melted_grouped$conc != "12uM",]
write_xlsx(tecan_melted_grouped,"Bri2NNK_ThT.xlsx")
