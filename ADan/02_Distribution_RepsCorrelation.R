library(tidyverse)
library(ggpubr)
library(ggrepel)

dir.create("02_Distribution_RepsCorrelation")
path="02_Distribution_RepsCorrelation"

load("nscore_df_ADan.RData")

#################################################################################
# distribution singles (Missense, Nonsense and Synonymous)

dist_singles<-singles_stops[,c("Mut", "nscore_c")]
dist_singles<-rbind(dist_singles, silent[silent$Nmut_codons==1,c("Mut", "nscore_c")])
dist_singles$type<-"Missense"
dist_singles[dist_singles$Mut=="*",]$type<-"Nonsense"
dist_singles[dist_singles$Mut=="silent",]$type<-"Synonymous"

p_hist<-ggplot(dist_singles, aes(x=factor(type, levels=c("Synonymous", "Nonsense", "Missense")), y=nscore_c))+
  geom_jitter(color="grey20", height = 0, width = 0.1, size=0.5)+
  geom_violin(fill="grey90", alpha=0.5, scale="width")+
  geom_boxplot(width = 0.02, outlier.shape=NA)+
  labs(x="", y="Nucleation Score", fill="")+
  theme_bw()+
  theme(legend.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        strip.text = element_text(size=14),
        plot.margin = margin(0,0,0,0, unit = 'cm'),
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(name="Nucleation Score", breaks=c(-8, -4, 0, 2), labels=c(-8, -4, 0, 2), limits=c(-8, 2))

p_hist

ggsave(p_hist, file="p_hist.jpg", path=path, width = 5, height = 4)

################################################################################
# correlation between replicates

#1 and 2
subset<-singles_stops[!is.na(singles_stops$nscore1_c),]
subset<-subset[!is.na(subset$nscore2_c),]
n<-length(subset$ID)

corr<-cor.test(singles_stops$nscore1_c, singles_stops$nscore2_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr_12<-ggplot(singles_stops, aes(x=nscore1_c, y=nscore2_c) )+
  stat_binhex()+
  theme_bw()+
  annotate("text", x = 1, y = -7.2, label = paste0("R=", round(R, 2)), size=5)+
  annotate("text", x = 1, y = -8, label = paste0("p=",format(p, digits = 2, scientific = T)), size=5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_fill_gradient(high="grey30", low="grey90")+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")+
  scale_x_continuous(name="Replicate 1", breaks=c(-8, -4, 0, 2), labels=c(-8, -4, 0, 2), limits=c(-8, 2))+
  scale_y_continuous(name="Replicate 2", breaks=c(-8, -4, 0, 2), labels=c(-8, -4, 0, 2), limits=c(-8, 2))

p_corr_12

ggsave(p_corr_12,path=path, file="p_corr_rep1vs2.jpg", width=5, height=4)

#1 and 3
subset<-singles_stops[!is.na(singles_stops$nscore1_c),]
subset<-subset[!is.na(subset$nscore3_c),]
n<-length(subset$ID)

corr<-cor.test(singles_stops$nscore1_c, singles_stops$nscore3_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr_13<-ggplot(singles_stops, aes(x=nscore1_c, y=nscore3_c) )+
  stat_binhex()+
  theme_bw()+
  annotate("text", x = 1, y = -7.2, label = paste0("R=", round(R, 2)), size=5)+
  annotate("text", x = 1, y = -8, label = paste0("p=",format(p, digits = 2, scientific = T)), size=5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_fill_gradient(high="grey30", low="grey90")+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")+
  scale_x_continuous(name="Replicate 1", breaks=c(-8, -4, 0, 4), labels=c(-8, -4, 0, 4), limits=c(-8, 4))+
  scale_y_continuous(name="Replicate 3", breaks=c(-8, -4, 0, 4), labels=c(-8, -4, 0, 4), limits=c(-8, 4))

p_corr_13

ggsave(p_corr_13, path=path, file="p_corr_reps1vs3.jpg", width=5, height=4)

#2 and 3
subset<-singles_stops[!is.na(singles_stops$nscore2_c),]
subset<-subset[!is.na(subset$nscore3_c),]
n<-length(subset$ID)

corr<-cor.test(singles_stops$nscore2_c, singles_stops$nscore3_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr_23<-ggplot(singles_stops, aes(x=nscore2_c, y=nscore3_c) )+
  stat_binhex()+
  theme_bw()+
  annotate("text", x = 1, y = -7.2, label = paste0("R=", round(R, 2)), size=5)+
  annotate("text", x = 1, y = -8, label = paste0("p=",format(p, digits = 2, scientific = T)), size=5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_fill_gradient(high="grey30", low="grey90")+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")+
  scale_x_continuous(name="Replicate 2", breaks=c(-8, -4, 0, 4), labels=c(-8, -4, 0, 4), limits=c(-8, 4))+
  scale_y_continuous(name="Replicate 3", breaks=c(-8, -4, 0, 4), labels=c(-8, -4, 0, 4), limits=c(-8, 4))

p_corr_23

ggsave(p_corr_23,path=path, file="p_corr_reps2vs3.jpg", width=5, height=4)

#
p_corr<-ggarrange(p_corr_12, p_corr_13, p_corr_23, ncol=3, common.legend = TRUE)
p_corr

ggsave(p_corr, path=path, file="p_corr.jpg", width=12, height=4)

### Individual Validation

individual_validation<-read_tsv("ADan_IndividualValidation.tsv")

individual_validation<-individual_validation %>% group_by(ID) %>% summarise(mean_growth=mean(growth_rate), 
                                                                            std_growth=sd(growth_rate))

# ADan individual validation: mean growth rate = 21.1, std = 3.05

individual_validation<-left_join(individual_validation, singles_stops[c("ID", "nscore_c", "sigma")], by="ID")

ADan<-data.frame("ID"="ADan", "mean_growth" = 21.1, "std_growth" = 3.05, "nscore_c"=0.011554036, "sigma"=0.05481394 )

individual_validation<-rbind(individual_validation, ADan)

corr<-cor.test(individual_validation$nscore_c, individual_validation$mean_growth, use="complete.obs")
R<-corr$estimate
p<-corr$p.value


p_small_large_scale<-ggplot(individual_validation, aes(x=mean_growth, y=nscore_c))+
  geom_smooth(method = "lm", se=F, color="lightgrey", linetype="dashed")+
  geom_pointrange(aes(ymin=nscore_c-sigma, ymax=nscore_c+sigma), color="darkgrey")+
  geom_pointrange(aes(xmin=mean_growth-std_growth, xmax=mean_growth+std_growth), color="darkgrey")+
  geom_point(size=2)+
  annotate("text", x = 25, y = -7.2, label=paste0("R=",round(R, 2)), size=5)+
  annotate("text", x = 25, y = -8, label =paste0("p=",format(p, digits = 2, scientific = T)), size=5)+
  geom_label_repel(aes(label=ID), seed=42, box.padding = 0.5)+
  scale_y_continuous(breaks=c(-8, -4, 0, 2), labels=c(-8, -4, 0, 2), limits=c(-8, 2))+
  theme_bw()+
  labs(x="Small scale" ,y="Large scale")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text = element_text(size=14))

p_small_large_scale

ggsave(p_small_large_scale, path=path, file="p_individualvalidation.jpg", width=5, height=4)
