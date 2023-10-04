library(tidyverse)
library(growthcurver)
library(DescTools)



file_name<-"BB_growth (Modified)_20230929_090759.csv"

tecan_results<-read_csv(file_name, locale=locale(encoding="latin1"))

tecan_results$Time<-tecan_results$`Time [s]`/60/60


dir.create("20230929")
path="20230929"


bri2_1_nocu<-c("B2", "C2", "D2")
bri2_1_cu<-c("E2", "F2", "G2")

mean_bri2_1_nocu<-rowMeans(tecan_results[bri2_1_nocu])
mean_bri2_1_cu<-rowMeans(tecan_results[bri2_1_cu])


bri2_2_nocu<-c("B3", "C3", "D3")
bri2_2_cu<-c("E3", "F3", "G3")

mean_bri2_2_nocu<-rowMeans(tecan_results[bri2_2_nocu])
mean_bri2_2_cu<-rowMeans(tecan_results[bri2_2_cu])


bri2_3_nocu<-c("B4", "C4", "D4")
bri2_3_cu<-c("E4", "F4", "G4")

mean_bri2_3_nocu<-rowMeans(tecan_results[bri2_3_nocu])
mean_bri2_3_cu<-rowMeans(tecan_results[bri2_3_cu])


ABri_1_nocu<-c("B5", "C5", "D5")
ABri_1_cu<-c("E5", "F5", "G5")

mean_ABri_1_nocu<-rowMeans(tecan_results[ABri_1_nocu])
mean_ABri_1_cu<-rowMeans(tecan_results[ABri_1_cu])


ABri_2_nocu<-c("B6", "C6", "D6")
ABri_2_cu<-c("E6", "F6", "G6")

mean_ABri_2_nocu<-rowMeans(tecan_results[ABri_2_nocu])
mean_ABri_2_cu<-rowMeans(tecan_results[ABri_2_cu])


ABri_3_nocu<-c("B7", "C7", "D7")
ABri_3_cu<-c("E7", "F7", "G7")

mean_ABri_3_nocu<-rowMeans(tecan_results[ABri_3_nocu])
mean_ABri_3_cu<-rowMeans(tecan_results[ABri_3_cu])


ADan_1_nocu<-c("B8", "C8", "D8")
ADan_1_cu<-c("E8", "F8", "G8")

mean_ADan_1_nocu<-rowMeans(tecan_results[ADan_1_nocu])
mean_ADan_1_cu<-rowMeans(tecan_results[ADan_1_cu])


ADan_2_nocu<-c("B9", "C9", "D9")
ADan_2_cu<-c("E9", "F9", "G9")

mean_ADan_2_nocu<-rowMeans(tecan_results[ADan_2_nocu])
mean_ADan_2_cu<-rowMeans(tecan_results[ADan_2_cu])


ADan_3_nocu<-c("B10", "C10", "D10")
ADan_3_cu<-c("E10", "F10", "G10")

mean_ADan_3_nocu<-rowMeans(tecan_results[ADan_3_nocu])
mean_ADan_3_cu<-rowMeans(tecan_results[ADan_3_cu])


SupN_nocu<-c("B11", "C11", "D11")
SupN_cu<-c("E11", "F11", "G11")

mean_SupN_nocu<-rowMeans(tecan_results[SupN_nocu])
mean_SupN_cu<-rowMeans(tecan_results[SupN_cu])


mean_growth<-data.frame(tecan_results["Time"], 
                        mean_bri2_1_cu, mean_bri2_1_nocu,
                        mean_bri2_2_cu, mean_bri2_2_nocu,
                        mean_bri2_3_cu, mean_bri2_3_nocu,
                        mean_ABri_1_cu, mean_ABri_1_nocu,
                        mean_ABri_2_cu, mean_ABri_2_nocu,
                        mean_ABri_3_cu, mean_ABri_3_nocu,
                        mean_ADan_1_cu, mean_ADan_1_nocu,
                        mean_ADan_2_cu, mean_ADan_2_nocu,
                        mean_ADan_3_cu, mean_ADan_3_nocu,
                        mean_SupN_cu, mean_SupN_nocu)


gc_bri2_1_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_bri2_1_nocu)
plot(gc_bri2_1_nocu)
gc_bri2_1_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_bri2_1_cu)
plot(gc_bri2_1_cu)
gc_bri2_2_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_bri2_2_nocu)
plot(gc_bri2_2_nocu)
gc_bri2_2_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_bri2_2_cu)
plot(gc_bri2_2_cu)
gc_bri2_3_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_bri2_3_nocu)
plot(gc_bri2_3_nocu)
gc_bri2_3_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_bri2_3_cu)
plot(gc_bri2_3_cu)

gc_ABri_1_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ABri_1_nocu)
plot(gc_ABri_1_nocu)
gc_ABri_1_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ABri_1_cu)
plot(gc_ABri_1_cu)
gc_ABri_2_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ABri_2_nocu)
plot(gc_ABri_2_nocu)
gc_ABri_2_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ABri_2_cu)
plot(gc_ABri_2_cu)
gc_ABri_3_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ABri_3_nocu)
plot(gc_ABri_3_nocu)
gc_ABri_3_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ABri_3_cu)
plot(gc_ABri_3_cu)

gc_ADan_1_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ADan_1_nocu)
plot(gc_ADan_1_nocu)
gc_ADan_1_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ADan_1_cu)
plot(gc_ADan_1_cu)
gc_ADan_2_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ADan_2_nocu)
plot(gc_ADan_2_nocu)
gc_ADan_2_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ADan_2_cu)
plot(gc_ADan_2_cu)
gc_ADan_3_nocu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ADan_3_nocu)
plot(gc_ADan_3_nocu)
gc_ADan_3_cu <- SummarizeGrowth(mean_growth$Time, mean_growth$mean_ADan_3_cu)
plot(gc_ADan_3_cu)

bri2_nocu_r<-c(gc_bri2_1_nocu$vals$r, gc_bri2_2_nocu$vals$r, gc_bri2_3_nocu$vals$r)
bri2_cu_r<-c(gc_bri2_1_cu$vals$r, gc_bri2_2_cu$vals$r, gc_bri2_3_cu$vals$r)

ABri_nocu_r<-c(gc_ABri_1_nocu$vals$r, gc_ABri_2_nocu$vals$r, gc_ABri_3_nocu$vals$r)
ABri_cu_r<-c(gc_ABri_1_cu$vals$r, gc_ABri_2_cu$vals$r, gc_ABri_3_cu$vals$r)

ADan_nocu_r<-c(gc_ADan_1_nocu$vals$r, gc_ADan_2_nocu$vals$r, gc_ADan_3_nocu$vals$r)
ADan_cu_r<-c(gc_ADan_1_cu$vals$r, gc_ADan_2_cu$vals$r, gc_ADan_3_cu$vals$r)

growth_rates<-data.frame("sample"=c("Bri2", "Bri2", "Bri2", "ABri", "ABri", "ABri", "ADan", "ADan", "ADan",
                                    "Bri2", "Bri2", "Bri2", "ABri", "ABri", "ABri", "ADan", "ADan", "ADan"),
                         "treatment"=c("non-induced", "non-induced", "non-induced", "non-induced", "non-induced", "non-induced", "non-induced", "non-induced", "non-induced",
                                       "induced", "induced", "induced", "induced", "induced", "induced", "induced", "induced", "induced"),
                         "growth_rates"=c(bri2_nocu_r, ABri_nocu_r, ADan_nocu_r,
                                          bri2_cu_r, ABri_cu_r, ADan_cu_r))

growth_rates_plot<-ggplot(growth_rates, aes(x=factor(sample, levels=c("Bri2", "ABri", "ADan")), y=growth_rates))+
  #geom_boxplot(position=position_dodge(0.8), aes(fill=factor(treatment, levels=c("non-induced", "induced"))))+
  geom_point(size=2, position=position_dodge(0.3), aes(color=factor(treatment, levels=c("non-induced", "induced"))))+
  labs(y="Growth Rate", x="", fill="", color="")+
  scale_color_manual(values=c("grey70", "black"))+
  ylim(0, 0.7)+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14),
        plot.margin = margin(0,0,0,0, unit = 'cm'),
        plot.title = element_text(hjust = 0.5))

growth_rates_plot

ggsave(growth_rates_plot, file="Bri2_ABri_ADan_toxicity.jpg", width=6, height=4)



growth_rates$ID_condition<-paste0(growth_rates$sample, "-", growth_rates$treatment)

oneway.test(growth_rates ~ ID_condition,
            data = growth_rates[growth_rates$treatment == "induced",],
            var.equal = TRUE # assuming equal variances
)

DunnettTest(growth_rates$growth_rates, factor(growth_rates$ID_condition, levels=c("Bri2-induced",
                                                                                  "ABri-induced",
                                                                                  "ADan-induced")), na.rm=T)



oneway.test(growth_rates ~ ID_condition,
            data = growth_rates[growth_rates$treatment == "non-induced",],
            var.equal = TRUE # assuming equal variances
)

DunnettTest(growth_rates$growth_rates, factor(growth_rates$ID_condition, levels=c("Bri2-non-induced",
                                                                                  "ABri-non-induced",
                                                                                  "ADan-non-induced")), na.rm=T)

