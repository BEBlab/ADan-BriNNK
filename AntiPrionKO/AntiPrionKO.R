library(tidyverse)
library(ggpubr)
library(reshape2)
library(DescTools)
library(ggbreak) 

dir.create("AntiPrionKO")
path="AntiPrionKO"

file_name<-"Nucleation_AntiPrionKO.txt"

df<-read_delim(file_name, locale=locale(encoding="latin1"), delim="\t")
df<-remove_missing(df)

proteins_interest<-c("Bri2", "ABri", "ABri A7P", "ABri R24C", "ADan", "ADan K14V", "ADan L20R")
df_<-df[df$Protein %in% proteins_interest,]

# Anova - Dunetts' test
df_bri2<-df_[df_$Protein == "Bri2",]

oneway.test(NormalizedGR ~ Background,
            data = df_bri2,
            var.equal = T # assuming equal variances
)


DunnettTest(df_bri2$NormalizedGR, factor(df_bri2$Background, levels=c("WT",
                                                                      "Cur1 KO",
                                                                      "Btn2 KO",
                                                                      "Hsp104 KO", 
                                                                      "Ssz1 KO", 
                                                                      "Upf1 KO")), na.rm=T)
###
df_abri<-df_[df_$Protein == "ABri",]
oneway.test(NormalizedGR ~ Background,
            data = df_abri,
            var.equal = T # assuming equal variances
)

DunnettTest(df_abri$NormalizedGR, factor(df_abri$Background, levels=c("WT",
                                                                      "Cur1 KO",
                                                                      "Btn2 KO",
                                                                      "Hsp104 KO", 
                                                                      "Ssz1 KO", 
                                                                      "Upf1 KO")), na.rm=T)
###
df_A7P<-df_[df_$Protein == "ABri A7P",]
oneway.test(NormalizedGR ~ Background,
            data = df_A7P,
            var.equal = T # assuming equal variances
)

DunnettTest(df_A7P$NormalizedGR, factor(df_A7P$Background, levels=c("WT",
                                                                    "Cur1 KO",
                                                                    "Btn2 KO",
                                                                    "Hsp104 KO", 
                                                                    "Ssz1 KO", 
                                                                    "Upf1 KO")), na.rm=T)

###
df_R24C<-df_[df_$Protein == "ABri R24C",]
oneway.test(NormalizedGR ~ Background,
            data = df_R24C,
            var.equal = T # assuming equal variances
)

DunnettTest(df_R24C$NormalizedGR, factor(df_R24C$Background, levels=c("WT",
                                                                      "Cur1 KO",
                                                                      "Btn2 KO",
                                                                      "Hsp104 KO", 
                                                                      "Ssz1 KO", 
                                                                      "Upf1 KO")), na.rm=T)

###
df_adan<-df_[df_$Protein == "ADan",]
oneway.test(NormalizedGR ~ Background,
            data = df_adan,
            var.equal = T # assuming equal variances
)

DunnettTest(df_adan$NormalizedGR, factor(df_adan$Background, levels=c("WT",
                                                                      "Cur1 KO",
                                                                      "Btn2 KO",
                                                                      "Hsp104 KO", 
                                                                      "Ssz1 KO", 
                                                                      "Upf1 KO")), na.rm=T)

###
df_K14V<-df_[df_$Protein == "ADan K14V",]
oneway.test(NormalizedGR ~ Background,
            data = df_K14V,
            var.equal = T # assuming equal variances
)

DunnettTest(df_K14V$NormalizedGR, factor(df_K14V$Background, levels=c("WT",
                                                                      "Cur1 KO",
                                                                      "Btn2 KO",
                                                                      "Hsp104 KO", 
                                                                      "Ssz1 KO", 
                                                                      "Upf1 KO")), na.rm=T)

###
df_L20R<-df_[df_$Protein == "ADan L20R",]
oneway.test(NormalizedGR ~ Background,
            data = df_L20R,
            var.equal = T # assuming equal variances
)

DunnettTest(df_L20R$NormalizedGR, factor(df_L20R$Background, levels=c("WT",
                                                                      "Cur1 KO",
                                                                      "Btn2 KO",
                                                                      "Hsp104 KO", 
                                                                      "Ssz1 KO", 
                                                                      "Upf1 KO")), na.rm=T)


# Calculate mean and standard error for each group
df_summary <- df_ %>%
  group_by(Protein, Background) %>%
  summarise(
    mean_growth = mean(NormalizedGR),
    sd = sd(NormalizedGR)  # Calculate standard error
  )


# Create the bar plot with error bars using the summarized data (without the break in y axis)
p_4 <- ggplot(df_summary, aes(x=factor(Protein, levels=proteins_interest), 
                              y=mean_growth, 
                              fill=factor(Background, levels=c("WT", "Cur1 KO", "Btn2 KO", "Hsp104 KO", "Ssz1 KO", "Upf1 KO"))))+
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color="black") +
  geom_errorbar(aes(ymin = mean_growth, ymax = mean_growth + sd), width = 0.2, position = position_dodge(0.5)) +
  scale_fill_manual(values=c("grey10", "grey20", "grey40", "grey60", "grey80", "grey90"))+
  #ylim(0, 100)+
  labs(x = '', y = ' -URA-Ade / -URA ', fill="Yeast Genotype") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

p_4

ggsave(p_4, file="AntiPrionKO_grouped_nobreak.jpg", width=10, height=4, path=path)
ggsave(p_4, file="AntiPrionKO_grouped_nobreak.pdf", width=10, height=4, path=path)


##

# Growth Rate ratio between wt background and KOs

wt_growth <- df_summary[df_summary$Background == "WT", c("Protein", "mean_growth")]
wt_growth <- rename(wt_growth, wt_growth = mean_growth)

df_summary <- left_join(df_summary, wt_growth, by="Protein")

df_summary <- df_summary %>% group_by(Protein) %>% mutate(ratio=(mean_growth/wt_growth))

ratio_plot<-ggplot(df_summary, aes(x=factor(Protein, levels=proteins_interest), 
                                   y=ratio, 
                                   fill=factor(Background, levels=c("WT", "Cur1 KO", "Btn2 KO", "Hsp104 KO", "Ssz1 KO", "Upf1 KO"))))+
  geom_hline(yintercept=1, alpha=0.5, linetype="dashed")+
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color="black") +
  scale_fill_manual(values=c("grey10", "grey20", "grey40", "grey60", "grey80", "grey90"))+
  scale_y_continuous(breaks=c(2, 4, 6, 8, 10))+
  labs(x = '', y = 'Ratio (KO/WT)', fill="Yeast Genotype") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ratio_plot

ggsave(ratio_plot, file="AntiPrionKO_ratio.jpg", width=10, height=4, path=path)
ggsave(ratio_plot, file="AntiPrionKO_ratio.pdf", width=10, height=4, path=path)
