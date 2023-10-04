library(tidyverse)
library(ggbreak) 
library(ggpubr)


#
df <- data.frame(sequence=c("SupN","AB42","Bri2", "ABri", "ADan" ,"ADan-F26S"),
                 rep1=c(0.006, 9.89, 0, 0.0049, 21.22, 22.8),
                 rep2=c(0.002, 10.53, 0.000076, 0.0110, 27.42, 29.21),
                 rep3=c(0.003, 9.67, 0, 0.0037, 26.79, 32.37))


df$mean_growth<-rowMeans(df[,c("rep1", "rep2", "rep3")], na.rm=T)
df$std_growth<-apply(df[,c("rep1", "rep2", "rep3")], 1, sd, na.rm=TRUE)

p1 <- ggplot(df[df$sequence %in% c("Bri2", "ABri", "ADan"),], 
             aes(x=factor(sequence, levels=c("Bri2", "ABri", "ADan")), y=mean_growth)) + 
  geom_col(color='black', fill='lightgrey') +
  geom_errorbar(aes(ymin=mean_growth-std_growth, ymax=mean_growth+std_growth), width=.2,
                position=position_dodge(.9))+
  scale_y_break(c(0.02, 10), scales=1.5)+
  labs(x='', y=' -URA-Ade / -URA ')+
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
p1


ggsave(p1, file="Bri2_ABri_ADan_growth.jpg", width=3, height=4)

###

df_1 <- data.frame(sequence=c("SupN", "SupN", "SupN", "AB42", "AB42", "AB42", 
                              "Bri2", "Bri2", "Bri2", "ABri", "ABri", "ABri",
                              "ADan", "ADan", "ADan"),
                   growth=c(0.006, 0.002, 0.003, 9.89, 10.53, 9.67,
                            0, 0.000076, 0, 0.0049, 0.0110, 0.0037,
                            21.22, 27.42, 26.79))

p_2 <- ggplot(df_1, aes(factor(sequence, levels=c("SupN", "AB42", "Bri2", "ABri", "ADan")), growth)) +
       geom_bar(stat="identity", position="dodge", width=.5, color=NaN, fill="lightgrey") +
       geom_signif(comparisons = list(c("Bri2", "ABri"), 
                                      c("SupN", "AB42"),
                                      c("SupN", "Bri2"), 
                                      c("SupN", "ABri"),
                                      c("SupN", "ADan"),
                                      c("AB42", "ADan")), 
                   test="t.test",
                   step_increase = .2,
                   map_signif_level = TRUE, tip_length = 0)+
    ylim(0, 60)+  
    scale_y_break(c(0.02, 10), scales=1)+
    labs(x='', y=' -URA-Ade / -URA ')+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 0),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12))
p_2

ggsave(p_2, file="SupN_AB42_Bri2_ABri_ADan_growth.jpg", width=5, height=4)






