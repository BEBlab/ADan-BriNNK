library(tidyverse)
library(ggpubr)


dir.create("01_Centering_Correction")
path="01_Centering_Correction"

name<-'ABri'
load(paste0(name, '_fitness_replicates.RData'))

all_variants<-rename(all_variants, nscore = fitness, nscore1 = fitness1_uncorr, nscore2 = fitness2_uncorr,
                     nscore3 = fitness3_uncorr)
                   
all_variants<-select(all_variants, aa_seq, starts_with('nscore') | starts_with('count'))
synonymous<-rename(synonymous, nscore = fitness, nscore1 = fitness1_uncorr, nscore2 = fitness2_uncorr, 
                   nscore3 = fitness3_uncorr)
                   

singles<-select(singles, !starts_with('fitness'))

silent<-synonymous
singles<-inner_join(singles, all_variants, by='aa_seq')

singles$Pos<-singles$Pos

#centering to the weighted mean of synonymous with 1 mut codon
mean_syn_1codon<-weighted.mean(silent[silent$Nmut_codons==1,]$nscore, silent[silent$Nmut_codons==1,]$sigma^-2, na.rm = T)

###silent
silent$nscore_c<-as.numeric(paste(silent$nscore-mean_syn_1codon))
silent$ID<-"silent"
silent$Mut<-"silent"

###singles
singles$nscore_c<-as.numeric(paste(as.numeric(singles$nscore)-mean_syn_1codon))
singles$nscore1_c<-as.numeric(paste(as.numeric(singles$nscore1)-mean_syn_1codon))
singles$nscore2_c<-as.numeric(paste(as.numeric(singles$nscore2)-mean_syn_1codon))
singles$nscore3_c<-as.numeric(paste(as.numeric(singles$nscore3)-mean_syn_1codon))

singles$ID<-paste(singles$WT_AA, singles$Pos, singles$Mut, sep = "-")

# FDR=0.1 correction and assignment into categories

singles$zscore<-singles$nscore_c/singles$sigma
singles$p.adjust<-p.adjust(2*pnorm(-abs(singles$zscore)), method = "BH")

singles$sig_10<-FALSE
singles[singles$p.adjust<0.1,]$sig_10<-TRUE

singles$category_10<-"WT-like"
singles[singles$sig_10==T & singles$nscore_c<0,]$category_10<-"NS_dec"
singles[singles$sig_10==T & singles$nscore_c>0,]$category_10<-"NS_inc"

singles$Bri2<-"Bri2-extension"
singles[singles$STOP == T & singles$Pos<23,]$Bri2<-"Bri2-truncation"

#with stops
singles_stops<-singles

#remove stops
singles<-singles[singles$Mut!="*",]

#mean input and ouput reads counts
singles_stops <- singles_stops %>%
  rowwise() %>%
  mutate(input_mean_count=mean(c(count_e1_s0, count_e2_s0, count_e3_s0)),
         output_mean_count=mean(c(count_e1_s1, count_e2_s1, count_e3_s1)))

#################################################################################
# interquartile range == fitness range

summary(singles_stops$nscore_c)
iqr<-IQR(singles_stops$nscore_c)

first_to_wt<-summary(singles_stops$nscore_c)[[2]]

p_iqr<-ggplot(singles_stops, aes(x=nscore_c))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  labs(x="Nucleation score", title="Fitness range")+
  
  geom_vline(xintercept=0)+
  annotate("text", label="WT fitness= 0", x=0.1, y=110, color="black")+
  
  geom_vline(xintercept=summary(singles_stops$nscore_c)[[2]], color="red")+
  annotate("text", label=paste0("1st Qu= ", round(summary(singles_stops$nscore_c)[[2]], 2)), 
           x=summary(singles_stops$nscore_c)[[2]]-0.8, y=140, color="red")+
  
  geom_vline(xintercept=summary(singles_stops$nscore_c)[[5]], color="red")+
  annotate("text", label=paste0("3rd Qu= ", round(summary(singles_stops$nscore_c)[[5]], 2)), 
           x=summary(singles_stops$nscore_c)[[5]]+0.8, y=140, color="red")+
  
  geom_vline(xintercept=summary(singles_stops$nscore_c)[[3]], color="blue")+
  annotate("text", label=paste0("median= ",round(summary(singles_stops$nscore_c)[[3]], 2) ), 
           x=summary(singles_stops$nscore_c)[[3]]+0.1, y=100, color="blue")
p_iqr

ggsave(p_iqr, file="p_iqr.jpg", path=path, width = 5, height = 3)


#sigma distribution
p_sigma_dist<-ggplot(singles_stops, aes(x=sigma))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  scale_x_continuous(limits = c(0,3))
p_sigma_dist

ggsave(p_sigma_dist, file="p_sigma_dist.jpg", path=path, width = 5, height = 3)

# normalise sigmas to fitness range - if fitness range is IQR

singles_stops$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(singles_stops$nscore_c))

# or if fitness range is 1rst to WT fitness

singles_stops$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(singles_stops$nscore_c)[[2]])

singles_stops$sigma_norm_iqr<-singles_stops$sigma / fitness_range_iqr
singles_stops$sigma_norm_first_toWT<-singles_stops$sigma / fitness_range_first_toWT

singles_stops$sigma_norm_iqr<-as.numeric(singles_stops$sigma_norm_iqr)
singles_stops$sigma_norm_first_toWT<-as.numeric(singles_stops$sigma_norm_first_toWT)

#sigma_norm_iqr

p1<-ggplot(singles_stops, aes(x=sigma_norm_iqr))+
  geom_histogram(bins=100, aes(fill=factor(category_10, levels=c("NS_inc", "NS_dec", "WT-like", "Bri2_stops22"))))+
  scale_fill_manual(values=c("#7979BE", "#DF9292", "#F2F2F2", "black"))+
  facet_wrap(~category_10, ncol=1)+
  labs(title="sigma normalised to fitness range (1st Qu to 3r Qu)", fill="Category_FDR")+
  scale_x_continuous(limits = c(0,1))+
  geom_vline(xintercept = 0.2)+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))
p1

ggsave(p1, file="p_sigma_normalised_1.jpg", path=path, width = 8, height = 5)


## exclude those that are over 20% of the fitness range
singles_stops$low_sigma<-FALSE
singles_stops[singles_stops$sigma_norm_iqr<=0.2,]$low_sigma<-TRUE

# Categories based on sigma normalized to fitness range
singles_stops$category_10_sigma<-"unclassified"

singles_stops[singles_stops$category_10=="NS_inc"& singles_stops$low_sigma==T,]$category_10_sigma<-"NS_inc"
singles_stops[singles_stops$category_10=="NS_dec" & singles_stops$low_sigma==T,]$category_10_sigma<-"NS_dec"
singles_stops[singles_stops$category_10=="WT-like" & singles_stops$low_sigma==T,]$category_10_sigma<-"WT-like"

#save(singles_stops, file="singles_stops.RData")

# Summary of singles variants
print(paste('N variants in input:', nrow(singles_stops)))
print(paste('N in output:', nrow(singles_stops[singles_stops$input_mean_count>1 & singles_stops$output_mean_count>=1,])))
print(paste('N FDR 10%:', nrow(singles_stops[singles_stops$sig_10==T,])))
print(paste('N WT-like (real):', nrow(singles_stops[singles_stops$category_10_sigma=="WT-like",])))
print(paste('N unclassified:', nrow(singles_stops[singles_stops$category_10_sigma=="unclassified",])))
print(paste('N deads:', nrow(singles_stops[singles_stops$input_mean_count>1 & singles_stops$output_mean_count<1,])))


save(silent, singles, singles_stops, file="nscore_df_ABri.RData")

#load("nscore_df_ABri.RData")

library("writexl")
write_xlsx(singles_stops, "ABri_singles.xlsx")
write_xlsx(silent, "ABri_synonymous.xlsx")

write.csv(singles_stops, "ABri_singles.csv", row.names=FALSE)
write.csv(silent, "ABri_synonymous.csv", row.names=FALSE)
