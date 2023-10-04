library(tidyverse)
require(reshape2)
library(ggpubr)

dir.create("02_Truncations")
path="02_Truncations"

#required data:
load("nscore_df_ADan.RData")

#
peptide_seq<-"EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY"
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], '\n', n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "*")


#
stops<-singles_stops[singles_stops$STOP == TRUE,]

#centering to the weighted mean of synonymous with 1 mut codon
stops$length<-stops$Pos-1
mean_bri2<-weighted.mean(stops[stops$length<23, "nscore"], stops[stops$length<23, "sigma"]^-2, na.rm = T)

###stops
stops$nscore_c<-as.numeric(paste(as.numeric(stops$nscore)-mean_bri2))
stops$nscore1_c<-as.numeric(paste(as.numeric(stops$nscore1)-mean_bri2))
stops$nscore2_c<-as.numeric(paste(as.numeric(stops$nscore2)-mean_bri2))
stops$nscore3_c<-as.numeric(paste(as.numeric(stops$nscore3)-mean_bri2))

silent$nscore_c<-as.numeric(paste(as.numeric(silent$nscore)-mean_bri2))
silent$length<-34

stops<-rbind(silent[silent$Nham_nt == 0, c("length", "nscore_c", "sigma")], stops[c("length", "nscore_c", "sigma")])

# FDR=0.1 correction and assignment into categories

stops$zscore<-stops$nscore_c/stops$sigma
stops$p.adjust<-p.adjust(2*pnorm(-abs(stops$zscore)), method = "BH")

stops$sig_1<-FALSE
stops[stops$p.adjust<0.01,]$sig_1<-TRUE

stops$category_1<-"Bri2-like"
#stops[stops$sig_1==T & stops$nscore_c<0,]$category_1<-"NS_dec"
stops[stops$sig_1==T & stops$nscore_c>0,]$category_1<-"NS_inc"

#

mean(as.numeric(stops[stops$length<23,]$nscore_c))

min<-min(stops$nscore_c)
max<-max(stops$nscore_c)
cols <- c(colorRampPalette(c( "brown3", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(-min+max)*100)-0.5))


#### single AA deletions

p_truncations<-ggplot(stops, aes(x=factor(length, levels=c(0:34)), y=nscore_c ))+
  geom_hline(yintercept = 0, size=0.5, linetype="dashed", color="grey70")+
  geom_errorbar(aes(ymin=nscore_c-1.96*sigma, ymax=nscore_c+1.96*sigma), width=0, size=0.1)+
  geom_point(aes(fill=nscore_c), size=4, shape=21)+
  scale_x_discrete(breaks=seq(0, 34))+
  labs(y="Nucleation score", x="Peptide length", fill="Nucleation score")+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color="black", size=0.25),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=14, color=c(rep("grey50", 23), rep("grey10", 12))),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), breaks=c(-2, 0, 2))
p_truncations

ggsave(p_truncations, file="p_truncations_effect.jpg", width = 14, height = 4, path=path)

