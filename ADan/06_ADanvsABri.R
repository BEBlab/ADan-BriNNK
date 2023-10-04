library(tidyverse)
library(ggpubr)
require("ggrepel")
library(stringr)


dir.create("06_ADanvsABri")
path="06_ADanvsABri"

#required data:


all_aa<-c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H","P")

AA_type<-data.frame("wt_aa" = all_aa,
                    "name_AA"=c("Glycine", "Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                                "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                                "Cysteine","Asparagine","Glutamine","Histidine","Proline"),
                    "type"=c("glycine",rep("aliphatic",5),rep("aromatic",3),rep("positive",2),rep("negative",2),rep("polar",6),"proline"))

color_AAtype<-c("aliphatic"="darkgrey",
                "aromatic"="#9A703EFF",
                "negative"="#EE0011FF",
                "positive"="#0C5BB0FF", 
                "polar"="#15983DFF", 
                "glycine"="grey30", 
                "proline"="#FEC10BFF")



###
# Load ADan data

name<-"ADan"

load(paste0("nscore_df_", name, ".RData"))

# Let's work only with those trustable sequences
singles_stops<-singles_stops[singles_stops$low_sigma == TRUE | singles_stops$sig_10 == TRUE,]
# Drop the nonsense
singles<-singles_stops[singles_stops$STOP == F,]

# Pre-processing of singles
singles_ADan<-singles
singles_ADan$protein<-"ADan"
singles_ADan$Pos_c<-singles_ADan$Pos
singles_ADan$nscore_c_ADan<-singles_ADan$nscore_c
singles_ADan$wt_ADan<-singles_ADan$WT_AA
singles_ADan$pos_ADan<-singles_ADan$Pos
singles_ADan$mut_ADan<-singles_ADan$Mut

singles_ADan<-left_join(singles_ADan, AA_type[,c("wt_aa", "type")], by=c("wt_ADan"="wt_aa"))
singles_ADan<-rename(singles_ADan, type_wt = type)
singles_ADan<-left_join(singles_ADan, AA_type[,c("wt_aa", "type")], by=c("mut_ADan"="wt_aa"))
singles_ADan<-rename(singles_ADan, type_mut = type)

singles_ADan$residue_ADan<-paste0(singles_ADan$WT_AA, singles_ADan$Pos)
singles_ADan$change<-paste0(singles_ADan$Pos, singles_ADan$Mut)
singles_ADan$change_type<-paste0(singles_ADan$Pos, singles_ADan$type_mut)

###
# Load ABri data
load("nscore_df_ABri.RData")

# Let's work only with those trustable sequences
singles_stops<-singles_stops[singles_stops$low_sigma == T ,]

# Drop the nonsense
singles<-singles_stops[singles_stops$STOP == F,]


# Pre-processing of singles
singles_ABri<-singles
singles_ABri$protein<-"ABri"
singles_ABri$Pos_c<-singles_ABri$Pos
singles_ABri$nscore_c_ABri<-singles_ABri$nscore_c
singles_ABri$wt_ABri<-singles_ABri$WT_AA
singles_ABri$pos_ABri<-singles_ABri$Pos
singles_ABri$mut_ABri<-singles_ABri$Mut

singles_ABri<-left_join(singles_ABri, AA_type[,c("wt_aa", "type")], by=c("wt_ABri"="wt_aa"))
singles_ABri<-rename(singles_ABri, type_wt = type)
singles_ABri<-left_join(singles_ABri, AA_type[,c("wt_aa", "type")], by=c("mut_ABri"="wt_aa"))
singles_ABri<-rename(singles_ABri, type_mut = type)

singles_ABri$residue_ABri<-paste0(singles_ABri$WT_AA, singles_ABri$Pos)
singles_ABri$change<-paste0(singles_ABri$Pos, singles_ABri$Mut)


################################################################################
# Singles correlation 

singles_ADan$change<-paste0(singles_ADan$Pos_c, singles_ADan$Mut)
singles_ADan$change_pos<-paste0(singles_ADan$WT_AA, singles_ADan$Pos_c, singles_ADan$Mut)
singles_ADan$change_aa<-paste0(singles_ADan$wt_ADan, singles_ADan$mut_ADan)
singles_ADan$change_type<-paste0(singles_ADan$Pos_c, singles_ADan$type_mut)
singles_ADan$change_pos_type<-paste0(singles_ADan$type_wt, singles_ADan$Pos_c, singles_ADan$type_mut)

singles_ABri$change<-paste0(singles_ABri$Pos_c, singles_ABri$Mut)
singles_ABri$change_pos<-paste0(singles_ABri$WT_AA, singles_ABri$Pos_c, singles_ABri$Mut)
singles_ABri$change_aa<-paste0(singles_ABri$wt_ABri, singles_ABri$mut_ABri)
singles_ABri$change_type<-paste0(singles_ABri$Pos_c, singles_ABri$type_mut)
singles_ABri$change_pos_type<-paste0(singles_ABri$type_wt, singles_ABri$Pos_c, singles_ABri$type_mut)

### Correlation per position - mut
singles_ADan_ABri<-inner_join(singles_ADan, singles_ABri, by="change")
singles_ADan_ABri<-rename(singles_ADan_ABri, Pos = Pos_c.x, change_ADan = change_pos.x, change_ABri = change_pos.y)
singles_ADan_ABri<-singles_ADan_ABri[c("change", "change_ADan", "change_ABri", "Pos", "residue_ADan", "residue_ABri", "nscore_c_ADan", "nscore_c_ABri")]
singles_ADan_ABri<-arrange(singles_ADan_ABri, Pos)
singles_ADan_ABri$shared_pos<-FALSE
singles_ADan_ABri[singles_ADan_ABri$change_ADan == singles_ADan_ABri$change_ABri,]$shared_pos<-TRUE

correlation_all<-cor.test(singles_ADan_ABri$nscore_c_ADan, singles_ADan_ABri$nscore_c_ABri, use="complete.obs")

LR<-lm(singles_ADan_ABri$nscore_c_ABri~singles_ADan_ABri$nscore_c_ADan)

singles_ADan_ABri$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_ABri$outliers<-FALSE
singles_ADan_ABri[singles_ADan_ABri$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_ABri$OutliersLabel<-singles_ADan_ABri$change
singles_ADan_ABri[singles_ADan_ABri$outliers == FALSE,]$OutliersLabel<-''


p_corr_all<-ggplot(singles_ADan_ABri, aes(x=nscore_c_ADan, y=nscore_c_ABri))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(color="grey", size=2)+
  xlim(-6, 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=5, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.3, vjust=3.5, size=5, color="black")+
  labs(x="Nucleation Score ADan", y="Nucleation Score ABri")
p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_pos_ADan_ABri.jpg", width = 5, height = 4, path=path)

# 
corr_vector=c()
for(i in singles_ADan_ABri$Pos){
  if (nrow(singles_ADan_ABri[singles_ADan_ABri$Pos == i,])>4){
    correlation<-cor.test(singles_ADan_ABri[singles_ADan_ABri$Pos==i,][["nscore_c_ADan"]],
                          singles_ADan_ABri[singles_ADan_ABri$Pos==i,][["nscore_c_ABri"]], use="complete.obs")
    corr<-correlation$estimate
    p_value<-correlation$p.value
    corr_vector=c(corr_vector, i, corr, p_value)
    }
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}


###
p_corr_pos<-ggplot(singles_ADan_ABri[singles_ADan_ABri$Pos %in% corr_text$Pos,], aes(x=nscore_c_ADan, y=nscore_c_ABri))+
  geom_rect(data = singles_ADan_ABri[singles_ADan_ABri$Pos %in% corr_text$Pos,], fill = "white", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(color="grey50", alpha=0.5, size=1)+
  facet_wrap(~factor(Pos), ncol=6)+
  scale_y_continuous(breaks=c(-4, 0, 4, 8), labels=c(-4, 0, 4, 8), lim=c(-4, 8))+
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4), labels=c(-4, -2, 0, 2, 4), lim=c(-4, 4))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.1, vjust=1.5, colour="black", size=3)+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.1, vjust=3, colour="black", size=3)+
  geom_text(data=singles_ADan_ABri[singles_ADan_ABri$Pos %in% corr_text$Pos,], aes(label=residue_ADan, x=-Inf, y=-Inf, color="orange"), hjust=-0.2, vjust=-0.5)+
  geom_text(data=singles_ADan_ABri[singles_ADan_ABri$Pos %in% corr_text$Pos,], aes(label=residue_ABri, x=Inf, y=-Inf, color="brown"), hjust=1.1, vjust=-0.5)+
  scale_color_manual(name = "Residue", labels = c("Mouse", "Human"), values=c("brown", "orange"))+
  labs(x="Nucleation Score ADan", y="Nucleation Score ABri")
p_corr_pos

ggsave(p_corr_pos, file="Corr_singles_nscore_poseach_ADan_ABri.jpg", width = 6, height = 4, path=path)


### Correlation per mean per position

singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(Pos) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))

singles_ABri_grouped<-as.data.frame(singles_ABri %>% group_by(Pos) %>% dplyr::summarise(mean_ns_ABri=mean(nscore_c_ABri)))

singles_ADan_ABri<-inner_join(singles_ADan_grouped, singles_ABri_grouped, by="Pos")

correlation_all<-cor.test(singles_ADan_ABri$mean_ns_ADan, singles_ADan_ABri$mean_ns_ABri, use="complete.obs")

LR<-lm(singles_ADan_ABri$mean_ns_ABri~singles_ADan_ABri$mean_ns_ADan)

singles_ADan_ABri$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_ABri$outliers<-FALSE
singles_ADan_ABri[singles_ADan_ABri$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_ABri$OutliersLabel<-singles_ADan_ABri$type_change
singles_ADan_ABri[singles_ADan_ABri$outliers == FALSE,]$OutliersLabel<-''

p_corr_pos<-ggplot(singles_ADan_ABri, aes(x=mean_ns_ADan, y=mean_ns_ABri))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(color="grey", size=2)+
  xlim(-3, 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=5, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.3, vjust=3.5, size=5, color="black")+
  labs(x="Nucleation Score ADan", y="Nucleation Score ABri", color="")
p_corr_pos

ggsave(p_corr_pos, file="Corr_singles_nscore_meanpos_ADan_ABri.jpg", width = 5, height = 4, path=path)




### Correlation per aa change
singles_ADan$type_change<-paste0(singles_ADan$type_wt, '-', singles_ADan$type_mut)
singles_ABri$type_change<-paste0(singles_ABri$type_wt, '-', singles_ABri$type_mut)

singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(change_aa, type_change) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))
singles_ABri_grouped<-as.data.frame(singles_ABri %>% group_by(change_aa, type_change) %>% dplyr::summarise(mean_ns_ABri=mean(nscore_c_ABri)))

singles_ADan_ABri<-inner_join(singles_ADan_grouped, singles_ABri_grouped, by=c("change_aa", "type_change"))

correlation_all<-cor.test(singles_ADan_ABri$mean_ns_ADan, singles_ADan_ABri$mean_ns_ABri, use="complete.obs")

LR<-lm(singles_ADan_ABri$mean_ns_ABri~singles_ADan_ABri$mean_ns_ADan)

singles_ADan_ABri$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_ABri$outliers<-FALSE
singles_ADan_ABri[singles_ADan_ABri$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_ABri$OutliersLabel<-singles_ADan_ABri$change_aa
singles_ADan_ABri[singles_ADan_ABri$outliers == FALSE,]$OutliersLabel<-''

changes_types<-data.frame(
                type_change=c("aliphatic-polar", "aliphatic-negative",  
                              "aliphatic-aromatic", "aliphatic-glycine", 
                              "aliphatic-aliphatic", "aliphatic-proline", 
                              "aliphatic-positive", "polar-aliphatic",
                              "polar-negative", "polar-aromatic", 
                              "polar-glycine", "polar-polar", 
                              "polar-positive", "polar-proline", 
                              "negative-aliphatic", "negative-polar", 
                              "negative-negative", "negative-aromatic", 
                              "negative-glycine", "negative-positive", 
                              "negative-proline", "aromatic-aliphatic", 
                              "aromatic-polar", "aromatic-negative", 
                              "aromatic-glycine", "aromatic-proline", 
                              "aromatic-positive", "aromatic-aromatic", 
                              "positive-aliphatic", "positive-polar", 
                              "positive-negative", "positive-aromatic", 
                              "positive-glycine", "positive-proline", 
                              "positive-positive"),
                type_change_resume=c("aliphatic-polar", "aliphatic-charged",  
                                     "aliphatic-aromatic", "aliphatic-G/P", 
                                     "aliphatic-aliphatic", "aliphatic-G/P", 
                                     "aliphatic-charged", "polar-aliphatic",
                                     "polar-charged", "polar-aromatic", 
                                     "polar-G/P", "polar-polar", 
                                     "polar-charged", "polar-G/P", 
                                     "charged-aliphatic", "charged-polar", 
                                     "charged-charged", "charged-aromatic", 
                                     "charged-G/P", "charged-charged", 
                                     "charged-G/P", "aromatic-aliphatic", 
                                     "aromatic-polar", "aromatic-charged", 
                                     "aromatic-G/P", "aromatic-G/P", 
                                     "aromatic-charged", "aromatic-aromatic", 
                                     "charged-aliphatic", "charged-polar", 
                                     "charged-charged", "charged-aromatic", 
                                     "charged-G/P", "charged-G/P", 
                                     "charged-charged"),
                type_change_simplified=c("aliphatic-polar / polar-aliphatic", "aliphatic-charged / charged-aliphatic",
                                         "aliphatic-aromatic / aromatic-aliphatic", "aliphatic-GP / GP-aliphatic", 
                                         "no type change", "aliphatic-GP / GP-aliphatic", 
                                         "aliphatic-charged / charged-aliphatic", "aliphatic-polar / polar-aliphatic",
                                         "polar-charged / charged-polar", "polar-aromatic / aromatic-polar", 
                                         "polar-GP / GP-polar", "no type change", 
                                         "polar-charged / charged-polar", "polar-GP / GP-polar", 
                                         "aliphatic-charged / charged-aliphatic", "polar-charged / charged-polar", 
                                         "no type change", "charged-aromatic / aromatic-charged", 
                                         "charged-GP / GP-charged", "no type change", 
                                         "charged-GP / GP-charged", "aliphatic-aromatic / aromatic-aliphatic", 
                                         "polar-aromatic / aromatic-polar", "charged-aromatic / aromatic-charged", 
                                         "aromatic-GP / GP-aromatic", "aromatic-GP / GP-aromatic", 
                                         "charged-aromatic / aromatic-charged", "no type change", 
                                         "aliphatic-charged / charged-aliphatic", "polar-charged / charged-polar", 
                                         "no type change", "charged-aromatic / aromatic-charged", 
                                         "charged-GP / GP-charged", "charged-GP / GP-charged", 
                                         "no type change"))

singles_ADan_ABri<-inner_join(singles_ADan_ABri, changes_types, by=c("type_change"))

changes_colours<-c("grey10", "brown2", "brown3", "brown4", "brown", "lightblue",
                   "green3", "red2", "gold1", "purple", "orange")

p_corr_all<-ggplot(singles_ADan_ABri, aes(x=mean_ns_ADan, y=mean_ns_ABri))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(aes(color=factor(type_change_simplified, 
                              levels=c("no type change", "aliphatic-GP / GP-aliphatic", 
                                       "aromatic-GP / GP-aromatic", "polar-GP / GP-polar", 
                                       "charged-GP / GP-charged", "aliphatic-aromatic / aromatic-aliphatic",
                                       "aliphatic-polar / polar-aliphatic", "aliphatic-charged / charged-aliphatic",
                                       "polar-aromatic / aromatic-polar", "charged-aromatic / aromatic-charged",
                                       "polar-charged / charged-polar"))), 
             size=2)+
  scale_color_manual(values=changes_colours)+
  xlim(-5, 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=5, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.3, vjust=3.5, size=5, color="black")+
  labs(x="Nucleation Score ADan", y="Nucleation Score ABri", color="")
p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_aachange_ADan_ABri.jpg", width = 8.5, height = 4, path=path)


### Correlation per type change

#singles_ADan$type_change<-paste0(singles_ADan$type_wt, '-', singles_ADan$type_mut)
singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(type_change) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))

#singles_ABri$type_change<-paste0(singles_ABri$type_wt, '-', singles_ABri$type_mut)
singles_ABri_grouped<-as.data.frame(singles_ABri %>% group_by(type_change) %>% dplyr::summarise(mean_ns_ABri=mean(nscore_c_ABri)))

singles_ADan_ABri<-inner_join(singles_ADan_grouped, singles_ABri_grouped, by="type_change")

correlation_all<-cor.test(singles_ADan_ABri$mean_ns_ADan, singles_ADan_ABri$mean_ns_ABri, use="complete.obs")

LR<-lm(singles_ADan_ABri$mean_ns_ABri~singles_ADan_ABri$mean_ns_ADan)

singles_ADan_ABri$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_ABri$outliers<-FALSE
singles_ADan_ABri[singles_ADan_ABri$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_ABri$OutliersLabel<-singles_ADan_ABri$type_change
singles_ADan_ABri[singles_ADan_ABri$outliers == FALSE,]$OutliersLabel<-''

singles_ADan_ABri<-inner_join(singles_ADan_ABri, changes_types, by=c("type_change"))

p_corr_all_type<-ggplot(singles_ADan_ABri, aes(x=mean_ns_ADan, y=mean_ns_ABri))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(aes(color=factor(type_change_simplified, 
                              levels=c("no type change", "aliphatic-GP / GP-aliphatic", 
                                       "aromatic-GP / GP-aromatic", "polar-GP / GP-polar", 
                                       "charged-GP / GP-charged", "aliphatic-aromatic / aromatic-aliphatic",
                                       "aliphatic-polar / polar-aliphatic", "aliphatic-charged / charged-aliphatic",
                                       "polar-aromatic / aromatic-polar", "charged-aromatic / aromatic-charged",
                                       "polar-charged / charged-polar"))), 
             size=2)+
  scale_color_manual(values=changes_colours)+
  xlim(-3, 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=5, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.3, vjust=3.5, size=5, color="black")+
  labs(x="Nucleation Score ADan", y="Nucleation Score ABri", color="")
p_corr_all_type

ggsave(p_corr_all_type, file="Corr_singles_nscore_typechange_ADan_ABri.jpg", width = 8.5, height = 4, path=path)


### Correlation per change to

singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(Mut) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))

singles_ABri_grouped<-as.data.frame(singles_ABri %>% group_by(Mut) %>% dplyr::summarise(mean_ns_ABri=mean(nscore_c_ABri)))

singles_ADan_ABri<-inner_join(singles_ADan_grouped, singles_ABri_grouped, by="Mut")

correlation_all<-cor.test(singles_ADan_ABri$mean_ns_ADan, singles_ADan_ABri$mean_ns_ABri, use="complete.obs")

LR<-lm(singles_ADan_ABri$mean_ns_ABri~singles_ADan_ABri$mean_ns_ADan)

singles_ADan_ABri$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))


singles_ADan_ABri<-inner_join(singles_ADan_ABri, AA_type, by=c("Mut"="wt_aa"))

p_corr_to<-ggplot(singles_ADan_ABri, aes(x=mean_ns_ADan, y=mean_ns_ABri))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(aes(color=type), size=2)+
  geom_text_repel(aes(label=Mut, color=type), size=5, min.segment.length=0.1)+
  scale_color_manual(values=color_AAtype)+
  xlim(-2, 1)+
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4), lim=c(0, 4))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=5, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.3, vjust=3.5, size=5, color="black")+
  labs(x="Nucleation Score ADan", y="Nucleation Score ABri", color="")
p_corr_to

ggsave(p_corr_to, file="Corr_singles_nscore_to_ADan_ABri.jpg", width = 6, height = 4, path=path)


### Correlation per change from

singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(WT_AA) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))

singles_ABri_grouped<-as.data.frame(singles_ABri %>% group_by(WT_AA) %>% dplyr::summarise(mean_ns_ABri=mean(nscore_c_ABri)))

singles_ADan_ABri<-inner_join(singles_ADan_grouped, singles_ABri_grouped, by="WT_AA")

correlation_all<-cor.test(singles_ADan_ABri$mean_ns_ADan, singles_ADan_ABri$mean_ns_ABri, use="complete.obs")

LR<-lm(singles_ADan_ABri$mean_ns_ABri~singles_ADan_ABri$mean_ns_ADan)

singles_ADan_ABri$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))


singles_ADan_ABri<-inner_join(singles_ADan_ABri, AA_type, by=c("WT_AA"="wt_aa"))

p_corr_from<-ggplot(singles_ADan_ABri, aes(x=mean_ns_ADan, y=mean_ns_ABri))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(aes(color=type), size=2)+
  geom_text_repel(aes(label=WT_AA, color=type), size=5, min.segment.length=0.1)+
  scale_color_manual(values=color_AAtype)+
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4), lim=c(0, 4))+
  xlim(-2, 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=16),
        legend.text= element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=5, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.3, vjust=3.5, size=5, color="black")+
  labs(x="Nucleation Score ADan", y="Nucleation Score ABri", color="")
p_corr_from

ggsave(p_corr_from, file="Corr_singles_nscore_from_ADan_ABri.jpg", width = 6, height = 4, path=path)




