library(tidyverse)
library(ggpubr)
require("ggrepel")
library(stringr)

dir.create("06_ADanvsAB42")
path="06_ADanvsAB42"

#required data:


all_aa<-c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

AA_type<-data.frame("wt_aa"= all_aa,
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
#

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

changes_colours<-c("grey10", "brown2", "brown3", "brown4", "brown", "lightblue",
                   "green3", "red2", "gold1", "purple", "orange")

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
# Load AB42 data
load("AB42_INDEL_datasets.RData")

singles_AB42<-singles.df

singles_AB42$protein<-"AB42"
singles_AB42$Pos_c<-singles_AB42$Pos
singles_AB42$nscore_c_AB42<-singles_AB42$nscore_c
singles_AB42$wt_AB42<-singles_AB42$WT_AA
singles_AB42$pos_AB42<-singles_AB42$Pos
singles_AB42$mut_AB42<-singles_AB42$Mut

singles_AB42<-left_join(singles_AB42, AA_type[,c("wt_aa", "type")], by=c("wt_AB42"="wt_aa"))
singles_AB42<-rename(singles_AB42, type_wt = type)
singles_AB42<-left_join(singles_AB42, AA_type[,c("wt_aa", "type")], by=c("mut_AB42"="wt_aa"))
singles_AB42<-rename(singles_AB42, type_mut = type)

singles_AB42$residue_AB42<-paste0(singles_AB42$WT_AA, singles_AB42$Pos)
singles_AB42$change<-paste0(singles_AB42$Pos, singles_AB42$Mut)


################################################################################
# Singles correlation - aligned at N-term
# ADan      EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY--------	34
# AB42      DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA	42

singles_ADan$change<-paste0(singles_ADan$Pos_c, singles_ADan$Mut)
singles_ADan$change_pos<-paste0(singles_ADan$WT_AA, singles_ADan$Pos_c, singles_ADan$Mut)
singles_ADan$change_aa<-paste0(singles_ADan$wt_ADan, singles_ADan$mut_ADan)
singles_ADan$change_type<-paste0(singles_ADan$Pos_c, singles_ADan$type_mut)
singles_ADan$change_pos_type<-paste0(singles_ADan$type_wt, singles_ADan$Pos_c, singles_ADan$type_mut)

singles_AB42$change<-paste0(singles_AB42$Pos_c, singles_AB42$Mut)
singles_AB42$change_pos<-paste0(singles_AB42$WT_AA, singles_AB42$Pos_c, singles_AB42$Mut)
singles_AB42$change_aa<-paste0(singles_AB42$wt_AB42, singles_AB42$mut_AB42)
singles_AB42$change_type<-paste0(singles_AB42$Pos_c, singles_AB42$type_mut)
singles_AB42$change_pos_type<-paste0(singles_AB42$type_wt, singles_AB42$Pos_c, singles_AB42$type_mut)

### Correlation per position - mut
singles_ADan_AB42<-inner_join(singles_ADan, singles_AB42, by="change")
singles_ADan_AB42<-rename(singles_ADan_AB42, Pos = Pos_c.x, change_ADan = change_pos.x, change_AB42 = change_pos.y)
singles_ADan_AB42<-singles_ADan_AB42[c("change", "change_ADan", "change_AB42", "Pos", "residue_ADan", "residue_AB42", "nscore_c_ADan", "nscore_c_AB42")]
singles_ADan_AB42<-arrange(singles_ADan_AB42, Pos)

singles_ADan_AB42$shared_pos<-FALSE
singles_ADan_AB42[singles_ADan_AB42$change_ADan == singles_ADan_AB42$change_AB42,]$shared_pos<-TRUE

correlation_all<-cor.test(singles_ADan_AB42$nscore_c_ADan, singles_ADan_AB42$nscore_c_AB42, use="complete.obs")

LR<-lm(singles_ADan_AB42$nscore_c_AB42~singles_ADan_AB42$nscore_c_ADan)

singles_ADan_AB42$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_AB42$outliers<-FALSE
singles_ADan_AB42[singles_ADan_AB42$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_AB42$OutliersLabel<-singles_ADan_AB42$change
singles_ADan_AB42[singles_ADan_AB42$outliers == FALSE,]$OutliersLabel<-''


p_corr_all_pos_mut<-ggplot(singles_ADan_AB42, aes(x=nscore_c_ADan, y=nscore_c_AB42))+
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
  labs(x="Nucleation Score ADan", y="Nucleation Score AB42")
p_corr_all_pos_mut

ggsave(p_corr_all_pos_mut, file="Corr_singles_nscore_pos_ADan_AB42.jpg", width = 5, height = 4, path=path)



### Correlation per aa change
singles_ADan$type_change<-paste0(singles_ADan$type_wt, '-', singles_ADan$type_mut)
singles_AB42$type_change<-paste0(singles_AB42$type_wt, '-', singles_AB42$type_mut)

singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(change_aa, type_change) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))
singles_AB42_grouped<-as.data.frame(singles_AB42 %>% group_by(change_aa, type_change) %>% dplyr::summarise(mean_ns_AB42=mean(nscore_c_AB42)))

singles_ADan_AB42<-inner_join(singles_ADan_grouped, singles_AB42_grouped, by=c("change_aa", "type_change"))

correlation_all<-cor.test(singles_ADan_AB42$mean_ns_ADan, singles_ADan_AB42$mean_ns_AB42, use="complete.obs")
LR<-lm(singles_ADan_AB42$mean_ns_AB42~singles_ADan_AB42$mean_ns_ADan)
singles_ADan_AB42$residual_abs<-as.numeric(abs(LR$residuals))
SD2<-2*sd(resid(LR))

singles_ADan_AB42<-inner_join(singles_ADan_AB42, changes_types, by=c("type_change"))

p_corr_all<-ggplot(singles_ADan_AB42, aes(x=mean_ns_ADan, y=mean_ns_AB42))+
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
  labs(x="Nucleation Score ADan", y="Nucleation Score AB42", color="")
p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_aachange_ADan_AB42.jpg", width = 8.5, height = 4, path=path)


### Correlation per type change

#singles_ADan$type_change<-paste0(singles_ADan$type_wt, '-', singles_ADan$type_mut)
singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(type_change) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))

#singles_AB42$type_change<-paste0(singles_AB42$type_wt, '-', singles_AB42$type_mut)
singles_AB42_grouped<-as.data.frame(singles_AB42 %>% group_by(type_change) %>% dplyr::summarise(mean_ns_AB42=mean(nscore_c_AB42)))

singles_ADan_AB42<-inner_join(singles_ADan_grouped, singles_AB42_grouped, by="type_change")

correlation_all<-cor.test(singles_ADan_AB42$mean_ns_ADan, singles_ADan_AB42$mean_ns_AB42, use="complete.obs")

LR<-lm(singles_ADan_AB42$mean_ns_AB42~singles_ADan_AB42$mean_ns_ADan)

singles_ADan_AB42$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_AB42$outliers<-FALSE
singles_ADan_AB42[singles_ADan_AB42$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_AB42$OutliersLabel<-singles_ADan_AB42$type_change
singles_ADan_AB42[singles_ADan_AB42$outliers == FALSE,]$OutliersLabel<-''

singles_ADan_AB42<-inner_join(singles_ADan_AB42, changes_types, by=c("type_change"))

p_corr_all_type<-ggplot(singles_ADan_AB42, aes(x=mean_ns_ADan, y=mean_ns_AB42))+
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
  labs(x="Nucleation Score ADan", y="Nucleation Score AB42", color="")
p_corr_all_type

ggsave(p_corr_all_type, file="Corr_singles_nscore_typechange_ADan_AB42.jpg", width = 8.5, height = 4, path=path)


### Correlation per change to

singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(Mut) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))

singles_AB42_grouped<-as.data.frame(singles_AB42 %>% group_by(Mut) %>% dplyr::summarise(mean_ns_AB42=mean(nscore_c_AB42)))

singles_ADan_AB42<-inner_join(singles_ADan_grouped, singles_AB42_grouped, by="Mut")

correlation_all<-cor.test(singles_ADan_AB42$mean_ns_ADan, singles_ADan_AB42$mean_ns_AB42, use="complete.obs")

LR<-lm(singles_ADan_AB42$mean_ns_AB42~singles_ADan_AB42$mean_ns_ADan)

singles_ADan_AB42$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_AB42$outliers<-FALSE
singles_ADan_AB42[singles_ADan_AB42$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_AB42$OutliersLabel<-singles_ADan_AB42$type_change
singles_ADan_AB42[singles_ADan_AB42$outliers == FALSE,]$OutliersLabel<-''

singles_ADan_AB42<-inner_join(singles_ADan_AB42, AA_type, by=c("Mut"="wt_aa"))

p_corr_to<-ggplot(singles_ADan_AB42, aes(x=mean_ns_ADan, y=mean_ns_AB42))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(aes(color=type), size=2)+
  geom_text_repel(aes(label=Mut, color=type), size=5, min.segment.length=0.1)+
  scale_color_manual(values=color_AAtype)+
  scale_y_continuous(breaks=c(-2, -1, 0, 1), lim=c(-2, 1))+
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
  labs(x="Nucleation Score ADan", y="Nucleation Score AB42", color="")
p_corr_to

ggsave(p_corr_to, file="Corr_singles_nscore_to_ADan_AB42.jpg", width = 6, height = 4, path=path)


### Correlation per change from

singles_ADan_grouped<-as.data.frame(singles_ADan %>% group_by(WT_AA) %>% dplyr::summarise(mean_ns_ADan=mean(nscore_c_ADan)))

singles_AB42_grouped<-as.data.frame(singles_AB42 %>% group_by(WT_AA) %>% dplyr::summarise(mean_ns_AB42=mean(nscore_c_AB42)))

singles_ADan_AB42<-inner_join(singles_ADan_grouped, singles_AB42_grouped, by="WT_AA")

correlation_all<-cor.test(singles_ADan_AB42$mean_ns_ADan, singles_ADan_AB42$mean_ns_AB42, use="complete.obs")

LR<-lm(singles_ADan_AB42$mean_ns_AB42~singles_ADan_AB42$mean_ns_ADan)

singles_ADan_AB42$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_AB42$outliers<-FALSE
singles_ADan_AB42[singles_ADan_AB42$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_AB42$OutliersLabel<-singles_ADan_AB42$type_change
singles_ADan_AB42[singles_ADan_AB42$outliers == FALSE,]$OutliersLabel<-''

singles_ADan_AB42<-inner_join(singles_ADan_AB42, AA_type, by=c("WT_AA"="wt_aa"))

p_corr_from<-ggplot(singles_ADan_AB42, aes(x=mean_ns_ADan, y=mean_ns_AB42))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(aes(color=type), size=2)+
  geom_text_repel(aes(label=WT_AA, color=type), size=5, min.segment.length=0.1)+
  scale_color_manual(values=color_AAtype)+
  scale_y_continuous(breaks=c(-2, -1, 0, 1), lim=c(-2, 1))+
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
  labs(x="Nucleation Score ADan", y="Nucleation Score AB42", color="")
p_corr_from

ggsave(p_corr_from, file="Corr_singles_nscore_from_ADan_AB42.jpg", width = 6, height = 4, path=path)

################################################################################
# Singles correlation - aligned at C-term
# ADan      --------EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY	34
# AB42      DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA	42

singles_AB42_cterm<-singles_AB42[(singles_AB42["pos_AB42"]>8),]
singles_AB42_cterm$pos_AB42<-singles_AB42_cterm$pos_AB42-8

singles_AB42_cterm$change<-paste0(singles_AB42_cterm$pos_AB42, singles_AB42_cterm$Mut)
singles_AB42_cterm$residue_AB42<-paste0(singles_AB42_cterm$WT_AA, singles_AB42_cterm$Pos_c)


### Correlation per position - mut
singles_ADan_AB42_cterm<-inner_join(singles_ADan, singles_AB42_cterm, by="change")
singles_ADan_AB42_cterm<-singles_ADan_AB42_cterm[c("change", "pos_AB42", "nscore_c_ADan", "nscore_c_AB42", "residue_AB42", "residue_ADan")]
singles_ADan_AB42_cterm<-as.data.frame(singles_ADan_AB42_cterm)
singles_ADan_AB42_cterm<-arrange(singles_ADan_AB42_cterm, pos_AB42)

correlation_all<-cor.test(singles_ADan_AB42_cterm$nscore_c_ADan, singles_ADan_AB42_cterm$nscore_c_AB42, use="complete.obs")

LR<-lm(singles_ADan_AB42_cterm$nscore_c_AB42~singles_ADan_AB42_cterm$nscore_c_ADan)

singles_ADan_AB42_cterm$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_AB42_cterm$outliers<-FALSE
singles_ADan_AB42_cterm[singles_ADan_AB42_cterm$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_AB42_cterm$OutliersLabel<-singles_ADan_AB42_cterm$change
singles_ADan_AB42_cterm[singles_ADan_AB42_cterm$outliers == FALSE,]$OutliersLabel<-''


p_corr_all_pos_mut_cterm<-ggplot(singles_ADan_AB42_cterm, aes(x=nscore_c_ADan, y=nscore_c_AB42))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(color="grey", size=2)+
  #geom_point(aes(color=outliers))+
  #scale_color_manual(values = c("grey", "red"))+
  #geom_text_repel(aes(label=OutliersLabel), size=3, min.segment.length=0.1)+
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
  labs(x="Nucleation Score ADan", y="Nucleation Score AB42")
p_corr_all_pos_mut_cterm

ggsave(p_corr_all_pos_mut_cterm, file="Corr_singles_nscore_pos_ADan_AB42_cterm.jpg", width = 5, height = 4, path=path)




################################################################################
# Singles correlation - aligned with ClustalO

##################################################################
# CLUSTAL O(1.2.4) multiple sequence alignment                   #
# ADan      ----EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY-----------	34 #
# AB42      DAEFRHDSGYEVHHQ-------KLVFFAEDVGSNKGAIIGLMVGGVVIA	42 #
#               . .. : ::*        .*: *   :.*::                  #
##################################################################

singles_AB42_align<-singles_AB42[(singles_AB42["pos_AB42"]>4) & (singles_AB42["pos_AB42"]<32),]

singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 5, 1)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 6, 2)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 7, 3)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 8, 4)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 9, 5)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 10, 6)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 11, 7)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 12, 8)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 13, 9)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 14, 10)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 15, 11)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 16, 19)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 17, 20)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 18, 21)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 19, 22)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 20, 23)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 21, 24)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 22, 25)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 23, 26)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 24, 27)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 25, 28)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 26, 29)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 27, 30)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 28, 31)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 29, 32)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 30, 33)
singles_AB42_align$pos_AB42<-replace(singles_AB42_align$pos_AB42, singles_AB42_align$pos_AB42 == 31, 34)


singles_AB42_align$change<-paste0(singles_AB42_align$pos_AB42, singles_AB42_align$Mut)
#singles_AB42_align$residue_AB42<-paste0(singles_AB42_align$WT_AA, singles_AB42_align$Pos_c)


### Correlation per position - mut
singles_ADan_AB42_align<-inner_join(singles_ADan, singles_AB42_align, by="change")
#singles_ADan_AB42_align<-rename(singles_ADan_AB42_align, change = change.x)
singles_ADan_AB42_align<-singles_ADan_AB42_align[c("change", "pos_AB42", "nscore_c_ADan", "nscore_c_AB42", "residue_AB42", "residue_ADan")]
singles_ADan_AB42_align<-as.data.frame(singles_ADan_AB42_align)
singles_ADan_AB42_align<-arrange(singles_ADan_AB42_align, pos_AB42)

correlation_all<-cor.test(singles_ADan_AB42_align$nscore_c_ADan, singles_ADan_AB42_align$nscore_c_AB42, use="complete.obs")

LR<-lm(singles_ADan_AB42_align$nscore_c_AB42~singles_ADan_AB42_align$nscore_c_ADan)

singles_ADan_AB42_align$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ADan_AB42_align$outliers<-FALSE
singles_ADan_AB42_align[singles_ADan_AB42_align$residual_abs>SD2,]$outliers<-TRUE

singles_ADan_AB42_align$OutliersLabel<-singles_ADan_AB42_align$change
singles_ADan_AB42_align[singles_ADan_AB42_align$outliers == FALSE,]$OutliersLabel<-''


p_corr_all_pos_mut_align<-ggplot(singles_ADan_AB42_align, aes(x=nscore_c_ADan, y=nscore_c_AB42))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point(color="grey", size=2)+
  #geom_point(aes(color=outliers))+
  #scale_color_manual(values = c("grey", "red"))+
  #geom_text_repel(aes(label=OutliersLabel), size=3, min.segment.length=0.1)+
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
  labs(x="Nucleation Score ADan", y="Nucleation Score AB42")
p_corr_all_pos_mut_align

ggsave(p_corr_all_pos_mut_align, file="Corr_singles_nscore_pos_ADan_AB42_aligned.jpg", width = 5, height = 4, path=path)





