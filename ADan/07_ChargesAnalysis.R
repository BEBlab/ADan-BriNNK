library(tidyverse)
require(ggplot2)
require(reshape2)

dir.create("07_ChargesAnalysis")
path="07_ChargesAnalysis"

##  import required data 

name<-"ADan"

load(paste0("nscore_df_", name, ".RData"))

# Let's work only with those trustable sequences
singles_stops<-singles_stops[singles_stops$low_sigma == TRUE | singles_stops$sig_10 == TRUE,]
# Drop the nonsense
singles<-singles_stops[singles_stops$STOP == F,]

#
peptide_seq<-'EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY'
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], '\n', n))
}

singles$Pos<-as.numeric(singles$Pos)

singles$split_seq<-strsplit(singles$aa_seq, "")

singles<-as.data.frame(singles)

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

#
count_table_all<-data.frame("aa_seq" = "", 
                            "G" = 0, "A" = 0, "V" = 0, "L" = 0, "M" = 0, "I" = 0, "F" = 0,
                            "Y" = 0, "W" = 0, "K" = 0, "R" = 0, "D" = 0, "E" = 0, "S" = 0,
                            "T" = 0, "C" = 0, "N" = 0, "Q" = 0, "H" = 0, "P" = 0)
for (i in 1:nrow(singles)){
  count_table<-as.data.frame(table(singles[i, "split_seq"]))
  count_df<-data.frame("aa_seq" = singles[i, "aa_seq"])
  for (j in vectorAA){
    count_df[j]<-count_table[count_table$Var1 == j, "Freq"]  
  }
  count_df[is.na(count_df)] <- 0
  count_table_all<-rbind(count_table_all, count_df)
}

singles<-left_join(singles, count_table_all, by="aa_seq")
singles$Neg_charge<-(singles$E + singles$D)*-1
singles$Pos_charge<-singles$R + singles$K
singles$Net_charge<-singles$Pos_charge + singles$Neg_charge
singles$Number_charge_residues<-singles$Pos_charge + (singles$Neg_charge)*-1

mean_ns<- singles %>% group_by(Net_charge, Number_charge_residues) %>% summarise(across(nscore_c, mean))
text_map<- singles %>% group_by(Net_charge, Number_charge_residues) %>% count(Pos_charge, Neg_charge)

text_map<-left_join(text_map, mean_ns, by=c("Net_charge", "Number_charge_residues"))
text_map$label<-paste0("+", text_map$Pos_charge, " ", "-", text_map$Neg_charge*-1)

# Considering all the residues
# NS - Net Charge - Number of Charged Residues
p1<-ggplot(singles, aes(x=Nham_aa, y=nscore_c))+
  geom_hline(yintercept = 0, color="black", linetype="dashed", size=0.3)+
  geom_violin(fill="grey80", scale="width")+
  geom_boxplot(width=0.1, outlier.shape = NA, size=0.3)+
  geom_text(data=text_map, aes(label=label, x=-Inf, y=-Inf), size=4, hjust=-0.1, vjust=-0.5)+
  scale_y_continuous("Nucleation Score", sec.axis = sec_axis(~ . + 10,  name="Number of charged residues"))+
  facet_grid(rows = vars(Number_charge_residues), cols = vars(Net_charge), switch = "x")+
  labs(x= "Net Charge", fill= "", title="")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y.right=element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.position="right")
p1
  
ggsave(p1, file="NetCharge_TotalChargedResidues_NS.jpg", width = 10, height = 4, path=path)


# NS - Net Charge - Number of Charged Residues


text_map<-singles %>% group_by(Pos_charge, Neg_charge) %>% summarise(across(nscore_c, mean))
text_map$Net_charge<-text_map$Pos_charge+text_map$Neg_charge

p2<-ggplot(singles, aes(x=Nham_aa, y=nscore_c))+
  geom_hline(yintercept = 0, color="black", linetype="dashed", size=0.3)+
  geom_violin(fill="grey80")+
  geom_boxplot(width=0.1, outlier.shape = NA, size=0.3)+
  geom_text(data=text_map, aes(label=Net_charge, x=-Inf, y=-Inf), size=4, hjust=-1, vjust=-0.5)+
  scale_y_continuous("Nucleation Score", sec.axis = sec_axis(~ . + 10,  name="Positively charged residues"))+
  facet_grid(rows = vars(Pos_charge), cols = vars(Neg_charge*-1), switch = "x")+
  labs(x= "Negatively charged residues", fill= "", title="")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y.right=element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.position="right")
p2

ggsave(p2, file="PosChargedResiudes_NegChargedResidues_NS.jpg", width = 3, height = 4, path=path)
