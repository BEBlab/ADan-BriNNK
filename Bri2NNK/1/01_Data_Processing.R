library(tidyverse)
library(ggpubr)
library(readxl)
library(caret)
library(modeest)


dir.create("01_Data_Processing")
path="01_Data_Processing"

name<-"Bri2_12NNK_01"

###
# Known Bri2 sequences

Bri2<-"EASNCFAIRHFENKFAVETLICS"
Bri2_cterm<-substr(Bri2, 23, 23)
Bri2_nterm<-substr(Bri2, 1, 22)
ABri<-"EASNCFAIRHFENKFAVETLICSRTVKKNIIEEN"
ABri_cterm<-substr(ABri, 23, 34)
ABri_FCD<-"EASNCFAIRHFENKFAVETLICSLTVKKNIIEEN"
ABri_FCD_cterm<-substr(ABri_FCD, 23, 34)
ABri_FKD<-"EASNCFAIRHFENKFAVETLICSSTVKKNIIEEN"
ABri_FKD_cterm<-substr(ABri_FKD, 23, 34)
ADan<-"EASNCFAIRHFENKFAVETLICFNLFLNSQEKHY"
ADan_cterm<-substr(ADan, 23, 34)

# Mode value
mode <- function(x, na.rm = FALSE) {
  
  if(na.rm){ #if na.rm is TRUE, remove NA values from input x
    x = x[!is.na(x)]
  }
  
  val <- unique(x)
  return(val[which.max(tabulate(match(x, val)))])
}

###
# Load dimsum results

load(paste0(name, '_fitness_replicates.RData'))

all_variants_df<-all_variants

# Get the length of the sequences
all_variants_df$aa_len<-regexpr('[:*:]', all_variants_df$aa_seq)
all_variants_df$aa_len<-as.numeric(all_variants_df$aa_len)-1

# replace values -1 for 12 aa len
all_variants_df$aa_len<-replace(all_variants_df$aa_len, all_variants_df$aa_len == -2, 12)


# Write the corrected aa seq - remove aa after stop codon (it's easier for matching and analyzing)
all_variants_df$aa_seq_cor<-all_variants_df$aa_seq
for (i in 1:nrow(all_variants_df)){
  if (all_variants_df[i, "aa_len"] != 12){
    all_variants_df[i, "aa_seq_cor"]<-substr(all_variants_df[i, "aa_seq"], 1, all_variants_df[i, "aa_len"])
  }
}

# rename fitness 
all_variants_df<-rename(all_variants_df, nscore = fitness)

# Divide the sequences in datasets, check if there are any known sequence
# Truncated sequences have different nscore (after the stop codon the sequence may be different)
# Let's keep all of them separately
all_variants_df$dataset<-"random_full"
all_variants_df[all_variants_df$STOP == TRUE, "dataset"]<-"random_stop"
all_variants_df[all_variants_df$aa_seq_cor == Bri2_cterm, "dataset"]<-"Bri2"
all_variants_df[all_variants_df$aa_seq_cor == ABri_cterm, "dataset"]<-"ABri"
all_variants_df[all_variants_df$aa_seq_cor == ABri_FCD_cterm, "dataset"]<-"ABri_FCD"
all_variants_df[all_variants_df$aa_seq_cor == ABri_FKD_cterm, "dataset"]<-"ABri_FKD"
all_variants_df[all_variants_df$aa_seq_cor == ADan_cterm, "dataset"]<-"ADan"

# Drop the seq after stop codon
all_variants_df$aa_len<-nchar(all_variants_df$aa_seq_cor)

# Build the full sequence cloned
all_variants_df$full_seq<-paste0(Bri2_nterm, all_variants_df$aa_seq_cor)

## Calculate the mean value of Bri2 and mode of dead nucleators
mean_bri2<-weighted.mean(all_variants_df[all_variants_df$dataset == "Bri2", "nscore"], all_variants_df[all_variants_df$dataset == "Bri2", "sigma"]^-2, na.rm = T)
# Mode of distribution
#dead_mode<-mlv(all_variants_df$nscore, method = "meanshift")[1]
# Mode value (most repeated value)
dead_mode<-mode(all_variants_df$nscore)


# Plot the variants distribution
p_hist_00<-ggplot(all_variants_df)+
  geom_histogram(bins=200, aes(x=nscore), alpha=0.7)+
  geom_vline(xintercept=mean_bri2, color="grey", linetype="dashed")+
  annotate("text", label="Bri2", x=mean_bri2+0.3, y=90, color="grey")+
  geom_vline(xintercept=dead_mode, color="red", linetype="dashed")+
  annotate("text", label="Dead Mode", x=dead_mode+0.6, y=100, color="red")+
  labs(x="Nucleation Score", y="Counts")+
  theme_classic()
p_hist_00

ggsave(p_hist_00, file="p_hist_00.jpg", path=path, width = 4, height = 4)

# Centering at the mode of deads
all_variants_df$nscore<-as.numeric(paste(all_variants_df$nscore-dead_mode))


# test variants against WT at FDR=0.05
all_variants_df$zscore<-all_variants_df$nscore/all_variants_df$sigma
all_variants_df$p_adjust_dead_bh<-p.adjust(pnorm(all_variants_df$zscore, lower.tail=F), method = "BH")
#all_variants_df$p_adjust_dead_bonf<-p.adjust(pnorm(all_variants_df$zscore,lower.tail=F), method = "bonferroni")
all_variants_df$sig_dead_bh<-FALSE
#all_variants_df$sig_dead_bonf<-FALSE
all_variants_df[all_variants_df$p_adjust_dead_bh<0.05,]$sig_dead_bh<-TRUE
#all_variants_df[all_variants_df$p_adjust_dead_bonf<0.05,]$sig_dead_bonf<-TRUE
all_variants_df$seed_bh<-0
#all_variants_df$seed_bonf<-0
all_variants_df[all_variants_df$sig_dead_bh == T & all_variants_df$nscore > 0, 
                "seed_bh"]<-1
#all_variants_df[all_variants_df$sig_dead_bonf == T & all_variants_df$nscore > 0, 
#                "seed_bonf"]<-1

#
# interquartile range == fitness range

summary(all_variants_df$nscore)
iqr<-IQR(all_variants_df$nscore)

first_to_wt<-summary(all_variants_df$nscore)[[2]]

p_iqr<-ggplot(all_variants_df, aes(x=nscore))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  labs(x="Nucleation score", title="Fitness range")+
  
  geom_vline(xintercept=summary(all_variants_df$nscore)[[2]], color="red")+
  annotate("text", label=paste0("1st Qu= ", round(summary(all_variants_df$nscore)[[2]], 2)), 
           x=summary(all_variants_df$nscore)[[2]]-0.8, y=150, color="red")+
  
  geom_vline(xintercept=summary(all_variants_df$nscore)[[5]], color="red")+
  annotate("text", label=paste0("3rd Qu= ", round(summary(all_variants_df$nscore)[[5]], 2)), 
           x=summary(all_variants_df$nscore)[[5]]+0.8, y=150, color="red")+
  
  geom_vline(xintercept=summary(all_variants_df$nscore)[[3]], color="blue")+
  annotate("text", label=paste0("median= ",round(summary(all_variants_df$nscore)[[3]], 2) ), 
           x=summary(all_variants_df$nscore)[[3]]-0.55, y=150, color="blue")
p_iqr

ggsave(p_iqr, file="p_iqr.jpg", path=path, width = 5, height = 3)

#sigma distribution
p_sigma_dist<-ggplot(all_variants_df, aes(x=sigma))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  scale_x_continuous(limits = c(0,3))
p_sigma_dist

ggsave(p_sigma_dist, file="p_sigma_dist.jpg", path=path, width = 4, height = 4)

# normalise sigmas to fitness range - if fitness range is IQR

all_variants_df$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(all_variants_df$nscore))

# or if fitness range is 1rst to WT fitness

all_variants_df$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(all_variants_df$nscore)[[2]])

all_variants_df$sigma_norm_iqr<-all_variants_df$sigma / fitness_range_iqr
all_variants_df$sigma_norm_first_toWT<-all_variants_df$sigma / fitness_range_first_toWT

all_variants_df$sigma_norm_iqr<-as.numeric(all_variants_df$sigma_norm_iqr)
all_variants_df$sigma_norm_first_toWT<-as.numeric(all_variants_df$sigma_norm_first_toWT)

#sigma_norm_iqr

p1_sigma_normalised<-ggplot(all_variants_df, aes(x=sigma_norm_iqr))+
  geom_histogram(bins=200)+
  labs(title="sigma normalised to fitness range (1st Qu to 3r Qu)")+
  theme_bw()+
  scale_x_continuous(limits = c(0,1))+
  geom_vline(xintercept = 1/3, linetype="dashed", color="grey")+
  geom_vline(xintercept=1/4, linetype="dashed", color="grey")+
  geom_vline(xintercept=0.1, linetype="dashed", color="grey")
p1_sigma_normalised

p2_sigma_normalised<-ggplot(all_variants_df, aes(x=sigma_norm_first_toWT))+
  geom_histogram(bins=200)+
  labs(title="sigma normalised to WT")+
  theme_bw()+
  scale_x_continuous(limits = c(0,1))
p2_sigma_normalised

p_sigma_normalised<-ggarrange(p1_sigma_normalised, p2_sigma_normalised, nrow=2)
p_sigma_normalised
ggsave(p_sigma_normalised, file="p_sigma_normalised.jpg", path=path, width = 5, height = 4)

## Calculate the mean value of Bri2 and mode of dead nucleators
mean_bri2<-weighted.mean(all_variants_df[all_variants_df$dataset == "Bri2", "nscore"], all_variants_df[all_variants_df$dataset == "Bri2", "sigma"]^-2, na.rm = T)
dead_mode<-mode(all_variants_df$nscore)

# Plot the variants distribution
p_hist_01<-ggplot(all_variants_df)+
  geom_histogram(bins=200, aes(x=nscore), alpha=0.7)+
  geom_vline(xintercept=mean_bri2, color="grey", linetype="dashed")+
  annotate("text", label="Bri2", x=mean_bri2+0.3, y=90, color="grey")+
  geom_vline(xintercept=dead_mode, color="red", linetype="dashed")+
  annotate("text", label="Dead Mode", x=dead_mode+0.6, y=100, color="red")+
  labs(x="Nucleation Score", y="Counts")+
  theme_classic()
p_hist_01

ggsave(p_hist_01, file="p_hist_01.jpg", path=path, width = 4, height = 4)



all_variants_df<-all_variants_df[,c("aa_seq_cor", "dataset", "mean_count", "nscore", "sigma", "aa_len", "full_seq", "zscore", 
                                    "p_adjust_dead_bh", "sig_dead_bh", "seed_bh")]

####

# Read dimsum output file
variant_data_merge_df<-read_tsv(paste0(name,"_variant_data_merge.tsv"))

#find non-nucleating variants (they have >=100 input reads but 0 output reads- NS not calculated)
#each AA variant is only resulting from one nt sequence by design (non nuc nt seq results in non nuc aa seq)
non_nuc_df<-variant_data_merge_df[variant_data_merge_df$input1_e1_s0_bNA_count >= 100 & 
                                  variant_data_merge_df$output1_e1_s1_b1_count==0,]

# Get the length of the sequences
non_nuc_df$aa_len<-regexpr('[:*:]', non_nuc_df$aa_seq)
non_nuc_df$aa_len<-as.numeric(non_nuc_df$aa_len)-1

# replace values -1 for 12 aa len
non_nuc_df$aa_len<-replace(non_nuc_df$aa_len, non_nuc_df$aa_len == -2, 12)

# Write the corrected aa seq - remove aa after stop codon (it's easier for matching and analyzing)
non_nuc_df$aa_seq_cor<-non_nuc_df$aa_seq
for (i in 1:nrow(non_nuc_df)){
  if (non_nuc_df[i, "aa_len"] != 12){
    non_nuc_df[i, "aa_seq_cor"]<-substr(non_nuc_df[i, "aa_seq"], 1, non_nuc_df[i, "aa_len"])
  }
}

# remove duplicated aa_seq
#non_nuc_df<-non_nuc_df[!duplicated(non_nuc_df$aa_seq_cor), ]

non_nuc_df$full_seq<-paste0(Bri2_nterm, non_nuc_df$aa_seq_cor)

non_nuc_df$dataset<-"random_non_nuc"
non_nuc_df[,c("mean_count", "nscore", "sigma", "zscore", "p_adjust_dead_bh")]<-NA
non_nuc_df[,c("sig_dead_bh")]<-F
non_nuc_df[,c("seed_bh")]<-0
non_nuc_df<-non_nuc_df[,c("aa_seq_cor", "dataset", "mean_count", "nscore", "sigma", "aa_len", "full_seq", "zscore", "p_adjust_dead_bh", "sig_dead_bh", "seed_bh")]

all_variants_df<-rbind(all_variants_df, non_nuc_df)

all_variants_df[all_variants_df$aa_seq_cor == Bri2_cterm, "dataset"]<-"Bri2"
all_variants_df[all_variants_df$aa_seq_cor == ABri_cterm, "dataset"]<-"ABri"
all_variants_df[all_variants_df$aa_seq_cor == ABri_FCD_cterm, "dataset"]<-"ABri_FCD"
all_variants_df[all_variants_df$aa_seq_cor == ABri_FKD_cterm, "dataset"]<-"ABri_FKD"
all_variants_df[all_variants_df$aa_seq_cor == ADan_cterm, "dataset"]<-"ADan"

print(distinct(all_variants_df, dataset))

# Bri2 full-length sequences
all_variants_df$pos_23<-substring(all_variants_df$aa_seq_cor, 1, 1)
all_variants_df$Bri2_nterm<-F
all_variants_df[all_variants_df$pos_23 == "S", "Bri2_nterm"]<-T

print("Bri2 full-length variants:")
print(nrow(all_variants_df[all_variants_df$Bri2_nterm == T,]))
print("Bri2 full-length variants non nucleating (no NS):")
print(nrow(all_variants_df[all_variants_df$Bri2_nterm == T & all_variants_df$dataset == "random_non_nuc",]))


p_bar_01<-ggplot(non_nuc_df)+
  geom_bar(aes(x=factor(aa_len, levels=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
                        labels=c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34")), 
               fill=factor(aa_len, levels=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
                           labels=c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34"))))+
  labs(x="Sequence length", fill="Sequence length", title="Non nucleating sequences counts")+
  theme_classic()
p_bar_01

ggsave(p_bar_01, file="p_bar_01.jpg", path=path, width = 6, height = 4)

# Keep only those sequences with NS that aren't present in the non nucleating sequences pool
# This are the only ones that are not contradictory


# Group them and get mode for seeds and the weighted mean for nscore
all_variants_df_grouped<-all_variants_df %>% 
  group_by(aa_seq_cor) %>% 
  mutate(mean_nscore=weighted.mean(nscore, sigma^-2, na.rm = T), 
         mode_seed_bh=mode(seed_bh))
# Drop duplicated sequences
all_variants_df_grouped<-all_variants_df_grouped[!duplicated(all_variants_df_grouped$aa_seq_cor),]
all_variants_df_grouped<-all_variants_df_grouped[c("aa_seq_cor", "full_seq", "dataset", "aa_len", 
                                                   "mean_nscore", "mode_seed_bh", "Bri2_nterm")]

# Think what to do with duplicated sequences...
print("Sequences:")
print(nrow(all_variants_df))
print("No duplicated sequences:")
print(nrow(all_variants_df[!duplicated(all_variants_df$aa_seq_cor), ]))
print("No duplicated sequences grouped:")
print(nrow(all_variants_df_grouped))
print("Bri2 23 aa long non-duplicated:")
print(nrow(all_variants_df_grouped[all_variants_df_grouped$Bri2_nterm == T,]))


# Only those sequences with NS
NNK_variants_df<-all_variants_df_grouped[!is.na(all_variants_df_grouped$mean_nscore),]
NNK_variants_df$STOP<-F
NNK_variants_df[NNK_variants_df$aa_len <12, "STOP"]<-T

print("Sequences with NS:")
print(nrow(NNK_variants_df))
print("Bri2 23 aa long NS sequences:")
print(nrow(NNK_variants_df[NNK_variants_df$Bri2_nterm == T,]))

NNK_variants_df$category_fdr<-"empty"
# Sequences that are significantly different from the deads are nucleators
NNK_variants_df[NNK_variants_df$mode_seed_bh==1 & NNK_variants_df$mean_nscore>0, "category_fdr"]<-"nucleators"
# As there are some sequences that have big NS but aren't significant due to big errors, 
# we re-defined "non-nucleators" as:
NNK_variants_df[NNK_variants_df$mode_seed_bh==0 & 
                NNK_variants_df$mean_nscore < min(NNK_variants_df[NNK_variants_df$category_fdr == "nucleators", "mean_nscore"]), 
                "category_fdr"]<-"non-nucleators"
# Let's drop those sequences that have big errors and big NS
NNK_variants_df<-NNK_variants_df[NNK_variants_df$category_fdr != "empty",]
# TOP10% are Fast nucleators
percent_10<-round(nrow(NNK_variants_df)/10, 1)
NNK_variants_df<-arrange(NNK_variants_df, desc(mean_nscore))
percent_10_nscore<-min(NNK_variants_df[1:percent_10, "mean_nscore"])
NNK_variants_df[NNK_variants_df$mode_seed_bh==T & NNK_variants_df$mean_nscore >= percent_10_nscore, "category_fdr"]<-"Top10 nucleators"

print("Non-nucleating sequences:")
print(nrow(NNK_variants_df[NNK_variants_df$category_fdr == "non-nucleators",]))
print("Nucleating sequences:")
print(nrow(NNK_variants_df[NNK_variants_df$category_fdr == "nucleators" | NNK_variants_df$category_fdr == "Top10 nucleators",]))

print("Bri2 23 aa long NS sequences:")
print("Non-nucleating sequences:") 
print(nrow(NNK_variants_df[NNK_variants_df$category_fdr == "non-nucleators" & NNK_variants_df$Bri2_nterm == T,]))
print("Nucleating sequences:")
print(nrow(NNK_variants_df[NNK_variants_df$category_fdr != "non-nucleators" & NNK_variants_df$Bri2_nterm == T,]))


# Plot the variants distribution, faceting by the presence or not of premature STOP codon
p_hist_02<-ggplot(NNK_variants_df)+
  geom_histogram(aes(x=mean_nscore, fill=category_fdr))+
  scale_fill_manual(values=c("#DF9292", "#7979BE", "darkblue"))+
  facet_wrap(~STOP)+
  labs(x="Nucleation Score", fill="Category")+
  theme_bw()
p_hist_02

ggsave(p_hist_02, file="p_hist_02.jpg", path=path, width = 9, height = 4)

# Plot the variants distribution of truncated sequences, faceting by the length of the sequence
p_hist_03<-ggplot(NNK_variants_df[NNK_variants_df$STOP == T,])+
  geom_histogram(aes(x=mean_nscore, fill=category_fdr))+
  scale_fill_manual(values=c("#DF9292", "#7979BE", "darkblue"))+
  labs(x="Nucleation Score", fill="Category")+
  facet_wrap(~factor(aa_len, levels=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels=c("22", "23", "24","25", "26", "27", "28", "29", "30", "31", "32", "33", "34")))+
  theme_bw()+
  theme(legend.position = c(.9, .15))
p_hist_03

ggsave(p_hist_03, file="p_hist_03.jpg", path=path, width = 8, height = 6)


###
print("Total number of unique aa sequence in the pool:")
print(nrow(all_variants_df[!duplicated(all_variants_df$aa_seq_cor), ]))
print("Total number of unique aa sequence with undoubtely NS:")
print(nrow(NNK_variants_df[!duplicated(NNK_variants_df$aa_seq_cor), ]))

###

# Plot the variants distribution
p_hist_04<-ggplot(NNK_variants_df)+
  geom_histogram(bins=200, aes(x=mean_nscore, fill=category_fdr))+
  geom_vline(xintercept=0, color="grey50", linetype="dashed")+
  scale_fill_manual(values=c("#DF9292", "#7979BE", "darkblue"))+
  labs(x="Nucleation Score", fill="Category")+
  theme_classic()+
  theme(legend.position = c(.8, .85))
p_hist_04

ggsave(p_hist_04, file="p_hist_04.jpg", path=path, width = 4, height = 4)


###
# Let'see the distribution for those sequences that have Bri2 23 aa long at the N-term.

# Plot the variants distribution of truncated sequences that have Bri2 full-length in the N-term, faceting by the length of the sequence
p_hist_05<-ggplot(NNK_variants_df[NNK_variants_df$STOP == T & NNK_variants_df$Bri2_nterm == T,])+
  geom_histogram(aes(x=mean_nscore, fill=category_fdr))+
  scale_fill_manual(values=c("#DF9292", "#7979BE", "darkblue"))+
  labs(x="Nucleation Score", fill="Category", title="Bri2 23 aa long + random")+
  facet_wrap(~factor(aa_len, levels=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels=c("22", "23", "24","25", "26", "27", "28", "29", "30", "31", "32", "33", "34")))+
  theme_bw()+
  theme(legend.position = c(.9, .15))
p_hist_05

ggsave(p_hist_05, file="p_hist_05.jpg", path=path, width = 8, height = 6.5)

###

# Plot the variants distribution for sequences that have Bri2 full-length in the N-term
p_hist_06<-ggplot(NNK_variants_df[NNK_variants_df$Bri2_nterm == T, ])+
  geom_histogram(bins=200, aes(x=mean_nscore, fill=category_fdr))+
  geom_vline(xintercept=0, color="grey50", linetype="dashed")+
  scale_fill_manual(values=c("#DF9292", "#7979BE", "darkblue"))+
  labs(x="Nucleation Score", fill="Category", title="Bri2 23 aa long + random")+
  theme_classic()+
  theme(legend.position = c(.8, .85))
p_hist_06

ggsave(p_hist_06, file="p_hist_06.jpg", path=path, width = 4, height = 4.5)


###

all_variants_df_grouped[is.na(all_variants_df_grouped$mean_nscore), "mean_nscore"]<-min(all_variants_df_grouped[!is.na(all_variants_df_grouped$mean_nscore), "mean_nscore"])

all_variants_df_grouped$category_fdr<-"empty"
# Sequences that are significantly different from the deads are nucleators
all_variants_df_grouped[all_variants_df_grouped$mode_seed_bh==1 & all_variants_df_grouped$mean_nscore>0, "category_fdr"]<-"nucleators"
# As there are some sequences that have big NS but aren't significant due to big errors, 
# we re-defined "non-nucleators" as:
all_variants_df_grouped[all_variants_df_grouped$mode_seed_bh==0 & 
                          all_variants_df_grouped$mean_nscore < min(all_variants_df_grouped[all_variants_df_grouped$category_fdr == "nucleators", "mean_nscore"]), 
                        "category_fdr"]<-"non-nucleators"
# Let's drop those sequences that have big errors and big NS
all_variants_df_grouped<-all_variants_df_grouped[all_variants_df_grouped$category_fdr != "empty",]

all_variants_df_grouped[all_variants_df_grouped$aa_seq_cor %in% NNK_variants_df[NNK_variants_df$category_fdr == "Top10 nucleators",]$aa_seq_cor, "category_fdr"]<-"Top10 nucleators"

# Plot the variants distribution of all the sequences
p_hist_04_all<-ggplot(all_variants_df_grouped)+
  geom_histogram(bins=200, aes(x=mean_nscore, fill=category_fdr))+
  geom_vline(xintercept=0, color="grey50", linetype="dashed")+
  scale_fill_manual(values=c("#DF9292", "#7979BE", "darkblue"))+
  labs(x="Nucleation Score", fill="Category")+
  theme_classic()+
  theme(legend.position = c(.8, .85))
p_hist_04_all

ggsave(p_hist_04_all, file="p_hist_04_all.jpg", path=path, width = 4, height = 4)


###


save(all_variants_df, all_variants_df_grouped, NNK_variants_df,  file = paste0(name, "_df.RData"))

save(non_nuc_df, file=paste0(name, "_non_nucleating_variants.RData"))

save(all_variants_df_grouped, file = paste0(name, "_all_variants_grouped.RData"))

write.table(all_variants_df_grouped, file= paste0(name, "_processed_data.tsv"), sep="\t", quote = F, row.names = F)

################################################################################
### Individual Validation
load(paste0(name, "_all_variants_grouped.RData"))

individual_validation<-read_tsv("Individual_variant_plating.txt")

individual_validation<-individual_validation %>% group_by(aa_seq_cor) %>% summarise(mean_growth=mean(growth_rate), 
                                                                                    std_growth=sd(growth_rate))

# Bri2 extensions individual validation:

individual_validation<-left_join(individual_validation, all_variants_df_grouped[c("aa_seq_cor", "mean_nscore")], by="aa_seq_cor")

corr<-cor.test(individual_validation$mean_nscore, individual_validation$mean_growth, use="complete.obs")
R<-corr$estimate
p<-corr$p.value


p_small_large_scale<-ggplot(individual_validation, aes(x=mean_growth, y=mean_nscore))+
  geom_smooth(method = "lm", se=F, color="grey90", linetype="dashed")+
  geom_point(size=2)+
  annotate("text", x = -Inf, y = Inf, vjust=1, hjust=-0.5, label=paste0("R=",round(R, 2)), size=5)+
  annotate("text", x = -Inf, y = Inf, vjust=2.5, hjust=-0.2,label =paste0("p=",format(p, digits = 2, scientific = T)), size=5)+
  geom_label_repel(aes(label=ID), seed=42, box.padding = 0.25, min.segment.length = 0.1)+
  ylim(-6, 2.5)+
  theme_classic()+
  labs(x="Small scale" ,y="Large scale")

p_small_large_scale

ggsave(p_small_large_scale, path=path, file="p_individualvalidation.jpg", width=4, height=4)
