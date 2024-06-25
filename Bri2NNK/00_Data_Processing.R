library(tidyverse)
library("ggrepel")


# Mode value
mode <- function(x, na.rm = FALSE) {
  
  if(na.rm){ #if na.rm is TRUE, remove NA values from input x
    x = x[!is.na(x)]
  }
  
  val <- unique(x)
  return(val[which.max(tabulate(match(x, val)))])
}

#
load("1/Bri2_12NNK_01_all_variants_grouped.RData")
print("Sequences in total Rep 1:")
nrow(all_variants_df_grouped[!duplicated(all_variants_df_grouped$aa_seq_cor), ])
print("Sequences with NS Rep 1:")
nrow(all_variants_df_grouped[!duplicated(all_variants_df_grouped$aa_seq_cor), ])
print("Sequences low nucleators Rep 1:")
nrow(all_variants_df_grouped[all_variants_df_grouped$category_fdr == "non-nucleators", ])
print("Sequences nucleating Rep 1:")
nrow(all_variants_df_grouped[all_variants_df_grouped$category_fdr != "non-nucleators", ])
all_variants_df_grouped$replicate<-1

NNK_all_df<-all_variants_df_grouped

################################################################################
#
load("2/Bri2_12NNK_02_all_variants_grouped.RData")
print("Sequences in total Rep 2:")
nrow(all_variants_df_grouped[!duplicated(all_variants_df_grouped$aa_seq_cor), ])
print("Sequences with NS Rep 2:")
nrow(all_variants_df_grouped[!duplicated(all_variants_df_grouped$aa_seq_cor), ])
print("Sequences low nucleators Rep 2:")
nrow(all_variants_df_grouped[all_variants_df_grouped$category_fdr == "non-nucleators", ])
print("Sequences nucleating Rep 2:")
nrow(all_variants_df_grouped[all_variants_df_grouped$category_fdr != "non-nucleators", ])
all_variants_df_grouped$replicate<-2

NNK_all_df<-rbind(NNK_all_df, all_variants_df_grouped)

### Individual Validation

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
  geom_label_repel(aes(label=aa_seq_cor), seed=42, box.padding = 0.25, min.segment.length = 0.1)+
  theme_classic()+
  labs(x="Small scale" ,y="Large scale")

p_small_large_scale

ggsave(p_small_large_scale, file="p_individualvalidation.jpg", width=4, height=4)

#
load("3/Bri2_12NNK_03_all_variants_grouped.RData")
print("Sequences in total Rep 3:")
nrow(all_variants_df_grouped[!duplicated(all_variants_df_grouped$aa_seq_cor), ])
print("Sequences with NS Rep 3:")
nrow(all_variants_df_grouped[!duplicated(all_variants_df_grouped$aa_seq_cor), ])
print("Sequences low nucleators Rep 3:")
nrow(all_variants_df_grouped[all_variants_df_grouped$category_fdr == "non-nucleators", ])
print("Sequences nucleating Rep 3:")
nrow(all_variants_df_grouped[all_variants_df_grouped$category_fdr != "non-nucleators", ])
all_variants_df_grouped$replicate<-3

NNK_all_df<-rbind(NNK_all_df, all_variants_df_grouped)

#
print("Duplicated sequences in the pool:")
print(NNK_all_df[duplicated(NNK_all_df$aa_seq_cor), "aa_seq_cor"])
print(NNK_all_df[duplicated(NNK_all_df$aa_seq_cor),])

# Group them and get mode for seeds and the weighted mean for nscore
NNK_all_df<-NNK_all_df %>% 
  group_by(aa_seq_cor) %>% 
  mutate(mean_nscore=mean(mean_nscore), 
         mode_seed_bh=mode(mode_seed_bh),
         category_fdr=mode(category_fdr))

# Drop duplicated sequences
NNK_all_df<-NNK_all_df[!duplicated(NNK_all_df$aa_seq_cor),]

# Plot the variants distribution
p_hist_01<-ggplot(NNK_all_df)+
  geom_histogram(bins=200, aes(x=mean_nscore, fill=category_fdr))+
  geom_vline(xintercept=0, color="grey", linetype="dashed")+
  theme_bw()
p_hist_01

###

save(NNK_all_df,  file = "NNK_all_df.RData")
write.table(NNK_all_df, file= paste0("NNK_all_df.tsv"), sep="\t", quote = F, row.names = F)

