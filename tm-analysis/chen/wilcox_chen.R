# Wilcox Test #

library(tidyverse); packageVersion("tidyverse")   #version: 2.0.0
library(phyloseq); packageVersion("phyloseq")    #version: 1.42.0
library(lmerTest); packageVersion("lmerTest")    #version: 3.1.3
library(modeest); packageVersion("modeest")      #version: 2.4.0
library(MicrobiomeStat); packageVersion("MicrobiomeStat")  #version: 1.1
library(topicmodels); packageVersion("topicmodels")      #version: 0.2.14
library(tidytext); packageVersion("tidytext")          #version: 0.4.1
library(ldatuning); packageVersion("ldatuning")      #version: 1.0.2
library(cowplot); packageVersion("cowplot")          #version: 1.1.1
library(ape); packageVersion("ape")                  #version: 5.7.1

#### chen wilcox ####
otu_mat<- read.csv('chen_asv.csv', row.names = 1)
tax_mat<- read.csv('chen_taxa_table.csv', row.names = 1)
samples_df <-read.csv('chen_metadata.csv', row.names = 1)

#### transform otu and taxa to matrices #### 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#### phyloseq object creation #### 
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
ps <- phyloseq(OTU, TAX, samples)
ps

# filter low prevalence
minTotRelAbun <- 1e-5           
x <- taxa_sums(ps)
keepTaxa <- (x / sum(x)) > minTotRelAbun
ps2 <- prune_taxa(keepTaxa, ps)

ps2_genus <- tax_glom(ps2, taxrank = 'Genus')
ps2_genus

# total sum scaling #
ps2_genus <- transform_sample_counts(ps2_genus, function(ASV) ASV/sum(ASV) *1000000)
sample_sums(ps2_genus)

uni_data <- ps2_genus

data_otu_filt <- data.frame(otu_table(uni_data))
asv_data <- data.frame(phyloseq::tax_table(uni_data))
asv_data <- subset(asv_data, select=Genus)
wilcox_data <- merge(data_otu_filt,asv_data, by=0)

wilcox_data <- wilcox_data[,order(names(wilcox_data))]

# hold all data
wilcox_results <- list()
p_values <- data.frame(p=NA)

# to confirm range is correct:
names(wilcox_data[1,2:31]) 
names(wilcox_data[1,33:68])

for (ASV in 1:nrow(wilcox_data)) {
  wilcox_results[ASV] <- wilcox.test(as.numeric(unlist(wilcox_data[ASV,2:31])), 
                                     as.numeric(unlist(wilcox_data[ASV,33:68])), paired = F)$p.val
}

# calculate p value
for (p in 1:length(wilcox_results)) {
  p_values[[p]] <- wilcox_results[[p]]
}
Pval <- pivot_longer(p_values,cols = 1:175)
wilcox_data$p <- Pval$value 

# calculate q value
wilcox_data$q <- p.adjust(wilcox_data$p, method = 'fdr') # false discovery rate

# filter for significant
# !! choose q-value for filtering:
qval <- .25
ps_keepers <- which(wilcox_data$q <= qval)
ps_sig <- wilcox_data[ps_keepers,]
ps_sig1 <- ps_sig[,c(32,69,70,71)]
ps_sig2 <- ps_sig[,-c(32,69,70,71)]

# to confirm range is correct:
names(ps_sig2[1,2:31]) # MS
names(ps_sig2[1,32:67]) # HC


for(i in 1:nrow(ps_sig2)){
  
  # mean
  ps_sig2$MS_mean[i] <- mean(as.numeric(unlist(ps_sig2[i,2:31])))
  # sd
  ps_sig2$MS_sd[i] <- sd(as.numeric(unlist(ps_sig2[i,2:31])))
  
  # HC mean
  ps_sig2$HC_mean[i] <- mean(as.numeric(unlist(ps_sig2[i,32:67])))
  # HC sd
  ps_sig2$HC_sd[i] <- sd(as.numeric(unlist(ps_sig2[i,32:67])))
  
  
}

ps_sig3 <- ps_sig2[,68:71]

ps_sig4 <- cbind(ps_sig1,ps_sig3)

# write.csv(x = ps_sig4, file = "wilcox_chen_summary.csv")


#### functional analysis ####

otu_mat<- read.csv('CHEN_path_abun_unstrat_translated.csv', row.names = 1)
tax_mat<- read.csv('chen_path_table.csv', row.names = 1)
samples_df <-read.csv('chen_metadata.csv', row.names = 1)

#### transform otu and taxa to matrices #### 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#### phyloseq object creation #### 
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
ps <- phyloseq(OTU, TAX, samples)
ps

# filter low prevalence
minTotRelAbun <- 1e-5           
x <- taxa_sums(ps)
keepTaxa <- (x / sum(x)) > minTotRelAbun
ps2 <- prune_taxa(keepTaxa, ps)

# total sum scaling #
ps2 <- transform_sample_counts(ps2, function(ASV) ASV/sum(ASV) *1000000)
sample_sums(ps2)

uni_data <- ps2

data_otu_filt <- data.frame(otu_table(uni_data))
asv_data <- data.frame(phyloseq::tax_table(uni_data))
asv_data <- subset(asv_data, select=pathway)
wilcox_data <- merge(data_otu_filt,asv_data, by=0)

wilcox_data <- wilcox_data[,order(names(wilcox_data))]

# hold all data
wilcox_results <- list()
p_values <- data.frame(p=NA)

# to confirm range is correct:
names(wilcox_data[1,2:31]) 
names(wilcox_data[1,32:67])

for (ASV in 1:nrow(wilcox_data)) {
  wilcox_results[ASV] <- wilcox.test(as.numeric(unlist(wilcox_data[ASV,2:31])), 
                                     as.numeric(unlist(wilcox_data[ASV,32:67])), paired = F)$p.val
}

# calculate p value
for (p in 1:length(wilcox_results)) {
  p_values[[p]] <- wilcox_results[[p]]
}
Pval <- pivot_longer(p_values,cols = 1:314)
wilcox_data$p <- Pval$value 

# calculate q value
wilcox_data$q <- p.adjust(wilcox_data$p, method = 'fdr') # false discovery rate

# filter for significant
# !! choose q-value for filtering:
qval <- .25
ps_keepers <- which(wilcox_data$q <= qval)
ps_sig <- wilcox_data[ps_keepers,]
ps_sig1 <- ps_sig[,c(69,70,71)]
ps_sig2 <- ps_sig[,-c(69,70,71)]

# to confirm range is correct:
names(ps_sig2[1,2:31]) # MS
names(ps_sig2[1,32:67]) # HC


for(i in 1:nrow(ps_sig2)){
  
  # mean
  ps_sig2$MS_mean[i] <- mean(as.numeric(unlist(ps_sig2[i,2:31])))
  # sd
  ps_sig2$MS_sd[i] <- sd(as.numeric(unlist(ps_sig2[i,2:31])))
  
  # HC mean
  ps_sig2$HC_mean[i] <- mean(as.numeric(unlist(ps_sig2[i,32:67])))
  # HC sd
  ps_sig2$HC_sd[i] <- sd(as.numeric(unlist(ps_sig2[i,32:67])))
  
  
}

ps_sig3 <- ps_sig2[,68:72]

ps_sig4 <- cbind(ps_sig1,ps_sig3)

# write.csv(x = ps_sig4, file = "wilcox_chen_pathway_summary.csv")
