# Yadav networks
# Rachel L Shrode
# June 26 2023

library(GGally); packageVersion("GGally")    #version: 2.1.2
library(ggnet); packageVersion("ggnet")      #version: 0.1.0
library(network); packageVersion("network")   #version: 1.18.1
library(sna); packageVersion("sna")           #version: 2.7.1
library(ggplot2); packageVersion("ggplot2")   #version: 3.4.2
library(tidyr); packageVersion("tidyr")      #version: 1.3.0
library(dplyr); packageVersion("dplyr")      #version: 1.1.2

# positive and negative correlations (8,9)
k=8
b_plot <- b_df %>%
  dplyr::filter(topic == k) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 20) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))

b_plot2 <- merge(tax, b_plot, by = "term")
head(b_plot2)

keep_g <- b_plot2$term

log2_df <- linda_res_g_df %>%
  filter(term %in% keep_g) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))

log2_df2 <- merge(b_plot2, log2_df, by = "term")

log2_df2 <- log2_df2 %>% 
  mutate(higher_in = ifelse(log2_df2$log2FoldChange>=0,"> MS","< MS"))
col = c("> MS" = "red", "< MS" = "blue")

log2_df2 <- log2_df2 %>% 
  mutate(color = ifelse(higher_in=="> MS","red","blue"))

# pull out correlation data that goes with topic k #
cor_net_data <- data.frame()

for (i in 1:nrow(cor_df_clean)) {
  if (cor_df_clean[i,3] %in% log2_df2$Genus & cor_df_clean[i,4] %in% log2_df2$Genus) {
    cor.add <- cor_df_clean[i,]
    cor_net_data <- rbind(cor_net_data,cor.add)
  }
}

head(cor_net_data)

# remove rows with NA
cor_net_data <- na.omit(cor_net_data)

# checking significance of correlations
# reducing all nonsignificant to correlation of 0 to remove from plot
for (i in 1:nrow(cor_net_data)) {
  if(cor_net_data[i,2]>0.05){
    cor_net_data[i,1]=0
  }
}

head(cor_net_data)

a <- pivot_wider(cor_net_data[-2], 
                 names_from = microbe1, 
                 values_from = R)
microbe2_names <- a$microbe2
a <- a[-1]
a <- as.data.frame(a)
a[is.na(a)] <- t(a)[is.na(a)]
rownames(a) <- microbe2_names
a <- as.matrix(round(a,2)) # cannot have negatives edges, will supply pos/neg later
a[is.na(a)] <- 0

# build network with correlations as edges
net3 = rgraph(nrow(a), mode = "graph", tprob=a)
net3 <- network(a)
net3 %v% "Abundance MS vs HC" = log2_df2$higher_in
net3 %v% "LFC MS vs HC" = (round(abs(log2_df2$log2FoldChange),1))
net3 %v% "Abun2" = log2_df2$color
net3 %v% "Probability of Assignment" = round((b_plot2$beta),3)
net3 %e% "Correlation" = a
net3 %e% "Correlation2" = abs(a)*2 # higher correlation thicker value
net3 %e% "Edge_Cor_Color" = ifelse(net3 %e% "Correlation" > 0, "red", "blue")

ggnet2(net3, label=T,
       label.color = "Abun2",
       label.size = 7,
       mode = "circle",
       node.color = "Abundance MS vs HC",
       node.size = "LFC MS vs HC",
       node.alpha = .5,
       palette = col,
       edge.alpha = .2,
       edge.size = "Correlation2",
       edge.color = "Edge_Cor_Color",
       edge.label.color = "black",
       legend.size = 15,
       layout.exp = .35)

ggsave(paste0("yadav_network_topic_lfc",k,".png"),plot = last_plot(),
       device = 'png',
       height = 12,
       width = 14,
       dpi = 600)