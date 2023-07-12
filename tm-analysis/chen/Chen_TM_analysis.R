# Chen Topic Modeling Analysis

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


#### fit lda model ####
ms_vs_control <- ps2_genus
count_matrix <- data.frame(t(data.frame(otu_table(ms_vs_control))))

result <- FindTopicsNumber(
  count_matrix,
  topics = seq(from = 3, to = 50, by = 3), # trying anywhere from 3 to 50 topics, by increments of 3
  metrics = c("CaoJuan2009", "Arun2010"),
  method = "VEM",
  control = list(seed = 243),
  mc.cores = 4,
  verbose = TRUE
)

my_plot <- function (values)
{
  if ("LDA_model" %in% names(values)) {
    values <- values[!names(values) %in% c("LDA_model")]
  }
  columns <- base::subset(values, select = 2:ncol(values))
  values <- base::data.frame(values["topics"], base::apply(columns, 
                                                           2, function(column) {
                                                             scales::rescale(column, to = c(0, 1), from = range(column))
                                                           }))
  values <- reshape2::melt(values, id.vars = "topics", na.rm = TRUE)
  values$group <- values$variable %in% c("Griffiths2004", "Deveaud2014")
  values$group <- base::factor(values$group, levels = c(FALSE, 
                                                        TRUE), labels = c("Minimize", "Maximize"))
  p <- ggplot(values, aes_string(x = "topics", y = "value", 
                                 group = "variable"))
  p <- p + geom_line()
  p <- p + geom_point(aes_string(shape = "variable"), size = 7)
  p <- p + guides(size = FALSE, shape = guide_legend(title = "Metrics:"))
  p <- p + scale_x_continuous(breaks = values$topics)
  p <- p + labs(x = "\nNumber of Topics", y = NULL)
  p <- p + facet_grid(group ~ .)
  p <- p + theme_bw() %+replace% theme(panel.grid.major.y = element_blank(), 
                                       panel.grid.minor.y = element_blank(), panel.grid.major.x = element_line(colour = "grey70"), 
                                       panel.grid.minor.x = element_blank(), legend.key = element_blank(), 
                                       strip.text.y = element_text(angle = 90, size = 25),
                                       legend.position = "bottom")
  return(p)
}

(p1 <- my_plot(result))

p1.5 <- p1 + 
              ggtitle("Chen et. al. 2016") +
              theme(plot.title = element_text(size=28,face='bold', hjust = 0.5),
                    axis.text = element_text(size=20, face="bold"),
                   axis.title = element_text(size=25,face='bold'),
                   legend.text = element_text(size=25, face ='bold'),
                   legend.title=element_text(size=25,face='bold'))
p1.5

# 30 topics
lda_k30 <- LDA(count_matrix, k = 30, method = "VEM", control = list(seed = 243))

# beta = per-topic-per-word probabilities
b_df <- data.frame(tidy(lda_k30, matrix = "beta"))

# gamma = per-document-per-topic probabilities
g_df <- data.frame(tidy(lda_k30, matrix = "gamma")) %>%
  arrange(document, topic)

# notice fractional membership 
head(b_df)
head(g_df)

# build topic model ps object
lib_size_df <- data.frame(sample_sums(ms_vs_control)) %>%
  dplyr::rename("read_count" = "sample_sums.ms_vs_control.") %>%
  rownames_to_column(var = "document")

tm_df <- left_join(lib_size_df, g_df) %>%
  mutate(topic_count = read_count * gamma,
         topic_count = round(topic_count, 0)) %>%
  dplyr::select(-read_count, -gamma) %>%
  pivot_wider(names_from = topic, values_from = topic_count) %>%
  dplyr::rename_with(~ paste0("Topic_", .), -document) %>%
  column_to_rownames(var = "document") %>%
  t(.) %>%
  data.frame(.)

write.csv(tm_df,"Chen_document_topic_mat.csv")

(ps_topic_g <- phyloseq(
  sample_data(ms_vs_control),
  otu_table(tm_df, taxa_are_rows = TRUE))) 


#### fitting DA model to the topics
tm_linda <- linda(phyloseq.obj = ps_topic_g,
                  formula = '~ Disease_state2', 
                  feature.dat.type = "count",
                  is.winsor = F, p.adj.method = "BH", alpha = 0.05)

linda_rrms <- data.frame(tm_linda$output$Disease_state2MS) %>%
  dplyr::arrange(pvalue) %>%
  rownames_to_column(var = "Topic")
head(linda_rrms)

q=0.25
p=0.05

fdr_linda <- linda_rrms %>%
  mutate(Reject = ifelse(padj < q & pvalue <= p, "Yes", "No")) %>%
  separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
  mutate(Topic = gsub("_", " ", Topic)) %>%
  arrange(desc(t_n)) %>%
  mutate(Topic = factor(Topic, levels = Topic))
head(fdr_linda)
sum(fdr_linda$Reject=="Yes")

write.csv(fdr_linda, "fdr_linda_rrmsVShc.csv")

(p2 <- ggplot(data = fdr_linda, aes(x = Topic, y = log2FoldChange, fill = Reject)) +
    geom_col(alpha = .8) +
    ggtitle("MS vs HC (Chen)", subtitle =  paste0("p ≤ ",p," & q ≤ ", q)) +
    labs(y = "\nLog2 Fold-Change", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 30),
          axis.text.y = element_text(color = "black", size = 20, face="bold"),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 35),
          plot.title = element_text(size = 40, face = "bold"),
          plot.subtitle = element_text(size = 35),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20)) +
    coord_flip() +
    scale_fill_manual(values = c("grey", "#5317D5")) +
    geom_hline(yintercept = 0, linetype="dotted") +
    theme(legend.position = "none"))

sig_topics <- fdr_linda[fdr_linda$Reject=="Yes",]$t_n
sig_topics

tax <- as.data.frame(tax_mat)
tax$term <- row.names(tax)

for (t in 1:length(sig_topics)) {
  k=sig_topics[t]
  ##
  
  b_plot <- b_df %>%
    dplyr::filter(topic == k) %>%
    arrange(desc(beta)) %>%
    slice_head(n = 20) %>%
    arrange(desc(term)) %>%
    mutate(term = factor(term, levels = term))
  
  b_plot2 <- merge(tax, b_plot, by = "term")
  
  # grab info about topic
  fdr_linda[fdr_linda$t_n == k,]
  lfc_k <- abs(round(fdr_linda[fdr_linda$t_n == k,]$log2FoldChange,3))
  
  if (fdr_linda[fdr_linda$t_n == k,]$log2FoldChange > 0) {
    r = "MS"
  } else {
    r = "HC"
  }
  
  p <- round(fdr_linda[fdr_linda$t_n == k,]$pvalue,5)
  q <- round(fdr_linda[fdr_linda$t_n == k,]$padj,5)
  
  
  
  (p3 <- ggplot(data = b_plot2, aes(x = beta, y = Genus, color = term)) +
      geom_point(aes(stroke= 8)) +
      ggtitle(paste0("Community Type ",k, " (Chen)"), subtitle = paste0("Log2FC ", lfc_k, " in ", r, " (p=",p," q=", q,")")) +
      labs(x = paste0("\nTopic ", k,": Genus-Topic Probabilities (beta)"), y = "") +
      theme(legend.position = "none",
            axis.text.x = element_text(size=23, hjust = 1),
            axis.text.y = element_text(size=20),
            axis.title = element_text(size=25),
            plot.title = element_text(size=30, face='bold'),
            plot.subtitle = element_text(size=25, face = 'bold'))
  )
  
  print(p3)
  ggsave(paste0("topic_",k,"_msVShc.png"), p3, device = "png", dpi = 600,
         width = 14, height = 12)
}

## differential abundance analysis within topics ##
linda_g_all <- linda(phyloseq.obj = ms_vs_control,
                     formula = '~ Disease_state2', 
                     feature.dat.type = "count",
                     is.winsor = T, p.adj.method = "BH", alpha = 0.05)

linda_res_g_df <- data.frame(linda_g_all$output$Disease_state2MS) %>%
  dplyr::arrange(padj) %>%
  rownames_to_column(var = "term")

## Print Differential Abundance and GGSAVE
for (t in 1:length(sig_topics)) {
  k=sig_topics[t]
  ##
  
  b_plot <- b_df %>%
    dplyr::filter(topic == k) %>%
    arrange(desc(beta)) %>%
    slice_head(n = 20) %>%
    arrange(desc(term)) %>%
    mutate(term = factor(term, levels = term))
  
  keep_g <- b_plot$term
  
  log2_df <- linda_res_g_df %>%
    filter(term %in% keep_g) %>%
    arrange(desc(term)) %>%
    mutate(term = factor(term, levels = term))
  
  log2_df2 <- merge(b_plot, log2_df, by = "term")
  
  sig_altered <- log2_df2 %>%
    mutate(reject = ifelse(padj < 0.25, "Yes", "No"))
  
  sig_altered <- merge(tax, sig_altered, by = "term")
  sig_altered <- sig_altered %>% 
    arrange(Genus)
  
  sig_altered
  sum(sig_altered$reject=="Yes")
  
  bold.genus <- sig_altered[sig_altered$reject=="Yes",]$term
  bold.genus
  
  bold.labels <- ifelse(sig_altered$term %in% bold.genus, yes = "bold", no = "plain")
  
  (p4 <- ggplot(data = sig_altered, aes(x = reorder(Genus, Genus), y = log2FoldChange, fill = term)) +
      geom_col(alpha = .8) +
      labs(y = "\nLog2 Fold-Change", x = "") +
      theme(legend.position = "none",
            axis.text.x = element_text(size=23, hjust = 1),
            axis.text.y = element_text(size=20, face = bold.labels),
            axis.title = element_text(size=25),
            plot.title = element_text(size=30, face = 'bold'),
            plot.subtitle = element_text(size=25, face = 'bold')) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype="dotted") +
      theme(legend.position = "none") +
      ggtitle("MS vs HC (Chen)", 
              subtitle = paste0("Topic: ",k, " (Bolded: p < .05 and q < 0.25)")))
  
  # checking the right ones were bolded
  bold.genus
  
  print(p4)
  
  ggsave(paste0("topic_",k,"_msVShc_l2fc.png"), p4, device = "png", dpi = 600,
         width = 14, height = 12)
}


## extract topic term_matrix 
topic_term_matrix <- b_df %>%  pivot_wider(names_from = topic,
                      values_from = beta)
colnames(topic_term_matrix) <- paste0("chen_",colnames(topic_term_matrix))
write.csv(topic_term_matrix,"chen_topic_term_matrix.csv", row.names = F)

# correlations for network #
ms_genus <- subset_samples(ps2_genus,Disease_state2=="MS")
cor_genus <- as.data.frame(otu_table(ms_genus))
cor_genus$term <- row.names(cor_genus)
cor_genus_2 <- merge(tax,cor_genus,by="term")
cor_genus_2 <- cor_genus_2[,c(7,9:length(cor_genus_2))]
rownames(cor_genus_2) <- cor_genus_2$Genus
cor_genus_2 <- cor_genus_2[,-1]
cor_genus_mat <- (as.matrix(cor_genus_2))

cor_df <- data.frame()
for (i in 1:nrow(cor_genus_mat)) {
  for (j in 1:nrow(cor_genus_mat)) {
    cor.tmp <- cor.test(cor_genus_mat[i,], cor_genus_mat[j,])
    cor_sum <- data.frame("R"=cor.tmp$estimate,"P"=cor.tmp$p.value,
                          "microbe1"=row.names(cor_genus_2[i,]),
                          "microbe2"=row.names(cor_genus_2[j,]))
    cor_df <- rbind(cor_df,cor_sum)
  }
  
}
head(cor_df)

cor_df_clean <- cor_df %>% 
  filter(!duplicated(paste0(pmax(microbe1, microbe2), pmin(microbe1, microbe2))))

write.csv(cor_df_clean,"chen_correlation.csv")
