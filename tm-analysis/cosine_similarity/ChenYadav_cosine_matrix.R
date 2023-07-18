# Cosing Similarity Matrix

library(tidyverse); packageVersion("tidyverse")    #version: 2.0.0
library(lsa); packageVersion("lsa")                #version: 0.73.3
library(pheatmap); packageVersion("pheatmap")      #version: 1.0.12

## creating combined term-topic-matrix (these matrices were made while performing TM on each dataset)
# read in topic-term matrices
yadav <- read.csv('yadav_topic_term_matrix.csv', row.names = 1)
chen <- read.csv('chen_topic_term_matrix.csv', row.names = 1)

yadav$term <- row.names(yadav)
chen$term <- row.names(chen)

# merge datasets, keeping all microbes (maybe compare with removal of zeros)
yadav_chen <- merge(yadav,chen, by = 'term', all = T)
yadav_chen[is.na(yadav_chen)] <- 0
head(yadav_chen)

# topic-term proability matrix
topics_df <- yadav_chen

# Get separate datasets for each study
yadav_df <- topics_df %>%
  dplyr::select(term, contains("yadav")) %>%
  arrange(term) %>%
  column_to_rownames(var = "term")

chen_df <- topics_df %>%
  dplyr::select(term, contains("chen")) %>%
  arrange(term) %>%
  column_to_rownames(var = "term")

# checking out the names of topics
names(yadav_df); names(chen_df)

# confirming the terms are the same
head(rownames(yadav_df), 5); head(rownames(chen_df), 5)
tail(rownames(yadav_df), 5); tail(rownames(chen_df), 5)

#Calculate cosine similarity matrix
cos_mat <- data.frame(outer(yadav_df, chen_df, Vectorize(cosine)))

# convert to tibble, add row identifier, and shape "long"
cos_mat2 <-
  cos_mat %>%
  as_tibble() %>%
  rownames_to_column("Yadav") %>%
  pivot_longer(-Yadav, names_to = "Chen", values_to = "value") %>%
  mutate(
    Yadav = factor(Yadav, levels = 1:30),
    Chen = factor(gsub("chen_", "", Chen), levels = 1:30)
  )

p1 <- ggplot(cos_mat2, aes(Yadav, Chen)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2)), size=5) +
  xlab("Yadav Topics")+
  ylab("Chen Topics")+
  scale_fill_gradient2(low = "blue",
                       mid="white", 
                       high = "red",
                       midpoint = 0.79, na.value = 'red',
                       limits=c(0.00,.8)) +
  theme(axis.title = element_text(face="bold",size=25),
        axis.text = element_text(size=22, face='bold'),
        legend.position = "none")

p1

ggsave(filename = "cosine_similarity_heatmap.png",p1,device="png",
       dpi=600, width=17,height=15)
