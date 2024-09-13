library(tidyverse)
library(readr)

#######load data##########
ortho_dro<-read.csv("~/Documents/Project/data/ortho_df_all_withallele.csv") #orthologous from drosophila genes to mosquito genes
ortho_mos<-read.csv("~/Documents/Project/data/ortho_df_all_withallele_mos.csv") #orthologous from  mosquito genes to drosophila genes

###################find unique orthologous pairs#########################
unique_pair_dro<- unique(ortho_dro[,c('id','Flybase_ID')])
unique_pair_mos <- unique(ortho_mos[,c('id','Flybase_ID')])

unique_pair_all <- rbind(unique_pair_dro,unique_pair_mos)
unique_pairs <-unique(unique_pair_all[,c('id','Flybase_ID')]) #22,696 unique orthologous pairs

write_csv(unique_pairs,"~/Documents/Project/data/all_ortho.csv")#df with all unique orthologous pairs

################Count occurrences of each Mosquito_ID and Fly_ID###########
mosquito_counts <- unique_pairs %>% # count the occurence of each mosquito gene
  group_by(id) %>%
  summarise(mosquito_count = n())

fly_counts <- unique_pairs %>% # count the occurence of each drosophila gene
  group_by(Flybase_ID) %>%
  summarise(fly_count = n())

unique_pairs <- unique_pairs %>% ## merge them
  left_join(mosquito_counts, by = "id") %>%
  left_join(fly_counts, by = "Flybase_ID")


######## Classify the orthologous relationships#############
unique_pairs <- unique_pairs %>%
  mutate(relationship = case_when(
    mosquito_count == 1 & fly_count == 1 ~ "1:1",
    mosquito_count == 1 & fly_count > 1 ~ "1:m", #one mosquito gene has multiple drosophila orthologs
    mosquito_count > 1 & fly_count == 1 ~ "m:1", #many mosquito genes have one drosophila orthologous
    mosquito_count > 1 & fly_count > 1 ~ "m:m"
  ))

write_csv(unique_pairs, "~/Documents/Project/data/ortholog_pair_wrelation.csv")
# Extract the specific relationships
one_to_one <- unique_pairs %>% filter(relationship == "1:1") #5457
one_to_many <- unique_pairs %>% filter(relationship == "1:m") #1161
many_to_many <- unique_pairs %>% filter(relationship == "m:m") #14144
many_to_one <- unique_pairs %>% filter(relationship == "m:1") #1932

#count how many unique mosuito id in many_one and many_many
unique_mosquito_ids_many_to_many <- unique(many_to_many$id) #2120
unique_mosquito_ids_many_to_one <- unique(many_to_one$id) #702
