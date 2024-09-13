library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)

##load data
essentialornot <- read.csv("../data/flybase/new_essential_not.csv") #new essentiality annotation
big_table <- read.csv("../data/FINAL_TABLE_FIXED_CODING_LENGTH.csv")
stop_table <- read.csv("../data/per_stop_gained_df_corrected.csv")
feature_df <- read.csv("../data/genome_features.csv")

####################calculate distance of PSC to the end of CDS####################
#####get transcript that only occurs once########
big_table$ID <- sub("-.*", "", big_table$ID)
big_table_unique <- big_table[!duplicated(big_table$ID), ]
big_table_unique$ID <- paste0(big_table_unique$ID , "-RA")
essentialornot$transcript <- paste0(essentialornot$id , "-RA")

# Subsetting those single transcripts
stop_table_single <- stop_table[stop_table$transcript %in% big_table_unique$ID, ]#2194941 transcripts, 12054 genes
stop_table_new <- stop_table[stop_table$transcript %in% essentialornot$transcript, ] #8972 unique out lf 9440 are single transcripts
write.csv(stop_table_new, "../data/unique_stop_gained.csv")

stop_table_gam <- subset(stop_table, stop_table$count_gam != 0)#433611 transcript, 12018 genes
stop_table_gam_new <- stop_table_gam[stop_table_gam$transcript %in% essentialornot$transcript, ] #8951

######get the distance from the end of the gene########
feature_df_unique <- feature_df[feature_df$ID %in% big_table_unique$ID, ]

#get CDS for each single copy unique transcript
cds_df<- subset(feature_df, feature_df$type =='CDS')
cds_unique <- cds_df[cds_df$Parent %in% big_table_unique$ID, ] #13048 unique genes in mosquito

# Calculate total length of each CDS chunk
cds_unique <- cds_unique %>%
  mutate(chunk_length = end - start + 1) ## cds end position - cds start position +1

# Sum the lengths to get the total CDS length for each gene (Parent)
total_cds_length <- cds_unique %>% #sum each chunk
  group_by(Parent) %>%
  summarise(total_length = sum(chunk_length))

colnames(cds_unique)[which(names(cds_unique) == "Parent")] <- "transcript"
colnames(total_cds_length)[which(names(total_cds_length) == "Parent")] <- "transcript"

###get the start, end, and total length for each transcript###
cds_summary <- cds_unique %>%
  group_by(transcript) %>%
  summarise(
    cds_start = min(start),
    total_length = sum(chunk_length)
  ) %>%
  mutate(cds_end = cds_start + total_length - 1) # -1 important

write.csv(cds_summary, "../data/cds_summary.csv")

#calculate the distance of position of stop codon to the end of cds
stop_table_cds_end <- merge(x = stop_table_gam, y = cds_summary, by = "transcript",
                            all.x = TRUE)
# Calculate the nucleotide position of the stop codon within the CDS
stop_table_cds_end$gene_pos <- stop_table_cds_end$cds_start + (stop_table_cds_end$aa_pos * 3) ## start of cds + aa position x3

# Calculate the distance from the stop codon to the end of the CDS
stop_table_cds_end$distance_to_cds <- stop_table_cds_end$cds_end - stop_table_cds_end$gene_pos


#########################compare it between essential vs non-essential#############################
###annotate cds df
stop_table_cds_end$id <- sub("-RA.*", "", stop_table_cds_end$transcript)
stop_table_cds_end <- merge(x = stop_table_cds_end, y = essentialornot, by = "id",
                            all.x = TRUE)
write.csv(stop_table_cds_end, "../data/unique_stop_gained_cds_distance_annotated_new.csv")

###Comparison between essential and non-essential genes
essential_stop_cds <- subset(stop_table_cds_end, stop_table_cds_end$Essentiality == "Essential") #1918 transcript
non_essential_stop_cds <-subset(stop_table_cds_end, stop_table_cds_end$Essentiality == "Non-essential") #248258, 7033 transcript
#get the relative length to total CDS length for each gene
essential_stop_cds$relative_dis <- (essential_stop_cds$distance_to_cds)/essential_stop_cds$total_length
non_essential_stop_cds$relative_dis <- (non_essential_stop_cds$distance_to_cds)/non_essential_stop_cds$total_length

#merge essential and non-essential genes
combined_df_cds <- rbind(essential_stop_cds, non_essential_stop_cds) #8951 genes in total

#mean&median
mean(essential_stop_cds$total_length) #4602
mean(non_essential_stop_cds$total_length) #2669.9
median(essential_stop_cds$total_length) #2979
median(non_essential_stop_cds$total_length) #1884

write.csv(combined_df_cds, "../data/cds_distance_annotated_with_start.csv")



#################################visualisations#########################################
#Cumulative density plot
#relative distance
ggplot(combined_df_cds, aes(x = relative_dis, color = Essentiality)) +
  stat_ecdf(geom = "point") +
  labs(title = "CDF of Stop Codon Distance for Essential and Non-Essential Genes",
       y = "Cumulative Density", x = "Distance from End of CDS") +
  theme_minimal() +
  scale_color_manual(values = c("Essential" = "blue", "Non-essential" = "red")) +
  theme(legend.title = element_blank())

#absolute distance
ggplot(combined_df_cds, aes(x = distance_to_cds, color = Essentiality)) +
  stat_ecdf(geom = "point") +
  labs(title = "CDF of Stop Codon Distance for Essential and Non-Essential Genes",
       y = "Cumulative Density", x = "Distance from End of CDS") +
  theme_minimal() +
  scale_color_manual(values = c("Essential" = "blue", "Non-essential" = "red")) +
  theme(legend.title = element_blank())


# Plot histogram for absolute distance with essential genes on top
combined_df_cds$Essentiality <- factor(combined_df_cds$Essentiality, levels = c("Non-essential", "Essential"))

combined_df_cds %>%
  ggplot(aes(x = distance_to_cds, fill = Essentiality)) +
  geom_histogram(position = "identity", alpha = 0.8, binwidth = 50) +  
  #geom_density(aes(y = ..count.. * 50), alpha = 0.5, adjust = 1.5, color = "black") + # Adjust bin width for better visualization# Limiting x-axis to relevant range
  labs(title = "Distribution of Stop Codon relative to CDS End",
       subtitle = "Analysis of stop codon distances from the end of coding sequences in Anopheles gambiae",
       x = "Distance from CDS End (bp)",
       y = "Frequency") +
  theme_minimal(base_size = 14) +  # Increase base font size for better readability
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#FFA500"), 
                    labels = c("Non-essential Genes", "Essential Genes")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma) + 
  scale_x_continuous(limits = c(0, 10000), expand = expansion(mult = c(0, 0.05)), 
                     labels = scales::comma) +
  theme(legend.position = "top",
        legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(face = "bold", hjust= 0.5,size = 16),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.caption = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(reverse = TRUE))  # Reverse the legend order if necessary

#ks test
ks.test(essential_stop_cds$distance_to_cds, non_essential_stop_cds$distance_to_cds)
ks.test(essential_stop_cds$relative_dis, non_essential_stop_cds$relative_dis)

# Calculate the bin with the highest frequency
#from end
bin_width <- 100
combined_df_cds <- combined_df_cds %>%
  mutate(bin = floor(distance_to_cds / bin_width) * bin_width)

# Summarize counts by bin and essentiality
bin_summary <- combined_df_cds %>%
  group_by(bin, Essentiality) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Display the bin with the highest frequency for essential and non-essential genes
highest_frequency_essential <- bin_summary %>%
  filter(Essentiality == "Essential") %>%
  arrange(desc(count)) %>%
  slice(1)

highest_frequency_non_essential <- bin_summary %>%
  filter(Essentiality == "Non-essential") %>%
  arrange(desc(count)) %>%
  slice(1)

print(highest_frequency_essential) #0-100 bp
print(highest_frequency_non_essential) #0-100 bp


######plot the distribution fo cds######
combined_df_cds %>%
  # filter(distance_to_end<30000 ) %>% 
  ggplot(aes(x = total_length, fill = Essentiality)) +
  geom_histogram(position = "identity", alpha = 0.8, binwidth = 100) + 
  xlim(0,20000)+
  labs(title = "Histogram of CDS length",
       x = "CDS length",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("Essential" = "#ffa500", "Non-essential" = "#404080")) +
  theme(legend.title = element_blank())

#ks test
ks.test(essential_stop_cds$total_length, non_essential_stop_cds$total_length)

##############Visualise the distributions for separate chromosome arms########
#get two separate datasets for X and autosomes
x_combined_cds <- subset(combined_df_cds, combined_df_cds$contig=='X')
aut_combined_cds <- subset(combined_df_cds, combined_df_cds$contig !='X') # chromosome arms that are not X

#for distance to cds end
x_combined_cds %>%
  # filter(distance_to_end<30000 ) %>% 
  ggplot(aes(x = distance_to_cds, fill = Essentiality)) +
  geom_histogram(position = "identity", alpha = 0.8, binwidth = 50) + 
  xlim(0,10000)+
  labs(title = "Histogram of Stop Codon Distances in X Chromosome",
       x = "Distance from End of CDS",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("Essential" = "#ffa500", "Non-essential" = "#404080")) +
  theme(legend.title = element_blank())

aut_combined_cds %>%
  # filter(distance_to_end<30000 ) %>% 
  ggplot(aes(x = distance_to_cds, fill = Essentiality)) +
  geom_histogram(position = "identity", alpha = 0.8, binwidth = 50) + 
  xlim(0,10000)+
  labs(title = "Histogram of Stop Codon Distances in Autosomes",
       x = "Distance from End of CDS",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("Essential" = "#ffa500", "Non-essential" = "#404080")) +
  theme(legend.title = element_blank())

#cdf
ggplot(x_combined_cds, aes(x = distance_to_cds, color = Essentiality)) +
  stat_ecdf(geom = "point") +
  labs(title = "CDF of Stop Codon Distance for Essential and Non-Essential Genes in X contig",
       y = "Cumulative Density", x = "Distance from End of CDS") +
  theme_minimal() +
  scale_color_manual(values = c("Essential" = "blue", "Non-essential" = "red")) +
  theme(legend.title = element_blank())

ggplot(aut_combined_cds, aes(x = distance_to_cds, color = Essentiality)) +
  stat_ecdf(geom = "point") +
  labs(title = "CDF of Stop Codon Distance for Essential and Non-Essential Genes in Autosomes",
       y = "Cumulative Density", x = "Distance from End of CDS") +
  theme_minimal() +
  scale_color_manual(values = c("Essential" = "blue", "Non-essential" = "red")) +
  theme(legend.title = element_blank())

#for relative length 
x_combined_cds %>%
  # filter(distance_to_end<30000 ) %>% 
  ggplot(aes(x = relative_dis, fill = Essentiality)) +
  geom_histogram(position = "identity", alpha = 0.8, binwidth = 0.01) + 
  xlim(0,1)+
  labs(title = "Histogram of Stop Codon Distances in X Contig",
       x = "Relative Distance from End of CDS",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("Essential" = "#ffa500", "Non-essential" = "#404080")) +
  theme(legend.title = element_blank())

aut_combined_cds %>%
  # filter(distance_to_end<30000 ) %>% 
  ggplot(aes(x = relative_dis, fill = Essentiality)) +
  geom_histogram(position = "identity", alpha = 0.8, binwidth = 0.01) + 
  xlim(0,1)+
  labs(title = "Histogram of Stop Codon Distances in Autosomes",
       x = "Relative Distance from End of CDS",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("Essential" = "#ffa500", "Non-essential" = "#404080")) +
  theme(legend.title = element_blank())


#statistical test to see the difference
ks.test(x_combined_cds$distance_to_cds[x_combined_cds$Essentiality == 'Essential'],
        aut_combined_cds$distance_to_cds[aut_combined_cds$Essentiality == 'Essential'])

ks.test(x_combined_cds$distance_to_cds[x_combined_cds$essentiality == 'Essential'],
        x_combined_cds$distance_to_cds[x_combined_cds$essentiality == 'Non-essential'])

ks.test(aut_combined_cds$distance_to_cds[x_combined_cds$essentiality == 'Essential'],
        aut_combined_cds$distance_to_cds[x_combined_cds$essentiality == 'Non-essential'])

##get the number of samples for essential/non-essential
table(combined_df_cds$Essentiality)
