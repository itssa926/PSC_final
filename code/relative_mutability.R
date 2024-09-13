library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

#setwd("~/Documents/Project/code")

##load data
essentialornot <- read.csv("../data/flybase/new_essential_not.csv") #new essentiality annotation
big_table <- read.csv("../data/FINAL_TABLE_FIXED_CODING_LENGTH.csv")
stop_table <- read.csv("../data/per_stop_gained_df_corrected.csv")
combined_df_cds <-read.csv("../data/cds_distance_annotated_with_start.csv") #df compiled when calculate CDS length, contain genes with only one single transcript
stop_table_annotated <- read.csv("../data/stop_table_annotated.csv")

########################calculate mutation probability for each potential PSC##################
# Define codons and stop codons5
bases <- c("A", "T", "C", "G")
codons <- as.vector(outer(outer(bases, bases, paste0), bases, paste0))
stop_codons <- c("TAA", "TAG", "TGA")

# Identify potential stop codons (not actual stop codons)
potential_stop_codons <-  c("TAT", "TAC", "TTA", "TCA", "TGG", 
                            "CAA", "CAG", "CGA", "GAA", "GAG", "GGA", 
                            "AAA", "AAG", "AGA", "TGT", "TGC", "TTG", "TCG")


# Function to generate all codons one mutation away
generate_one_mutation_away <- function(codon) {
  mutations <- character()
  for (i in 1:3) {
    for (base in bases) {
      if (substr(codon, i, i) != base) {
        mutated_codon <- codon
        substr(mutated_codon, i, i) <- base
        mutations <- c(mutations, mutated_codon)
      }
    }
  }
  return(mutations)
}


# Define the mutation rates #from josh's proportion to the reference
mutation_rates <- list(
  "A>T" = 0.022854504, "A>C" = 0.009066047, "A>G" = 0.019325919,
  "T>A" =0.022854574, "T>C" = 0.019247999, "T>G" = 0.008989717,
  "C>A" = 0.018248178, "C>T" = 0.043948747, "C>G" = 0.009234053,
  "G>A" = 0.043813703 , "G>T" = 0.018228395, "G>C" = 0.009192993
) 

# Calculate the number of stop codons one or two mutations away and the probabilities
calculate_probabilities <- function(codon) {
  one_step_mutations <- generate_one_mutation_away(codon)
  
  # Count unique stop codons one mutation away
  stop_codons_one_step <- unique(one_step_mutations[one_step_mutations %in% stop_codons])
  stop_count_one_step <- length(stop_codons_one_step)
  
  # Calculate probability based on mutation rates
  prob_one_step <- 0
  
  # Calculate for one-step mutations
  for (i in 1:3) {
    original_base <- substr(codon, i, i)
    for (mutated_base in bases) {
      if (original_base != mutated_base) {
        mutated_codon <- codon
        substr(mutated_codon, i, i) <- mutated_base
        if (mutated_codon %in% stop_codons) {
          prob_one_step <- prob_one_step + mutation_rates[[paste(original_base, ">", mutated_base, sep = "")]] #sum of all possibility
        }
      }
    }
  }
  
  return(list(stop_count_one_step = stop_count_one_step, 
              probability_one_step = prob_one_step))
}

# Apply the function to each potential stop codon
results <- lapply(potential_stop_codons, calculate_probabilities)

# Extract the results into separate vectors
stop_counts <- sapply(results, `[[`, "stop_count_one_step")
probabilities <- sapply(results, `[[`, "probability_one_step")

# Create the data frame
results_df <- data.frame(
  codon = potential_stop_codons,
  stop_count_one_step = stop_counts,
  probability_one_step = probabilities
)
print(results_df)

#visualise it
# Create bar plot for probability_one_step
ggplot(results_df, aes(x = reorder(codon, -probability_one_step), y = probability_one_step)) +
  geom_bar(stat="identity", fill = "#69b3a2") +
  labs(title = "Probability of Mutating to Stop Codon (One Mutation Away)",
       x = "Codon",
       y = "Probability") +
  theme_minimal()

#################visualise mutation rate#######################
# Convert to a data frame
mutation_df <- data.frame(
  Mutation = names(mutation_rates),
  Rate = unlist(mutation_rates)
)

# Sort the data frame by Rate in descending order
mutation_df <- mutation_df[order(-mutation_df$Rate), ]

# Create the bar chart
ggplot(mutation_df, aes(x = reorder(Mutation, -Rate), y = Rate)) +
  geom_bar(stat = "identity", fill = "#1f78b4", color = "black", width = 0.7) +
  labs(x = "Mutation Type", y = "Singleton Mutation Propotion", title = " Singleton Mutation Propotion to Reference Genome by Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")) +
  theme_minimal(base_size = 14)+
  geom_text(aes(label = round(Rate, 4)), vjust = -0.5, size = 3.5)+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = c(0.25, 0.8))# Add data labels on top of the bars


############################empirical count########################################################
################plot of no. potential stop codon mutate to stop codon######
#annotate stop table df
stop_table$id <-substring(stop_table$transcript, 1, nchar(stop_table$transcript) - 3)
stop_table_annotated <- merge(stop_table, essentialornot, by='id') #9440 genes

# Convert ref_codon and potential_codon to uppercase
stop_table_annotated$ref_codon <- toupper(stop_table_annotated$ref_codon)
stop_table_annotated$alt_codon <- toupper(stop_table_annotated$alt_codon)

#separate essential/non-essential
essential_stop_table <- subset(stop_table_annotated, stop_table_annotated$Essentiality=="Essential") #1920
Non_essential_stop_table <-subset(stop_table_annotated, stop_table_annotated$Essentiality=="Non-essential")#7059

#find unique potential stop
stop_table_unique <- stop_table_annotated %>%
  distinct(transcript, POS, ref_codon, .keep_all = TRUE) ### ensure no repeat potential stop at the same position of same gene


############### Group by ref_codon and alt_codon and count the occurrences of each type of PSC mutation############
stop_table_gam <- subset(stop_table_annotated, stop_table_annotated$count_gam != 0)
#PSC mutation count
codon_combinations <- stop_table_annotated %>%
  dplyr::group_by(ref_codon, alt_codon, Essentiality) %>%
  dplyr::summarise(count = n(), .groups = 'drop')

# Summing the counts for each reference codon
combined_codon_data <- codon_combinations %>%
  group_by(ref_codon,Essentiality) %>%
  summarise(total_count = sum(count))

## Calculate the overall total count#
#that is all the count of potential stop codon in essential/non-essential
overall_total_count <- combined_codon_data %>%
  group_by(Essentiality) %>%
  summarise(overall_total_count = sum(total_count))
# Merge the overall total count back into the combined_codon_data
combined_codon_data <- combined_codon_data %>%
  left_join(overall_total_count, by = "Essentiality")

# Calculate the relative counts
combined_codon_data <- combined_codon_data %>%
  group_by(Essentiality) %>%
  mutate(relative_count = total_count / overall_total_count)

# Print the updated data frame
print(combined_codon_data) # this give the proportion of potential PSC actually mutated in essential & non-essential genes

######################do the same for calculated probability#####################
# Calculate the overall total probability
overall_probability <- sum(results_df$stop_probability)
# Calculate the relative probability
results_df <- results_df %>%
  mutate(relative_prob = stop_probability / overall_probability) %>%
  arrange(desc(relative_prob))

print(results_df)

# Ensure the codon columns have the same name#
results_df <- rename(results_df, ref_codon = codon)

# Merge the data frames
merged_prob_df <- merge(results_df, combined_codon_data, by = "ref_codon")

# Print the merged data frame
print(merged_prob_df)


###########################get the no. of potential stop for each gene (include not mutated)#####################
#count the occurence for each potential stop in each transcript
potential_count_trans <- stop_table_unique %>%
  group_by(ref_codon, transcript, Essentiality) %>%
  summarise(count = n(), .groups = 'drop')

#put mutational weight (mutation probability) on the potential - get relative mutativity to stop
potential_count_weight <- merge(potential_count_trans, results_df, by = "ref_codon")
# Calculate the weighted counts
potential_count_weight$rela_weighted_count <- potential_count_weight$count * potential_count_weight$relative_prob
# View the adjusted dataframe
print(potential_count_weight)
write_csv(potential_count_weight, "~/Documents/Project/data/potential_psc_count_with_mutation.csv")

# Sum the weighted counts for each transcript
weighted_count_trans <- potential_count_weight %>%
  group_by(transcript) %>%
  summarize(rela_total_weighted_count = sum(rela_weighted_count))

#merge with the big df
combined_df_cds_w_mutation <- merge(combined_df_cds_w_mutation, weighted_count_trans, by.x = "transcript", by.y = "transcript")
write_csv(combined_df_cds_w_mutation, "~/Documents/Project/data/cds_distance_annotated_with_mutation.csv")

#####################get the observed count for each gene's potential stop codon mutation##################
stop_table_gam <- subset(stop_table_annotated, stop_table_annotated$count_gam != 0)
codon_combinations <- stop_table_gam %>%
  group_by(ref_codon, alt_codon) %>%
  summarise(count = n(), .groups = 'drop')

# Summing the counts for each reference codon
combined_codon_data <- codon_combinations %>%
  dplyr::group_by(ref_codon) %>%
  dplyr::summarise(total_count = sum(count))

# Calculate the overall total count
#that is all the count of potential stop codon in essential/non-essential
overall_total_count <- stop_table_annotated %>%
  dplyr::group_by(ref_codon) %>%
  dplyr::summarise(overall_total_count = n(), .groups = 'drop')
# Merge the overall total count back into the combined_codon_data
combined_codon_data<-combined_codon_data %>%
  left_join(overall_total_count, by = "ref_codon")

# Calculate the relative counts
combined_codon_data <- combined_codon_data %>%
  mutate(relative_count = ifelse(overall_total_count > 0, total_count / overall_total_count, 0))

# Calculate normalized relative counts across all data
total_relative_count <- sum(combined_codon_data$relative_count)
combined_codon_data <- combined_codon_data %>%
  mutate(normalized_relative_count = relative_count / total_relative_count)


#####################get the mosquito count for each potential stop codon mutation (Empirical mutation probability)##################
#group by ref and alt to get the mosquito count in total for the overall genome
codon_combinations_mos <- stop_table_annotated %>%
  dplyr::group_by(ref_codon, alt_codon) %>%
  dplyr::summarise(count = sum(count_gam), .groups = 'drop')

#get overall mosquito count
total_mut_count <-codon_combinations_mos %>%
  dplyr::summarise(total_count = sum(count), .groups = 'drop')
# 
# #get mutation rate -- number of mosquito has the mutation / overall mosquito count
codon_combinations_mos$total_mut_count <- total_mut_count$total_count
codon_combinations_mos$relative_prob <- codon_combinations_mos$count/codon_combinations_mos$total_mut_count 
# 
# # Sum the counts for each reference codon
summarized_data_mos <- codon_combinations_mos %>%
  dplyr::group_by(ref_codon) %>%
  dplyr::summarise(
    total_count = sum(count),
    total_mut_count = first(total_mut_count),
    relative_prob = sum(count) / first(total_mut_count)
  )

# Merge the data frames
results_df <- dplyr::rename(results_df, ref_codon=codon)
summarized_data_mos <- merge(results_df,codon_combinations_mos, by ='ref_codon') 
summarized_data_mos <- dplyr::rename(summarized_data_mos,  relative_sim_prob=relative_prob.x)
summarized_data_mos <- dplyr::rename(summarized_data_mos,  relative_prob=relative_prob.y)
summarized_data_mos

#summarise it
summarized_data_mos <- summarized_data_mos %>%
  dplyr::group_by(ref_codon) %>%
  dplyr::summarise(
    stop_count_one_step = mean(stop_count_one_step),  # or you might want to use sum, depending on context
    relative_prob = sum(relative_prob),
    relative_sim_prob = mean(relative_sim_prob),  # Assuming averaging probability makes sense
    total_count = sum(count),  # Summing the count if it makes sense
    .groups = "drop"  # Drops the grouping structure from the result
  )

#get long data
long_prob_mos <- summarized_data_mos %>%
  dplyr::select(ref_codon,relative_prob, relative_sim_prob) %>%
  pivot_longer(cols = c("relative_prob", "relative_sim_prob"), names_to = "source", values_to = "value")

# Filter out the empirical probabilities
empirical_probs <- long_prob_mos %>%
  filter(source == "relative_prob")

#################################just the stop count##############
summarized_data <- merge(results_df,combined_codon_data, by ='ref_codon') 

# Order codons by empirical probability in descending order, ensure they are unique
ordered_codons <- empirical_probs %>%
  arrange(desc(value)) %>%
  pull(ref_codon)%>%
  unique()

# Reorder ref_codon levels based on empirical count
long_prob_mos$ref_codon <- factor(long_prob_mos$ref_codon, levels =ordered_codons)

rownames(summarized_data_mos) <- summarized_data_mos$ref_codon

##############scatter plot#######################
ggplot(summarized_data_mos, aes(x = relative_prob, y = relative_sim_prob)) +
  geom_point(color="red") +
  labs(#title = "Relative Probabilities of Potnetial Stop Codon Mutating to Actual Stop Codon",
    #subtitle =  "Emprical Probability vs. Calculated Probability",
    x = "Relative Emperical Mutation Probability",
    y = "Calculated Relative Mutability") +
  geom_text_repel(
    aes(label = rownames(summarized_data_mos)),
    family = "Poppins",
    size = 4,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .004,
    nudge_y = .009,
    color = "grey50"
  )+
  theme_minimal()+
  theme(# Axis lines are now lighter than default
    axis.line = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 15),
    plot.subtitle = element_text(size = 15))+
  geom_abline(slope =1, intercept = 0, linetype = "dashed", color = "red")+
  ylim(0,NA)+
  xlim(0,NA)

#correlation tests
cor.test(summarized_data_mos$relative_prob,summarized_data_mos$relative_sim_prob, method = "spearman")
