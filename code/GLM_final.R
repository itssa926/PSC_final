#####################load require packages#############
library(dplyr)
library(ggplot2)
library(MASS)
#effect plot
library(ggeffects)
library(patchwork)
library(ggsignif)
library(grid)
#visualise interactions
library(interactions)
library(effects)
require(ggiraph)
require(ggiraphExtra)
require(plyr)

###################load data#################################
big_df <- read.csv("../data/FINAL_TABLE_FIXED_CODING_LENGTH.csv")
combined_df_cds_w_mutation <- read.csv("../data/cds_distance_annotated_with_mutation.csv")
essentialornot <- read.csv("../data/flybase/new_essential_not.csv") #new essentiality annotation

############log-transformed##################
# Log-Transform data -- all +1 to avoid 0
combined_df_cds_w_mutation <- combined_df_cds_w_mutation %>%
  mutate(log_potential_stop_gained_mut = log1p(potential_stop_gained_mut),
         log_gam_sequences_stop_gained_count = log1p(gam_sequences_stop_gained_count),
         log_rela_potential_weighted = log1p(rela_total_weighted_count),
         log_coding_length = log1p(coding_length))

#############Classified chromosomes into autosome and X chromosome#################
combined_df_cds_w_mutation <- combined_df_cds_w_mutation %>%
  mutate(chromosome_group = case_when(
    contig.y %in% c("2L", "2R", "3R", "3L") ~ "Autosome",
    contig.y == "X" ~ "X",
    TRUE ~ contig.y  # Catch-all in case there are other categories
  ))

######################################GLM ON Actual PSC vs coding length (population-level)#######################
reduced_stop_model <-glm.nb(gam_sequences_stop_gained_count ~ Essentiality,data = combined_df_cds_w_mutation)
reduced_stop_model_2 <- glm.nb(gam_sequences_stop_gained_count ~ Essentiality + log_coding_length,data = combined_df_cds_w_mutation)
full_stop_model <- glm.nb(gam_sequences_stop_gained_count ~ Essentiality +log_coding_length+ contig.y, data = combined_df_cds_w_mutation)


# Perform the likelihood ratio test
anova(reduced_stop_model,reduced_stop_model_2,full_stop_model, test = "Chisq")

#diagnostic plot
plot(full_stop_model)

#see p value
summary(full_stop_model)

#####################visualise actual stop efect plot#############################
# Generate partial effect data
# Effect data for Essentiality
effect_data_essentiality <- ggpredict(full_stop_model, terms = "Essentiality")
# Effect data for coding sequence length
effect_data_coding_length <- ggpredict(full_stop_model, terms = "log_coding_length")
# Effect data for contig
effect_data_contig <- ggpredict(full_stop_model, terms = "contig.y")
# Reorder the contig levels
effect_data_contig$x <- factor(effect_data_contig$x, levels = c("2L", "2R", "3R", "3L", "X"))
# Significance labels based on the provided coefficients
significance_labels <- data.frame(
  x = c("2R", "3L", "3R", "X"),
  y = c(315, 315, 315, 315),  # Adjust y position as necessary
  label = c("***", "***", "***", "***")
)
# Calculate positions for the significance label in the coding length plot
x_mid <- mean(range(effect_data_coding_length$x))
y_top <- max(effect_data_coding_length$conf.high)*0.9  # Adjust multiplier for appropriate padding

# Add significance labels for coding length plot
significance_labels_coding_length <- data.frame(
  x = c(x_mid),
  y = c(y_top),
  label = c("***")
)

###### Generate the effect plots##
plot_essentiality <- ggplot(effect_data_essentiality, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.5)) +
  labs(#title = "Effect of Gene Essentiality on PSC Counts per Gene",
    x = "Essentiality",
    y = "PSC Counts per Gene") +
  theme_classic() +
  theme(legend.position = "none",
        #plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16.5)) +
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +
  geom_signif(data=effect_data_essentiality,stat="signif",position="identity",
              comparisons=list(c("Essential","Non-essential")),map_signif_level = TRUE,annotations="***", textsize = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  # Reducing the gap by setting no expansion at the lower end-


plot_coding_length <- ggplot(effect_data_coding_length, aes(x = x, y = predicted)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "red") +
  labs(#title = "Effect of Log(Coding Length) on PSC Counts per Gene",
    x = "Log(Coding Length) (bp)",
    y = "PSC Count per Gene") +
  geom_text(data = significance_labels_coding_length, aes(x = x, y = y, label = label), vjust = -0.5, size = 10) +  # Adjust size as neede
  theme_classic() +
  theme(#plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16.5))

plot_contig <- ggplot(effect_data_contig, aes(x = x, y = predicted)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "red", position = position_dodge(0.5)) +
  geom_point(color = "black", size = 1.5, position = position_dodge(0.5)) +
  labs(#title = "Effect of Gene Chromosomal Location on PSC Counts per Gene",
    x = "Gene Chromosomal Location",
    y = "PSC Count per Gene") +
  theme_classic() +
  geom_text(data = significance_labels, aes(x = x, y = y, label = label), vjust = -0.5, size = 10) +
  theme(legend.position = "none",
        #plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16.5)) +
  coord_cartesian(ylim = c(170, 330)) 
# Combine the plots
combined_plot <- (plot_essentiality | (plot_coding_length/ plot_contig)) +
  plot_annotation(#title = "Marginal Effects of Gene Essentiality, Coding Length, and Chromosomal Location \non Predicted Population-level PSC Counts per Gene",
    tag_levels = 'a') &
  theme(plot.tag=element_text(size = 16, face = "bold"))+
  plot_layout(guides = "collect") 

# Print the combined plot
print(combined_plot)

##########all separated plot#################
ggplot(effect_data_stop_all, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 0.5) +  # Adjust line thickness
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color=NA) +  # Add confidence interval
  facet_wrap(~ facet, scales = "free_y") +  # Facet by chromosome, allowing free scaling of y-axis
  scale_color_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +# Custom color palette
  labs(
    title = "Predicted PSC Counts Per Gene  (Population-level) vs Coding Length\nby Chromosomal Location, Gene Essentiality",
    x = "Log(Coding Length) (bp)",
    y = "Predicted PSC Count Per Gene (Population-level)",
    color = "Gene Essentiality",  # Combine legends by using the same title
    fill = "Gene Essentiality"  # Combine legends by using the same title
  ) +
  theme_bw(base_size = 14) +  # Clean black and white theme with larger base font size
  theme(
    strip.text = element_text(face = "bold", size = 12),  # Bold and larger facet labels
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Center and bold plot title
    legend.position = "bottom",  # Move legend to the bottom
    legend.title = element_text(face = "bold"),  # Bold legend title
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)  # Adjust margins to make room for annotation
  )

# Use grid to add text to the bottom right corner of the whole plot
print(p, vp = viewport(layout = grid.layout(1, 1)))  # Print the plot in a single layout
grid.text(
  "Data Source: Ag1000G (3.0-3.8)", 
  x = unit(1, "npc") - unit(2, "mm"), 
  y = unit(2, "mm"), 
  just = c("right", "bottom"),
  gp = gpar(fontsize = 12)
)

#######################################GLM for potential PSC count vs coding length#########################
reduced_model <- glm(potential_stop_gained_mut ~ Essentiality + log_coding_length, data= combined_df_cds_w_mutation, family = quasipoisson)
reduced_model_more <- glm(potential_stop_gained_mut ~ Essentiality, data= combined_df_cds_w_mutation, family = quasipoisson)
full_model<-glm(potential_stop_gained_mut ~ Essentiality + log_coding_length + contig.y, data= combined_df_cds_w_mutation, family = quasipoisson)

# Compare the models using ANOVA
anova(reduced_model_more, reduced_model, full_model, test = "Chisq")
summary(full_model)

##################################visualise potential PSC effect plot######################################
# Generate partial effect data
# Effect data for Essentiality
effect_data_essentiality_potential <- ggpredict(full_model, terms = "Essentiality")
# Effect data for coding sequence length
effect_data_coding_length_potential <- ggpredict(full_model, terms = "log_coding_length")
# Effect data for contig
effect_data_contig_potential <- ggpredict(full_model, terms = "contig.y")
# Reorder the contig levels
effect_data_contig_potential$x <- factor(effect_data_contig_potential$x, levels = c("2L", "2R", "3R", "3L", "X"))
# Significance labels based on the provided coefficients
significance_labels_potential <- data.frame(
  x = c("2R", "3L", "3R", "X"),
  y = c(263, 263, 263, 263),  # Adjust y position as necessary
  label = c("***", "", "", "***")
)

# Add significance labels for coding length plot
significance_labels_coding_length_potential <- data.frame(
  x = c(mean(range(effect_data_coding_length_potential$x))),
  y = c( max(effect_data_coding_length_potential$conf.high, na.rm = TRUE)*0.8),
  label = c("***")
)

# Generate the effect plots
plot_essentiality_potential <- ggplot(effect_data_essentiality_potential, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.5)) +
  labs(#title = "Effect of Gene Essentiality on Potential PSC Counts",
    x = "Essentiality",
    y = "Potential PSC Counts per Gene") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 16.5),
        axis.text = element_text(size = 16)) +
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +
  geom_signif(data=effect_data_essentiality,stat="signif",position="identity",
              comparisons=list(c("Essential","Non-essential")),map_signif_level = TRUE,y_position = 271, annotations="***", textsize = 10)+
  coord_cartesian(ylim = c(150, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  # Reducing the gap by setting no expansion at the lower end

plot_coding_length_potential <- ggplot(effect_data_coding_length_potential, aes(x = x, y = predicted)) +
  geom_line(color = "black", size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha =0.8, fill = "red", size=1) +
  labs(#title = "Effect of Log(Coding Length) on Potential PSC Counts",
    x = "Log(Coding Length) ",
    y = "Potential PSC Counts per Gene") +
  geom_text(data = significance_labels_coding_length_potential, aes(x = x, y = y, label = label), vjust = -0.5, size = 10) +  # Adjust size as neede
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 16.5),
        axis.text = element_text(size = 16))

plot_contig_potential <- ggplot(effect_data_contig_potential, aes(x = x, y = predicted)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "red", position = position_dodge(0.5)) +
  geom_point(color = "black", size = 1.5, position = position_dodge(0.5)) +
  labs(#title = "Effect of Chromosomal Location on Potential PSC Counts",
    x = "Chromosomal Location",
    y = "Potential PSC Counts per Gene") +
  theme_classic() +
  geom_text(data = significance_labels_potential, aes(x = x, y = y, label = label), vjust = -0.5, size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 16.5),
        axis.text = element_text(size = 16)) +
  coord_cartesian(ylim = c(190, 280))
# Combine the plots
combined_plot_potential <- (plot_essentiality_potential | (plot_coding_length_potential / plot_contig_potential)) +
  plot_annotation(#title = "Marginal Effects of Gene Essentiality, Coding Length, and Chromosomal Location \non Predicted Population-level PSC Counts per Gene",
    tag_levels = 'a') &
  theme(plot.tag=element_text(size = 16, face = "bold"))+
  plot_layout(guides = "collect") 

# Print the combined plot
print(combined_plot_potential)

#############All plots##############
effect_data_all <-ggpredict(full_model, terms = c("log_coding_length","Essentiality", "contig.y"))

ggplot(effect_data_all, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 0.5) +  # Adjust line thickness
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color=NA) +  # Add confidence interval
  facet_wrap(~ facet, scales = "free_y") +  # Facet by chromosome, allowing free scaling of y-axis
  scale_color_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +# Custom color palette
  labs(
    title = "Predicted Potential Premature Stop Codon (PSC) Counts Per Gene vs. Coding Length \nby Chromosomal Location and Gene Essentiality",
    x = "Log(Coding Length) ",
    y = "Predicted Potential PSC Count Per Gene (Genome-level)",
    color = "Gene Essentiality",  # Combine legends by using the same title
    fill = "Gene Essentiality"  # Combine legends by using the same title
  ) +
  theme_bw(base_size = 14) +  # Clean black and white theme with larger base font size
  theme(
    strip.text = element_text(face = "bold", size = 12),  # Bold and larger facet labels
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Center and bold plot title
    legend.position = "bottom",  # Move legend to the bottom
    legend.title = element_text(face = "bold")  # Bold legend title
  )

grid.text(
  "Data Source: AgamP4 Reference Genome", 
  x = unit(1, "npc") - unit(2, "mm"), 
  y = unit(2, "mm"), 
  just = c("right", "bottom"),
  gp = gpar(fontsize = 12)
)

#############################################GLMs on actual PSC count vs. relative mutability#########################################
# Fit Poisson regression model with log-transformed data
full_model_weighted <- glm.nb(gam_sequences_stop_gained_count ~ log_rela_potential_weighted+Essentiality +  chromosome_group, data = combined_df_cds_w_mutation)
reduced_model_weighted <- glm.nb(gam_sequences_stop_gained_count ~ Essentiality+log_rela_potential_weighted, data = combined_df_cds_w_mutation)
reduced_model_weighted_more<-glm.nb(gam_sequences_stop_gained_count ~ log_rela_potential_weighted, data = combined_df_cds_w_mutation)
full_model_weighted_more <- glm.nb( gam_sequences_stop_gained_count~ log_rela_potential_weighted+Essentiality + log_rela_potential_weighted*Essentiality+ contig.y, data = combined_df_cds_w_mutation)

most<-glm.nb( gam_sequences_stop_gained_count~ log_rela_potential_weighted+Essentiality + log_rela_potential_weighted*Essentiality+ chromosome_group +log_rela_potential_weighted*chromosome_group, data = combined_df_cds_w_mutation)

#model selction
anova(reduced_model_weighted_more,reduced_model_weighted,full_model_weighted,full_model_weighted_more, most,test = "Chisq")
summary(most)

#############add effect plot##############
effect_data_essentiality_most <- ggpredict(most, terms = "Essentiality")

# Effect data for interaction between log_rela_potential_weighted and Essentiality
effect_data_interaction_most <- ggpredict(most, terms = c("log_rela_potential_weighted", "Essentiality"))

#interaction between mutability &chromosome
effect_data_interaction_chromosme_most <- ggpredict(most, terms = c("log_rela_potential_weighted", "chromosome_group"))
#Effect data for contig
effect_data_contig_most <- ggpredict(most, terms = "chromosome_group")
# Reorder the contig levels
effect_data_contig_most$x <- factor(effect_data_contig_most$x, levels = c("Autosome", "X"))

# Significance labels based on the provided coefficients
significance_labels_most<- data.frame(
  x =  "X",
  y = 260,  # Adjust y position as necessary
  label =  "***"
)

# Add significance labels for coding length plot
significance_labels_coding_length_most <- data.frame(
  x = mean(range(effect_data_interaction_most$x)),
  y = max(effect_data_interaction_most$conf.high, na.rm = TRUE) * 0.9,
  label = "***"
)

significance_labels_contig_most <- data.frame(
  x = mean(range(effect_data_interaction_chromosme_most$x)),
  y = max(effect_data_interaction_chromosme_most$conf.high, na.rm = TRUE) * 0.9,
  label = "***"
)

###plot them###
plot_essentiality_potential_actual_weighted <- ggplot(effect_data_essentiality_most, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.5)) +
  labs(#title = "Predicted Observed PSC Counts by Gene Essentiality",
    x = "Gene Essentiality",
    y = "PSC Counts per Gene") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.2, size = 10.5),
        axis.title = element_text(size = 16.5),
        axis.text = element_text(size = 16)) +
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +
  geom_signif(comparisons = list(c("Essential", "Non-essential")),
              annotations = "***", 
              y_position = 265, 
              tip_length = 0.05, 
              textsize = 10) +
  coord_cartesian(ylim = c(80, 290)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
#ylim(40,NA)# Adjust the expansion here

plot_interaction_potential_actual_weighted <- ggplot(effect_data_interaction_most, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color=NA) +
  labs(#title = "Interaction Effect of Log(Relative Mutability) and Gene Essentiality",
    x = "Log(Relative Mutability of Gene)",
    y = "PSC Counts per Gene",
    color = "Gene Essentiality",
    fill = "Gene Essentiality") +
  scale_color_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +
  theme_classic() +
  geom_text(data = significance_labels_coding_length_most, aes(x = x, y = y, label = label), vjust = -0.5, size = 10, inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5),
        axis.title = element_text(size = 16.5),
        axis.text = element_text(size = 16),
        legend.position= c(0.25,0.8),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))


plot_contig_potential_actual_weighted <- ggplot(effect_data_contig_most, aes(x = x, y = predicted)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "red", position = position_dodge(0.5)) +
  geom_point(color = "black", size = 1.5, position = position_dodge(0.5)) +
  labs(#title = "Predicted Observed PSC Counts by Chromosomes",
    x = "Chromosomal Location",
    y = "PSC Counts per Gene") +
  theme_classic() +
  geom_signif(comparisons = list(c("Autosome", "X")),
              annotations = "***", 
              y_position = 265, 
              tip_length = 0.05, 
              textsize = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5),
        axis.title = element_text(size = 16.5),
        axis.text = element_text(size = 16),
        legend.position= c(0.2,0.8))+
  coord_cartesian(ylim = c(175, 280)) 

plot_interaction_chromosome <- ggplot(effect_data_interaction_chromosme_most, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  labs(#title = "Interaction Effect of Log(Relative Mutability) and Chromosomal Location",
    x = "Log(Relative Mutability of Gene)",
    y = "PSC Counts per Gene",
    color = "Chromosomal Location",
    fill = "Chromosomal Location") +
  scale_color_manual(values = c("X" = "red", "Autosome" = "blue")) +
  scale_fill_manual(values = c("X" = "red", "Autosome" = "blue")) +
  theme_classic() +
  geom_text(data = significance_labels_contig_most, aes(x = x, y = y, label = label), vjust = -0.5, size = 10, inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5),
        axis.title = element_text(size = 16.5),
        axis.text = element_text(size = 16),
        legend.position = c(0.35, 0.8),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

# Combine the plots using patchwork
combined_plot_potential_weighted <- ((plot_essentiality_potential_actual_weighted/plot_contig_potential_actual_weighted) | (plot_interaction_potential_actual_weighted / plot_interaction_chromosome)) +
  plot_annotation(#title = "Marginal Effects of Gene Essentiality, Relative Mutability, and Chromosomal Location \non Predicted Population-level PSC Counts Per Gene",
    tag_levels = 'a')&
  theme(plot.tag=element_text(size = 16, face = "bold"))+
  plot_layout(guides = "collect")*
  theme(legend.position= c(0.8,0.8), legend.text = element_text(size = 16), legend.title = element_text(size = 16))
#theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

print(combined_plot_potential_weighted)

############plot all full model##########
# Generate predictions for interaction effects
predicted_most <- ggpredict(most, terms = c("log_rela_potential_weighted [all]", "Essentiality", "chromosome_group"))

ggplot(predicted_most, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 0.8) +  # Thicker lines for visibility
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3, color = NA) +  # Confidence intervals
  facet_wrap(~ facet, scales = "free_y") +  # Facet by chromosome group
  scale_color_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +  # Custom color palette for lines
  scale_fill_manual(values = c("Essential" = "#404080", "Non-essential" = "#ffa500")) +  # Custom color palette for ribbons
  labs(
    title = "Predicted Population-Level PSC Counts Per Gene vs. Genome-Level Relative Mutability\nby Chromosome Groups and Gene Essentiality",
    x = "Log(Relative Mutability of Potential PSC Per Gene) (Genome-level)",
    y = "Predicted PSC Count Per Gene (Population-level)",
    color = "Gene Essentiality",
    fill = "Gene Essentiality"
  ) +
  theme_bw(base_size = 14) +  # Clean black and white theme with larger base font size
  theme(
    strip.text = element_text(face = "bold", size = 12),  # Bold and larger facet labels
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Center and bold plot title
    legend.position = "bottom",  # Move legend to the bottom
    legend.title = element_text(face = "bold"),  # Bold legend title
    plot.margin = margin(10, 10, 40, 10)  # Adjust margins to make room for annotation
  )
grid.text(
  "Data Source: Ag1000G (3.0-3.8) & AgamP4 reference genome", 
  x = unit(1, "npc") - unit(2, "mm"), 
  y = unit(2, "mm"), 
  just = c("right", "bottom"),
  gp = gpar(fontsize = 12)
)
