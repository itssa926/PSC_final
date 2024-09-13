#####load flybase gene data########
sterile_dominant <- read.delim("../data/flybase/FlyBase_Fields_sterile_dominant.txt",sep="\t")
sterile_recessive <- read.delim("../data/flybase/FlyBase_Fields_sterile_recessive.txt",sep="\t")
lethal_dominant <- read.delim("../data/flybase/FlyBase_Fields_lethal_dominant.txt",sep="\t")
lethal_recessive <- read.delim("../data/flybase/FlyBase_lethal_recessive.txt",sep="\t")

#ortholog data
all_ortho <- read.csv("../data/ortholog_pair_wrelation.csv") ### all unique orthologs with annotation

###########combine them into a big df###############
sterile_dominant$phenotype <- "sterile_dominant"
sterile_recessive$phenotype <- "sterile_recessive"
lethal_dominant$phenotype <- "lethal_dominant"
lethal_recessive$phenotype <- "lethal_recessive"

gene_list_new <-rbind(sterile_dominant,sterile_recessive,lethal_dominant,lethal_recessive) #### contains all essential genes in drosopjila
write.csv(gene_list_new,"../data/flybase/gene_list_new.csv")

##############annotate essentiality############
colnames(gene_list_new)[which(names(gene_list_new) == "FBID_KEY")] <- "Flybase_ID"
gene_list_new$Essentiality <- "Essential" #7690 genes

#get a df with all essential and non-essential genes
essential_not_df <- merge(x=gene_list_new, y=all_ortho, by = "Flybase_ID",
                          all.y = TRUE)

#remove duplicate genes
unique_essential_not_df <- essential_not_df[!duplicated(essential_not_df$id), ]#9,440 unique An.gambiae genes.
#classify genes without annotation as non-essential
unique_essential_not_df$Essentiality[is.na(unique_essential_not_df$Essentiality)] <- "Non-essential"
write.csv(unique_essential_not_df, "../data/flybase/new_essential_not.csv")

#count essential genes and non-essential genes
essential_genes <- subset(unique_essential_not_df, unique_essential_not_df$Essentiality=="Essential" ) #1973 essential genes, other don't have orthologus in mosquito
non_essential_genes <- subset(unique_essential_not_df, unique_essential_not_df$Essentiality=="Non-essential" ) #7467

#######################visualise in bar chart#####################
#count each relationship
total_counts <- table(unique_essential_not_df$relationship)
essential_counts <-table(essential_genes$relationship)
non_essential_genes <-subset(unique_essential_not_df, unique_essential_not_df$Essentiality == "Non-essential")
non_essential_counts <- table(non_essential_genes$relationship)

# Combine counts into a dataframe
counts_df <- data.frame(
  Relationship = names(total_counts),
  Total = as.numeric(total_counts),
  Essential = as.numeric(essential_counts),
  Non_Essential = as.numeric(non_essential_counts)
)

relationship_long <- pivot_longer(counts_df, cols = c(Essential, Non_Essential, Total), 
                                  names_to = "Category", values_to = "Count") #convert a long dataframe for plotting


#visualise
my_bar <- ggplot(relationship_long, aes(x = Relationship, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  #labs(title = "Anopheles Gene Counts by Orthologous Relationship to Drosophila and Essentiality",
       #x = "Orthologous Relationship to Drosophila", y = "Count") +
  geom_text(aes(label = paste("n=", Count)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.4, 
            size = 3.9) +
  scale_x_discrete(labels = c('1:1', '1:Many', 'Many:1', 'Many:Many')) +
  scale_fill_manual(values = c("Essential" = "#404080", "Non_Essential" = "#ffa500", "Total" = "forestgreen"), 
                    labels = c("Essential Genes", "Non-essential Genes", "Total")) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme(legend.position = "top",
        legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(face = "bold", hjust= 0.5,size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.caption = element_text(size = 14),
        legend.text = element_text(size = 16)) 
my_bar


