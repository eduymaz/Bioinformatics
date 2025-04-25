data_path <- "your_data_path"
#install.packages("readxl")
library(readxl)

all_files <- list.files(data_path, pattern = "^data_05.*\\.xlsx$", full.names = TRUE)

#make a list
data_list <- list()

# for loop 
for (i in seq_along(all_files)) {
  data_list[[i]] <- read_excel(all_files[i], sheet = 1)
}

names(data_list) <- paste0("data", 1:length(all_files))

####### FOR Bar plot figures | Fig 1 ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

for (i in seq_along(all_files)) {
  
  # Read file
  data <- read_excel(all_files[i], sheet = 1)
  
  # Remove NA 
  cleaned_data <- data %>%
    filter(!is.na(Ratio_PSI))

  #Filtration non-numeric values
  cleaned_data <- cleaned_data %>%
    mutate(Ratio_PSI = as.numeric(Ratio_PSI)) %>%
    filter(!is.na(Ratio_PSI))
  
  #Computituonal AS events
  as_events <- cleaned_data %>%
    group_by(gp1, gp2) %>%
    summarise(avg_ratio_psi = mean(Ratio_PSI))
  
  # Bar Plot
  p <- ggplot(as_events, aes(x = interaction(gp1, gp2), y = avg_ratio_psi, fill = gp2)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = paste0("Average Alternative Splicing Events by Stage: ", basename(all_files[i])), 
         x = "Stage Comparison", y = "Average Ratio PSI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save a plot
  ggsave(filename = paste0("AS_Barplot_", i, ".png"), plot = p, width = 10, height = 6)
  
  # List of results
  data_list[[i]] <- list(file_name = basename(all_files[i]), plot = p)
}

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
# Figure 2 | Box plots 
data$FDR <- p.adjust(data$p_value, method = "BH")


selected_genes <- c("AMD1", "ARMC8", "ATP6V1C2", "TDP1", "THOC2", "UBAP2L") # genes are identified from pre-analysis results
data <- data{i}
filtered_data <- data %>% filter(gene_name %in% selected_genes)

data_tidy <- filtered_data %>% 
  gather(key = "Group", value = "PSI", PSI_gp1, PSI_gp2)

# Grup info is merged
data_tidy <- data_tidy %>% 
  mutate(Group_Label = paste(gp1, "vs", gp2))

# Box plot  - with facet_wrap for all gp1 ve gp2 ccombination:
ggplot(data_tidy, aes(x = gene_name, y = PSI, fill = Group)) +
  geom_boxplot() +
  labs(title = "PSI Changes for Selected Genes", 
       x = "Genes", y = "PSI Values") +
  theme_minimal() +
  facet_wrap(~ Group_Label, scales = "free_y", ncol = 2)  # Multi box plot 


####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
# Figure 4 | Heatmap
merged_data <- rbind(data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11)
# FDR hesaplama
merged_data$FDR <- p.adjust(merged_data$p_value, method = "BH")
selected_genes <- c("AMD1", "ARMC8", "ATP6V1C2", "TDP1", "THOC2", "UBAP2L")
scats_results <- read_excel("/Users/edyilmaz/Desktop/embryo/results/smaller_than_0.05_pvalues_new2.xlsx", sheet=1)

####################### Fig 4| Heatmap with diagonal 

#FILTERED_DATA
filtered_data <- merged_data %>%
  filter(gene_name %in% selected_genes)
#write.csv(filtered_data, "filtered_data.csv", row.names = FALSE)
#filtered_data <- read.csv("filtered_data.csv")

# PLOT CODES: 
library(corrplot)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

AMD1_data <- filtered_data %>% filter(gene_name == "AMD1")
ARMC8_data <- filtered_data %>% filter(gene_name == "ARMC8")
ATP6V1C2_data <- filtered_data %>% filter(gene_name == "ATP6V1C2")
TDP1_data <- filtered_data %>% filter(gene_name == "TDP1")
THOC2_data <- filtered_data %>% filter(gene_name == "THOC2")
UBAP2L_data <- filtered_data %>% filter(gene_name == "UBAP2L")

##########

#data: AMD1_data | ARMC8_data | ATP6V1C2_data | TDP1_data | THOC2_data | UBAP2L_data
data <- UBAP2L_data
  
# calculate average values
average_values <- data %>%
  group_by(gp1, gp2) %>%
  summarize(
    avg_PSI_gp1 = mean(PSI_gp1, na.rm = TRUE),
    avg_PSI_gp2 = mean(PSI_gp2, na.rm = TRUE)
  ) %>%
  ungroup()

# 
unique_gp <- unique(c(average_values$gp1, average_values$gp2))
unique_gp <- sort(unique_gp) 

#make a matrix
heatmap_matrix <- matrix(NA, nrow = length(unique_gp), ncol = length(unique_gp))
rownames(heatmap_matrix) <- unique_gp
colnames(heatmap_matrix) <- unique_gp

# put the average values in map
for (i in unique_gp) {
  for (j in unique_gp) {
    avg1 <- average_values %>%
      filter(gp1 == i & gp2 == j) %>%
      summarize(avg = mean(avg_PSI_gp2, na.rm = TRUE)) %>%
      pull(avg)
    avg2 <- average_values %>%
      filter(gp1 == j & gp2 == i) %>%
      summarize(avg = mean(avg_PSI_gp1, na.rm = TRUE)) %>%
      pull(avg)
    if (is.na(avg1) && is.na(avg2)) {
      heatmap_matrix[i, j] <- NA
    } else {
      heatmap_matrix[i, j] <- mean(c(avg1, avg2), na.rm = TRUE)
    }
  }
}

# Make a correlation matrx
numeric_data <- heatmap_matrix
correlation_matrix <- cor(numeric_data, use = "pairwise.complete.obs", method = "pearson")

corrplot(correlation_matrix, method = "color", type = "full", tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", number.cex = 0.1)
