library(tidyverse)

load("C:/Users/tanma/Desktop/ALL/New project(Rupam da)/cancer_data (File)(1).rda")

## View(data_clinical_sample)
## View(data_expression)

## data_expression_polya <- data_expression[data_expression$platform == "polya", ]
## View(data_expression_polya)

data_expression_polya <- data_expression %>%
  filter(platform == "polya")
data_clinical_sample_polya <- data_clinical_sample %>%
  filter(sample_id %in% data_expression_polya$sample_id)
##At first, filtered the data to keep only samples related to the "polya" platform for further consistent analysis.##At first, filtered the data to keep only samples related to the "polya" platform for further consistent analysis.


## View(data_clinical_sample_polya)

sort(unique(data_clinical_sample_polya$pathology_classification))
data_clinical_sample_polya <- data_clinical_sample_polya %>%
  mutate(
    pathology_classification_simpler = ifelse(
      pathology_classification %in%
        c("Inadequate for diagnosis", "Not available"),
      "Unknown",
      pathology_classification
    )
  )
#####A simpler pathology classification variable, grouping  missing values into a new category called "Unknown".
table(data_clinical_sample_polya$pathology_classification)
table(data_clinical_sample_polya$pathology_classification_simpler)

summary(
  lm(
    data_expression_polya$a1bg ~ data_clinical_sample_polya$pathology_classification_simpler
  )
)
####A quick test linear regression showed how a1bg gene expression varies across cancer types.
#Here p-value is 0.03948  less than 0.05 

################################################################################
######################### New Codes ############################################
################################################################################

genes <- sort(unique(colnames(data_expression_polya)[-(1:2)]))

## dependent_vars <- c("a1bg", "a1cf", "a2m", "a2ml1")
dependent_vars <- genes

x <- data_clinical_sample_polya$pathology_classification_simpler
summary_tables <- list()

for (var in dependent_vars) {
  y <- data_expression_polya[, var, drop = T]
  model <- lm(y ~ x)
  ##Regression analysis for 19,272 genes against simplified cancer types, 
  #extracting coefficients (estimate, SE, p-value) and confidence intervals for each gene-cancer pair.
  
  summary_tables[[length(summary_tables) + 1]] <- data.frame(
    gene = var,
    data.frame(summary(model)$coefficient)[-1, ] %>%
      rownames_to_column(var = "cancer_type") %>%
      mutate(cancer_type = gsub("x", "", cancer_type)) %>%
      set_names(c(
        "cancer_type", "estimate", "se", "t", "p"
      )) %>%
      mutate(
        ci_low = estimate - qnorm(0.975) * se,
        ci_up = estimate + qnorm(0.975) * se,
        sig = p < 0.05
      )
  )

 # print(length(summary_tables))
}

summary_tables <- bind_rows(summary_tables)

################################################################################
######################### Tasks ################################################
################################################################################

# Run this loop for all genes.
# Which cancer types are significant for which groups?
# Make nice plots out of this data.

################################################################################
################################################################################
################################################################################

table(summary_tables$cancer_type[summary_tables$sig == TRUE])
#### here,identified which cancer types had the most significant associations.such as,
#Small cell carcinoma had the highest number of significant gene associations.
#Other notable types are Adenocarcinoma, Squamous cell carcinoma.

summary_tables %>%
  group_by(gene) %>%
  summarize(sig = sum(sig)) %>%
  arrange(sig)
summary_tables %>%
  group_by(gene) %>%
  summarize(sig = sum(sig))
summary_tables %>%
  group_by(gene) %>%
  summarize(sig = sum(sig)) %>%
  arrange(-sig)

# Filter for significant cancer types
summary_tables$q <- p.adjust(summary_tables$p, method = "fdr")
summary_tables$sig_q <- summary_tables$q < 0.05

significant_results <- summary_tables %>%
  filter(sig_q == TRUE) %>%
  select(gene, cancer_type, estimate, ci_low, ci_up, p, q) %>%
  arrange(gene, cancer_type)

# View the significant results
summary(significant_results)
head(significant_results)

length(unique(significant_results$gene))

################################################################################
################################################################################
################################################################################

summary_tables %>%
  group_by(gene) %>%
  summarize(sig_count = sum(sig)) %>%
  arrange(-sig_count) %>%
  ggplot(aes(x = reorder(gene, -sig_count), y = sig_count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme(axis.text.x = element_blank()) +
  labs(x = "Gene", y = "Number of Significant Associations", title = "Significant Associations per Gene")
##Most genes have few significant associations, but some (e.g., A1BG) are associated with multiple cancer types.
################################################################################
################################################################################
################################################################################

summary_tables <- summary_tables %>%
  mutate(color = ifelse(sig_q, ifelse(estimate >= 0, "darkred", "darkblue"), "grey"))

ggplot(
  data = summary_tables %>%
    filter(cancer_type == "Small cell" & abs(estimate) < 4000),
  mapping = aes(x = estimate, y = -log10(q), color = color)
) +
  geom_point() +
  theme_classic() +
  theme(text = element_text(face = "bold")) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = 2,
    size = 1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "black",
    linetype = 2,
    size = 1
  ) +
  scale_color_identity()

################################################################################
################################################################################
################################################################################

cancer_type_significance <- significant_results %>%
  group_by(cancer_type) %>%
  summarize(significant_count = n()) %>%
  arrange(desc(significant_count))

ggplot(
  cancer_type_significance,
  aes(x = reorder(cancer_type, -significant_count), y = significant_count)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Significant Results by Cancer Type", x = "Cancer Type", y = "Number of Significant Associations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
###"Small cell" has the most significant gene association between other cancer type.
################################################################################
################################################################################
################################################################################

library(reshape2)

heatmap_data <- significant_results %>%
  select(gene, cancer_type, estimate) %>%
  spread(key = cancer_type, value = estimate)

heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$gene
heatmap_matrix[is.na(heatmap_matrix)] <- 0

heatmap(
  heatmap_matrix,
  main = "Heatmap of Effect Sizes",
  col = colorRampPalette(c("blue", "white", "red"))(50),
  Rowv = T
)
###Clear patterns of up/down-regulation in specific cancer types (e.g., "Small cell" vs. "Adenocarcinoma").
#########################################
#########################################


