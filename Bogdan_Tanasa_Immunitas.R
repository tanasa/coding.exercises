# the files were downloaded from https://github.com/immunitastx/coding_exercises
# I would have to assume that we work with normalized values (as provided by Kalisto, Salmon, Sailfish)
# If the gene expression values reflect the raw counts, we will have to use limma, edgeR, or DESeq2.

library("ggplot2")
library("reshape2")
library("dplyr")
library("tidyr")
library("nortest")
library("pheatmap")
library("car")
library("nortest")



# Approach Outline:

# 1. Data Reading: import the dataset for analysis.

# 2. Exploratory Data Analysis (EDA):

# 2.1. Dimensionality Reduction: PCA and MDS to identify clusters in the data.

# 2.2. Hierarchical Clustering to detect clusters, complemented by PCA and MDS.

# 2.3. Make scatter plots to compare average expression values between Condition 1 and Condition 2.

# 3. Normality Check: Assess whether gene expression data for each gene follows a normal distribution in both Condition 1 and Condition 2.

# 4. Variance Equality Check: Determine if gene expression data exhibits equal variance between Condition 1 and Condition 2.

# Observations: A large number of genes show Gaussian distribution and equal variance in their expression levels.

# 5. Statistical Testing: Apply a t-test assuming equal variance, or more robustly, use Welch’s t-test to account for any variance differences.

# 6. Correction for Multiple Testing: use FDR correction to adjust p-values for multiple comparisons.        



# 1. Data Reading: import the dataset for analysis.



condition = read.table("condition.tsv", header = T, sep ="\t", stringsAsFactors = FALSE)
head(condition, 2)

omics = read.table("omics.tsv", header = T, sep ="\t", stringsAsFactors = FALSE)
head(omics, 2)

# Counting the number of cells in Condition1 and in Condition2
print("the number of cells in condition1 and in condition2")
count_condition_cells <- table(condition$condition)
print(count_condition_cells)

x = merge(condition, omics, by = "X")
head(x,2)

y = x %>% select(-X)
head(y,2)
dim(y)

# z = t(y)
# head(z,4)
# dim(z)
# we can work with either z or y data frame

# checking whether we have any NA : 

total_na_count <- sum(is.na(y))

cat("Total number of NA values in the data frame:\n")
print(total_na_count)

# using dplyr / tidyr / tidyverse functions to reformat the data

y_long <- y %>%
  pivot_longer(cols = -condition, names_to = "gene", values_to = "expression")

head(y_long, 2)
tail(y_long, 2)
dim(y_long)



# 2. Exploratory Data Analysis (EDA):



# 2.1. Dimensionality Reduction: PCA and MDS to identify clusters in the data.

# PCA

y_pca <- y %>%
  select(-condition) 

# Perform PCA and display the results
pca_result <- prcomp(y_pca, scale. = TRUE)
# summary(pca_result)

pca_scores <- as.data.frame(pca_result$x)
pca_scores$condition <- y$condition


ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition)) +
  geom_point() +
  labs(title = "PCA of Gene Expression Data",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

print("PCA shows good separability of the expression data")

# MDS

y_mds <- y %>%
  select(-condition) 

# Compute the distance matrix and perform MDS

distance_matrix <- dist(y_mds)
mds_result <- cmdscale(distance_matrix, k = 2)
mds_scores <- as.data.frame(mds_result)
mds_scores$condition <- y$condition

ggplot(mds_scores, aes(x = V1, y = V2, color = condition)) +
  geom_point() +
  labs(title = "MDS of Gene Expression Data",
       x = "dim 1",
       y = "dim 2") +
  theme_minimal()

print("MDS shows good separability of the expression data")



# 2.3. Make scatter plots to compare average expression values between Condition 1 and Condition 2.

# Calculate the average expression for each gene in each condition

average_expression <- y_long %>%
  group_by(condition, gene) %>%
  summarise(average_expression = mean(expression, na.rm = TRUE)) %>%
  pivot_wider(names_from = condition, values_from = average_expression, names_prefix = "average_") 

head(average_expression, 2)
tail(average_expression, 2)

# Make the scatter plot
ggplot(average_expression, aes(x = average_condition_1, y = average_condition_2)) +
  geom_point() +
  labs(title = "Scatter Plot of Average Gene Expression Levels",
       x = "Average Expression in Condition 1",
       y = "Average Expression in Condition 2") +
  theme_minimal()



# 2.2. Hierarchical Clustering to detect clusters

z_long <- y %>%
  select(-condition) %>%
  t() %>%
  as.data.frame()

# Set column names to the first row and remove the first row
head(z_long, 2)
colnames(z_long) <- z_long[1,]

z_long <- as.data.frame(lapply(z_long, as.numeric))

# Scale the data
z_scaled <- scale(z_long)
dist_matrix <- dist(z_scaled)

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "complete")

# Plot the dendrogram
# plot(hc, main = "Hierarchical Clustering of Genes")

print("Hierarchical clustering of the genes :")

options(repr.plot.width=16, repr.plot.height=8)
pheatmap(z_scaled, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         treeheight_row = 0,          # Set tree height for rows to 0 (hide dendrogram)
         treeheight_col = 0, 
         main = "Heatmap of Gene Clustering")      



# Exploratory data analysis

# The density plots of the genes in Condition_1 and in _Condition_2

# y_long <- y %>%
#  pivot_longer(cols = -condition, names_to = "gene", values_to = "expression")
# head(y_long, 2)
# dim(y_long)

# The Density Plots of the expression of ALl the genes in Condition 1 and in Condition 2
# ggplot(y_long, aes(x = expression, fill = condition)) +
#  geom_density(alpha = 0.5) +
#  labs(title = "Density Plots of Gene Expressions",
#       x = "Expression Values",
#       y = "Density") +
#  theme_minimal()

# expression_condition_1 <- y_long %>%
#  filter(condition == 'condition_1') %>%
#  pull(expression)

# expression_condition_2 <- y_long %>%
#  filter(condition == 'condition_2') %>%
#  pull(expression)

# using Anderson-Darling Test to check the Gaussian Distribution of the expression of all the genes in those two conditions
# ad.test(expression_condition_1)
# ad.test(expression_condition_2)

# print("the results of AD test indicate that there is strong evidence that the data deviates significantly from a normal distribution")



# 3. Normality Check: Assess whether gene expression data for each gene follows a normal distribution in both Condition 1 and Condition 2.



# A function to apply Anderson-Darling test and check normality of the gene expression distribution

check_normality <- function(expression_values) {
  result <- ad.test(expression_values)
  return(result$p.value)
}

# Apply Anderson-Darling test to each gene and each condition

normality_results <- y_long %>%
  group_by(condition, gene) %>%
  summarise(p_value = check_normality(expression)) %>%
  ungroup()

dim(normality_results)

# Define a significance level
alpha <- 0.05

# Count how many genes follow a normal distribution in each condition

normal_genes_count <- normality_results %>%
  group_by(condition) %>%
  summarise(normal_genes = sum(p_value > alpha)) %>%
  ungroup()

# Print the results
# print(normality_results)
print("the number of the genes that follow a normal distribution :")
print(as.data.frame(normal_genes_count))

# Count how many genes do not follow a normal distribution in each condition

non_normal_genes_count <- normality_results %>%
  group_by(condition) %>%
  summarise(non_normal_genes = sum(p_value <= alpha)) %>%
  ungroup()

# Print the results
# print(normality_results)
print("the number of the genes that do not follow a normal distribution :")
print(as.data.frame(non_normal_genes_count))

print("as we can see a large majority of the genes follow a normal distribution")



# 4. Variance Equality Check: Determine if gene expression data exhibits equal variance between Condition 1 and Condition 2.



# I would like to double check the equality of VARIANCES : we use VAR.TEST (F.TEST)
# Testing the equality of Variances for each gene in order to assess Homoscedasticity or Heteroscedasticity 

test_equal_variance <- function(expression_condition_1, expression_condition_2) {
  var_test <- var.test(expression_condition_1, expression_condition_2)
  return(var_test$p.value)
}

# Pivot the data to get separate columns for each condition
# Gene are on the row, and the expression data in Condition 1 and in Condition 2 is contained in LISTS

variance_results <- y_long %>%
  pivot_wider(names_from = condition, values_from = expression) %>%
  group_by(gene) %>%
  rowwise() %>%
  mutate(p_value = test_equal_variance(condition_1, condition_2)) %>%
  ungroup()

# Define a significance level
alpha <- 0.05

# Count how many genes have equal variance
equal_variance_genes_count <- variance_results %>%
  summarise(equal_variance_genes = sum(p_value > alpha)) %>%
  pull(equal_variance_genes)  # pull value to get a single number

# Count how many genes do not have equal variance
total_genes_count <- nrow(variance_results)
different_variance_genes_count <- total_genes_count - equal_variance_genes_count

# Print the results
print("Number of gene whose expressions have equal or different variances in Condition 1 and in Condition 2 :")
print(paste("Number of genes with equal variance:", equal_variance_genes_count))
print(paste("Number of genes with different variance:", different_variance_genes_count))

# Observations: A large number of genes show Gaussian distribution and equal variance in their expression levels.

# Given the gaussian distribution of the data and its homoscedasticity, 
# we can apply T test to identify the genes with the differential expression.
# we can use the standard t.test with equal variances ; 
# t_test <- t.test(expression_condition_1, expression_condition_2, var.equal = TRUE)
# or Welch T-test, given the fact that several genes do not have equal variance : 
# t_test <- t.test(expression_condition_1, expression_condition_2



# 5. Statistical Testing: Apply a t-test assuming equal variance, or more robustly, use Welch’s t-test to account for any variance differences.



test_equal_expression  <- function(expression_condition_1, expression_condition_2) {
  t_test <- t.test(expression_condition_1, expression_condition_2)
  return(t_test$p.value)
}

# Pivot the data to get separate columns for each condition
# Genes are on the row, and the expression data in Condition 1 and in Condition 2 is contained in LISTS

expression  <- y_long %>%
  pivot_wider(names_from = condition, values_from = expression) %>%
  group_by(gene) %>%
  rowwise() %>%
  mutate(p_value = test_equal_expression (condition_1, condition_2)) %>%
  ungroup()

# Define a significance level
alpha <- 0.05

# Count how many genes have non-differential expression
equal_expression_genes_count <- expression  %>%
  summarise(equal_expression_genes = sum(p_value > alpha)) %>%
  pull(equal_expression_genes)  # pull value to get a single number

# Count how many genes are differentially expressed 
total_genes_count <- nrow(expression )
different_expression_genes_count <- total_genes_count - equal_expression_genes_count

# Print the results
print("Number of genes whose expressions have equal or different expression s in Condition 1 and in Condition 2:")
print(paste("Number of genes with equal expression :", equal_expression_genes_count))
print(paste("Number of genes with different expression :", different_expression_genes_count))

# Printing the genes that are differentially expressed 

genes_differential_expression <- expression %>%
                                filter(p_value < alpha)  

print("The number of the genes with differential expression before FDR-correction :")
dim(genes_differential_expression)

print("Genes with differential expression:")
print(genes_differential_expression)



# 6. Correction for Multiple Testing: use FDR correction to adjust p-values for multiple comparisons.  



# Applying FDR correction : we use p.adjust

expression <- expression %>%
              mutate(adj_p_value = p.adjust(p_value, method = "fdr"))

# Count how many genes have non-differential expression
equal_expression_genes_count <- expression %>%
  summarise(equal_expression_genes = sum(adj_p_value > alpha)) %>%
  pull(equal_expression_genes) 

# Count how many genes are differentially expressed
total_genes_count <- nrow(expression)
different_expression_genes_count <- total_genes_count - equal_expression_genes_count

# Print the results
print("Number of genes whose expressions have equal or different expressions in Condition 1 and Condition 2:")
print(paste("Number of genes with equal expression:", equal_expression_genes_count))
print(paste("Number of genes with different expression:", different_expression_genes_count))

# Filter and print the genes with different expressions (FDR-corrected)
genes_different_expression <- expression %>%
                filter(adj_p_value < alpha)  # Genes with significant differences after FDR correction

# print("Genes with differential expression:")
# print(as.data.frame(genes_differential_expression))


