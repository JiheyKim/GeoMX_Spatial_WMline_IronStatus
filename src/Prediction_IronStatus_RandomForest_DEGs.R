library(dplyr)
library(tibble)
library(randomForest)

# Step 1: Read the expression matrix (ensure Gene column is treated properly)
#expr_data <- read.csv("train_expression_matrix.csv", check.names = FALSE)
expr_data <- read.csv("train_expression_matrix_CORRECTED.csv", check.names = FALSE)
#expr_data <- read.csv("train_expression_matrix_allGenes.csv", check.names = FALSE)

# Step 2: Read metadata
labels <- read.csv("train_labels_CORRECTED.csv")

# Step 3: Ensure metadata sample names match expression data row names
labels <- labels %>% rename(SampleID = ROI)

# Step 4: Remove rows in 'labels' where 'Iron_Status' is empty or NA
labels <- labels %>% filter(!is.na(Iron_Status) & Iron_Status != "")

# Step 5: Transpose the expression data to have samples as rows and genes as columns
expr_data_t <- as.data.frame(t(expr_data[-1]))  # Transpose, removing the first column (Gene names)
colnames(expr_data_t) <- expr_data$Gene  # Set genes as column names
rownames(expr_data_t) <- colnames(expr_data)[-1]  # Set sample IDs as row names

# Step 6: Merge on SampleID (now with cleaned labels)
train_data <- labels %>% inner_join(expr_data_t %>% rownames_to_column("SampleID"), by = "SampleID")

# Step 7: Check the structure of the merged data
#str(train_data)

# Step 8: Prepare the training data
X_train <- train_data %>% select(-SampleID, -Iron_Status)  # Features (remove SampleID and Iron_Status)
y_train <- train_data$Iron_Status  # Target variable


# Ensure the response variable y_train is a factor (since you're likely performing classification)
y_train <- factor(y_train)

# Step 9: Train the Random Forest Model
set.seed(123)  # For reproducibility
rf_model <- randomForest(x = X_train, y = y_train, ntree = 1000, mtry = sqrt(ncol(X_train)), importance = TRUE)

# Step 10: Read new data (new_ROI_expression.csv)

file_name <- "MS177-10-6_ROI_expression.csv"  # Load new dataset
#file_name <- "MS160-4-2_ROI_expression.csv"  # Load new dataset
#file_name <- "MS160-4-4_ROI_expression.csv"  # Load new dataset
#file_name <- "MS177-6-5_1_24_23_ROI_expression.csv"   # Load new dataset
#file_name <- "MS177-6-5_1_24_23_022to028_ROI_expression.csv"  # Load new dataset
#file_name <- "MS184-6-4_ROI_expression.csv"  # Load new dataset
#file_name <- "MS184-6-3_ROI_expression.csv"  # Load new dataset
#file_name <- "MS160-8-8_ROI_expression.csv"  # Load new dataset
#file_name <- "MS160-8-8_016to018_ROI_expression.csv"  # Load new dataset
#file_name <- "MS160-8-7_ROI_expression.csv"  # Load new dataset
#file_name <- "MS_160-8-8_ROI_expression.csv"  # Load new dataset

base_name <- tools::file_path_sans_ext(basename(file_name))

new_data <- read.csv(file_name, row.names = 1)

# Step 11: Transpose new data so genes are columns
new_data_t <- as.data.frame(t(new_data))


# Step 12: Find common genes between the training and new data
common_genes <- intersect(colnames(X_train), colnames(new_data_t))

# Step 13: Subset new data to match the training genes
new_data_selected <- new_data_t[, common_genes]

# Step 14: Predict on new data
predictions <- predict(rf_model, newdata = new_data_selected)

# Step 15: Output predictions
predictions
# Step 14: Predict probabilities on new data
prob_predictions <- predict(rf_model, newdata = new_data_selected, type = "prob")

# Step 15: View probabilities
print(prob_predictions)
new_data_selected$Predicted_Iron_Status <- predictions

# Save probability predictions to a CSV file
write.csv(prob_predictions, paste0(base_name, "_probabilities_RF.csv"), row.names = TRUE)
write.csv(new_data_selected, paste0(base_name, "_Predicted_Iron_Status_RF.csv"), row.names = TRUE)

#### OOB error rate
oob_error_rate <- rf_model$err.rate[500, "OOB"]
print(paste("OOB Error Rate:", oob_error_rate))
accuracy <- 1 - oob_error_rate
print(paste("Accuracy:", accuracy))

