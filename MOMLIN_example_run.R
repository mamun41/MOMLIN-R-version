"Author: Md. Mamunur Rashid
< mamun.stat92@gmail.com > 
MOMLIN software R version; June 14, 2024

"


# Clear environment and set working directory
rm(list = ls())
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)


# load MONLIN functions
source("momlin_func.R")


# Load required libraries
library(caret)
library(pheatmap)
library(reshape2)
library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(matrixStats)
library(gridExtra)



# Set Output directory

ResultsFile <- "momlin_out"
if (!dir.exists(ResultsFile)) {
  dir.create(ResultsFile)
}


# Load BRCA multi-Omics data
mRNA <- fread(paste0(current_dir,'/data/01_rna_tpm.csv'))
X1 <- as.matrix(mRNA[, -1, with=FALSE])
mRNA_name <- colnames(mRNA)[2:ncol(mRNA)]


mutation <- fread(paste0(current_dir,'/data/02_dna_count.csv'))
Y1 <- as.matrix(mutation[, -1, with=FALSE])
mutation_name <- colnames(mutation)[2:ncol(mutation)]


clinical <- fread(paste0(current_dir,'/data/03_clinical_pheno.csv'))
Y2 <- as.matrix(clinical[, -1, with=FALSE])
clinical_name <- colnames(clinical)[2:ncol(clinical)]


SubtypeL <- fread(paste0(current_dir,'/data/00_class_metadata.csv'))
hormon_status <- as.vector(SubtypeL$ERHER2.status)
ER_HER2_status <- as.numeric(hormon_status == 'ER- HER2-')
ER_HER2_status[ER_HER2_status == 0] <- -1


Y2 <- cbind(Y2, ER_HER2_status)
clinical_name <- c(clinical_name, "ER-HER2-")


TiME <- fread(paste0(current_dir,'/data/04_TiME_score.csv'))
Y3 <- as.matrix(TiME[, -1, with=FALSE])
TiME_name <- colnames(TiME)[2:ncol(TiME)]


path <- fread(paste0(current_dir,'/data/05_p.acivity.gsva.csv'))
Y4 <- as.matrix(path[, -1, with=FALSE])
path_name <- colnames(path)[2:ncol(path)]


Slabel <- as.vector(SubtypeL$RCB.category)
n_class <- length(unique(Slabel))


clasID <- as.numeric(factor(Slabel, levels = c('pCR', 'RCB-I', 'RCB-II', 'RCB-III')))
Class_b <- one_hot_encode(labels = clasID,  num_classes = n_class)
Z <- as.matrix(Class_b)


# Variable Feature selection in high dim data
X1raw <- X1
Y1raw <- Y1
Y2raw <- Y2
Y3raw <- Y3
Y4raw <- Y4


# rank & find variable features based on MAD values
X1var <- sort(colMads(X1raw), decreasing = TRUE, index.return = TRUE)$ix
Y1var <- sort(colMads(Y1raw), decreasing = TRUE, index.return = TRUE)$ix
Y2var <- sort(colMads(Y2raw), decreasing = TRUE, index.return = TRUE)$ix
Y3var <- sort(colMads(Y3raw), decreasing = TRUE, index.return = TRUE)$ix
Y4var <- sort(colMads(Y4raw), decreasing = TRUE, index.return = TRUE)$ix

# selected features from each modality (user choice)
X1indices2keep <- floor(0.1 * ncol(X1raw))  # 10% RNA since large feature set
Y1indices2keep <- floor(1.0 * ncol(Y1raw))  # 100% DNA since small feature set
Y2indices2keep <- floor(1.0 * ncol(Y2raw))  # 100% Clinical since small feature set
Y3indices2keep <- floor(1.0 * ncol(Y3raw))  # 100% TiME since small feature set
Y4indices2keep <- floor(1.0 * ncol(Y4raw))  # 100% Pathways since small feature set


# selected datasets 
XX1 <- X1raw[, X1var[1:X1indices2keep]]
YY1 <- Y1raw[, Y1var[1:Y1indices2keep]]
YY2 <- Y2raw[, Y2var[1:Y2indices2keep]]
YY3 <- Y3raw[, Y3var[1:Y3indices2keep]]
YY4 <- Y4raw[, Y4var[1:Y4indices2keep]]


# 
mRNA_name <- mRNA_name[X1var[1:X1indices2keep]]
mutation_name <- mutation_name[Y1var[1:Y1indices2keep]]
clinical_name <- clinical_name[Y2var[1:Y2indices2keep]]
TiME_name <- TiME_name[Y3var[1:Y3indices2keep]]
path_name <- path_name[Y4var[1:Y4indices2keep]]


# aggregate Y's
YY <- list(YY1, YY2, YY3, YY4)
YY_names <- list(mutation_name, clinical_name, TiME_name, path_name)


# Set tuned parameters
opts <- list(
  alpha_u = 0.80, 
  beta = 0.5,
  
  alpha_v1 = 0.60,
  alpha_v2 = 0.40,
  alpha_v3 = 0.60,
  alpha_v4 = 0.70,
  
  scale_type = "std"
)


# Number of features in X, and Y
p <- ncol(XX1)
q <- sapply(YY, ncol)


# Main Algorithm MOMLIN
nModels <- 5
k_fold <- 3


# For reproducibility
rCVF_train <- array(FALSE, dim = c(length(clasID), k_fold, nModels))
rCVF_test <- array(FALSE, dim = c(length(clasID), k_fold, nModels))


# Initialize arrays to store the means over K-fold 
U_mean_i <- array(0, dim = c(p, n_class, nModels))
V_mean_i <- lapply(q, function(r) array(0, dim = c(r, n_class, nModels)))


# To save means CCCs over the K-fold 
CCCs_train_mean_i <- array(0, dim = c(length(YY), n_class, nModels))
CCCs_test_mean_i <- array(0, dim = c(length(YY), n_class, nModels))


# Setup parallel backend
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

results <- foreach(ii = 1:nModels, .packages = c('caret', 'data.table')) %dopar% {
  # Create stratified folds
  # set.seed(ii)  # Set seed for reproducibility within each model
  indices <- createFolds(clasID, k = k_fold, list = TRUE, returnTrain = TRUE)
  
  # Initialize U and V for this model
  U_model <- array(0, dim = c(p, n_class, k_fold))
  V_model <- lapply(q, function(r) array(0, dim = c(r, n_class, k_fold)))
  
  CCCs_train_model <- array(0, dim = c(length(YY), n_class, k_fold))
  CCCs_test_model <- array(0, dim = c(length(YY), n_class, k_fold))
  
  fold_indices <- list()
  
  for (k in seq_along(indices)) {
    cat(sprintf("[conduct fold %d ", k))
    
    opts$k <- k
    opts$ii <- ii
    
    idx_train <- indices[[k]]
    idx_test <- setdiff(seq_along(clasID), idx_train)
    
    fold_indices[[k]] <- list(train = idx_train, test = idx_test)
    
    # rCVF_train[idx_train, k, ii] <- TRUE
    # rCVF_test[idx_test, k, ii] <- TRUE
    
    trainData <- list(
      X = list(
        XX1[idx_train, ],
        XX1[idx_train, ],
        XX1[idx_train, ],
        XX1[idx_train, ]
      ),
      Y = list(
        YY1[idx_train, ],
        YY2[idx_train, ],
        YY3[idx_train, ],
        YY4[idx_train, ]
      ),
      Z = Z[idx_train, ],
      n_class = n_class
    )
    
    testData <- list(
      X = list(
        XX1[idx_test, ],
        XX1[idx_test, ],
        XX1[idx_test, ],
        XX1[idx_test, ]
      ),
      Y = list(
        YY1[idx_test, ],
        YY2[idx_test, ],
        YY3[idx_test, ],
        YY4[idx_test, ]
      ),
      Z = Z[idx_test, ],
      n_class = n_class
    )
    
    # Train momlin
    tic <- Sys.time()
    result <- momlin_func(trainData, testData, opts)
    
    # Save U and V with k-fold dimensions
    U_model[, , k] <- result[[1]]
    for (i in seq_along(YY)) {
      V_model[[i]][, , k] <- result[[2]][[i]]
    }
    
    if (k == 1) {
      lo <- result[[3]]
    }
    
    # Correlation values 
    CCCs_train <- calcCCC(trainData, result[[1]], result[[2]], scale_type = 'std')
    CCCs_train_model[, , k] <- CCCs_train
    
    CCCs_test <- calcCCC(testData, result[[1]], result[[2]], scale_type = 'std')
    CCCs_test_model[, , k] <- CCCs_test
    
    toc <- Sys.time()
    cat(sprintf("completed in %f seconds]\n", as.numeric(difftime(toc, Sys.time(), units = "secs"))))
  }
  
  
  # # convergence loss 
  # loss <- lo
  # 
  # # k-fold raw
  # U_model_kF <- U_model
  # V_model_kF <- V_model
  # 
  # CCCs_train_model_kF <- CCCs_train_model
  # CCCs_test_model_kF <- CCCs_test_model
  
  # Aggregate results across folds
  U_mean <- apply(U_model, c(1, 2), mean)
  V_mean <- lapply(seq_along(V_model), function(i) apply(V_model[[i]], c(1, 2), mean))
  
  CCCs_train_mean <- apply(CCCs_train_model, c(1, 2), mean)
  CCCs_test_mean <- apply(CCCs_test_model, c(1, 2), mean)
  
  
  return(list(U_mean = U_mean, V_mean = V_mean, CCCs_train_mean = CCCs_train_mean, CCCs_test_mean = CCCs_test_mean,
              U_model_kF = U_model, V_model_kF = V_model, CCCs_train_model_kF = CCCs_train_model, 
              CCCs_test_model_kF = CCCs_test_model, loss = lo, fold_indices = fold_indices))
}

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()


# Extract and aggregate results from parallel computation
for (ii in seq_along(results)) {
  U_mean_i[, , ii] <- results[[ii]]$U_mean
  for (i in seq_along(YY)) {
    V_mean_i[[i]][, , ii] <- results[[ii]]$V_mean[[i]]
  }
  
  CCCs_train_mean_i[, , ii] <- results[[ii]]$CCCs_train_mean
  CCCs_test_mean_i[, , ii] <- results[[ii]]$CCCs_test_mean
  
  for (f in 1:k_fold) {
    train_id <- results[[ii]]$fold_indices[[f]]$train
    test_id <- results[[ii]]$fold_indices[[f]]$test
    
    rCVF_train[train_id, f, ii] <- TRUE
    rCVF_test[test_id, f, ii] <- TRUE
    
  }
  
}


# # Save the results
# saveRDS(list(U_mean_i = U_mean_i, V_mean_i = V_mean_i, CCCs_train_mean_i = CCCs_train_mean_i, 
#              CCCs_test_mean_i = CCCs_test_mean_i), file = paste0(ResultsFile, "/momlin_results.rds"))

# readRDS(paste0(current_dir, "/momlin_out/momlin_results.rds"))


#=====================================================================

# Calculate Canonical Weights
# U_mean <- sweep(apply(U_mean_i, c(1, 2), sum), 2, sqrt(colSums(apply(U_mean_i, c(1, 2), sum)^2)), "/")
U_mean <- apply(U_mean_i, c(1, 2), mean)
rownames(U_mean) <- mRNA_name

V_mean <- lapply(seq_along(V_mean_i), function(i) {
          data <- apply(V_mean_i[[i]], c(1, 2), mean)
          rownames(data) <- YY_names[[i]]
          return(data)
          })


# mean Correlation over iterations
CCCs_mean_train <- apply(CCCs_train_mean_i, c(1, 2), mean)
rownames(CCCs_mean_train) <- c("rna_dna", "rna_clinic","rna_TiME", "rna_path")
colnames(CCCs_mean_train) <- c("pCR", "RCB-i","RCB-ii", "RCB-iii")

CCCs_mean_test <- apply(CCCs_test_mean_i, c(1, 2), mean)
rownames(CCCs_mean_test) <- c("rna_dna", "rna_clinic","rna_TiME", "rna_path")
colnames(CCCs_mean_test) <- c("pCR", "RCB-i","RCB-ii", "RCB-iii")


# Load required libraries
library(reshape2)
library(ggplot2)

plot_heatmap <- function(data, cutof = 0, y_names = NULL, c_limit = c(-0.03, 0.03), title, xlab, ylab) {
  
  s.Features <- reshape2::melt(data[, ncol(data):1])
  keep <- abs(s.Features$value) >= cutof
  
  heatmap <- ggplot(s.Features[keep, ], aes(Var1, Var2, fill = value)) + 
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c_limit) +
    labs(title = title, x = xlab, y = ylab) + 
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "lines")  # Adjust the bottom margin
    ) 
  
  # If y_labels are provided, customize the Y-axis labels
  if (!is.null(y_names)) {
    heatmap <- heatmap + scale_y_continuous(breaks = 1:length(y_names), labels = y_names)
  }
  
  return(heatmap)
}


U_mean_plot <- plot_heatmap(data = U_mean, cutof = .0019,
                            y_names = c("RCB-III", "RCB-II", "RCB-I", "pCR"),
                            c_limit = c(-.009, 0.009), title = "Weights for mRNA", 
                            xlab = "mRNA", ylab = " ")


V1_mean_plot <- plot_heatmap(data = V_mean[[1]], cutof = .0019,
                            y_names = c("RCB-III", "RCB-II", "RCB-I", "pCR"),
                            c_limit = c(-.05, 0.05), title = "Weights for DNA", 
                            xlab = "DNA", ylab = " ")


V2_mean_plot <- plot_heatmap(data = V_mean[[2]], cutof = .0019,
                             y_names = c("RCB-III", "RCB-II", "RCB-I", "pCR"),
                             c_limit = c(-.07, 0.07), title = "Weights for Clinical attibutes", 
                             xlab = "Clinical", ylab = " ")


V3_mean_plot <- plot_heatmap(data = V_mean[[3]], cutof = .008,
                             y_names = c("RCB-III", "RCB-II", "RCB-I", "pCR"),
                             c_limit = c(-.04, 0.04), title = "Weights for TiME", 
                             xlab = "TiME", ylab = " ")



V4_mean_plot <- plot_heatmap(data = V_mean[[4]], cutof = .004,
                             y_names = c("RCB-III", "RCB-II", "RCB-I", "pCR"),
                             c_limit = c(-.04, 0.04), title = "Weights for Pathways", 
                             xlab = "Pathways", ylab = " ")



# Arrange plots in a single figure
grid.arrange(U_mean_plot, V1_mean_plot, V2_mean_plot, V3_mean_plot, V4_mean_plot, ncol = 3)


# Plotting train/test loss for each iteration during optimization
# Assuming lossC is a list with train, test, and iter elements
loss_data <- data.frame(
  Iteration = results[[1]]$loss$iter[-1],
  Train_Loss = results[[1]]$loss$train[-1],
  Test_Loss = results[[1]]$loss$test[-1]
)

loss_plot <- ggplot(loss_data, aes(x = Iteration)) + 
  geom_line(aes(y = Train_Loss, color = "Train"), size = 1.2) + 
  geom_line(aes(y = Test_Loss, color = "Test"), linetype = "dashed", size = 1.2) +
  scale_color_manual(values = c("Train" = "black", "Test" = "red")) +
  labs(title = "Train-test convergence curve", x = "Iteration", y = "Loss") + 
  theme_minimal() + 
  theme(legend.title = element_blank(), text = element_text(size = 16)) 

print(loss_plot)

# Entire work space
# save.image(file = paste0(current_dir,"/momlin_out/momlin_work_space_fromR_macbook.RData"))


# ================ Plot sample distributions =======================
# 
X_cov <- XX1 %*% U_mean
Y1_cov <- YY1 %*% V_mean[[1]]
Y2_cov <- YY2 %*% V_mean[[2]]
Y3_cov <- YY3 %*% V_mean[[3]]
Y4_cov <- YY4 %*% V_mean[[4]]

plot(X_cov[,1], X_cov[,3], col = clasID)

plot(Y1_cov[,1], Y1_cov[,2], col = clasID)

plot(Y2_cov[,1], Y2_cov[,2], col = clasID)

plot(Y3_cov[,1], Y3_cov[,2], col = clasID)

plot(Y4_cov[,1], Y4_cov[,2], col = clasID)



# TSNE visua;ization

library(Rtsne)
library(ggplot2)

# Running t-SNE
set.seed(12)
# tsne_result <- Rtsne(cbind(X_cov, Y4_cov), dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
tsne_result <- Rtsne(cbind(X_cov, Y3_cov, Y4_cov), dims = 2, perplexity = 20, verbose = TRUE, max_iter = 500)

# Extract t-SNE coordinates
tsne_coords <- as.data.frame(tsne_result$Y)
colnames(tsne_coords) <- c("Dim1", "Dim2")
plot(tsne_coords$Dim1, tsne_coords$Dim2, col = clasID, pch=20, cex=1.6)


# XX_pca <- prcomp(XX1)
# tsne_result <- Rtsne(XX_pca$x, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
# 
# # Extract t-SNE coordinates
# tsne_coords <- as.data.frame(tsne_result$Y)
# colnames(tsne_coords) <- c("Dim1", "Dim2")
# plot(tsne_coords$Dim1, tsne_coords$Dim2, col = clasID)









