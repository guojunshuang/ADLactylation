
rm(list = ls())
gc()
library(glmnet)
library(foreign)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(broom)

lactylation_matrix <- read.table("GSE132903 as the training cohort.txt", header = TRUE, row.names = NULL, sep = "\t")
lactylation_matrix <- lactylation_matrix[, -c(1, 2)]
data <- lactylation_matrix
set.seed(123)
data <- na.omit(data)
y<-as.matrix(data[,1])
print(y)
x<-as.matrix(data[,c(2:9)])
print(x)
f1 =glmnet(x, y, family="binomial", nlambda=100, alpha=1) 
print(f1)
cvfit=cv.glmnet(x,y,family="binomial")
tidy_df <- broom::tidy(f1)
tidy_cvdf <- broom::tidy(cvfit)
head(tidy_df)
head(tidy_cvdf)

tmp <- as_tibble(as.matrix(coef(f1)), rownames = "coef") %>% 
  pivot_longer(cols = -coef, 
               names_to = "variable", 
               names_transform = list(variable = parse_number), 
               values_to = "value") %>% 
  filter(coef != "(Intercept)") %>% 
  group_by(variable) %>% 
  mutate(lambda = f1$lambda[variable + 1], 
         norm = sum(abs(value)))

ggplot(tmp, aes(norm, value, color = coef, group = coef)) +
  geom_line(size = 1.2) +
  labs(x = "Log Lambda", y = "Coefficients") +
  theme_bw()

cvfit=cv.glmnet(x,y,family="binomial")
plot(cvfit)

cvfit$lambda.min
cvfit$lambda.1se
l.coef2<-coef(cvfit$glmnet.fit,s=0.01819986,exact= F)
l.coef1<-coef(cvfit$glmnet.fit,s=0.1065969,exact= F)
l.coef1
l.coef2


library(pROC)
library(openxlsx)
library(ggplot2)
library(lattice)
library(caret)
library(class)
library(kknn)
library(e1071)
library(kernlab)
library(reshape2)
library(pROC)
library(bnlearn)
library(MASS)
library(klaR)
library(rpart) 
library(grid)
library(libcoin)
library(mvtnorm)
library(partykit) 
library(ElemStatLearn) 
library(randomForest) 
library(xgboost) 
library(lattice)
library(caret) 
library(mlbench)
library(neuralnet)
library(nnet)
library(vcd)
library(Rcpp) 
library(RSNNS)
library(keras)
library(rpart.plot)
library(pROC)

AE <- data

print(AE)
set.seed(123)

features <- c("LDHA", "ACSF2", "MECP2", "ARF1")
target <- "Group"

n_bootstrap <- 1000
bootstrap_results <- vector("list", n_bootstrap)
for (i in 1:n_bootstrap) {

  bootstrap_sample <- AE[sample(nrow(AE), replace = TRUE), ]
  
  model <- rpart(formula = Group ~ ., data = bootstrap_sample[, c(features, target)], method = "class")
  
  predictions <- predict(model, newdata = AE[, features], type = "class")
  
  bootstrap_results[[i]] <- predictions
}
model
print(model)
rpart.plot(model)
library(dplyr)
importances <- varImp(model)
importances <- importances %>%
  arrange(desc(Overall))
importance_plot <- ggplot(importances, aes(x = reorder(rownames(importances), Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "#FFCB5B") +
  labs(x = "Features", y = "Overall Importance", title = "Variable Importance") +
  geom_text(aes(label = round(Overall, 2)), hjust = -0.1, color = "black") + 
  theme_minimal() +
  coord_flip()
print(importance_plot)
plot(importance_plot)

library(pROC)
probabilities <- predict(model, newdata = AE[, features], type = "prob")
positive_probs <- probabilities[, "1"]
roc_objDT <- roc(AE$Group, positive_probs)
auc <- auc(roc_objDT)
auc_ci <- ci.auc(roc_objDT)
plot(roc_objDT, main = "ROC Curve")
confusion_matrix <- table(AE$Group, ifelse(positive_probs >= 0.5, "1", "0"))
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
negative_precision <- sum(ifelse(positive_probs < 0.5, "0", "1") == AE$Group & AE$Group == "0") / sum(AE$Group == "0")
sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])
f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)

####################################################################

library(randomForest)
library(pROC)
set.seed(123)
AE$Group <- factor(AE$Group)

n_bootstrap <- 1000
bootstrap_results <- vector("list", n_bootstrap)
for (i in 1:n_bootstrap) {
  bootstrap_sample <- AE[sample(nrow(AE), replace = TRUE), ]
  model_rf <- randomForest(formula = Group ~ ., data = bootstrap_sample[, c(features, target)], ntree = 100)

  predictions_rf <- predict(model_rf, newdata = AE[, features], type = "response")

  bootstrap_results[[i]] <- predictions_rf
}
model_rf

plot(model_rf, main = "Random Forest", cex.main = 0.9)

var_importance <- importance(model_rf)
var_names <- rownames(var_importance)
var_importance <- data.frame(Variables = var_names, Importance = var_importance[, "MeanDecreaseGini"])
var_importance <- var_importance[order(var_importance$Importance, decreasing = F), ]

barplot(var_importance$Importance, names.arg = var_importance$Variables, horiz = TRUE,
        main = "Variable Importance", xlab = "Mean Decrease Gini", xlim = c(0, 15), col = "#F2A1A7")

ggplot(var_importance, aes(x = reorder(Variables, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#F2A1A7") +
  geom_text(aes(label = round(Importance, 2)), hjust = -0.1, color = "black") +  
  labs(x = "Variables", y = "Mean Decrease Gini", title = "Variable Importance") +
  theme_minimal() +
  coord_flip()

probabilities_rf <- predict(model_rf, newdata = AE[, features], type = "prob")
positive_probs_rf <- probabilities_rf[, "1"]
roc_obj_rf <- roc(AE$Group, positive_probs_rf)
auc_rf <- auc(roc_obj_rf)
auc_rf
auc_ci_rf <- ci.auc(roc_obj_rf)
auc_ci_rf

confusion_matrix_rf <- table(AE$Group, ifelse(positive_probs_rf >= 0.5, "1", "0"))
AE$RF_predict <- ifelse(positive_probs_rf >= 0.5, "1", "0")
confusion_matrix_RF <- table(AE$Group, AE$RF_predict)

accuracy_rf <- sum(diag(confusion_matrix_rf)) / sum(confusion_matrix_rf)
precision_rf <- confusion_matrix_rf[2, 2] / sum(confusion_matrix_rf[, 2])
negative_precision_rf <- sum(ifelse(positive_probs_rf < 0.5, "0", "1") == AE$Group & AE$Group == "0") / sum(AE$Group == "0")
sensitivity_rf <- confusion_matrix_rf[2, 2] / sum(confusion_matrix_rf[2, ])
specificity_rf <- confusion_matrix_rf[1, 1] / sum(confusion_matrix_rf[1, ])
f1_score_rf <- 2 * (precision_rf * sensitivity_rf) / (precision_rf + sensitivity_rf)

####################################################################
library(xgboost)
library(pROC)

AE$Group <- ifelse(AE$Group == "0", 0, 1)
n_bootstrap <- 1000
bootstrap_results <- vector("list", n_bootstrap)
for (i in 1:n_bootstrap) {
  bootstrap_sample <- AE[sample(nrow(AE), replace = TRUE), ]
  
  dtrain <- xgb.DMatrix(data = as.matrix(bootstrap_sample[, features]), label = bootstrap_sample[, target], missing = NA)
  model_xgb <- xgboost(data = dtrain, nrounds = 100, objective = "binary:logistic")
  
  dtest <- xgb.DMatrix(data = as.matrix(AE[, features]), missing = NA)
  predictions_xgb <- predict(model_xgb, dtest)
  
  bootstrap_results[[i]] <- predictions_xgb
}

model_xgb

var_importance <- xgb.importance(feature_names = features, model = model_xgb)
var_importance
var_names <- var_importance$Feature
importance_scores <- var_importance$Gain
importance_scores

var_importance_df <- data.frame(Variables = var_names, Importance = importance_scores)
library(ggplot2)
var_importance_df <- var_importance_df[order(var_importance_df$Importance, decreasing = TRUE), ]
ggplot(var_importance_df, aes(x = reorder(Variables, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#F8B072") +
  labs(x = "Variables", y = "Importance", title = "Variable Importance") +
  geom_text(aes(label = round(Importance, 3)), hjust = -0.1, color = "black") + 
  coord_flip()+theme_minimal() 

probabilities_xgb <- predict(model_xgb, dtest)
positive_probs_xgb <- probabilities_xgb
roc_obj_xgb <- roc(AE$Group, positive_probs_xgb)
auc_xgb <- auc(roc_obj_xgb)
auc_xgb
auc_ci_xgb <- ci.auc(roc_obj_xgb)
auc_ci_xgb

confusion_matrix_xgb <- table(AE$Group, ifelse(positive_probs_xgb >= 0.5, 1, 0))
AE$XGB_predict <- ifelse(positive_probs_xgb >= 0.5, "1", "0")
confusion_matrix_XGB <- table(AE$Group, AE$XGB_predict)
accuracy_xgb <- sum(diag(confusion_matrix_xgb)) / sum(confusion_matrix_xgb)
precision_xgb <- confusion_matrix_xgb[2, 2] / sum(confusion_matrix_xgb[, 2])
negative_precision_xgb <- sum(ifelse(positive_probs_xgb < 0.5, 0, 1) == AE$Group & AE$Group == 0) / sum(AE$Group == 0)
sensitivity_xgb <- confusion_matrix_xgb[2, 2] / sum(confusion_matrix_xgb[2, ])
specificity_xgb <- confusion_matrix_xgb[1, 1] / sum(confusion_matrix_xgb[1, ])
f1_score_xgb <- 2 * (precision_xgb * sensitivity_xgb) / (precision_xgb + sensitivity_xgb)

####################################################################
library(class)
library(pROC)
AE$Group <- ifelse(AE$Group == "0", 0, 1)
n_bootstrap <- 1000
bootstrap_results <- vector("list", n_bootstrap)
for (i in 1:n_bootstrap) {
  bootstrap_sample <- AE[sample(nrow(AE), replace = TRUE), ]
  model_knn <- knn(train = bootstrap_sample[, features], test = AE[, features], cl = bootstrap_sample[, target], k = 5)

  bootstrap_results[[i]] <- model_knn
}
model_knn

confusion_matrix_knn <- table(AE$Group, bootstrap_results[[n_bootstrap]])
accuracy_knn <- sum(diag(confusion_matrix_knn)) / sum(confusion_matrix_knn)
precision_knn <- confusion_matrix_knn[2, 2] / sum(confusion_matrix_knn[, 2])
negative_precision_knn <- sum(ifelse(bootstrap_results[[n_bootstrap]] == 0 & AE$Group == 0, 1, 0)) / sum(AE$Group == 0)
sensitivity_knn <- confusion_matrix_knn[2, 2] / sum(confusion_matrix_knn[2, ])
specificity_knn <- confusion_matrix_knn[1, 1] / sum(confusion_matrix_knn[1, ])
f1_score_knn <- 2 * (precision_knn * sensitivity_knn) / (precision_knn + sensitivity_knn)

probabilities_knn <- as.numeric(bootstrap_results[[n_bootstrap]])

roc_obj_knn <- roc(AE$Group, probabilities_knn)
auc_knn <- auc(roc_obj_knn)
auc_ci_knn <- ci.auc(roc_obj_knn)
auc_ci_knn
auc_knn

####################################################################
library(pROC)
AE$Group <- ifelse(AE$Group == "0", 0, 1)
n_bootstrap <- 1000
nb_bootstrap_results <- vector("list", n_bootstrap)
model_nb <- naiveBayes(formula = Group ~ ., data = AE[, c(features, target)])
for (i in 1:n_bootstrap) {
  bootstrap_sample <- AE[sample(nrow(AE), replace = TRUE), ]

  model_nb_bootstrap <- naiveBayes(formula = Group ~ ., data = bootstrap_sample[, c(features, target)])

  predictions_nb <- predict(model_nb_bootstrap, newdata = AE[, features])

  nb_bootstrap_results[[i]] <- predictions_nb
}

probabilities_nb <- predict(model_nb, newdata = AE[, features], type = "raw")
positive_probs_nb <- probabilities_nb[, "1"]
roc_obj_nb <- roc(AE$Group, positive_probs_nb)
auc_nb <- auc(roc_obj_nb)
auc_nb
auc_ci_nb <- ci.auc(roc_obj_nb)
auc_ci_nb

plot(roc_obj_nb, main = "ROC Curve ", print.auc = TRUE, auc.polygon = FALSE)

confusion_matrix_nb <- table(AE$Group, ifelse(positive_probs_nb >= 0.5, "1", "0"))
accuracy_nb <- sum(diag(confusion_matrix_nb)) / sum(confusion_matrix_nb)
precision_nb <- confusion_matrix_nb[2, 2] / sum(confusion_matrix_nb[, 2])
negative_precision_nb <- sum(ifelse(positive_probs_nb < 0.5, "0", "1") == AE$Group & AE$Group == 0) / sum(AE$Group == 0)
sensitivity_nb <- confusion_matrix_nb[2, 2] / sum(confusion_matrix_nb[2, ])
specificity_nb <- confusion_matrix_nb[1, 1] / sum(confusion_matrix_nb[1, ])
f1_score_nb <- 2 * (precision_nb * sensitivity_nb) / (precision_nb + sensitivity_nb)

#######################################################################
library(pROC)
library(nnet)
AE$Group <- as.factor(AE$Group)
n_bootstrap <- 1000

roc_data_nn <- vector("list", n_bootstrap)
for (i in 1:n_bootstrap) {
  bootstrap_sample <- AE[sample(nrow(AE), replace = TRUE), ]

  model_nn_bootstrap <- nnet(formula = Group ~ ., data = bootstrap_sample[, c(features, target)], size = 5, decay = 0.1)

  predictions_nn <- predict(model_nn_bootstrap, newdata = AE[, features], type = "raw")

  predictions_nn <- as.numeric(predictions_nn)

  roc_data_nn[[i]] <- roc(AE$Group, predictions_nn, ci = TRUE)  
}

auc_ci <- sapply(roc_data_nn, function(x) x$ci) 
plot(roc_data_nn[[1]], print.auc = T, auc.polygon = F, max.auc.polygon = F)

confusion_matrix <- table(AE$Group, predictions_nn > 0.5)
tp <- confusion_matrix[2, 2]  
tn <- confusion_matrix[1, 1]  
fp <- confusion_matrix[1, 2]  
fn <- confusion_matrix[2, 1] 

accuracy <- (tp + tn) / sum(confusion_matrix)
precision <- tp / (tp + fp)
negative_predictive_value <- tn / (tn + fn)
sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
f1_score <- 2 * precision * sensitivity / (precision + sensitivity)

################   lightgbm   
library(PRROC)
library(lightgbm)
AE$Group <- ifelse(AE$Group == "0", 0, 1)
n_bootstrap <- 1000
bootstrap_results <- vector("list", n_bootstrap)
for (i in 1:n_bootstrap) {

  bootstrap_sample <- AE[sample(nrow(AE), replace = TRUE), ]
  
  lgb_data <- lgb.Dataset(data = as.matrix(bootstrap_sample[, features]),
                          label = bootstrap_sample[, target])

  params <- list(objective = "binary", nrounds = 100)  
  model_lgb <- lgb.train(params = params, data = lgb_data)

  predictions_lgb <- predict(model_lgb, data = as.matrix(AE[, features]))  

  bootstrap_results[[i]] <- predictions_lgb
}

library(ROCR) 
var_importance <- lgb.importance(model_lgb)
var_names <- var_importance$Feature
importance_scores <- var_importance$Gain
barplot(importance_scores, names.arg = var_names, horiz = TRUE,
        main = "Variable Importance", xlab = "Gain", xlim = c(0, 0.8))
text(importance_scores + 0.02, 1:length(var_names), round(importance_scores, 2), pos = 4)

var_importance_df <- data.frame(Variables = var_names, Importance = importance_scores)
var_importance_df <- var_importance_df[order(var_importance_df$Importance, decreasing = TRUE), ]

ggplot(var_importance_df, aes(x = reorder(Variables, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#FA8072") +
  geom_text(aes(label = round(Importance, 2)), hjust = -0.1, color = "black") + 
  labs(x = "Variables", y = "Gain", title = "Variable Importance") +
  theme_minimal() +
  coord_flip()

probabilities_lgb <- predict(model_lgb, data = as.matrix(AE[, features]))
positive_probs_lgb <- 1 / (1 + exp(-probabilities_lgb))
roc_obj_lgb <- roc(AE$Group, positive_probs_lgb)
auc_lgb <- auc(roc_obj_lgb)
auc_lgb
auc_ci_lgb <- ci.auc(roc_obj_lgb)
auc_ci_lgb

confusion_matrix_lgb <- table(AE$Group, ifelse(probabilities_lgb >= 0.5, 1, 0))
accuracy_lgb <- sum(diag(confusion_matrix_lgb)) / sum(confusion_matrix_lgb)
precision_lgb <- confusion_matrix_lgb[2, 2] / sum(confusion_matrix_lgb[, 2])
negative_precision_lgb <- sum(ifelse(probabilities_lgb < 0.5, 0, 1) == AE$Group & AE$Group == 0) / sum(AE$Group == 0)
sensitivity_lgb <- confusion_matrix_lgb[2, 2] / sum(confusion_matrix_lgb[2, ])
specificity_lgb <- confusion_matrix_lgb[1, 1] / sum(confusion_matrix_lgb[1, ])
f1_score_lgb <- 2 * (precision_lgb * sensitivity_lgb) / (precision_lgb + sensitivity_lgb)




####################################################################
#logistic
library(boot)

boot_func <- function(data, index) {
  boot_sample <- data[index, ]

  model <- glm(as.formula(paste(target, "~", paste(features, collapse = "+"))), 
               data = boot_sample, family = binomial)

  predicted_probs <- predict(model, data, type = "response")

  return(predicted_probs)
}

bootstrap_results <- boot(data = AE, statistic = boot_func, R = 1000)
average_probs <- colMeans(bootstrap_results$t)
roc_objLR <- roc(AE$Group, average_probs)
auc <- auc(roc_objLR)
auc
auc_ciLR <- ci.auc(roc_objLR)
auc_ciLR

confusion_matrix <- table(AE$Group, ifelse(average_probs >= 0.5, "Positive", "Negative"))
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

precision <- confusion_matrix[2, 2] / (confusion_matrix[2, 2] + confusion_matrix[1, 2])
npv <- confusion_matrix[1, 1] / (confusion_matrix[1, 1] + confusion_matrix[2, 1])
sensitivity <- confusion_matrix[2, 2] / (confusion_matrix[2, 2] + confusion_matrix[2, 1])
specificity <- confusion_matrix[1, 1] / (confusion_matrix[1, 1] + confusion_matrix[1, 2])
f1_score <- 2 * precision * sensitivity / (precision + sensitivity)


#roc
####################################################################
plot(roc_objDT, main = "ROC Curve",auc.polygon = F, 
     legacy.axes = TRUE,col = "#FFCB5B")
lines(roc_obj_rf,col = "#F2A1A7")
lines(roc_obj_xgb,col = "#F8B072")
lines(roc_obj_knn,col = "#A8D3A0")
lines(roc_obj_nb,col = "#7DC69B")
lines(roc_data_nn[[1]],col = "#C9A1CA")
lines(roc_obj_lgb[[1]],col = "#FA8072")
lines(roc_objLR,col = "#9BD7F3")
legend(0.48,0.50, legend = c("DT AUC: 0.851", "RF AUC: 0.954",
                             "XGB AUC: 0.928", "KNN AUC: 0.764", 
                              "NB AUC: 0.819",
                             "NN AUC: 0.822", "LGBM AUC: 0.919","LR AUC: 0.815"),
       col = c("#FFCB5B","#F2A1A7","#F8B072","#A8D3A0","#7DC69B","#C9A1CA","#FA8072","#9BD7F3"), 
       lwd = 2, bg = "white",bty = "n")