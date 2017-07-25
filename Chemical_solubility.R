source("my_functions.R")
library(glmnet)
install.packages("corrplot")
library(corrplot)
library(car)
library(pls)
#install.packages('lmtest')
#install.packages('gvlma')
library(lmtest)
library(gvlma)
library(MASS)


#for reproducibility
set.seed(12345)

all_data = read.csv("solubility_data.csv",  header = TRUE)
response_name = "solubility"

#data description
# factors: 208 binary "fingerprints" that indicate the presence or absence of a particular chemical sub-structure
# counts: 16 count descriptors (such as the number of bonds or the number of Bromine atoms)
# continuous: 4 continuous descriptors (such as molecular weight or surface area).
# response: solubility

#see all the data
options(max.print=10000)
summary(all_data)

## Find the predictors that are fingerprints
fpVars <- names(all_data)[grepl("FP", names(all_data))]

#look at corrplots:
#don't look at factors
all_data_no_factors <- all_data[ , !(names(all_data) %in% fpVars)]
cor_data = cor(all_data_no_factors)#don't look at factors

corrplot(cor_data, type = "upper")

#combine linearly related predictors
all_data$Avg_surfaceArea = (all_data$SurfaceArea1 + all_data$SurfaceArea2)/2
all_data$Avg_bonds = (all_data$NumAromaticBonds + all_data$NumMultBonds)/2
drops <- c("SurfaceArea1", "SurfaceArea2", "NumAromaticBonds", "NumMultBonds")
all_data <- all_data[ , !(names(all_data) %in% drops)]

#create dummy variables for categorical data
all_data <- cbind(all_data, model.matrix( ~ as.factor(NumRings) - 1, data = all_data))
#all_data <- cbind(all_data, model.matrix( ~ as.factor(NumHalogen) - 1, data = all_data))
#all_data <- cbind(all_data, model.matrix( ~ as.factor(NumChlorine) - 1, data = all_data))
all_data <- cbind(all_data, model.matrix( ~ as.factor(NumSulfer) - 1, data = all_data))
#all_data <- cbind(all_data, model.matrix( ~ as.factor(NumOxygen) - 1, data = all_data))
all_data <- cbind(all_data, model.matrix( ~ as.factor(NumNitrogen) - 1, data = all_data))
#all_data <- cbind(all_data, model.matrix( ~ as.factor(NumDblBonds) - 1, data = all_data))
#all_data <- cbind(all_data, model.matrix( ~ as.factor(NumRotBonds) - 1, data = all_data))

#now drop the original predictors
#predictors_with_factors <- c("NumRings", "NumHalogen", "NumChlorine", "NumSulfer", "NumOxygen", "NumNitrogen", "NumAromaticBonds", "NumDblBonds", "NumRotBonds")
predictors_with_factors <- c("NumRings", "NumSulfer","NumNitrogen", "NumAromaticBonds")
all_data <- all_data[ , !(names(all_data) %in% predictors_with_factors)]

#remove or combine highly-correlated predictors
all_data <- drop_combine_correlated_predictors(all_data, response_name)

#run it again to see if newly-created predictors are also correlated
all_data <-  drop_combine_correlated_predictors(all_data, response_name)

#transform continuous predictors
continous_predictors <- c("Avg_bonds", "Avg_surfaceArea", "HydrophilicFactor", "MolWeight")
hetero_results <- list()
for (i in 1:length(continous_predictors))
{
  
  model2 = lm(
      as.formula(paste("solubility", "~",
                       paste(continous_predictors[i]),
                       sep = "" )),data=all_data)
    
  hetero_results[[i]] <- gvlma(model2)
}

#Avg_bonds and MolWeight are heteroscedastic
plot(all_data$solubility, all_data$MolWeight)
#log transform molweight
all_data$MolWeight <- log(all_data$MolWeight)
plot(all_data$solubility, all_data$MolWeight)

boxplot.matrix(as.matrix(all_data[ , (names(all_data) %in% continous_predictors)]))

all_data[ , (names(all_data) %in% continous_predictors)] <- scale(all_data[ , (names(all_data) %in% continous_predictors)], FALSE, TRUE)
boxplot.matrix(as.matrix(all_data[ , (names(all_data) %in% continous_predictors)]))


#look for outliers
basic_fit <- lm(solubility ~., data=all_data)
outliers <- influencePlot(basic_fit)
n <- 3
std_res <- as.numeric(tail(row.names(outliers[order(outliers$StudRes), ]), n))
outlierRows <- all_data[which.names(names=std_res, rownames(all_data)),]
outTest <- outlierTest(basic_fit)

#drop high-leverage points
all_data <- all_data[-which.names(names=std_res, rownames(all_data)),]


##MODELING
#OLS
no.of.folds = 5
index.values = sample(1:no.of.folds, size = dim(all_data)[1], replace = TRUE)
test.mse = rep(0, no.of.folds)
for (i in 1:no.of.folds)
{
  index.out            = which(index.values == i)                             ### These are the indices of the rows that will be left out.
  left.out.data        = all_data[  index.out, ]                                  ### This subset of the data is left out. (about 1/10)
  left.in.data         = all_data[ -index.out, ]                                  ### This subset of the data is used to get our regression model. (about 9/10)
  tmp.lm               = lm(solubility ~ ., data = left.in.data)                 ### Perform regression using the data that is left in.
  tmp.predicted.values = predict(tmp.lm, newdata = left.out.data)             ### Predict the y values for the data that was left out
  test.mse[i]          = mean((left.out.data$solubility - tmp.predicted.values)^2)  ### Get one of the test.mse's
}
#drops 5 predictors
test_error_ols = sqrt(mean(test.mse))

#Lasso
x_data <- as.matrix(all_data[, !(names(all_data) %in% response_name)])
lasso_model = cv.glmnet(x=x_data, y=all_data$solubility, alpha = 1, nfolds = 5)
test_error_lasso = sqrt(min(lasso_model$cvm))
coef_lasso <- coef(lasso_model)
#Ridge
ridge_model = cv.glmnet(x = x_data, y = all_data$solubility, alpha = 0, nfolds = 5)
test_error_ridge = sqrt(min(ridge_model$cvm))

#PCA
drops_unit_var <- get_unit_var_predictors(all_data,0)
n_col <- ncol(all_data)

pcr_model = pcr(solubility ~ ., data = all_data, scale = TRUE, validation = "CV", ncomp = n_col - 1, segments=5)
pcr_cv = RMSEP(pcr_model, estimate = "CV")

###
### Let's plot it.  The [-1] leaves out the intercept only model.
###

plot(pcr_cv$val[-1], pch = 19, type = "b", ylab = "Test RMSE", xlab = "Number of Components")

best_comp = which.min(pcr_cv$val[-1])
abline(v = best_comp, col = "red")
test_error_pcr <- pcr_cv$val[ best_comp ]

coef_pcr=abs (coef(pcr_model, ncomp = 1, intercept = FALSE) )
coef_pcr = data.frame(coef_pcr)
colnames(coef_pcr)
newdata <- coef_pcr[order( coef_pcr$solubility,decreasing = TRUE),,drop=FALSE] 
top_5 <- rownames(newdata)[1:5]
