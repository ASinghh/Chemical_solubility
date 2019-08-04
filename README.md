# Comparision of Regression methods for modelling Chemical solubility
The goal of this project is to predict a response variable using a variety of techniques. As a first
step, an exploratory analysis of the data is conducted. After pre-processing the data, each of the
following models are trained and tested using k-fold cross-validation: linear regression, Lasso, Ridge
and Principal Components Regression (PCR). The final section analyzes the model results and the
important predictors.


## Modelling Software
R

### Packages Used
* glmnet
* corrplot
* corrplot
* car
* pls
* lmtest
* gvlma
* lmtest
* gvlma
* MASS


## Hardware
* i7 Quadcore Processor
* 16 GBs RAM
 


# Results
We Observed the following results from the modelling exercise

## Test RMSE
* OLS     0.722622
* Lasso   0.669559
* Ridge   0.679769
* PCR     0.694552

## Importance of predictiors
The 5 most important predictors, based on the P-score are given below, in their order of importance,

1. Avg_surf aceArea
2. FP067
3. NumNonHBonds
4. FP068
5. FP091

# Conclusion
 For our use-case, we found Ordinary Least Square method to be the most effective Linear method to model chemical solubility,
 for the given dataset.For detailed project report, please read <strong>predicting-chemical-solubility.pdf</strong>.







