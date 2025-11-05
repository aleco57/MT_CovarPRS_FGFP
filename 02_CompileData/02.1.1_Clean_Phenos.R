#### Lets generate some R code to QC the variables we are looking at

# Load necessary libraries
library(ggplot2)
library(performance)

# Fit the linear model
i = 9
#Extract name of pheno and prs
pheno <- data4prs$matchedvars[i, "pheno"]
prs <- paste0(data4prs$matchedvars[i, "prs"], ".Pt_5e.08")
m <- lm(reformulate(c(confounders, prs), response = pheno), data = data4prs$pheno_covariate_prs)

# Summary of the model
summary(model)

# 1. Distribution of x and y (with outlier detection via boxplots)
par(mfrow = c(1, 2))  # Set plot area to show 2 plots side-by-side
boxplot(data4prs$pheno_covariate_prs[[pheno]], main = "Boxplot of x", col = "lightblue", horizontal = TRUE)
boxplot(data4prs$pheno_covariate_prs[[prs]], main = "Boxplot of y", col = "lightgreen", horizontal = TRUE)

# 2. Histograms with density plots
par(mfrow = c(1, 2))  # Reset to 2 plots side-by-side
hist(data4prs$pheno_covariate_prs[[pheno]], breaks = 20, col = "lightblue", main = "Histogram of x", xlab = "x", freq = FALSE)
lines(density(na.omit(data4prs$pheno_covariate_prs[[pheno]])), col = "blue", lwd = 2)

hist(data4prs$pheno_covariate_prs[[prs]], breaks = 20, col = "lightgreen", main = "Histogram of y", xlab = "y", freq = FALSE)
lines(density(na.omit(data4prs$pheno_covariate_prs[[prs]])), col = "darkgreen", lwd = 2)

# 3. Scatter plot with regression line
ggplotRegression <- function (fit, var) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = var, y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(m)

# Optional: Diagnostic plots of the model
check_model(m)


#I want to look at some of the models of GGT PRS -> Microbiome
mt <- "G_unclassified_P_Firmicutes_RNTRes"
lm <- lm(reformulate(c(pheno, age_sex), response = mt), data = micropheno_df_scale)

plot(x = micropheno_df[[pheno]], y = micropheno_df[[mt]], 
     xlab = pheno, ylab = mt, main = paste("Scatter plot of", pheno, "vs", mt))

check_model(lm)
ggplotRegression(lm, var = pheno)

#Lets filter out outliers in the model to see how much this effects the model
lm1 <- lm(reformulate(c(pheno, age_sex), response = mt), data = filter(micropheno_df_scale, Gamma.GT_UL < 10))
ggplotRegression(lm1, var = pheno)
check_model(lm)


lm2 <- lm(reformulate(c(prs, age_sex), response = mt), data = prsmicro_df_scale)
check_model(lm2)
ggplotRegression(lm2, var = prs)
