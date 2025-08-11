#install pacakges

library(tidyverse)
library(ape)
#install.packages("DHARMa")
library(DHARMa)
library(phylolm)
library(car)
library(phytools)

#data
d_pca = read_csv("Data/d_pca.csv") #pca scores for native and non-native urban plant species
d_pca
phy = read.tree("Data/phy.tre")#phylogeny for species in the PCA data set
plot(phy, type = "fan", show.tip.label = F)



####PC1-All taxa####

#make data frame
d_pca = as.data.frame(d_pca)
row.names(d_pca) = d_pca$tip

#lm for provenance only
lm1<-lm(PC1 ~ provenance, data = d_pca)
summary(lm1)
#plot(lm1)

#phylm for provenance only
phym1 = phylolm(PC1 ~ provenance, data = d_pca, phy)
summary(phym1)


# Compare Q-Q plots 
par(mfrow = c(1, 2)) 

# Q-Q plot for lm model
qqnorm(residuals(lm1), main = "Q-Q Plot: lm")
qqline(residuals(lm1), col = "blue")

# Q-Q plot for phylolm model
qqnorm(residuals(phym1), main = "Q-Q Plot: phylolm")
qqline(residuals(phym1), col = "red")

# Reset plotting area
par(mfrow = c(1, 1))

#Check for phylogenetic signal (yes) TAKES A WHILE TO RUN

#order and match tree tips and rownames
#trait_vector <- d_pca$PC2
#names(trait_vector) <- rownames(d_pca)
#trait_vector <- trait_vector[phy$tip.label]
#head(trait_vector)

#phylosig(phy, trait_vector, method = "lambda", test = TRUE)



#models with other predictors

plot(d_pca$ave_tmean, d_pca$PC1)
plot(d_pca$ave_precip, d_pca$PC1)
plot(d_pca$ave_tmean, d_pca$PC2)
plot(d_pca$ave_precip, d_pca$PC2)

# scale variables
d_pca$ave_tmean = scale(d_pca$ave_tmean)[,1]
d_pca$ave_precip = scale(d_pca$ave_precip)[,1]


#all subsets (almost)
results <- data.frame(model_name=character(), formula=character(), R2=numeric(), AIC=numeric(), deltaAIC=numeric(), stringsAsFactors=FALSE)

model_list <- list()

model_1 <- phylolm(PC1 ~ provenance, data = d_pca, phy = phy)
model_list[['model_1']] <- model_1

model_2 <- phylolm(PC1 ~ provenance + ave_tmean, data = d_pca, phy = phy)
model_list[['model_2']] <- model_2

model_3 <- phylolm(PC1 ~ provenance + ave_precip, data = d_pca, phy = phy)
model_list[['model_3']] <- model_3

model_4 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip, data = d_pca, phy = phy)
model_list[['model_4']] <- model_4

model_5 <- phylolm(PC1 ~ provenance + ave_tmean + ave_tmean:provenance, data = d_pca, phy = phy)
model_list[['model_5']] <- model_5

model_6 <- phylolm(PC1 ~ provenance + ave_tmean + I(ave_tmean^2), data = d_pca, phy = phy)
model_list[['model_6']] <- model_6

model_7 <- phylolm(PC1 ~ provenance + ave_precip + ave_precip:provenance, data = d_pca, phy = phy)
model_list[['model_7']] <- model_7

model_8 <- phylolm(PC1 ~ provenance + ave_precip + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_8']] <- model_8

model_9 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance, data = d_pca, phy = phy)
model_list[['model_9']] <- model_9

model_10 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_precip:provenance, data = d_pca, phy = phy)
model_list[['model_10']] <- model_10

model_11 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + I(ave_tmean^2), data = d_pca, phy = phy)
model_list[['model_11']] <- model_11

model_12 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_12']] <- model_12

model_13 <- phylolm(PC1 ~ provenance + ave_tmean + ave_tmean:provenance + I(ave_tmean^2), data = d_pca, phy = phy)
model_list[['model_13']] <- model_13

model_14 <- phylolm(PC1 ~ provenance + ave_precip + ave_precip:provenance + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_14']] <- model_14

model_15 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance, data = d_pca, phy = phy)
model_list[['model_15']] <- model_15

model_16 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + I(ave_tmean^2), data = d_pca, phy = phy)
model_list[['model_16']] <- model_16

model_17 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_precip:provenance + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_17']] <- model_17

model_18 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_18']] <- model_18

model_19 <- phylolm(PC1 ~ provenance + ave_tmean + ave_tmean:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_19']] <- model_19

model_20 <- phylolm(PC1 ~ provenance + ave_precip + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_20']] <- model_20

model_21 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2), data = d_pca, phy = phy)
model_list[['model_21']] <- model_21

model_22 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_22']] <- model_22

model_23 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_23']] <- model_23

model_24 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_24']] <- model_24

model_25 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)
model_list[['model_25']] <- model_25

aic_values <- sapply(model_list, AIC)
r2_values <- sapply(model_list, function(m) summary(m)$r.squared)
delta_aic <- aic_values - min(aic_values)

results <- data.frame(
          model_name = names(model_list),
          formula = sapply(model_list, function(m) as.character(formula(m))),
          R2 = r2_values,
          AIC = aic_values,
          deltaAIC = delta_aic
           )

results <- results[order(results$deltaAIC), ]
print(results)

# get results from top 5 models
top_models <- head(results$model_name, 5)

# Print summaries for top 5 models
for (model_name in top_models) {
  cat("\n--- Summary of", model_name, "---\n")
  print(summary(get(model_name)))
}

#Variance inflation

vif(lm(PC1 ~ provenance + ave_tmean + ave_precip + ave_precip:provenance + I(ave_precip^2), data = d_pca), type="predictor")

#check residuals with dharma
simulationOutput <- simulateResiduals(fittedModel = model_17, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

#residuals for reg qqplot
qqnorm(residuals(model_17), main = "Q-Q Plot: phylolm")
qqline(residuals(model_17), col = "red")

#PLOTS FOR PC1
# Get means for centering
mean_tmean <- mean(d_pca$ave_tmean, na.rm = TRUE)
mean_precip <- mean(d_pca$ave_precip, na.rm = TRUE)

# 1. Effect of ave_precip (with quadratic term), by provenance
precip_seq <- seq(min(d_pca$ave_precip), max(d_pca$ave_precip), length.out = 100)
new_precip <- expand.grid(
  ave_precip = precip_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_tmean = mean_tmean,
    `I(ave_precip^2)` = ave_precip^2
  )

new_precip$PC1_pred <- predict(model_22, newdata = new_precip)

# Define custom colors
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Plot for ave_precip(scled to full pc1 breadth)
ggplot(new_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 vs Precipitation - All Taxa",
    x = "Average Precipitation",
    y = "Predicted PC1",
    color = "Provenance")+
  coord_cartesian(ylim = c(-3, 3)) +
  theme_minimal()  
#add raw data  
#geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),
                              #alpha = 0.4, size = 1.5)


# 2. Effect of ave_tmean, by provenance
tmean_seq <- seq(min(d_pca$ave_tmean), max(d_pca$ave_tmean), length.out = 100)
new_tmean <- expand.grid(
  ave_tmean = tmean_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_precip = mean_precip,
    `I(ave_precip^2)` = mean_precip^2
  )

new_tmean$PC1_pred <- predict(model_17, newdata = new_tmean)

# 2. Plot for ave_tmean
ggplot(new_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 vs Temperature - All Taxa",
    x = "Average Temperature",
    y = "Predicted PC1",
    color = "Provenance") +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_minimal()

#3 Create a simple prediction dataset with only provenance varying
box_data <- data.frame(
  ave_tmean = mean(d_pca$ave_tmean, na.rm = TRUE),
  ave_precip = mean(d_pca$ave_precip, na.rm = TRUE),
  provenance = c("native", "non_native")
) %>%
  mutate(`I(ave_precip^2)` = ave_precip^2)

# Predict PC1
box_data$PC1_pred <- predict(model_17, newdata = box_data)


# Duplicate rows to simulate a distribution (optional)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]  # simulate 20 values per group

# Add small noise to simulate variability (optional)
set.seed(123)
box_data_expanded$PC1_pred <- box_data_expanded$PC1_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

# Predicted Plot
ggplot(box_data_expanded, aes(x = provenance, y = PC1_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("native" = "blue", "non_native" = "orange")) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Predicted PC1 by Provenance",
    x = "Provenance",
    y = "Predicted PC1"
  ) +
  theme_minimal()

#observed box plot
geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("native" = "blue", "non_native" = "orange")) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(
    title = "Observed PC1 by Provenance",
    x = "Provenance",
    y = "Observed PC1"
  ) +
  theme_minimal()



####PC2####

#lm
lm2<-lm(PC2 ~ provenance, data = d_pca)
summary(lm2)
#plot(lm2)

#phylm
phym2 = phylolm(PC2 ~ provenance, data = d_pca, phy)
summary(phym2)


# Compare Q-Q plots 
par(mfrow = c(1, 2))  
# Q-Q plot for lm model
qqnorm(residuals(lm2), main = "Q-Q Plot: lm")
qqline(residuals(lm2), col = "blue")

# Q-Q plot for phylolm model
qqnorm(residuals(phym2), main = "Q-Q Plot: phylolm")
qqline(residuals(phym2), col = "red")

# Reset plotting area
par(mfrow = c(1, 1))

#Check for phylogenetic signal (yes)

#phylosig(phy, trait_vector, method = "lambda", test = TRUE)



#compare models 
AIC(lm2,phym2)# lm us waaaaaay better(but only significant because of sample size)

# Extract residuals
predicted_values <- predict(lm2, newdata = d_pca)
print(predicted_values)
residuals <- residuals(lm2)

# Plot residuals
plot(predicted_values, residuals, main = "Residuals of phylolm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for phylolm Model")
qqline(residuals, col = "red")


#PC2 models with other predictors
phym3a = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = d_pca, phy) # no polynomial
summary(phym3a)

lm3a = phylolm(PC2 ~ ave_tmean+ provenance + ave_precip+I(ave_precip^2), data = d_pca, phy) # no polynomial # lm no polynomial
summary(lm3a)

#best model per aic below
phym3b = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca, phy) # polynomial and interaction
summary(phym3b)

phym3c = phylolm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = d_pca, phy) # polynomial and no interaction
summary(phym3c) 

#this is the top model via AIC below
phym3d = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = d_pca, phy) # interaction+ polynomial for only precip (based on plot)
summary(phym3d)

#check vif best model
vif(lm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = d_pca), type='predictor')

# Calculate AIC for each model
aic3a <- AIC(phym3a)
aic3b <- AIC(phym3b)#best mode per aic
aic3c <- AIC(phym3c)
aic3d <- AIC(phym3d)#lowes AIC
aic3e <- AIC(lm3a)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"),
  AIC = c(aic3a, aic3b, aic3c,aic3d,aic3e)
)


# Order by AIC
aic_values <- aic_values[order(aic_values$AIC), ]

# Calculate delta AIC
aic_values$Delta_AIC <- aic_values$AIC - min(aic_values$AIC)
aic_values

#check residuals with dharma
simulationOutput <- simulateResiduals(fittedModel = phym3b, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)


#plot residuals by hand
residuals <- residuals(lm3b)
predicted<-predict.lm(lm3b)

# Plot residuals
plot(residuals, main = "Residuals of phylolm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for phylolm Model")
qqline(residuals, col = "red")


##plot interactions

# Means for holding other variables constant
mean_tmean <- mean(d_pca$ave_tmean, na.rm = TRUE)
mean_precip <- mean(d_pca$ave_precip, na.rm = TRUE)

# Custom colors
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Effect of ave_tmean (with quadratic term)
tmean_seq <- seq(min(d_pca$ave_tmean), max(d_pca$ave_tmean), length.out = 100)
new_tmean <- expand.grid(
  ave_tmean = tmean_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_precip = mean_precip,
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = mean_precip^2
  )

new_tmean$PC2_pred <- predict(phym3b, newdata = new_tmean)

ggplot(new_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 vs Temperature",
    x = "Average Temperature",
    y = "Predicted PC2"
  ) +
  theme_minimal()

# 2. Effect of ave_precip (with quadratic term)
precip_seq <- seq(min(d_pca$ave_precip), max(d_pca$ave_precip), length.out = 100)
new_precip <- expand.grid(
  ave_precip = precip_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_tmean = mean_tmean,
    `I(ave_tmean^2)` = mean_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

new_precip$PC2_pred <- predict(phym3b, newdata = new_precip)

ggplot(new_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 vs Precipitation",
    x = "Average Precipitation",
    y = "Predicted PC2"
  ) +
  theme_minimal()


#3 box plot with only provenance varying


#  Create prediction grid
box_data <- expand.grid(
  ave_tmean = mean(d_pca$ave_tmean, na.rm = TRUE),
  ave_precip = mean(d_pca$ave_precip, na.rm = TRUE),
  provenance = c("native", "non_native")
) %>%
  mutate(
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

#  Predict PC1 using your phylolm model
box_data$PC2_pred <- predict(phym3b, newdata = box_data)

# (optional): Expand to simulate distribution for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC2_pred <- box_data_expanded$PC2_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

# Plot
ggplot(box_data_expanded, aes(x = provenance, y = PC2_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("native" = "blue", "non_native" = "orange")) +
  labs(
    title = "Predicted PC2 by Provenance",
    x = "Provenance",
    y = "Predicted PC2"
  ) +
  theme_minimal()
####PC1-herb####

herb<-d_pca%>%
      filter(growth_form=="herb") 

herb = as.data.frame(herb)
row.names(herb) = herb$tip



#lm
lm1<-lm(PC1 ~ provenance, data = herb)
summary(lm1)
#plot(lm2)

#phylm
phym1 = phylolm(PC1 ~ provenance, data = herb, phy)
summary(phym1)



# Compare Q-Q plots 
par(mfrow = c(1, 2))  
# Q-Q plot for lm model
qqnorm(residuals(lm1), main = "Q-Q Plot: lm")
qqline(residuals(lm1), col = "blue")

# Q-Q plot for phylolm model
qqnorm(residuals(phym1), main = "Q-Q Plot: phylolm")
qqline(residuals(phym1), col = "red")

# Reset plotting area
par(mfrow = c(1, 1))

#compare models 
AIC(lm1)#
AIC(phym1)


#env variables
plot(herb$ave_tmean, herb$PC1)
plot(herb$ave_precip, herb$PC1)
plot(herb$ave_tmean, herb$PC2)
plot(herb$ave_precip, herb$PC2)



#models with predictors 
phylm3a = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance, data = herb,phy)# no polynomial
summary(phylm3a)

#best fit model based on AIC below
phylm3b = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = herb,phy) # polynomial and interaction
summary(phylm3b)

phylm3c = phylolm(PC1 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = herb,phy) # polynomial and no interaction
summary(phylm3c) 


phylm3d = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data =herb,phy) # interaction+ polynomial for only precip (based on plot)
summary(phylm3d)



# Calculate AIC for each model
aic3a <- AIC(phylm3a)
aic3b <- AIC(phylm3b)#lowes AIC
aic3c <- AIC(phylm3c)
aic3d <- AIC(phylm3d)
#aic3e <- AIC(lm3a)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4"),
  AIC = c(aic3a, aic3b, aic3c,aic3d)
)


# Order by AIC
aic_values <- aic_values[order(aic_values$AIC), ]

# Calculate delta AIC
aic_values$Delta_AIC <- aic_values$AIC - min(aic_values$AIC)
aic_values


# Q-Q plot for phylolm model
qqnorm(residuals(phylm3b), main = "Q-Q Plot: phylolm")
qqline(residuals(phylm3b), col = "red")


#plot residuals by hand
residuals <- residuals(phylm3b)
predicted<-predict(phylm3b,newdata = herb)

# Plot residuals
plot(predicted,residuals, main = "Residuals of lm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for lm Model")
qqline(residuals, col = "red")


##plot interactions
# Set up means
mean_tmean <- mean(herb$ave_tmean, na.rm = TRUE)
mean_precip <- mean(herb$ave_precip, na.rm = TRUE)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean (with squared term)
tmean_seq <- seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)
pred_tmean <- expand.grid(
  ave_tmean = tmean_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_precip = mean_precip,
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = mean_precip^2
  )

pred_tmean$PC1_pred <- predict(phylm3b, newdata = pred_tmean)

ggplot(pred_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 vs Temperature - Herbs",
    x = "Average Temperature",
    y = "Predicted PC1"
  ) +
  theme_minimal()

# 2. Prediction over ave_precip (with squared term + interaction)
precip_seq <- seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100)
pred_precip <- expand.grid(
  ave_precip = precip_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_tmean = mean_tmean,
    `I(ave_tmean^2)` = mean_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

pred_precip$PC1_pred <- predict(phylm3b, newdata = pred_precip)

ggplot(pred_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 vs Precipitation- Herbs",
    x = "Average Precipitation",
    y = "Predicted PC1"
  ) +
  theme_minimal()



# 3. Boxplot of predicted PC1 by provenance
box_data <- expand.grid(
  ave_tmean = mean_tmean,
  ave_precip = mean_precip,
  provenance = c("native", "non_native")
) %>%
  mutate(
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

box_data$PC1_pred <- predict(phylm3b, newdata = box_data)

# Optional: simulate variability for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC1_pred <- box_data_expanded$PC1_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

ggplot(box_data_expanded, aes(x = provenance, y = PC1_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 by Provenance - Herbs",
    x = "Provenance",
    y = "Predicted PC1"
  ) +
  theme_minimal()






####PC2-herbs####

phymh2 = phylolm(PC2 ~ provenance, data = herb, phy)
summary(phymh2) # still significantly different after considering phylogenetic relationship

lmh2 = lm(PC2 ~ provenance, data = herb)
summary(lmh2) 

#better
AIC(phymh2)
AIC(lmh2)


#env variables
plot(herb$ave_tmean, herb$PC1)
plot(herb$ave_precip, herb$PC1)
plot(herb$ave_tmean, herb$PC2)
plot(herb$ave_precip, herb$PC2)



#models with predictors 
phylm3a = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = herb,phy)# no polynomial
summary(phylm3a)

#best fit model based on AIC below
phylm3b = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = herb,phy) # polynomial and interaction
summary(phylm3b)

phylm3c = phylolm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = herb,phy) # polynomial and no interaction
summary(phylm3c) 


phylm3d = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data =herb,phy) # interaction+ polynomial for only precip (based on plot)
summary(phylm3d)


# Calculate AIC for each model
aic3a <- AIC(phylm3a)
aic3b <- AIC(phylm3b)#lowes AIC
aic3c <- AIC(phylm3c)
aic3d <- AIC(phylm3d)
#aic3e <- AIC(lm3a)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4"),
  AIC = c(aic3a, aic3b, aic3c,aic3d)
)


# Order by AIC
aic_values <- aic_values[order(aic_values$AIC), ]

# Calculate delta AIC
aic_values$Delta_AIC <- aic_values$AIC - min(aic_values$AIC)
aic_values


# Q-Q plot for phylolm model
qqnorm(residuals(phylm3b), main = "Q-Q Plot: phylolm")
qqline(residuals(phylm3b), col = "red")


#plot residuals by hand
residuals <- residuals(phylm3b)
predicted<-predict(phylm3b,newdata = herb)

# Plot residuals
plot(predicted,residuals, main = "Residuals of lm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for lm Model")
qqline(residuals, col = "red")


##plot interactions

# Set up means
mean_tmean <- mean(herb$ave_tmean, na.rm = TRUE)
mean_precip <- mean(herb$ave_precip, na.rm = TRUE)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean (with squared term + interaction)
tmean_seq <- seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)
pred_tmean <- expand.grid(
  ave_tmean = tmean_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_precip = mean_precip,
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = mean_precip^2
  )

pred_tmean$PC2_pred <- predict(phylm3b, newdata = pred_tmean)

ggplot(pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 vs Temperature",
    x = "Average Temperature",
    y = "Predicted PC2"
  ) +
  theme_minimal()

# 2. Prediction over ave_precip (with squared term + interaction)
precip_seq <- seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100)
pred_precip <- expand.grid(
  ave_precip = precip_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_tmean = mean_tmean,
    `I(ave_tmean^2)` = mean_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

pred_precip$PC2_pred <- predict(phylm3b, newdata = pred_precip)

ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 vs Precipitation - Herbs",
    x = "Average Precipitation",
    y = "Predicted PC2"
  ) +
  theme_minimal()

# 3. Boxplot of predicted PC2 by provenance
box_data <- expand.grid(
  ave_tmean = mean_tmean,
  ave_precip = mean_precip,
  provenance = c("native", "non_native")
) %>%
  mutate(
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

box_data$PC2_pred <- predict(phylm3b, newdata = box_data)

# Optional: simulate variability for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC2_pred <- box_data_expanded$PC2_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

ggplot(box_data_expanded, aes(x = provenance, y = PC2_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 by Provenance - Herbs",
    x = "Provenance",
    y = "Predicted PC2"
  ) +
  theme_minimal()

####PC1-woody####


woody<-d_pca%>%
  filter(growth_form=="woody")

woody = as.data.frame(woody)
row.names(woody) = woody$tip
phywm1 = phylolm(PC1 ~ provenance, data = woody, phy)
summary(phywm1) # still significantly different after considering phylogenetic relationship

lwm1 = lm(PC1 ~ provenance, data = woody)
summary(lwm1)



plot(woody$ave_tmean, woody$PC1)
plot(woody$ave_precip, woody$PC1)
plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)


# with other variables
mha = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance, data = woody, phy) # no interaction
summary(mha)

mhb = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = woody, phy) # no interaction
summary(mhb)

mhd = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = woody, phy) # no interaction
summary(mhd)

mhc = phylolm(PC1 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = woody, phy) # no interaction
summary(mhc) # still significantly different after considering phylogenetic relationship and temp and precip


# Calculate AIC for each model
aicha <- AIC(mha)
aichb <- AIC(mhb)#lowes AIC
aichc <- AIC(mhc)
aichd <- AIC(mhd)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4"),
  AIC = c(aicha, aichb, aichc,aichd)
)

# Order by AIC
aic_values <- aic_values[order(aic_values$AIC), ]

# Calculate delta AIC
aic_values$Delta_AIC <- aic_values$AIC - min(aic_values$AIC)
aic_values


# Q-Q plot for phylolm model
qqnorm(residuals(mhb), main = "Q-Q Plot: phylolm")
qqline(residuals(mhb), col = "red")


#plot residuals by hand
residuals <- residuals(mhb)
predicted<-predict(mhb,newdata = woody)

# Plot residuals
plot(predicted,residuals, main = "Residuals of lm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for lm Model")
qqline(residuals, col = "red")


##plot interactions

# Set up means and colors
mean_tmean <- mean(woody$ave_tmean, na.rm = TRUE)
mean_precip <- mean(woody$ave_precip, na.rm = TRUE)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean
tmean_seq <- seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)
pred_tmean <- expand.grid(
  ave_tmean = tmean_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_precip = mean_precip,
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = mean_precip^2
  )

pred_tmean$PC1_pred <- predict(mhb, newdata = pred_tmean)

ggplot(pred_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 vs Temperature-Woody",
    x = "Average Temperature",
    y = "Predicted PC1"
  ) +
  ylim(-1, 1) +
  theme_minimal()

# 2. Prediction over ave_precip
precip_seq <- seq(min(woody$ave_precip, na.rm = TRUE), max(woody$ave_precip, na.rm = TRUE), length.out = 100)
pred_precip <- expand.grid(
  ave_precip = precip_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_tmean = mean_tmean,
    `I(ave_tmean^2)` = mean_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

pred_precip$PC1_pred <- predict(mhb, newdata = pred_precip)

ggplot(pred_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 vs Precipitation -Woody",
    x = "Average Precipitation",
    y = "Predicted PC1"
  ) +
  theme_minimal()

# 3. Boxplot of predicted PC1 by provenance


# Create prediction data
box_data <- expand.grid(
  ave_tmean = mean_tmean,
  ave_precip = mean_precip,
  provenance = c("native", "non_native")
) %>%
  mutate(
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

# Predict PC1
box_data$PC1_pred <- predict(mhb, newdata = box_data)

# Simulate variability for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC1_pred <- box_data_expanded$PC1_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

# Plot
ggplot(box_data_expanded, aes(x = provenance, y = PC1_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Predicted PC1 by Provenance - Woody",
    x = "Provenance",
    y = "Predicted PC1"
  ) +
  ylim(-.5, .5) +  # Set y-axis range here
  theme_minimal()




####PC2-Woody####
#
phywm1 = phylolm(PC2 ~ provenance, data = woody, phy)
summary(phywm1) # still significantly different after considering phylogenetic relationship

lwm1 = lm(PC2 ~ provenance, data = woody)
summary(lwm1)



plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)
plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)


#PC2 models with other variables
mha = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = woody, phy) # no interaction
summary(mha)

#best model
mhb = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = woody, phy) # no interaction
summary(mhb)

vif(lm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = woody),type = 'predictor') # no interaction
summary(mhb)


mhd = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = woody, phy) # no interaction
summary(mhd)

mhc = phylolm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = woody, phy) # no interaction
summary(mhc) 


# Calculate AIC for each model
aicha <- AIC(mha)
aichb <- AIC(mhb)#lowes AIC
aichc <- AIC(mhc)
aichd <- AIC(mhd)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4"),
  AIC = c(aicha, aichb, aichc,aichd)
)

# Order by AIC
aic_values <- aic_values[order(aic_values$AIC), ]

# Calculate delta AIC
aic_values$Delta_AIC <- aic_values$AIC - min(aic_values$AIC)
aic_values


# Q-Q plot for phylolm model
qqnorm(residuals(mhb), main = "Q-Q Plot: phylolm")
qqline(residuals(mhb), col = "red")


#plot residuals by hand
residuals <- residuals(mhb)
predicted<-predict(mhb,newdata = woody)

# Plot residuals
plot(predicted,residuals, main = "Residuals of lm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for lm Model")
qqline(residuals, col = "red")


##plot interactions

# Set up means and colors
mean_tmean <- mean(woody$ave_tmean, na.rm = TRUE)
mean_precip <- mean(woody$ave_precip, na.rm = TRUE)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean
tmean_seq <- seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)
pred_tmean <- expand.grid(
  ave_tmean = tmean_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_precip = mean_precip,
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = mean_precip^2
  )

pred_tmean$PC2_pred <- predict(mhb, newdata = pred_tmean)

ggplot(pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 vs Temperature-Woody",
    x = "Average Temperature",
    y = "Predicted PC2"
  ) +
  ylim(-1, 1) +
  theme_minimal()

# 2. Prediction over ave_precip
precip_seq <- seq(min(woody$ave_precip, na.rm = TRUE), max(woody$ave_precip, na.rm = TRUE), length.out = 100)
pred_precip <- expand.grid(
  ave_precip = precip_seq,
  provenance = c("native", "non_native")
) %>%
  mutate(
    ave_tmean = mean_tmean,
    `I(ave_tmean^2)` = mean_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

pred_precip$PC2_pred <- predict(mhb, newdata = pred_precip)

ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 vs Precipitation -Woody",
    x = "Average Precipitation",
    y = "Predicted PC2"
  ) +
  theme_minimal()

# 3. Boxplot of predicted PC2 by provenance


# Create prediction data
box_data <- expand.grid(
  ave_tmean = mean_tmean,
  ave_precip = mean_precip,
  provenance = c("native", "non_native")
) %>%
  mutate(
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

# Predict PC2
box_data$PC2_pred <- predict(mhb, newdata = box_data)

# Simulate variability for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC2_pred <- box_data_expanded$PC2_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

# Plot
ggplot(box_data_expanded, aes(x = provenance, y = PC2_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 by Provenance - Woody",
    x = "Provenance",
    y = "Predicted PC2"
  ) +
  #ylim(-.5, .5) +  # Set y-axis range here
  theme_minimal()


#add raw data to plots

ggplot() +
  # Raw data points
  geom_point(data = woody, aes(x = ave_tmean, y = PC2, color = provenance),
             alpha = 0.4, size = 1.5) +
  
  # Predicted lines
  geom_line(data = pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance),
            size = 1.2) +
  
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Predicted PC2 vs Temperature - Woody",
    x = "Average Temperature",
    y = "PC2"
  ) +
  ylim(-4, 4) +
  theme_minimal()


















      
