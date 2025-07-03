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



####All Taxa####

#make data frame
d_pca = as.data.frame(d_pca)
row.names(d_pca) = d_pca$tip

#PC1

#lm
lm1<-lm(PC1 ~ provenance, data = d_pca)
summary(lm1)
#plot(lm1)

#phylm
phym1 = phylolm(PC1 ~ provenance, data = d_pca, phy)
summary(phym1)

phym1 = phylolm(yj$x.t ~ provenance, data = d_pca, phy)
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

#Check for phylogenetic signal (yes)

#order and match tree tips and rownames
trait_vector <- d_pca$PC2
names(trait_vector) <- rownames(d_pca)
trait_vector <- trait_vector[phy$tip.label]
head(trait_vector)

phylosig(phy, trait_vector, method = "lambda", test = TRUE)

#compare models, yes we need the phylolm
AIC(lm1)
AIC(phym1)


# Extract residuals
predicted_values <- predict(phym1, newdata = d_pca)
print(predicted_values)
residuals <- residuals(phym1)

# Plot residuals
plot(predicted_values, residuals, main = "Residuals of phylolm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for phylolm Model")
qqline(residuals, col = "red")


#models with other predictors

plot(d_pca$ave_tmean, d_pca$PC1)
plot(d_pca$ave_precip, d_pca$PC1)
plot(d_pca$ave_tmean, d_pca$PC2)
plot(d_pca$ave_precip, d_pca$PC2)

# scale variables
d_pca$ave_tmean = scale(d_pca$ave_tmean)[,1]
d_pca$ave_precip = scale(d_pca$ave_precip)[,1]


#PC1 models with other predictors
phym3a = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance, data = d_pca, phy) # no polynomial
summary(phym3a)

lm3a = lm(PC1 ~ ave_tmean*provenance + ave_precip*provenance , data = d_pca) # lm no polynomial
summary(lm3a)

phym3b = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca, phy) # polynomial and interaction
summary(phym3b)

phym3c = phylolm(PC1 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = d_pca, phy) # polynomial and no interaction
summary(phym3c) 

#this is the top model via AIC below
phym3d = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = d_pca, phy) # interaction+ polynomial for only precip (based on plot)
summary(phym3d)

#check vif best model
vif(lm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = d_pca), type='predictor')

# Calculate AIC for each model
aic3a <- AIC(phym3a)
aic3b <- AIC(phym3b)
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
simulationOutput <- simulateResiduals(fittedModel = phym3d, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)


#plot residuals by hand
residuals <- residuals(phym3b)
predicted<-predict.lm(phym3b)#not working because phylolm

# Plot residuals
plot(residuals, main = "Residuals of phylolm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals, main = "QQ Plot of Residuals for phylolm Model")
qqline(residuals, col = "red")



#PC2

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

phylosig(phy, trait_vector, method = "lambda", test = TRUE)



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


#models with predictors (only using LM here becasue it is such a better fit re: aic)
lm3a = lm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = d_pca)# no polynomial
summary(lm3a)

#best fit model based on AIC below
lm3b = lm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca) # polynomial and interaction
summary(lm3b)

lm3c = lm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = d_pca) # polynomial and no interaction
summary(lm3c) 

#this is the top model via AIC below
lm3d = lm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = d_pca) # interaction+ polynomial for only precip (based on plot)
summary(lm3d)



# Calculate AIC for each model
aic3a <- AIC(lm3a)
aic3b <- AIC(lm3b)
aic3c <- AIC(lm3c)
aic3d <- AIC(lm3d)#lowes AIC
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

#check residuals with dharma
simulationOutput <- simulateResiduals(fittedModel = phym3d, plot = F)
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

# Create prediction grid for both interactions
grid_tmean <- expand.grid(
  ave_tmean = seq(min(d_pca$ave_tmean), max(d_pca$ave_tmean), length.out = 100),
  ave_precip = mean(d_pca$ave_precip),
  provenance = levels(factor(d_pca$provenance))
)




grid_tmean$interaction <- "ave_tmean × provenance"
grid_tmean$I.ave_tmean.2 <- grid_tmean$ave_tmean^2
grid_tmean$I.ave_precip.2 <- grid_tmean$ave_precip^2

grid_precip <- expand.grid(
  ave_tmean = mean(d_pca$ave_tmean),
  ave_precip = seq(min(d_pca$ave_precip), max(d_pca$ave_precip), length.out = 100),
  provenance = levels(factor(d_pca$provenance))
)

grid_precip$interaction <- "ave_precip × provenance"
grid_precip$I.ave_tmean.2 <- grid_precip$ave_tmean^2
grid_precip$I.ave_precip.2 <- grid_precip$ave_precip^2

# Combine both grids
new_data <- bind_rows(grid_tmean, grid_precip)

# Predict from the model
model <- lm(PC2 ~ ave_tmean * provenance + ave_precip * provenance + 
              I(ave_tmean^2) + I(ave_precip^2), data = d_pca)

new_data$PC2_pred <- predict(model, newdata = new_data)

# Plot both interactions
ggplot(new_data, aes(x = ifelse(interaction == "ave_tmean × provenance", ave_tmean, ave_precip),
                     y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  facet_wrap(~interaction, scales = "free_x") +
  labs(x = "Predictor", y = "Predicted PC2", color = "Provenance",
       title = "Interaction Effects on PC2") +
  theme_minimal()




####herb####

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



#models with predictors (only using LM here becasue it is such a better fit re: aic)
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

# Create prediction grid for both interactions
grid_tmean <- expand.grid(
  ave_tmean = seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100),
  ave_precip = mean(herb$ave_precip),
  provenance = levels(factor(herb$provenance))
)




grid_tmean$interaction <- "ave_tmean × provenance"
grid_tmean$I.ave_tmean.2 <- grid_tmean$ave_tmean^2
grid_tmean$I.ave_precip.2 <- grid_tmean$ave_precip^2

grid_precip <- expand.grid(
  ave_tmean = mean(herb$ave_tmean),
  ave_precip = seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100),
  provenance = levels(factor(herb$provenance))
)

grid_precip$interaction <- "ave_precip × provenance"
grid_precip$I.ave_tmean.2 <- grid_precip$ave_tmean^2
grid_precip$I.ave_precip.2 <- grid_precip$ave_precip^2

# Combine both grids
new_data <- bind_rows(grid_tmean, grid_precip)

# Predict from the model
model <- phylm3b

new_data$PC1_pred <- predict(model, newdata = new_data)

# Plot both interactions
ggplot(new_data, aes(x = ifelse(interaction == "ave_tmean × provenance", ave_tmean, ave_precip),
                     y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  facet_wrap(~interaction, scales = "free_x") +
  labs(x = "Predictor", y = "Predicted PC1", color = "Provenance",
       title = "Interaction Effects on PC1") +
  coord_cartesian(ylim = c(-4, 2)) + 
  theme_minimal()





#pc2

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

# Create prediction grid for both interactions
grid_tmean <- expand.grid(
  ave_tmean = seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100),
  ave_precip = mean(herb$ave_precip),
  provenance = levels(factor(herb$provenance))
)




grid_tmean$interaction <- "ave_tmean × provenance"
grid_tmean$I.ave_tmean.2 <- grid_tmean$ave_tmean^2
grid_tmean$I.ave_precip.2 <- grid_tmean$ave_precip^2

grid_precip <- expand.grid(
  ave_tmean = mean(herb$ave_tmean),
  ave_precip = seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100),
  provenance = levels(factor(herb$provenance))
)

grid_precip$interaction <- "ave_precip × provenance"
grid_precip$I.ave_tmean.2 <- grid_precip$ave_tmean^2
grid_precip$I.ave_precip.2 <- grid_precip$ave_precip^2

# Combine both grids
new_data <- bind_rows(grid_tmean, grid_precip)

# Predict from the model
model <- phylm3b

new_data$PC2_pred <- predict(model, newdata = new_data)

# Plot both interactions
ggplot(new_data, aes(x = ifelse(interaction == "ave_tmean × provenance", ave_tmean, ave_precip),
                     y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  facet_wrap(~interaction, scales = "free_x") +
  labs(x = "Predictor", y = "Predicted PC2", color = "Provenance",
       title = "Interaction Effects on PC2") +
  coord_cartesian(ylim = c(-3, 3)) + 
  theme_minimal()

phylm3b

####woody####

woody<-d_pca%>%
       filter(growth_form=="woody")

woody = as.data.frame(woody)
row.names(woody) = woody$tip

#pc1
phywm1 = phylolm(PC1 ~ provenance, data = woody, phy)
summary(phywm1) # still significantly different after considering phylogenetic relationship

lwm1 = lm(PC1 ~ provenance, data = woody)
summary(lwm1)



plot(woody$ave_tmean, woody$PC1)
plot(woody$ave_precip, woody$PC1)
plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)


#PC1
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

# Create prediction grid for both interactions
grid_tmean <- expand.grid(
  ave_tmean = seq(min(woody$ave_tmean), max(woody$ave_tmean), length.out = 100),
  ave_precip = mean(woody$ave_precip),
  provenance = levels(factor(woody$provenance))
)




grid_tmean$interaction <- "ave_tmean × provenance"
grid_tmean$I.ave_tmean.2 <- grid_tmean$ave_tmean^2
grid_tmean$I.ave_precip.2 <- grid_tmean$ave_precip^2

grid_precip <- expand.grid(
  ave_tmean = mean(woody$ave_tmean),
  ave_precip = seq(min(woody$ave_precip), max(woody$ave_precip), length.out = 100),
  provenance = levels(factor(woody$provenance))
)

grid_precip$interaction <- "ave_precip × provenance"
grid_precip$I.ave_tmean.2 <- grid_precip$ave_tmean^2
grid_precip$I.ave_precip.2 <- grid_precip$ave_precip^2

# Combine both grids
new_data <- bind_rows(grid_tmean, grid_precip)

# Predict from the model
model <- mhb

new_data$PC1_pred <- predict(model, newdata = new_data)

# Plot both interactions
ggplot(new_data, aes(x = ifelse(interaction == "ave_tmean × provenance", ave_tmean, ave_precip),
                     y = PC1_pred, color = provenance)) +
  geom_line(size = 1.2) +
  facet_wrap(~interaction, scales = "free_x") +
  labs(x = "Predictor", y = "Predicted PC1", color = "Provenance",
       title = "Interaction Effects on PC1") +
  coord_cartesian(ylim = c(-4, 2)) + 
  theme_minimal()







#pc2
#
phywm1 = phylolm(PC2 ~ provenance, data = woody, phy)
summary(phywm1) # still significantly different after considering phylogenetic relationship

lwm1 = lm(PC2 ~ provenance, data = woody)
summary(lwm1)



plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)
plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)


#PC2
mha = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = woody, phy) # no interaction
summary(mha)

mhb = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = woody, phy) # no interaction
summary(mhb)

vif(lm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = woody),type = 'predictor') # no interaction
summary(mhb)


mhd = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = woody, phy) # no interaction
summary(mhd)

mhc = phylolm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = woody, phy) # no interaction
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

# Create prediction grid for both interactions
grid_tmean <- expand.grid(
  ave_tmean = seq(min(woody$ave_tmean), max(woody$ave_tmean), length.out = 100),
  ave_precip = mean(woody$ave_precip),
  provenance = levels(factor(woody$provenance))
)




grid_tmean$interaction <- "ave_tmean × provenance"
grid_tmean$I.ave_tmean.2 <- grid_tmean$ave_tmean^2
grid_tmean$I.ave_precip.2 <- grid_tmean$ave_precip^2

grid_precip <- expand.grid(
  ave_tmean = mean(woody$ave_tmean),
  ave_precip = seq(min(woody$ave_precip), max(woody$ave_precip), length.out = 100),
  provenance = levels(factor(woody$provenance))
)

grid_precip$interaction <- "ave_precip × provenance"
grid_precip$I.ave_tmean.2 <- grid_precip$ave_tmean^2
grid_precip$I.ave_precip.2 <- grid_precip$ave_precip^2

# Combine both grids
new_data <- bind_rows(grid_tmean, grid_precip)

# Predict from the model
model <- mhb

new_data$PC2_pred <- predict(model, newdata = new_data)

# Plot both interactions
ggplot(new_data, aes(x = ifelse(interaction == "ave_tmean × provenance", ave_tmean, ave_precip),
                     y = PC2_pred, color = provenance)) +
  geom_line(size = 1.2) +
  facet_wrap(~interaction, scales = "free_x") +
  labs(x = "Predictor", y = "Predicted PC2", color = "Provenance",
       title = "Interaction Effects on PC2") +
  coord_cartesian(ylim = c(-4, 2)) + 
  theme_minimal()



















      
