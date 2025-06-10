library(tidyverse)
library(ape)
install.packages("DHARMa")
library(DHARMa)
library(phylolm)

#data
d_pca = read_csv("Data/d_pca.csv")
d_pca
phy = read.tree("Data/phy.tre")
plot(phy, type = "fan", show.tip.label = F)



#All Taxa

#make data frame
d_pca = as.data.frame(d_pca)
row.names(d_pca) = d_pca$tip

#do we need to include Phylogeny? yes
lm1<-lm(PC1 ~ provenance, data = d_pca)
summary(lm1)
plot(lm1)

phym1 = phylolm(PC1 ~ provenance, data = d_pca, phy)
summary(phym1)

lm1<-lm(PC1 ~ provenance, data = d_pca)

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



phym2 = phylolm(PC2 ~ provenance, data = d_pca, phy)
summary(phym2) # still significantly different after considering phylogenetic relationship

lm2 = lm(PC2 ~ provenance, data = d_pca)
summary(lm2)

AIC(phym2)
AIC(lm2)

plot(d_pca$ave_tmean, d_pca$PC1)
plot(d_pca$ave_precip, d_pca$PC1)
plot(d_pca$ave_tmean, d_pca$PC2)
plot(d_pca$ave_precip, d_pca$PC2)

# scale variables
d_pca$ave_tmean = scale(d_pca$ave_tmean)[,1]
d_pca$ave_precip = scale(d_pca$ave_precip)[,1]


#PC1
m3a = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance, data = d_pca, phy) # no interaction
summary(m3)

phym3b = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca, phy) # no interaction
summary(phym3b)
plot(m3b)

lm3b = lm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca) # no interaction
summary(lm3b)
plot(m3b)

AIC(phym3b)
AIC(lm3b)
simulationOutput <- simulateResiduals(fittedModel = phym3b, plot = F)

plot(simulationOutput)

testDispersion(simulationOutput)

residuals <- residuals(phym3b)


predicted<-predict.lm(phym3b)

# Plot residuals
plot(residuals, main = "Residuals of phylolm Model", xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")



qqnorm(residuals, main = "QQ Plot of Residuals for phylolm Model")
qqline(residuals, col = "red")


m3d = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = d_pca, phy) # no interaction
summary(m3d)

m3c = phylolm(PC1 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = d_pca, phy) # no interaction
summary(m3c) # still significantly different after considering phylogenetic relationship and temp and precip


# Calculate AIC for each model
aic3a <- AIC(m3a)
aic3b <- AIC(m3b)#lowes AIC
aic3c <- AIC(m3c)
aic3d <- AIC(m3d)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4"),
  AIC = c(aic3a, aic3b, aic3c,aic3d)
)


#PC2
m4a = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = d_pca, phy) # no interaction
summary(m4a)

phym4b = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca, phy) # no interaction
summary(m4b)

lm4b = lm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca) # no interaction
summary(m4b)

#lm is better for second
AIC(phym4b)
AIC(lm4b)


simulationOutput <- simulateResiduals(fittedModel = m4b, plot = F)

plot(simulationOutput)

testDispersion(simulationOutput)

m4c = phylolm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = d_pca, phy) # no interaction
summary(m4c) # still significantly different after considering phylogenetic relationship and temp and precip


# Calculate AIC for each model
aic4a <- AIC(m4a)
aic4b <- AIC(m4b)#lowes AIC
aic4c <- AIC(m4c)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3"),
  AIC = c(aic4a, aic4b, aic4c)
)


####herb####

herb<-d_pca%>%
      filter(growth_form=="herb") 

herb = as.data.frame(herb)
row.names(herb) = herb$tip

phymh1 = phylolm(PC1 ~ provenance, data = herb, phy)
summary(phymh1) # still significantly different after considering phylogenetic relationship
lmh1 = lm(PC1 ~ provenance, data = herb)
summary(lmh1) 

#lm better
AIC(phymh1)

AIC(lmh1)

#pc2

phymh2 = phylolm(PC2 ~ provenance, data = herb, phy)
summary(phymh2) # still significantly different after considering phylogenetic relationship

lmh2 = lm(PC2 ~ provenance, data = herb)
summary(lmh2) 

#lm better
AIC(phymh2)
AIC(lmh2)


#env variables
plot(herb$ave_tmean, herb$PC1)
plot(herb$ave_precip, herb$PC1)
plot(herb$ave_tmean, herb$PC2)
plot(herb$ave_precip, herb$PC2)


#PC1
mha = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance, data = herb, phy) # no interaction
summary(mha)

phymhb2 = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = herb, phy) # no interaction
summary(phymhb2)

lmhb2 = lm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = herb) # no interaction
summary(lmhb2)

#lm is better
AIC(phymhb2)
AIC(lmhb2)

mhd = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_precip^2) , data = herb, phy) # no interaction
summary(mhd)

mhc = phylolm(PC1 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = herb, phy) # no interaction
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


#PC2
m4a = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = d_pca, phy) # no interaction
summary(m4a)

m4b = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca, phy) # no interaction
summary(m4b)

m4c = phylolm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = d_pca, phy) # no interaction
summary(m4c) # still significantly different after considering phylogenetic relationship and temp and precip


# Calculate AIC for each model
aic4a <- AIC(m4a)
aic4b <- AIC(m4b)#lowes AIC
aic4c <- AIC(m4c)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3"),
  AIC = c(aic4a, aic4b, aic4c)
)




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

AIC(phywm1)
AIC(lwm1)

#pc2
m2 = phylolm(PC2 ~ provenance, data = woody, phy)
summary(m2) # still significantly different after considering phylogenetic relationship

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


#PC2
m4a = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance, data = d_pca, phy) # no interaction
summary(m4a)

m4b = phylolm(PC2 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca, phy) # no interaction
summary(m4b)

m4c = phylolm(PC2 ~ ave_tmean + ave_precip + provenance+I(ave_tmean^2)+I(ave_precip^2), data = d_pca, phy) # no interaction
summary(m4c) # still significantly different after considering phylogenetic relationship and temp and precip


# Calculate AIC for each model
aic4a <- AIC(m4a)
aic4b <- AIC(m4b)#lowes AIC
aic4c <- AIC(m4c)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3"),
  AIC = c(aic4a, aic4b, aic4c)
)


      
