library(tidyverse)
library(ape)
d_pca = read_csv("Data/d_pca.csv")
d_pca
phy = read.tree("Data/phy.tre")
plot(phy, type = "fan", show.tip.label = F)



#All Taxa
library(phylolm)
d_pca = as.data.frame(d_pca)
row.names(d_pca) = d_pca$tip

m1 = phylolm(PC1 ~ provenance, data = d_pca, phy)
summary(m1) # still significantly different after considering phylogenetic relationship

m2 = phylolm(PC2 ~ provenance, data = d_pca, phy)
summary(m2) # still significantly different after considering phylogenetic relationship

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

m3b = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = d_pca, phy) # no interaction
summary(m3b)

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


####herb####

herb<-d_pca%>%
      filter(growth_form=="herb") 

herb = as.data.frame(herb)
row.names(herb) = herb$tip

m1 = phylolm(PC1 ~ provenance, data = herb, phy)
summary(m1) # still significantly different after considering phylogenetic relationship

m2 = phylolm(PC2 ~ provenance, data = herb, phy)
summary(m2) # still significantly different after considering phylogenetic relationship

plot(herb$ave_tmean, herb$PC1)
plot(herb$ave_precip, herb$PC1)
plot(herb$ave_tmean, herb$PC2)
plot(herb$ave_precip, herb$PC2)


#PC1
mha = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance, data = herb, phy) # no interaction
summary(mha)

mhb = phylolm(PC1 ~ ave_tmean*provenance + ave_precip*provenance+I(ave_tmean^2)+I(ave_precip^2) , data = herb, phy) # no interaction
summary(mhb)

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

m1 = phylolm(PC1 ~ provenance, data = woody, phy)
summary(m1) # still significantly different after considering phylogenetic relationship

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


      
