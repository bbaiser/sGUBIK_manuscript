#Phylogenetic Linear Models and Plots

#install packages
library(phylolm.hp)
library(tidyverse)
library(ape)
library(DHARMa)
library(phylolm)
library(car)
library(phytools)
library(rr2)
library(patchwork)

####data####
#pca scores for native and non-native urban plant species including climate data
d_pca <- read_csv("Data/d_pca.csv") 

#phylogeny for species in the PCA data set
phy <- read.tree("Data/phy.tre") 

#Plot Phylogeny
plot(phy, type = "fan", show.tip.label = F)


####PC1-herbs####


#data

#filter for herbs (1174 species)
herb<-d_pca%>%
      filter(growth_form=="herb") 
      herb <- as.data.frame(herb)
      row.names(herb) = herb$tip


#scale climate variables
herb$ave_tmean<-scale(herb$ave_tmean)[,1]
herb$ave_precip<-scale(herb$ave_precip)[,1]

#Phylolm model
model_pc1_Herb <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb, phy = phy, model="lambda")

summary(model_pc1_Herb)

#VIF with lm
vif(lm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb),type= 'predictor')


# Q-Q plot for phylolm model
qqnorm(residuals(model_pc1_Herb), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc1_Herb), col = "red")


##plots
# Calculate means for plots
mean_tmean <- round(mean(herb$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(herb$ave_precip, na.rm = TRUE),2)


#Prediction over temperature

#make prediction grid
tmean_seq <-seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                    `I(ave_tmean^2)` = ave_tmean^2,
                    `I(ave_precip^2)` = mean_precip^2)

#predict onto grid
pred_tmean$PC1_pred <- predict(model_pc1_Herb, newdata = pred_tmean)

#plot
H_PC1_temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
            geom_line(size = 1.2) +
            geom_rug(data = herb, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
            scale_color_manual(values =c("native" = "blue", "non_native" = "orange")) +
            labs(title = "Predicted PC1 vs Temperature - Herbs",
                 x = " Temperature",
                 y = "PC1") +
            coord_cartesian(ylim = c(-3, 3)) +
            theme_minimal()

#Prediction over precipitation 

#make prediction grid
precip_seq <-seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)

#predict onto grid
pred_precip$PC1_pred <- predict(model_pc1_Herb, newdata = pred_precip)



H_PC1_precip<-ggplot(pred_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
              geom_line(size = 1.2) +
              geom_rug(data = herb, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
              scale_color_manual(values =c("native" = "blue", "non_native" = "orange")) +
              labs(title = "Predicted PC1 vs Precipitation- Herbs",
                   x = "Average Precipitation",
                   y = "PC1") +
              coord_cartesian(ylim = c(-3, 3)) +
              theme_minimal()




####PC2-herbs####


#Phylolm model
model_pc2_Herb <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb, phy = phy, model="lambda")

summary(model_pc2_Herb)


#VIF with lm
vif(lm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb),type= 'predictor')


# Q-Q plot for phylolm model
qqnorm(residuals(model_pc2_Herb ), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc2_Herb ), col = "red")


##plots
# Calculate means
mean_tmean <- round(mean(herb$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(herb$ave_precip, na.rm = TRUE),2)


# Prediction over temperature

#make prediction grid
tmean_seq <-seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                    `I(ave_tmean^2)` = ave_tmean^2,
                    `I(ave_precip^2)` = mean_precip^2)

#predict onto grid
pred_tmean$PC2_pred <- predict(model_pc2_Herb, newdata = pred_tmean)


#plot
H_PC2_temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
            geom_line(size = 1.2) +
            geom_rug(data = herb, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
            scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
            labs(title = "Predicted PC2 vs Temperature - Herbs",
                 x = "Average Temperature",
                 y = "PC2") +
            coord_cartesian(ylim = c(-3, 3)) +
            theme_minimal()

#Prediction over precipitation 

#make prediction grid
precip_seq <-seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)
#predict onto grid
pred_precip$PC2_pred <- predict(model_pc2_Herb, newdata = pred_precip)

H_PC2_precip<-ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
              geom_line(size = 1.2) +
              geom_rug(data = herb, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
              scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
              labs(title = "Predicted PC2 vs Precipitation- Herbs",
                   x = "Average Precipitation",
                   y = "PC2") +
              coord_cartesian(ylim = c(-3, 3)) +
              theme_minimal()



####PC1-woody####

#data
#filter for woody species
woody<-d_pca%>%
       filter(growth_form=="woody")
woody <-as.data.frame(woody)
row.names(woody) = woody$tip


#scale climate variables
woody$ave_tmean<-scale(woody$ave_tmean)[,1]
woody$ave_precip<-scale(woody$ave_precip)[,1]



#Phlyolm model
model_pc1_woody <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody, phy = phy, model="lambda")

summary(model_pc1_woody)

#VIF with lm
vif(lm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody),type= 'predictor')

# Q-Q plot for phylolm model
qqnorm(residuals(model_pc1_woody), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc1_woody), col = "red")


##Plots

# Calculate means
mean_tmean <- round(mean(woody$ave_tmean, na.rm = TRUE),2)
mean_precip <-round(mean(woody$ave_precip, na.rm = TRUE),2)


#Prediction over temeprature

#make prediction grid
tmean_seq <-seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                     `I(ave_tmean^2)` = ave_tmean^2,
                     `I(ave_precip^2)` = mean_precip^2)

#predict over grid
pred_tmean$PC1_pred <- predict(model_pc1_woody, newdata = pred_tmean)

#Plot
W_Pc1_Temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
            geom_line(size = 1.2) +
            geom_rug(data = woody, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
            scale_color_manual(values =  c("native" = "blue", "non_native" = "orange")) +
            labs(title = "Predicted PC1 vs Temperature-Woody",
                 x = "Temperature",
                 y = "PC1") +
            coord_cartesian(ylim=c(-3,3))+
            theme_minimal()



#Prediction over Precipitation
#make prediction grid
precip_seq <-seq(min(woody$ave_precip), max(woody$ave_precip), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)

#Predict over grid
pred_precip$PC1_pred <- predict(model_pc1_woody, newdata = pred_precip)

W_Pc1_Precip<-ggplot(pred_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
              geom_line(size = 1.2) +
              geom_rug(data = woody, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
              scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
              labs(title = "Predicted PC1 vs Precipitation -Woody",
                   x = "Average Precipitation",
                   y = "Predicted PC1") +
              coord_cartesian(ylim=c(-3,3))+
              theme_minimal()






####PC2-woody####

#Phylolm model
model_pc2_woody <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody, phy = phy, model="lambda")

summary(model_pc2_woody)

#VIF with lm
vif(lm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody),type= 'predictor')

# Q-Q plot for phylolm model
qqnorm(residuals(model_pc2_woody), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc2_woody), col = "red")



##Plots

# Calculate means 
mean_tmean <- mean(woody$ave_tmean, na.rm = TRUE)
mean_precip <- mean(woody$ave_precip, na.rm = TRUE)

#Prediction over temperature

#make prediction grid
tmean_seq <-seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                     `I(ave_tmean^2)` = ave_tmean^2,
                     `I(ave_precip^2)` = mean_precip^2)

#predict over grid
pred_tmean$PC2_pred <- predict(model_pc2_woody, newdata = pred_tmean)

W_Pc2_Temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
            geom_line(size = 1.2) +
            geom_rug(data = woody, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
            scale_color_manual(values =c("native" = "blue", "non_native" = "orange")) +
            labs(title = "Predicted PC2 vs Temperature-Woody",
                 x = "Temperature",
                 y = "PC2") +
            ylim(-3, 3) +
            theme_minimal()

#Prediction over precipitation

#make prediction grid
precip_seq <-seq(min(woody$ave_precip, na.rm = TRUE), max(woody$ave_precip, na.rm = TRUE), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)

#predict over grid
pred_precip$PC2_pred <- predict(model_pc2_woody, newdata = pred_precip)

#Plot
W_Pc2_Precip<-ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
              geom_line(size = 1.2) +
              geom_rug(data = woody, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
              scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
              labs(title = "Predicted PC2 vs Precipitation -Woody",
                   x = " Precipitation",
                   y = "PC2") +
              ylim(-3, 3) +
              theme_minimal()




####Make panel figure####


# Remove titles or legends from individual plots:
p1 <- H_PC1_temp + 
      labs(y = "PC 1 (Plant Structure)")+
      theme(plot.title = element_blank(), 
            legend.position = "none", 
            axis.title.x = element_blank())



p2 <- H_PC1_precip + 
      labs(y = "PC 1 (Plant Stucture)")+
       theme(plot.title = element_blank(), 
            axis.title.x = element_blank(),
            legend.position = c(0.99, 0.99),   # inside top right
            legend.justification = c("right", "top"),
            legend.background = element_rect(fill = alpha('white', 0.7), color = NA))

p3 <- H_PC2_temp + 
      labs(y = "PC 2 (Leaf Economics)") +
      theme(plot.title = element_blank(), 
           axis.title.x = element_blank(),
           legend.position = "none")



p4 <- H_PC2_precip + 
      labs(y = "PC 2 (Leaf Economics)") +
      theme(plot.title = element_blank(), 
            axis.title.x = element_blank(),
            legend.position = c(0.99, 0.99),   # inside top right
            legend.justification = c("right", "top"),
            legend.background = element_rect(fill = alpha('white', 0.7), color = NA))

p5 <- W_Pc1_Temp +   
      labs(y = "PC 1 (Plant Structure)")+
      theme(plot.title = element_blank(), 
            legend.position = "none")

p6 <- W_Pc1_Precip + 
      labs(y = "PC 1 (Plant Structure)")+
      theme(plot.title = element_blank(), 
            legend.position = "none")

p7 <- W_Pc2_Temp +  
      labs(y = "PC 2 (Leaf Economics)")+
      theme(plot.title = element_blank(), 
            legend.position = "none")

p8 <- W_Pc2_Precip +  
      labs(y = "PC 2 (Leaf Economics)")+
      theme(plot.title = element_blank(), 
            legend.position = "none")


#Figure 2 in the manuscript
# Add manual titles to each plot
p1_labeled <- p1 + labs(title = "a) Herbaceous") +
  theme(plot.title = element_text(hjust = 0))

p2_labeled <- p2 + labs(title = "b) Herbaceous") +
  theme(plot.title = element_text(hjust = 0))

p5_labeled <- p5 + labs(title = "c) Woody") +
  theme(plot.title = element_text(hjust = 0))

p6_labeled <- p6 + labs(title = "d) Woody") +
  theme(plot.title = element_text(hjust = 0))

# Combine manually labeled plots WITHOUT automatic tags
combined <- (p1_labeled + p2_labeled) / (p5_labeled + p6_labeled)



#Figure 3 in the manuscript
# Add manual titles to each plot
p3_labeled <- p3 + labs(title = "a) Herbaceous") +
  theme(plot.title = element_text(hjust = 0))

p4_labeled <- p4 + labs(title = "b) Herbaceous") +
  theme(plot.title = element_text(hjust = 0))

p7_labeled <- p7 + labs(title = "c) Woody") +
  theme(plot.title = element_text(hjust = 0))

p8_labeled <- p8 + labs(title = "d) Woody") +
  theme(plot.title = element_text(hjust = 0))

# Combine plots 
combined2 <- (p3_labeled + p4_labeled) / (p7_labeled + p8_labeled)

