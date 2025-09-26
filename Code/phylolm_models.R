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
d_pca = read_csv("Data/d_pca.csv") 

#phylogeny for species in the PCA data set
phy = read.tree("Data/phy.tre") 

#Plot Phylogeny
plot(phy, type = "fan", show.tip.label = F)


####PC1-All taxa (not in manuscript)#### 

#make data frame 
d_pca <- as.data.frame(d_pca)
row.names(d_pca) <- d_pca$tip

#scale 
d_pca$ave_tmean<-scale(d_pca$ave_tmean)[,1]
d_pca$ave_precip<-scale(d_pca$ave_precip)[,1]

#phylolm for provenance only
phym1 = phylolm(PC1 ~ provenance, data = d_pca, phy, model="lambda")
summary(phym1)

# Q-Q plot for phylolm model
qqnorm(residuals(phym1), main = "Q-Q Plot: phylolm")
qqline(residuals(phym1), col = "red")


#Full model (no model selection)
model_pc1_All <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy, model= "lambda")
summary(model_pc1_All)


#chi square test comparing native/non-native across woody/herbaceous
t<-table(d_pca$growth_form, d_pca$provenance)
chisq.test(t,correct= F )


#Variance inflation (using lm)
vif(lm(PC1 ~ provenance + ave_tmean + ave_precip + ave_precip:provenance + I(ave_precip^2)+I(ave_tmean^2), data = d_pca), type="predictor")

#residuals for reg qqplot
qqnorm(residuals(model_pc1_All), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc1_All), col = "red")

#PLOTS FOR PC1
#Get means for centering
mean_tmean <- round(mean(d_pca$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(d_pca$ave_precip, na.rm = TRUE),2)

# Plot Effect of ave_precip (with quadratic term), by provenance
#precip_seq <- seq(quantile(d_pca$ave_precip, 0.025), quantile(d_pca$ave_precip, 0.975), length.out = 100)

#Make prediction grid
precip_seq <-seq(min(d_pca$ave_precip), max(d_pca$ave_precip), length.out = 100)

new_precip <- expand.grid(
                          ave_precip = precip_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(
                     ave_tmean = mean_tmean,
                    `I(ave_precip^2)` = ave_precip^2,
                    `I(ave_tmean^2)` = ave_tmean^2)

#predict on the grid
new_precip$PC1_pred <- predict(model_pc1_All, newdata = new_precip)



# 1. Plot for ave_precip(scaled to full pc1 breadth)
ggplot(new_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = d_pca, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      #geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
      labs(title = "Predicted PC1 vs Precipitation - All Taxa",
           x = "Average Precipitation",
           y = "Predicted PC1",
           color = "Provenance")+
      coord_cartesian(ylim = c(-3, 3)) +
      theme_minimal() 


# 2. Effect of ave_tmean, by provenance

#make prediction grid
#tmean_seq <- seq(quantile(d_pca$ave_tmean, 0.025), quantile(d_pca$ave_tmean, 0.975), length.out = 100)
tmean_seq <- seq(min(d_pca$ave_tmean), max(d_pca$ave_tmean), length.out = 100)

new_tmean <- expand.grid(
                        ave_tmean = tmean_seq,
                        provenance = c("native", "non_native")) %>%
            mutate(
                   ave_precip = mean_precip,
                  `I(ave_precip^2)` = ave_precip^2,
                  `I(ave_tmean^2)` = ave_tmean^2)
#predict onto grid
new_tmean$PC1_pred <- predict(model_pc1_All, newdata = new_tmean)

# 2. Plot for ave_tmean
ggplot(new_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = d_pca, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      #geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
      labs(title = "Predicted PC1 vs Temperature - All Taxa",
           x = "Average Temperature",
           y = "Predicted PC1",
           color = "Provenance") +
      coord_cartesian(ylim = c(-3, 3)) +
      theme_minimal()

summary(model_pc1_All)

#3 Create a simple prediction dataset with only provenance varying
box_data <- data.frame(ave_tmean = mean(d_pca$ave_tmean, na.rm = TRUE),
                       ave_precip = mean(d_pca$ave_precip, na.rm = TRUE),
                       provenance = c("native", "non_native")) %>%
  mutate(`I(ave_precip^2)` = ave_precip^2,
         `I(ave_tmean^2)` = ave_tmean^2)

# Predict PC1
box_data$PC1_pred <- predict(model_25, newdata = box_data)


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
  labs(title = "Predicted PC1 vs Provenance - All Taxa",
       x = "Provenance",
       y = "Predicted PC1") +
  theme_minimal()

####PC2- All taxa (not in manuscript)####
#phylm
phym2 = phylolm(PC2 ~ provenance, data = d_pca, phy, model= "lambda")
summary(phym2)

# Q-Q plot for phylolm model
qqnorm(residuals(phym2), main = "Q-Q Plot: phylolm")
qqline(residuals(phym2), col = "red")

#Full model (no model selection)
model_pc2_All <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy, model= "lambda")

summary(model_pc2_All)

#residuals for reg qqplot
qqnorm(residuals(model_pc2_All), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc2_All), col = "red")


##Plots

# Means for holding other variables constant
mean_tmean <- round(mean(d_pca$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(d_pca$ave_precip, na.rm = TRUE),2)

# 1. Effect of ave_tmean (with quadratic term)

#make prediction grid
#tmean_seq <- seq(quantile(d_pca$ave_precip, 0.025), quantile(d_pca$ave_precip, 0.975), length.out = 100)
tmean_seq <- seq(min(d_pca$ave_tmean), max(d_pca$ave_tmean), length.out = 100)

new_tmean <- expand.grid(ave_tmean = tmean_seq,
                         provenance = c("native", "non_native")) %>%
             mutate(ave_precip = mean_precip,
                   `I(ave_tmean^2)` = ave_tmean^2,
                   `I(ave_precip^2)` = mean_precip^2)

#predict onto grid
new_tmean$PC2_pred <- predict(model_pc2_All, newdata = new_tmean)

ggplot(new_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = d_pca, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      #geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
      labs(
        title = "Predicted PC2 vs Temperature - All Taxa",
        x = "Average Temperature",
        y = "Predicted PC2"
      ) +
      coord_cartesian(ylim = c(-3, 3)) +
      theme_minimal()

# 2. Effect of ave_precip (with quadratic term)
#make prediction grid
precip_seq <- seq(min(d_pca$ave_precip), max(d_pca$ave_precip), length.out = 100)

new_precip <- expand.grid(ave_precip = precip_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)

#predict onto grid
new_precip$PC2_pred <- predict(model_pc2_All, newdata = new_precip)

ggplot(new_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = d_pca, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      #geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = custom_colors) +
      labs(title = "Predicted PC2 vs Precipitation - All Taxa",
           x = "Average Precipitation",
           y = "Predicted PC2") +
      coord_cartesian(ylim = c(-3, 3)) +
      theme_minimal()


#3 box plot with only provenance varying


#  Create prediction grid
box_data <- expand.grid(ave_tmean = mean(d_pca$ave_tmean, na.rm = TRUE),
                        ave_precip = mean(d_pca$ave_precip, na.rm = TRUE),
                        provenance = c("native", "non_native")) %>%
  mutate(`I(ave_tmean^2)` = ave_tmean^2,
         `I(ave_precip^2)` = ave_precip^2)

#  Predict PC2 using your phylolm model
box_data$PC2_pred <- predict(model_25, newdata = box_data)

# (optional): Expand to simulate distribution for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC2_pred <- box_data_expanded$PC2_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

# Plot
ggplot(box_data_expanded, aes(x = provenance, y = PC2_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("native" = "blue", "non_native" = "orange")) +
  labs(title = "Predicted PC2 by Provenance - All Taxa",
       x = "Provenance",
       y = "Predicted PC2") +
  coord_cartesian(ylim = c(-1.25, 1)) +
  theme_minimal()
####PC1-herb####


####data
d_pca = read_csv("Data/d_pca.csv") #pca scores for native and non-native urban plant species

phy = read.tree("Data/phy.tre") #phylogeny for species in the PCA data set
#plot(phy, type = "fan", show.tip.label = F)

#filter for herbs
herb<-d_pca%>%
      filter(growth_form=="herb") 
herb <- as.data.frame(herb)
row.names(herb) = herb$tip


#scale
herb$ave_tmean<-scale(herb$ave_tmean)[,1]
herb$ave_precip<-scale(herb$ave_precip)[,1]

#phylm
phym1 = phylolm(PC1 ~ provenance, data = herb, phy, model="lambda")
summary(phym1)

# Q-Q plot for phylolm model
qqnorm(residuals(phym1), main = "Q-Q Plot: phylolm")
qqline(residuals(phym1), col = "red")



#Full mode (no model selection)
model_pc1_Herb <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb, phy = phy, model="lambda")

summary(model_pc1_Herb)
# Q-Q plot for phylolm model
qqnorm(residuals(model_pc1_Herb), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc1_Herb), col = "red")


##plots
# Set up means
mean_tmean <- round(mean(herb$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(herb$ave_precip, na.rm = TRUE),2)


# 1. Prediction over ave_tmean (with squared term)
#make prediction grid
#tmean_seq <- seq(quantile(herb$ave_tmean, 0.025), quantile(herb$ave_tmean, 0.975), length.out = 100)
tmean_seq <-seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                     `I(ave_tmean^2)` = ave_tmean^2,
                     `I(ave_precip^2)` = mean_precip^2)

#predict onto grid
pred_tmean$PC1_pred <- predict(model_pc1_Herb, newdata = pred_tmean)

H_PC1_temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
                  geom_line(size = 1.2) +
                  #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
                  geom_rug(data = herb, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
                  #geom_point(data = herb, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
                  scale_color_manual(values =c("native" = "blue", "non_native" = "orange")) +
                  labs(title = "Predicted PC1 vs Temperature - Herbs",
                       x = "Average Temperature",
                       y = "Predicted PC1") +
                  coord_cartesian(ylim = c(-3, 3)) +
                  theme_minimal()

# 2. Prediction over ave_precip (with squared term + interaction)

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
                    #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
                    geom_rug(data = herb, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
                    #geom_point(data = herb, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
                    scale_color_manual(values =c("native" = "blue", "non_native" = "orange")) +
                    labs(title = "Predicted PC1 vs Precipitation- Herbs",
                         x = "Average Precipitation",
                         y = "Predicted PC1") +
                    coord_cartesian(ylim = c(-3, 3)) +
                    theme_minimal()



# 3. Boxplot of predicted PC1 by provenance
box_data <- expand.grid(ave_tmean = mean_tmean,
                        ave_precip = mean_precip,
                        provenance = c("native", "non_native")) %>%
                mutate(`I(ave_tmean^2)` = ave_tmean^2,
                       `I(ave_precip^2)` = ave_precip^2)

box_data$PC1_pred <- predict(model_25, newdata = box_data)

# Optional: simulate variability for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC1_pred <- box_data_expanded$PC1_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

ggplot(box_data_expanded, aes(x = provenance, y = PC1_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Predicted PC1 by Provenance - Herbs",
       x = "Provenance",
       y = "Predicted PC1") +
  coord_cartesian(ylim = c(-3, 1)) +
  theme_minimal()

phylolm
####PC2-herbs####

phymh2 = phylolm(PC2 ~ provenance, data = herb, phy, model="lambda")

summary(phymh2) # still significantly different after considering phylogenetic relationship



#Full model(no model selection)

model_pc2_Herb <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb, phy = phy, model="lambda")

summary(model_pc2_Herb)


#VIF with lm
vif(lm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb),type= 'predictor')


# Q-Q plot for phylolm model
qqnorm(residuals(model_pc2_Herb ), main = "Q-Q Plot: phylolm")
qqline(residuals(model_pc2_Herb ), col = "red")


##plot interactions
# Set up means
mean_tmean <- round(mean(herb$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(herb$ave_precip, na.rm = TRUE),2)


# 1. Prediction over ave_tmean (with squared term)
#tmean_seq <- seq(quantile(herb$ave_tmean, 0.025), quantile(herb$ave_tmean, 0.975), length.out = 100)

#make prediction grid
tmean_seq <-seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
               mutate(ave_precip = mean_precip,
                     `I(ave_tmean^2)` = ave_tmean^2,
                     `I(ave_precip^2)` = mean_precip^2)

#predict onto grid
pred_tmean$PC2_pred <- predict(model_pc2_Herb, newdata = pred_tmean)



H_PC2_temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
        geom_line(size = 1.2) +
        geom_rug(data = herb, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
        #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA)+
        #geom_point(data = herb, aes(x = ave_tmean, y = PC2, color = provenance),alpha = 0.4, size = 1.5)+
        scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
        labs(title = "Predicted PC2 vs Temperature - Herbs",
             x = "Average Temperature",
             y = "Predicted PC2") +
        coord_cartesian(ylim = c(-3, 3)) +
        theme_minimal()

# 2. Prediction over ave_precip (with squared term + interaction)
#precip_seq <- seq(quantile(herb$ave_precip, 0.025), quantile(herb$ave_precip, 0.975), length.out = 100)

#make prediction grid
precip_seq <-seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)

pred_precip$PC2_pred <- predict(model_pc2_Herb, newdata = pred_precip)

H_PC2_precip<-ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
                    geom_line(size = 1.2) +
                    geom_rug(data = herb, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
                    #geom_point(data = herb, aes(x = ave_precip, y = PC2, color = provenance),alpha = 0.4, size = 1.5)+
                    scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
                    labs(title = "Predicted PC2 vs Precipitation- Herbs",
                         x = "Average Precipitation",
                         y = "Predicted PC2") +
                    coord_cartesian(ylim = c(-3, 3)) +
                    theme_minimal()



# 3. Boxplot of predicted PC2 by provenance
box_data <- expand.grid(ave_tmean = mean_tmean,
                        ave_precip = mean_precip,
                        provenance = c("native", "non_native")) %>%
  mutate(`I(ave_tmean^2)` = ave_tmean^2,
         `I(ave_precip^2)` = ave_precip^2)

box_data$PC2_pred <- predict(model_25, newdata = box_data)

# Optional: simulate variability for boxplot
set.seed(123)
box_data_expanded <- box_data[rep(1:nrow(box_data), each = 20), ]
box_data_expanded$PC2_pred <- box_data_expanded$PC2_pred + rnorm(nrow(box_data_expanded), sd = 0.1)

ggplot(box_data_expanded, aes(x = provenance, y = PC2_pred, fill = provenance)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Predicted PC2 by Provenance - Herbs",
       x = "Provenance",
       y = "Predicted PC2") +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal()


####PC1-woody####

woody<-d_pca%>%
       filter(growth_form=="woody")

woody = as.data.frame(woody)
row.names(woody) = woody$tip


#scale
woody$ave_tmean<-scale(woody$ave_tmean)[,1]
woody$ave_precip<-scale(woody$ave_precip)[,1]



phywm1 = phylolm(PC1 ~ provenance, data = woody, phy=phy, model="lambda")
summary(phywm1) # still significantly different after considering phylogenetic relationship


#Full model (no model selection)

model_pc1_woody <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody, phy = phy, model="lambda")

summary(model_pc1_woody)
# Q-Q plot for phylolm model
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")



##plot interactions

# Set up means and colors
mean_tmean <- round(mean(woody$ave_tmean, na.rm = TRUE),2)
mean_precip <-round(mean(woody$ave_precip, na.rm = TRUE),2)


# 1. Prediction over ave_tmean
#make prediction grid
#tmean_seq <- seq(quantile(woody$ave_tmean, 0.025), quantile(woody$ave_tmean, 0.975), length.out = 100)
tmean_seq <-seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                    `I(ave_tmean^2)` = ave_tmean^2,
                    `I(ave_precip^2)` = mean_precip^2)

pred_tmean$PC1_pred <- predict(model_pc1_woody, newdata = pred_tmean)

W_Pc1_Temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
                  geom_line(size = 1.2) +
                  geom_rug(data = woody, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
                  #geom_point(data = woody, aes(x = ave_tmean, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
                  scale_color_manual(values =  c("native" = "blue", "non_native" = "orange")) +
                  labs(
                    title = "Predicted PC1 vs Temperature-Woody",
                    x = "Average Temperature",
                    y = "Predicted PC1") +
                  coord_cartesian(ylim=c(-3,3))+
                  theme_minimal()



# 2. Prediction over ave_precip
#make prediction grid
#precip_seq <- seq(quantile(woody$ave_precip, 0.025), quantile(woody$ave_precip, 0.975), length.out = 100)
precip_seq <-seq(min(woody$ave_precip), max(woody$ave_precip), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
                mutate(ave_tmean = mean_tmean,
                      `I(ave_tmean^2)` = mean_tmean^2,
                      `I(ave_precip^2)` = ave_precip^2)

pred_precip$PC1_pred <- predict(model_pc1_woody, newdata = pred_precip)

W_Pc1_Precip<-ggplot(pred_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
                    geom_line(size = 1.2) +
                    geom_rug(data = woody, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
                   # geom_point(data = woody, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
                    scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
                    labs(
                      title = "Predicted PC1 vs Precipitation -Woody",
                      x = "Average Precipitation",
                      y = "Predicted PC1") +
                    coord_cartesian(ylim=c(-3,3))+
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
box_data$PC1_pred <- predict(model_25, newdata = box_data)

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
  ylim(-1, 1) +  # Set y-axis range here
  theme_minimal()




####PC2-Woody####
#model with provenance only
phywm1 = phylolm(PC2 ~ provenance, data = woody, phy, model="lambda")
summary(phywm1) # still significantly different after considering phylogenetic relationship


#full model (no model selection)
model_pc2_woody <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody, phy = phy, model="lambda")

summary(model_pc2_woody)

# Q-Q plot for phylolm model
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")



##plot interactions

# Set up means and colors
mean_tmean <- mean(woody$ave_tmean, na.rm = TRUE)
mean_precip <- mean(woody$ave_precip, na.rm = TRUE)

# 1. Prediction over ave_tmean

#make prediction grid
#tmean_seq <- seq(quantile(woody$ave_tmean, 0.025), quantile(woody$ave_tmean, 0.975), length.out = 100)

tmean_seq <-seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                    `I(ave_tmean^2)` = ave_tmean^2,
                    `I(ave_precip^2)` = mean_precip^2)

pred_tmean$PC2_pred <- predict(model_pc2_woody, newdata = pred_tmean)

W_Pc2_Temp<-ggplot(pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
                  geom_line(size = 1.2) +
                  geom_rug(data = woody, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
                  #geom_point(data = woody, aes(x = ave_tmean, y = PC2, color = provenance),alpha = 0.4, size = 1.5)+
                  scale_color_manual(values =c("native" = "blue", "non_native" = "orange")) +
                  labs(
                    title = "Predicted PC2 vs Temperature-Woody",
                    x = "Average Temperature",
                    y = "Predicted PC2"
                  ) +
                  ylim(-3, 3) +
                  theme_minimal()

# 2. Prediction over ave_precip

#make prediction grid
#precip_seq <- seq(quantile(woody$ave_precip, 0.025), quantile(woody$ave_precip, 0.975), length.out = 100)
precip_seq <-seq(min(woody$ave_precip, na.rm = TRUE), max(woody$ave_precip, na.rm = TRUE), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                      `I(ave_tmean^2)` = mean_tmean^2,
                      `I(ave_precip^2)` = ave_precip^2)

pred_precip$PC2_pred <- predict(model_pc2_woody, newdata = pred_precip)

W_Pc2_Precip<-ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
              geom_line(size = 1.2) +
              geom_rug(data = woody, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
              #geom_point(data = woody, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
              scale_color_manual(values = c("native" = "blue", "non_native" = "orange")) +
              labs(
                title = "Predicted PC2 vs Precipitation -Woody",
                x = "Average Precipitation",
                y = "Predicted PC2") +
              ylim(-3, 3) +
              theme_minimal()

# 3. Boxplot of predicted PC2 by provenance


# Create prediction data
box_data <- expand.grid(
  ave_tmean = mean_tmean,
  ave_precip = mean_precip,
  provenance = c("native", "non_native")) %>%
  mutate(
    `I(ave_tmean^2)` = ave_tmean^2,
    `I(ave_precip^2)` = ave_precip^2
  )

# Predict PC2
box_data$PC2_pred <- predict(model_25, newdata = box_data)

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
  ylim(-1, 1) +  # Set y-axis range here
  theme_minimal()


#make panel figure
library(ggplot2)
library(patchwork)

# Suppose you have 8 plots: p1, p2, ..., p8

# Remove titles or legends from individual plots:
p1 <- H_PC1_temp + 
      labs(y = "PC 1 (Herbaceous)")+
      theme(plot.title = element_blank(), 
            legend.position = "none", 
            axis.title.x = element_blank(),
            )



p2 <- H_PC1_precip + 
          labs(y = "PC 1 (Herbaceous)")+
          theme(plot.title = element_blank(), 
                axis.title.x = element_blank(),
                legend.position = c(0.99, 0.99),   # inside top right
                legend.justification = c("right", "top"),
                legend.background = element_rect(fill = alpha('white', 0.7), color = NA)
                
          )

p3 <- H_PC2_temp + 
            labs(y = "PC 2 (Herbaceous)") +
            theme(
              plot.title = element_blank(), 
              axis.title.x = element_blank(),
              legend.position = "none",   # inside top right
           )



p4 <- H_PC2_precip + 
          labs(y = "PC 2 (Herbaceous)") +
          theme(
            plot.title = element_blank(), 
            axis.title.x = element_blank(),
            legend.position = c(0.99, 0.99),   # inside top right
            legend.justification = c("right", "top"),
            legend.background = element_rect(fill = alpha('white', 0.7), color = NA)
            
            
            )

p5 <- W_Pc1_Temp +   
      labs(y = "PC 1 (Woody)")+
      theme(plot.title = element_blank(), 
        legend.position = "none")

p6 <- W_Pc1_Precip + 
      labs(y = "PC 1 (Woody)")+
      theme(plot.title = element_blank(), 
            legend.position = "none")
p7 <- W_Pc2_Temp +  
      labs(y = "PC 2 (Woody)")+
      theme(plot.title = element_blank(), 
        legend.position = "none")

p8 <- W_Pc2_Precip +  
      labs(y = "PC 2 (Woody)")+
      theme(plot.title = element_blank(), 
            legend.position = "none")

# Combine in 2 rows x 4 columns grid:
combined <- (p1 + p2)/ (p5 + p6)+
             plot_annotation(tag_levels = 'a', 
                              tag_prefix = "", 
                              tag_suffix = ")")

combined2<-(p3 + p4)  /(p7 + p8)+
            plot_annotation(tag_levels = 'a', 
                            tag_prefix = "", 
                            tag_suffix = ")")


print(combined)
print(combined2)





