#install packages
install.packages("phylolm.hp")
library(phylolm.hp)
library(tidyverse)
library(ape)
#install.packages("DHARMa")
library(DHARMa)
library(phylolm)
library(car)
library(phytools)
library(rr2)

####data####
#use if not removing outliers
d_pca = read_csv("Data/d_pca.csv") #pca scores for native and non-native urban plant species

phy = read.tree("Data/phy.tre")#phylogeny for species in the PCA data set
#plot(phy, type = "fan", show.tip.label = F)



####PC1-All taxa (not in maniscript)#### 

#make data frame 
d_pca = as.data.frame(d_pca_filtered)
row.names(d_pca) = d_pca$tip


#phylolm for provenance only
phym1 = phylolm(PC1 ~ provenance, data = d_pca, phy)
summary(phym1)

# Q-Q plot for phylolm model
qqnorm(residuals(phym1), main = "Q-Q Plot: phylolm")
qqline(residuals(phym1), col = "red")

#Check for phylogenetic signal (yes) TAKES A WHILE TO RUN

#order and match tree tips and rownames
#trait_vector <- d_pca$PC2
#names(trait_vector) <- rownames(d_pca)
#trait_vector <- trait_vector[phy$tip.label]
#head(trait_vector)

#phylosig(phy, trait_vector, method = "lambda", test = TRUE)


#plot relationships with climate varaibles for each pc axis
plot(d_pca$ave_tmean, d_pca$PC1)
plot(d_pca$ave_precip, d_pca$PC1)
plot(d_pca$ave_tmean, d_pca$PC2)
plot(d_pca$ave_precip, d_pca$PC2)



#Full model (no model selection)
model_25 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)
summary(model_25)


#chi square test comparing native/non-native across woddy/herbaceous

t<-table(d_pca_raw$growth_form, d_pca_raw$provenance)
chisq.test(t,correct= F )


#Variance inflation (using lm)
vif(lm(PC1 ~ provenance + ave_tmean + ave_precip + ave_precip:provenance + I(ave_precip^2)+I(ave_tmean^2), data = d_pca), type="predictor")

#check residuals with dharma
simulationOutput <- simulateResiduals(fittedModel = model_25, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

#residuals for reg qqplot
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")

#PLOTS FOR PC1
#Get means for centering
mean_tmean <- round(mean(d_pca$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(d_pca$ave_precip, na.rm = TRUE),2)

# 1. Effect of ave_precip (with quadratic term), by provenance
#precip_seq <- seq(quantile(d_pca$ave_precip, 0.025), quantile(d_pca$ave_precip, 0.975), length.out = 100)
precip_seq <-seq(min(d_pca$ave_precip), max(d_pca$ave_precip), length.out = 100)

new_precip <- expand.grid(
                 ave_precip = precip_seq,
                 provenance = c("native", "non_native")) %>%
              mutate(
                 ave_tmean = mean_tmean,
                `I(ave_precip^2)` = ave_precip^2,
                `I(ave_tmean^2)` = ave_tmean^2)

new_precip$PC1_pred <- predict(model_25, newdata = new_precip)

# Define custom colors
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Plot for ave_precip(scled to full pc1 breadth)
ggplot(new_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = d_pca_filtered, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      #geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = custom_colors) +
      labs(title = "Predicted PC1 vs Precipitation - All Taxa",
           x = "Average Precipitation",
           y = "Predicted PC1",
           color = "Provenance")+
     coord_cartesian(ylim = c(-3, 3)) +
     theme_minimal() 


# 2. Effect of ave_tmean, by provenance
#tmean_seq <- seq(quantile(d_pca$ave_tmean, 0.025), quantile(d_pca$ave_tmean, 0.975), length.out = 100)
tmean_seq <- seq(min(d_pca$ave_tmean), max(d_pca$ave_tmean), length.out = 100)
new_tmean <- expand.grid(
                ave_tmean = tmean_seq,
                provenance = c("native", "non_native")) %>%
             mutate(
               ave_precip = mean_precip,
              `I(ave_precip^2)` = ave_precip^2,
              `I(ave_tmean^2)` = ave_tmean^2)

new_tmean$PC1_pred <- predict(model_25, newdata = new_tmean)

# 2. Plot for ave_tmean
ggplot(new_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = d_pca_filtered, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = custom_colors) +
      labs(title = "Predicted PC1 vs Temperature - All Taxa",
          x = "Average Temperature",
          y = "Predicted PC1",
          color = "Provenance") +
     coord_cartesian(ylim = c(-3, 3)) +
     theme_minimal()

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
phym2 = phylolm(PC2 ~ provenance, data = d_pca, phy)
summary(phym2)

# Q-Q plot for phylolm model
qqnorm(residuals(phym2), main = "Q-Q Plot: phylolm")
qqline(residuals(phym2), col = "red")


#Check for phylogenetic signal (yes)

#phylosig(phy, trait_vector, method = "lambda", test = TRUE)


#Full model (no model selection)
model_25 <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = d_pca, phy = phy)

summary(model_25)

#residuals for reg qqplot
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")


#check residuals with dharma
simulationOutput <- simulateResiduals(fittedModel = model_25, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

##Plots

# Means for holding other variables constant
mean_tmean <- round(mean(d_pca$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(d_pca$ave_precip, na.rm = TRUE),2)

# Custom colors
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Effect of ave_tmean (with quadratic term)
#tmean_seq <- seq(quantile(d_pca$ave_precip, 0.025), quantile(d_pca$ave_precip, 0.975), length.out = 100)
tmean_seq <- seq(min(d_pca$ave_tmean), max(d_pca$ave_tmean), length.out = 100)
new_tmean <- expand.grid(ave_tmean = tmean_seq,
                         provenance = c("native", "non_native")) %>%
             mutate(ave_precip = mean_precip,
                  `I(ave_tmean^2)` = ave_tmean^2,
                  `I(ave_precip^2)` = mean_precip^2)

new_tmean$PC2_pred <- predict(model_25, newdata = new_tmean)

ggplot(new_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
            geom_line(size = 1.2) +
            geom_rug(data = d_pca, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
            #geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
            scale_color_manual(values = custom_colors) +
            labs(
              title = "Predicted PC2 vs Temperature - All Taxa",
              x = "Average Temperature",
              y = "Predicted PC2"
            ) +
            coord_cartesian(ylim = c(-3, 3)) +
            theme_minimal()

# 2. Effect of ave_precip (with quadratic term)
precip_seq <- seq(min(d_pca$ave_precip), max(d_pca$ave_precip), length.out = 100)
new_precip <- expand.grid(ave_precip = precip_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_tmean = mean_tmean,
                    `I(ave_tmean^2)` = mean_tmean^2,
                    `I(ave_precip^2)` = ave_precip^2)

new_precip$PC2_pred <- predict(model_25, newdata = new_precip)

ggplot(new_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
       geom_line(size = 1.2) +
       geom_rug(data = d_pca, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
       geom_point(data = d_pca, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
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

#herb<-d_pca%>%
      #filter(growth_form=="herb") 

herb = as.data.frame(herb_filtered)
row.names(herb) = herb$tip

#phylm
phym1 = phylolm(PC1 ~ provenance, data = herb, phy)
summary(phym1)

# Q-Q plot for phylolm model
qqnorm(residuals(phym1), main = "Q-Q Plot: phylolm")
qqline(residuals(phym1), col = "red")


#env variables
plot(herb$ave_tmean, herb$PC1)
plot(herb$ave_precip, herb$PC1)
plot(herb$ave_tmean, herb$PC2)
plot(herb$ave_precip, herb$PC2)


#Full mode (no model selection)
model_25 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb, phy = phy)

summary(model_25)

# Q-Q plot for phylolm model
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")


##plots
# Set up means
mean_tmean <- round(mean(herb$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(herb$ave_precip, na.rm = TRUE),2)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean (with squared term)
#tmean_seq <- seq(quantile(herb$ave_tmean, 0.025), quantile(herb$ave_tmean, 0.975), length.out = 100)
tmean_seq <-seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)
pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(ave_precip = mean_precip,
                    `I(ave_tmean^2)` = ave_tmean^2,
                    `I(ave_precip^2)` = mean_precip^2)

#pred_tmean$PC1_pred <- predict(model_25, newdata = pred_tmean)

preds <- predict(model_25, newdata = pred_tmean,boot=TRUE, nboot=1000)
pred_tmean$PC1_pred <- as.vector(preds$fit)
pred_tmean$PC1_se <- preds$se.fit[1:200]


# 3. Calculate 95% confidence intervals
pred_tmean2 <- pred_tmean %>%
  mutate(
    lower = PC1_pred - 1.96 * PC1_se,
    upper = PC1_pred + 1.96 * PC1_se
  )

H_PC1_temp<-ggplot(pred_tmean2, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
           geom_line(size = 1.2) +
           #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
           geom_rug(data = herb, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
           #geom_point(data = herb, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
           scale_color_manual(values = custom_colors) +
           labs(title = "Predicted PC1 vs Temperature - Herbs",
                x = "Average Temperature",
                y = "Predicted PC1") +
           coord_cartesian(ylim = c(-3, 3)) +
           theme_minimal()

# 2. Prediction over ave_precip (with squared term + interaction)
precip_seq <-seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)

#pred_precip$PC1_pred <- predict(model_25, newdata = pred_precip)




preds <- predict(model_25, newdata = pred_precip, se.fit = TRUE)
pred_precip$PC1_pred <- as.vector(preds$fit)
pred_precip$PC1_se <- preds$se.fit[1:200]


# 3. Calculate 95% confidence intervals
pred_precip2 <- pred_precip %>%
  mutate(
    lower = PC1_pred - 1.96 * PC1_se,
    upper = PC1_pred + 1.96 * PC1_se
  )


ggplot(pred_precip2, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
       geom_line(size = 1.2) +
       #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
       geom_rug(data = herb, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
       geom_point(data = herb, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
       scale_color_manual(values = custom_colors) +
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


####PC2-herbs####

phymh2 = phylolm(PC2 ~ provenance, data = herb, phy)

summary(phymh2) # still significantly different after considering phylogenetic relationship



#Full model(no model selection)

model_25 <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb, phy = phy, model="lambda")

summary(model_25)

#VIF with lm
vif(lm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb),type= 'predictor')
model_25<-lm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = herb)


# Q-Q plot for phylolm model
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")


##plot interactions
# Set up means
mean_tmean <- round(mean(herb$ave_tmean, na.rm = TRUE),2)
mean_precip <- round(mean(herb$ave_precip, na.rm = TRUE),2)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean (with squared term)
#tmean_seq <- seq(quantile(herb$ave_tmean, 0.025), quantile(herb$ave_tmean, 0.975), length.out = 100)
tmean_seq <-seq(min(herb$ave_tmean), max(herb$ave_tmean), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
                mutate(ave_precip = mean_precip,
                    `I(ave_tmean^2)` = ave_tmean^2,
                    `I(ave_precip^2)` = mean_precip^2)

#pred_tmean$PC2_pred <- predict(model_25, newdata = pred_tmean)


#pred_precip <- expand.grid(provenance = c("native", "non_native"),ave_precip = precip_seq) 


preds <- predict(model_25, newdata = pred_tmean, interval = "confidence")


pred_tmean$PC2_pred <- as.vector(preds$fit)
pred_tmean$PC2_se <- preds_list$se.fit[1:200]
plot(pred_tmean$PC2_pred , pred_tmean$PC2_se)
pred_tmean2 <- pred_tmean %>%
                mutate(
                  lower = PC2_pred - 1.96 * PC2_se,
                  upper = PC2_pred + 1.96 * PC2_se
                )


ggplot(pred_df, aes(x = ave_tmean, y = fit, color = provenance)) +
       geom_line(size = 1.2) +
       geom_rug(data = herb, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
       geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA)+
       #geom_point(data = herb, aes(x = ave_tmean, y = PC2, color = provenance),alpha = 0.4, size = 1.5)+
       scale_color_manual(values = custom_colors) +
       labs(title = "Predicted PC2 vs Temperature - Herbs",
            x = "Average Temperature",
            y = "Predicted PC2") +
       coord_cartesian(ylim = c(-5, 5)) +
       theme_minimal()

# 2. Prediction over ave_precip (with squared term + interaction)
#precip_seq <- seq(quantile(herb$ave_precip, 0.025), quantile(herb$ave_precip, 0.975), length.out = 100)
precip_seq <-seq(min(herb$ave_precip), max(herb$ave_precip), length.out = 100)

pred_precip <- expand.grid(ave_precip = precip_seq,
                           provenance = c("native", "non_native")) %>%
               mutate(ave_tmean = mean_tmean,
                     `I(ave_tmean^2)` = mean_tmean^2,
                     `I(ave_precip^2)` = ave_precip^2)

pred_precip$PC2_pred <- predict(model_25, newdata = pred_precip)

ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
       geom_line(size = 1.2) +
       geom_rug(data = herb, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
       #geom_point(data = herb, aes(x = ave_precip, y = PC2, color = provenance),alpha = 0.4, size = 1.5)+
       scale_color_manual(values = custom_colors) +
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

#woody<-d_pca%>%
      # filter(growth_form=="woody")

woody = as.data.frame(woody_filtered)
row.names(woody) = woody$tip

phywm1 = phylolm(PC1 ~ provenance, data = woody, phy)
summary(phywm1) # still significantly different after considering phylogenetic relationship


#plot pc vs climate
plot(woody$ave_tmean, woody$PC1)
plot(woody$ave_precip, woody$PC1)
plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)


#Full model (no model selection)

model_25 <- phylolm(PC1 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody, phy = phy)

summary(model_25)


# Q-Q plot for phylolm model
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")



##plot interactions

# Set up means and colors
mean_tmean <- round(mean(woody$ave_tmean, na.rm = TRUE),2)
mean_precip <-round(mean(woody$ave_precip, na.rm = TRUE),2)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean
#tmean_seq <- seq(quantile(woody$ave_tmean, 0.025), quantile(woody$ave_tmean, 0.975), length.out = 100)
tmean_seq <-seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)

pred_tmean <- expand.grid(ave_tmean = tmean_seq,
                          provenance = c("native", "non_native")) %>%
              mutate(
                ave_precip = mean_precip,
                `I(ave_tmean^2)` = ave_tmean^2,
                `I(ave_precip^2)` = mean_precip^2
              )

pred_tmean$PC1_pred <- predict(model_25, newdata = pred_tmean)

ggplot(pred_tmean, aes(x = ave_tmean, y = PC1_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = woody, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      geom_point(data = woody, aes(x = ave_tmean, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = custom_colors) +
      labs(
        title = "Predicted PC1 vs Temperature-Woody",
        x = "Average Temperature",
        y = "Predicted PC1") +
      coord_cartesian(ylim=c(-3,3))+
      theme_minimal()



# 2. Prediction over ave_precip
#precip_seq <- seq(quantile(woody$ave_precip, 0.025), quantile(woody$ave_precip, 0.975), length.out = 100)
precip_seq <-seq(min(woody$ave_precip), max(woody$ave_precip), length.out = 100)

pred_precip <- expand.grid(
              ave_precip = precip_seq,
              provenance = c("native", "non_native")
            ) %>%
              mutate(
                ave_tmean = mean_tmean,
                `I(ave_tmean^2)` = mean_tmean^2,
                `I(ave_precip^2)` = ave_precip^2
              )

pred_precip$PC1_pred <- predict(model_25, newdata = pred_precip)

ggplot(pred_precip, aes(x = ave_precip, y = PC1_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = woody, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      geom_point(data = woody, aes(x = ave_precip, y = PC1, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = custom_colors) +
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
#
phywm1 = phylolm(PC2 ~ provenance, data = woody, phy)
summary(phywm1) # still significantly different after considering phylogenetic relationship



plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)
plot(woody$ave_tmean, woody$PC2)
plot(woody$ave_precip, woody$PC2)



model_25 <- phylolm(PC2 ~ provenance + ave_tmean + ave_precip + ave_tmean:provenance + ave_precip:provenance + I(ave_tmean^2) + I(ave_precip^2), data = woody, phy = phy)


summary(model_25)

# Q-Q plot for phylolm model
qqnorm(residuals(model_25), main = "Q-Q Plot: phylolm")
qqline(residuals(model_25), col = "red")



##plot interactions

# Set up means and colors
mean_tmean <- mean(woody$ave_tmean, na.rm = TRUE)
mean_precip <- mean(woody$ave_precip, na.rm = TRUE)
custom_colors <- c("native" = "blue", "non_native" = "orange")

# 1. Prediction over ave_tmean
#tmean_seq <- seq(quantile(woody$ave_tmean, 0.025), quantile(woody$ave_tmean, 0.975), length.out = 100)

tmean_seq <-seq(min(woody$ave_tmean, na.rm = TRUE), max(woody$ave_tmean, na.rm = TRUE), length.out = 100)

pred_tmean <- expand.grid(
              ave_tmean = tmean_seq,
              provenance = c("native", "non_native")
            ) %>%
              mutate(
                ave_precip = mean_precip,
                `I(ave_tmean^2)` = ave_tmean^2,
                `I(ave_precip^2)` = mean_precip^2
              )

pred_tmean$PC2_pred <- predict(model_25, newdata = pred_tmean)

ggplot(pred_tmean, aes(x = ave_tmean, y = PC2_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = woody, aes(x = ave_tmean), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      #geom_point(data = woody, aes(x = ave_tmean, y = PC2, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = custom_colors) +
      labs(
        title = "Predicted PC2 vs Temperature-Woody",
        x = "Average Temperature",
        y = "Predicted PC2"
      ) +
      ylim(-3, 3) +
      theme_minimal()

# 2. Prediction over ave_precip
#precip_seq <- seq(quantile(woody$ave_precip, 0.025), quantile(woody$ave_precip, 0.975), length.out = 100)
precip_seq <-seq(min(woody$ave_precip, na.rm = TRUE), max(woody$ave_precip, na.rm = TRUE), length.out = 100)

pred_precip <- expand.grid(
              ave_precip = precip_seq,
              provenance = c("native", "non_native")
            ) %>%
              mutate(
                ave_tmean = mean_tmean,
                `I(ave_tmean^2)` = mean_tmean^2,
                `I(ave_precip^2)` = ave_precip^2
              )

pred_precip$PC2_pred <- predict(model_25, newdata = pred_precip)

ggplot(pred_precip, aes(x = ave_precip, y = PC2_pred, color = provenance)) +
      geom_line(size = 1.2) +
      geom_rug(data = woody, aes(x = ave_precip), inherit.aes = FALSE, sides = "b", alpha = 0.3)+
      geom_point(data = woody, aes(x = ave_precip, y = PC2, color = provenance),alpha = 0.4, size = 1.5)+
      scale_color_manual(values = custom_colors) +
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


#add raw data to plots

ggplot() +
  # Raw data points
  geom_point(data = woody, aes(x = ave_tmean, y = PC2, color = provenance),
             alpha = 0.2, size = 1.5) +
  
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


















      
