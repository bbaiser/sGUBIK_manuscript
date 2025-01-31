#PCA AND KS tests of pca score distributions for All, native, and non-native urban plants
#figure 2

####Set up####

#load packages
library(plyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(devtools)
library(cowplot)
library(hypervolume)
library(dplyr)
library(visdat)
library(purrr)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(kSamples)


#Load data sets 

#Data from Diaz et al. 2022 and GBIF for urban for speices that have all 6 diaz traits. taxonomically matched to WCVP 
Diaz_Final <- read.csv("Data/Diaz_Gubic_final_tot2.csv",row=1)


#counts of native (n=1564) and non-native(n=1213)
prov_counts<-Diaz_Final %>%
            group_by(provenance_glonaf) %>%
            tally 
prov_counts

#### PCA for all species####
#scale data
z_Diaz_Final<-scale(Diaz_Final[, c("Leaf_area","LMA","Leaf_N","Seed_mass","Stem_density","Height")])

#run PCA
PCA <- prcomp(z_Diaz_Final)

#look at PCA results and extract loadings and scores
summary(PCA) 
PCAvalues   <- data.frame(Species = Diaz_Final$Species_WCVP, provenance = Diaz_Final$provenance_glonaf, growth_form= Diaz_Final$growth_form, PCA$x)# Extract PC axes for plotting
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation) # Extract loading of the variables


#change axis direction to match Diaz et al. if necessary
#Axis 1
#PCAvalues$PC1 <- (PCAvalues$PC1*-1)
#PCAloadings$PC1 <- (PCAloadings$PC1*-1)

#Axis2
PCAvalues$PC2 <- (PCAvalues$PC2*-1)
PCAloadings$PC2 <- (PCAloadings$PC2*-1)


###plot pca and pca score distributions

#set themes
p_my_theme <-  theme( legend.text = element_text(size=9), axis.title= element_text(size=9),
                      axis.title.x= element_text(colour="black", size=9),
                      axis.title.y= element_text(colour="black", size=9),
                      axis.text.x= element_text(colour= "black", size=9),
                      axis.text.y= element_text(colour= "black", size=9),
                      panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank(),
                      panel.border=element_rect(fill=NA,colour="grey50"))

p_my_theme2 <- theme( legend.text = element_text(size=9), axis.title= element_text(size=9),
                      axis.text.x= element_text(colour= "black", size=9), axis.text.y= element_text(colour= "black", size=9))

All_Sp_PCA <- ggplot(PCAvalues) +   
              geom_point(size = .7, alpha=0.5, 
                         aes(x = PC1, y = PC2, group = provenance, colour = provenance, size = provenance),  # label = T
                         show.legend = FALSE) +
              # geom_text( aes(x = PC1, y = PC2, label = Species, colour = source, size = source), hjust = 0)+  # add species names
              coord_fixed() + 
              scale_colour_manual(values=c("blue1","orange")) +
              geom_segment(data = PCAloadings, size = 1,
                           aes(x = 0, xend = PC1*4.5, y = 0, yend = PC2*4.5),
                           arrow = arrow(length = unit(.3, "cm")),colour = "black")   +
              geom_text(data = PCAloadings, aes(x = PC1*5, y = PC2*5, label = Variables), size = 5,
                        hjust=c(0, 0, 0, 0, 0, 0) , vjust=c(0, 0, 0, 0, 0, 1))    + 
              xlab("PC1 (47%) ") + ylab("PC2 (26%)")  +
              # ggtitle("All Species")+
              xlim(c(-5.5,5.5))+
              ylim(c(-5,5))+
              p_my_theme


# Marginal density distribution exotic and native PCA scores
#axis 1
xplot <- ggdensity(PCAvalues, "PC1", fill = "provenance", color = "provenance", palette = c("blue","orange")) +
          p_my_theme2+
          theme(legend.position = "none")



#axis2 (can get rid of rotate)
yplot <- ggdensity(PCAvalues, "PC2", fill = "provenance", color = "provenance", palette = c("blue","orange"))  + 
          p_my_theme2+
          theme(legend.position = "none")



#Boxplots
boxX<-ggplot(PCAvalues, aes(x=provenance, y=PC1, fill=provenance)) +
      geom_boxplot(fill = c("blue","orange"), width=0.5)+
      p_my_theme+
      xlab(" ") +
      coord_flip()+
      theme_void() +
      theme(legend.position = "none")+
      scale_x_discrete(expand = c(2, 2))


boxY<-ggplot(PCAvalues, aes(x=provenance, y=PC2, fill=provenance)) +
      geom_boxplot(fill = c("blue","orange"), width=0.5, outlier.shape = NA)+
      p_my_theme+
      xlab(" ") +
      coord_flip()+
      theme_void() +
      theme(legend.position = "none")+
      scale_x_discrete(expand = c(2, 2))


####Test distributions of pc scores for native and non-native species on the first two axes####

# extract pca values for the first 2 axes for native and non-native
nm_native <- filter(PCAvalues, provenance == "native") %>% 
             select(PC1,PC2)
nm_non_native <- filter(PCAvalues, provenance == "non_native")%>% 
                 select(PC1,PC2)


#KS-TEST
ks.test(nm_non_native$PC1, nm_native$PC1)#D = 0.28431, p-value < 2.2e-16

ks.test(nm_non_native$PC2, nm_native$PC2)#D = 0.07241, p-value = 0.001548

#t-TEST (data violates the assumptions of this test)
t.test(nm_non_native$PC1, nm_native$PC1)#t = -15.564, df = 2598.7, p-value < 2.2e-16

t.test(nm_non_native$PC2, nm_native$PC2)#t = 3.326, df = 2688.4, p-value = 0.0008929

#AD-TEST
ad.test(nm_non_native$PC1, nm_native$PC1)#AD =111.5,  T.AD= 145.29, P= 1.5802e-61

ad.test(nm_non_native$PC2, nm_native$PC2)#AD =5.8953, T.AD= 6.4338, P=0.0010553



####PCA score distributions and tests for herbaceous species only####
#Filter PCA for herb species
PCAvalues_herb<- PCAvalues%>%
                filter(growth_form=="herb") 

#counts of native (n=467) and non-native(n=707) herbs
prov_herb_counts<-PCAvalues_herb %>%
                  group_by(provenance) %>%
                  tally 
prov_herb_counts

#Urban NPPS versus urban natives for Herbs only
# Marginal density distribution   
xplot_herb <- ggdensity(PCAvalues_herb, "PC1", fill = "provenance", color = "provenance", palette = c("blue","orange"))+
              ggtitle("Herbaceous Species")+
              p_my_theme2+
              theme(legend.position = "none")

yplot_herb <- ggdensity(PCAvalues_herb, "PC2", fill = "provenance", color = "provenance", palette = c("blue","orange"))+
              p_my_theme2+
              ggtitle(" ")+
              theme(legend.position = "none")

#Boxplots
boxX_herb<-ggplot(PCAvalues_herb, aes(x=provenance, y=PC1, fill=provenance)) +
            geom_boxplot(fill = c("blue","orange"), width=0.5, outlier.shape = NA)+
            p_my_theme+
            xlab(" ") +
            coord_flip()+
            theme_void() +
            theme(legend.position = "none")+
            scale_x_discrete(expand = c(2, 2))

boxY_herb<-ggplot(PCAvalues_herb, aes(x=provenance, y=PC2, fill=provenance)) +
            geom_boxplot(fill = c("blue","orange"), width=0.5, outlier.shape = NA)+
            p_my_theme+
            xlab(" ") +
            coord_flip()+
            theme_void() +
            theme(legend.position = "none")+
            scale_x_discrete(expand = c(2, 2))



#Test distributions of pca scores for native and non-native herbaceous species
# extract pca values for the first two axes for herbaceous species

nm_native_herb <- filter(PCAvalues_herb, provenance == "native") %>% 
                  select(PC1,PC2)
nm_non_native_herb <- filter(PCAvalues_herb, provenance == "non_native")%>% 
                      select(PC1,PC2)

#KS-TEST

ks.test(nm_native_herb$PC1, nm_non_native_herb$PC1)#D = 0.059848, p-value = 0.2661

ks.test(nm_non_native_herb$PC2, nm_native_herb$PC2)#D = 0.21329, p-value = 1.541e-11

#ttest
t.test(nm_native_herb$PC1, nm_non_native_herb$PC1)#t = 1.1501, df = 964.08, p-value = 0.2504

t.test(nm_non_native_herb$PC2, nm_native_herb$PC2)#t = 8.1786, df = 975.45, p-value = 8.89e-16

#ad test
ad.test(nm_non_native_herb$PC1, nm_native_herb$PC1)# AD =1.5127, T.AD= 0.67441, P=0.17323

ad.test(nm_non_native_herb$PC2, nm_native_herb$PC2)#AD= 32.305, T.AD= 41.182, P=5.4227e-18

####PCA score distributions and tests for WOODY species only####
#filter woody/shrub pca scores
PCAvalues_woody<- PCAvalues%>%
                  filter(growth_form=="tree"|growth_form=="shrub") 

#counts of native (n=1097) and non-native(n=506) woody
prov_woody_counts<-PCAvalues_woody %>%
                  group_by(provenance) %>%
                  tally 
prov_woody_counts

# Marginal density distribution woddy native and non_native PCA scores  
xplot_wood <- ggdensity(PCAvalues_woody, "PC1", fill = "provenance", color = "provenance", palette = c("blue","orange")) + 
              p_my_theme2+
              ggtitle("Woody Species")+
              theme(legend.position = "none")

yplot_wood <- ggdensity(PCAvalues_woody, "PC2", fill = "provenance", color = "provenance", palette = c("blue","orange")) +
              p_my_theme2+
              ggtitle(" ")+
              theme(legend.position = "none")
#Boxplots
boxX_wood<-ggplot(PCAvalues_woody, aes(x=provenance, y=PC1, fill=provenance)) +
          geom_boxplot(fill = c("blue","orange"), width=0.5, outlier.shape = NA)+
          p_my_theme+
          xlab(" ") +
          coord_flip()+
          theme_void() +
          theme(legend.position = "none")+
          scale_x_discrete(expand = c(2, 2))

boxY_wood<-ggplot(PCAvalues_woody, aes(x=provenance, y=PC2, fill=provenance)) +
          geom_boxplot(fill = c("blue","orange"), width=0.5, outlier.shape = NA)+
          p_my_theme+
          xlab(" ") +
          coord_flip()+
          theme_void() +
          theme(legend.position = "none")+
          scale_x_discrete(expand = c(2, 2))

#Test distributions of pc scores for native and non-native woody species on the first two axes
# extract pca values for the first 2 axes

nm_native_wood <- filter(PCAvalues_woody, provenance == "native") %>% 
                  select(PC1,PC2,PC3)
nm_non_native_wood <- filter(PCAvalues_woody, provenance == "non_native")%>% 
                      select(PC1,PC2,PC3)

#KS-TEST

ks.test(nm_non_native_wood$PC1, nm_native_wood$PC1)#D = 0.11531, p-value = 0.0002005

ks.test(nm_non_native_wood$PC2, nm_native_wood$PC2)#D = 0.056557, p-value = 0.218

#t-TEST

t.test(nm_non_native_wood$PC1, nm_native_wood$PC1)#t = -4.4872, df = 917, p-value = 8.134e-06

t.test(nm_non_native_wood$PC2, nm_native_wood$PC2)#t = -0.294, df = 966.97, p-value = 0.7688

#ad-TEST

ad.test(nm_non_native_wood$PC1, nm_native_wood$PC1)#AD= 10.386, T.AD= 12.342, P=6.2306e-06

ad.test(nm_non_native_wood$PC2, nm_native_wood$PC2)#AD= 0.72532, T.AD -0.36117, P=0.53781





####Combine plots for manuscript####
#density plots


#herbaceous species denstiy plots and box plots for Fig. 2 (these plots were combined outside of R to produce figure 2)

#density plots
figure2herb <- ggarrange( xplot_herb,yplot_herb,
                      labels = c("A)", "B)"),
                      ncol = 2, nrow = 1, common.legend = F)
figure2herb

#boxplots
figure2herb_box<-boxX_herb|boxY_herb
figure2herb_box



#woody denstiy plots and box plots for Fig. 2
#density plots
figure2wood <- ggarrange( xplot_wood,yplot_wood,
                       labels = c("C)", "D)"),
                       ncol = 2, nrow = 1, common.legend = F)
figure2wood

#boxplots
figure2wood_box<-boxX_wood|boxY_wood
figure2wood_box



