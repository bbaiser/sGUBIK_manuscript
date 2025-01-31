#PCA with marginal distributions 


#PCA 
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

#Data from Diaz et al. 2022 and GBIF/Gubic for urban for species that have all 6 diaz traits. taxonomically matched to WCVP 
Diaz_Final <- read.csv("Data/Diaz_Gubic_final_tot2.csv",row=1) %>%
              rename("Leaf area"="Leaf_area","Leaf N"="Leaf_N",
                     "Seed mass"="Seed_mass","Stem density"="Stem_density" )%>%
              mutate(provenance_glonaf = replace(provenance_glonaf, str_detect(provenance_glonaf, "non_native"), "non-native"))


#counts of native (n=1564) and non-native(n=1213)
prov_counts<-Diaz_Final %>%
            group_by(provenance_glonaf) %>%
            tally 
prov_counts

#### PCA for all species####
#scale data
z_Diaz_Final<-scale(Diaz_Final[, c("Leaf area","LMA","Leaf N","Seed mass","Stem density","Height")])

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
                      panel.border=element_rect(fill=NA,colour="grey50"),
                      legend.title=element_blank())

p_my_theme2 <- theme( legend.text = element_text(size=9), axis.title= element_text(size=9),
                      axis.text.x= element_text(colour= "black", size=9), axis.text.y= element_text(colour= "black", size=9))

All_Sp_PCA <- ggplot(PCAvalues) +   
              geom_point(size = .7, alpha=0.5, 
                         aes(x = PC1, y = PC2, group = provenance, colour = provenance, size = provenance),  # label = T
                         show.legend = FALSE) +
              coord_fixed() + 
              scale_colour_manual(values=c("blue1","orange")) +
              geom_segment(data = PCAloadings, size = 1,
                           aes(x = 0, xend = PC1*5.5, y = 0, yend = PC2*5.5),
                           arrow = arrow(length = unit(.3, "cm")),colour = "black")   +
              geom_text(data = PCAloadings, aes(x = PC1*6, y = PC2*6, label = Variables), size = 5,
                        hjust=c(0, 0, 0, 0, 0, 0) , vjust=c(0, 0, 0, 0, 0, 1))    + 
              xlab("PC1 (47%) ") + ylab("PC2 (26%)")  +
              # ggtitle("All Species")+
              xlim(c(-6,6))+
              ylim(c(-5.7,5.7))+
              p_my_theme



#Add denisty plots to the PCA

dens_PCA<-ggplot(data=PCAvalues) +
          stat_density_2d(aes(x=PC1, y=PC2, color=provenance, fill=provenance), geom="polygon",contour = TRUE, bins=5, alpha=.35,show.legend = T) +
          geom_point(aes(x=PC1, y=PC2, color=provenance), size =.011, show.legend = FALSE) +
          scale_color_manual( values = c("blue","orange"))+
          scale_fill_manual( values = c("blue","orange"))+
          coord_fixed() + 
          geom_segment(data = PCAloadings, size = 1,
                     aes(x = 0, xend = PC1*5.5, y = 0, yend = PC2*5.5),
                     arrow = arrow(length = unit(.3, "cm")),colour = "black")   +
          geom_text(data = PCAloadings, aes(x = PC1*6, y = PC2*6, label = Variables), size = 5,
                  hjust=c(0, 0, 0, 0, 0, 0) , vjust=c(0, 0, 0, 0, 0, 1))    + 
          xlab("PC1 (47%)") + ylab("PC2 (26%)")  +
          # ggtitle("All Species")+
          xlim(c(-6,6))+
          ylim(c(-5.7,5.7))+
          theme(legend.position = c(.2,.85))+
          p_my_theme


# Add density curves for PCA scores to y and x axis
xdens <- axis_canvas(All_Sp_PCA, axis = "x") + 
        geom_density(data = PCAvalues, aes(x = PC1, fill = provenance, colour = provenance), alpha = 0.3)+
        scale_color_manual( values = c("blue","orange"))+
        scale_fill_manual( values = c("blue","orange"))
      

ydens <-axis_canvas(All_Sp_PCA, axis = "y", coord_flip = TRUE) + 
        geom_density(data = PCAvalues, aes(x = PC2, fill = provenance, colour = provenance), alpha = 0.3) +
        scale_color_manual( values = c("blue","orange"))+
        scale_fill_manual( values = c("blue","orange"))+
        coord_flip()

#Make figure 1
Fig_1<-dens_PCA %>%
      insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
      insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
      ggdraw()

Fig_1