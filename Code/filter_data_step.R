#Filter climate outliers >3 sd 

#load adata
#Species PCA Sscores and climate data matrix
d_pca_raw <- read_csv("Data/d_pca.csv") #pca scores for native and non-native urban plant species

####All Taxa ####
d_pca<-d_pca_raw
# scale variables temp and precip
d_pca$ave_tmean = scale(d_pca$ave_tmean)[,1]
d_pca$ave_precip = scale(d_pca$ave_precip)[,1]


# Plot temp- appears that there are no temp outliers
ggplot(d_pca, aes(x = ave_tmean)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

#remove outliers
d_pca_tfilt <- d_pca %>%
            filter(abs(ave_tmean) <= 3)

#confirmed, no outliers removed
dim (d_pca)
dim(d_pca_tfilt)
           

#Plot for precipitation- outliers present
ggplot(d_pca, aes(x = ave_precip)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

#remove outliers
d_pca_pfilt <- d_pca %>%
              filter(abs(ave_precip) <= 3)



#confirmed, 43 outliers removed (1.5% of species removed)
dim (d_pca)
dim(d_pca_pfilt)

#Remove both precip and temp outliers
d_pca_filtered <- d_pca %>%
                 filter(abs(ave_tmean) <= 3 & abs(ave_precip) <= 3)

####Herb Taxa####
#Filter Herb species
herb<-d_pca_raw%>%
      filter(growth_form=="herb") 

# scale variables
herb$ave_tmean = scale(herb$ave_tmean)[,1]
herb$ave_precip = scale(herb$ave_precip)[,1]


# Plot temp- outliers observed
ggplot(herb, aes(x = ave_tmean)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

#Filter temp outliers
herb_tfilt <- herb %>%
               filter(abs(ave_tmean) <= 3)

#confirmed - 29 ouliers removed (2.5% of herb species)

dim(herb)
dim(herb_tfilt)

#Plot precip- otliers observed
ggplot(herb, aes(x = ave_precip)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

#remove outliers
herb_pfilt<- herb %>%
             filter(abs(ave_precip) <= 3)

#confirmed- 22 species removed (~2% of species)
dim(herb)
dim(herb_pfilt)

#Final filtered data Herbs removing temp and precip outliers

herb_filtered <- herb %>%
                 filter(abs(ave_tmean) <= 3 & abs(ave_precip) <= 3)


#total outliers removed =37 (~3% of species)
dim(herb)
dim(herb_filtered)

####Woody Taxa####
#Filter for woody species
woody<-d_pca_raw%>%
       filter(growth_form=="woody") 


# scale variables
woody$ave_tmean = scale(woody$ave_tmean)[,1]
woody$ave_precip = scale(woody$ave_precip)[,1]


# Plot temp-no oultiers opbserved
ggplot(woody, aes(x = ave_tmean)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

#remove ouliers
woody_tfilt <- woody %>%
               filter(abs(ave_tmean) <= 3)

#confirmed- no temp outliers
dim(woody)
dim(woody_tfilt)

#Plot precip-outliers observed
ggplot(woody, aes(x = ave_precip)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

#remove precip outliers
woody_pfilt<- woody %>%
             filter(abs(ave_precip) <= 3)

#confirmed 18 species removed (~1% of the data)
dim(woody)
dim(woody_pfilt)

#Final filtered data Herbs
woody_filtered <- woody %>%
                 filter(abs(ave_tmean) <= 3 & abs(ave_precip) <= 3)



