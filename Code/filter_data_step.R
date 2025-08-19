
#data
d_pca = read_csv("Data/d_pca.csv") #pca scores for native and non-native urban plant species

####All Taxa####
# scale variables
d_pca$ave_tmean = scale(d_pca$ave_tmean)[,1]
d_pca$ave_precip = scale(d_pca$ave_precip)[,1]


# Ctemp
ggplot(d_pca, aes(x = ave_tmean)) +
    geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
    geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Histogram with ±3 Z-score Highlighted",
         x = "Z-score",
         y = "Frequency") +
    theme_minimal()

d_pca_tfilt <- d_pca %>%
            filter(abs(ave_tmean) <= 3)
           

#precip
ggplot(d_pca, aes(x = ave_precip)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Histogram with ±3 Z-score Highlighted",
       x = "Z-score",
       y = "Frequency") +
  theme_minimal()

d_pca_pfilt <- d_pca %>%
              filter(abs(ave_precip) <= 3)


#All Taxa final filtered data set


d_pca_filtered <- d_pca %>%
                 filter(abs(ave_tmean) <= 3 & abs(ave_precip) <= 3)

#herb
herb<-d_pca%>%
      filter(growth_form=="herb") 




# scale variables
herb$ave_tmean = scale(herb$ave_tmean)[,1]
herb$ave_precip = scale(herb$ave_precip)[,1]


# temp
ggplot(herb, aes(x = ave_tmean)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

herb_tfilt <- herb %>%
               filter(abs(ave_tmean) <= 3)


#precip
ggplot(herb, aes(x = ave_precip)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

herb_pfilt<- herb %>%
             filter(abs(ave_precip) <= 3)

#Final filtered data Herbs

herb_filtered <- herb %>%
                 filter(abs(ave_tmean) <= 3 & abs(ave_precip) <= 3)

#woody
woody<-d_pca%>%
       filter(growth_form=="woody") 


# scale variables
woody$ave_tmean = scale(woody$ave_tmean)[,1]
woody$ave_precip = scale(woody$ave_precip)[,1]


# temp
ggplot(woody, aes(x = ave_tmean)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()

woody_tfilt <- woody %>%
               filter(abs(ave_tmean) <= 3)


#precip
ggplot(woody, aes(x = ave_precip)) +
      geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
      geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = -3, color = "red", linetype = "dashed", size = 1) +
      labs(title = "Histogram with ±3 Z-score Highlighted",
           x = "Z-score",
           y = "Frequency") +
      theme_minimal()


woody_pfilt<- woody %>%
             filter(abs(ave_precip) <= 3)



#Final filtered data Herbs

woody_filtered <- woody %>%
                 filter(abs(ave_tmean) <= 3 & abs(ave_precip) <= 3)



