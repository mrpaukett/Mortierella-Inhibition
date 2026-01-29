#Low temperature growth advantage and inhibition of Fusarium solani by Mortierella spp.  
#Authors: Dr. Michelle Paukett, and Dr. Paul Esker
#Department of Plant Pathology & Environmental Microbiology, The Pennsylvania State University, University Park, PA 16802, U.S.A.

#####Citations----

#Enio G. Jelihovschi, Jose C. Faria, and Ivan Bezerra Allaman (2014). ScottKnott: A Package for Performing
#  the Scott-Knott Clustering Algorithm in R. Trends in Applied and Computational Mathematics 15(1), 3-17.
#  <https://tema.sbmac.org.br/tema/article/view/646/643>.

#Müller K (2020). _here: A Simpler Way to Find Your Files_. R package version 1.0.1.
#  <https://CRAN.R-project.org/package=here>.

#Pedersen T (2025). _patchwork: The Composer of Plots_. R package version 1.3.2.
#  <https://CRAN.R-project.org/package=patchwork>.

#R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for
#  Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.

#Signorell A (2024). _DescTools: Tools for Descriptive Statistics_. R package version 0.99.56.
#  <https://CRAN.R-project.org/package=DescTools>.

#Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J,
#  Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K,
#  Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_,
#  4(43), 1686. <https://doi.org/10.21105/joss.01686>.

#Wickham H, Hester J, Bryan J (2023). _readr: Read Rectangular Text Data_. R package version 2.1.4.
#  <https://CRAN.R-project.org/package=readr>.

#####Notes----

#See Mortierella Inhibition CSV File Descriptions text file for details on input files
#Additional formatting of figures was completed in Microsoft PowerPoint, with image resolutions adjusted using Krita 5.2.2.

#Sections include----

###Intro sections----

#####Load necessary packages----

library(here)                                                                    #File referencing/path-making for projects
library(tidyverse)                                                               #Data manipulation, includes dyplr, ggplot2, & readr
library(ScottKnott)                                                              #For mean clustering in dual culture experiment
library(DescTools)                                                               #To run Dunnett's test
library(patchwork)                                                               #To make grouped box plot in fungicide study

#####Load necessary files----

obs_data_inhibit <- read.csv(file = "InhibitionData.csv")                        #Contains dual culture experiment data 
obs_data_temp <- read.csv(file = "TempData.csv")                                 #Contains growth data at 5*C, 15*C, & 25*C (Control temp)
obs_data_chem <- read.csv(file = "FungicideData.csv")                            #Contains fungicide growth rate data at 0 (Control), 10x, 100x, & 1000x

file_list <- list.files(path = ".", pattern = "^M.*\\.csv$", 
                        full.names = TRUE)                                       #Get the individual isolate fungicide data files for analysis

print(file_list)                                                                 #View the selected files to confirm

#####Create data frames for analysis ----

tibbleIN <- as_tibble(obs_data_inhibit)                                          #Convert inhibition data to a tibble 
head(tibbleIN)                                                                   #Show results in console

tibbletemp <- as_tibble(obs_data_temp)                                           #Convert temp data to a tibble 
head(tibbletemp)                                                                 #Show results in console

tibblechem <- as_tibble(obs_data_chem)                                           #Convert fungicide data of all isolates to a tibble 
head(tibblechem)                                                                 #Show results in console

#####Summarize growth data----

pgrowthsum <- tibbleIN %>% group_by(catorg) %>% summarize(min = min(pgrow), 
         q1 = quantile(pgrow, 0.25), median = median(pgrow),mean = mean(pgrow),
         q3 = quantile(pgrow, 0.75), max = max(pgrow))                           #Summarize pathogen growth in dual culture & control plates                                                                             #Note: 0 represents single culture controls 
head(pgrowthsum)                                                                 #Show results in console

mgrowthsum <- tibbleIN %>% group_by(catorg) %>% summarize(min = min(mgrow), 
         q1 = quantile(mgrow, 0.25), median = median(mgrow),mean = mean(mgrow),
         q3 = quantile(mgrow, 0.75), max = max(mgrow))                           #Summarize Mortierella spp. growth in dual culture & control plates 
head(mgrowthsum)                                                                 #Show results in console

growthsum <- tibbletemp %>% group_by(catorg) %>% 
            summarize(min = min(measurement), q1 = quantile(measurement, 0.25),
            median = median(measurement), mean = mean(measurement),
            q3 = quantile(measurement, 0.75), max = max(measurement))            #Summarize growth by organism & temp
head(growthsum)                                                                  #Show results in console

growthchem <- tibblechem %>% group_by(Catorg) %>% summarize(min = min(GRate), 
           q1 = quantile(GRate, 0.25), median = median(GRate), 
           mean = mean(GRate), q3 = quantile(GRate, 0.75), max = max(GRate))     #Summarize growth by organism & growth on fungicide treated agar
head(growthchem)                                                                 #Show results in console

###Part 1: Dual culture inhibition experiment----

#####Calculate average growth & identify control growth----

pavegrowth <- tibbleIN %>% group_by(treat, org) %>% 
  mutate(Pave = mean(pgrow)) %>% distinct(Pave, .keep_all = TRUE)                #Calculate average pathogen growth 
head(pavegrowth)                                                                 #Show results in console

Pcon <- pavegrowth$Pave[21]                                                      #Save pathogen control average as a value for later use
                                                                                 #Important: Verify correct value before proceeding!! (currently pathogen control average found in row 21)
mavegrowth <- tibbleIN %>% group_by(treat, org) %>%
  mutate(Mave = mean(mgrow)) %>% distinct(Mave, .keep_all = TRUE)                #Get average Mortierella spp. growth
head(mavegrowth)                                                                 #Show results in console

FTgrowth <- pavegrowth[!(pavegrowth$treat=="C" | pavegrowth$org=="M"),]          #Create table with only Fusarium dual culture values
Mcon <- mavegrowth[!(mavegrowth$treat=="F" | mavegrowth$org=="F"),]              #Create table with only Mortierella spp. control growth values 
Mtreat <- mavegrowth[!(mavegrowth$treat=="C" | mavegrowth$org=="F"),]            #Create table with only Mortierella spp. dual culture growth values 
n <- Mcon %>% select(-c(rep, round, catorg, pgrow, mgrow))                       #Create table with only Mortierella spp. control values
h <- cbind(n, Mtreat[8], FTgrowth[8])                                            #Add to table dual culture growth values
Growthsummary <- h %>% group_by(org) %>% select(-c(treat)) %>% 
  mutate_if(is.numeric, ~round(., 2))                                            #Summarize averaged growth, rounding values
colnames(Growthsummary)[2] ="MCon"                                               #Rename Mortierella spp. control ave growth column  
colnames(Growthsummary)[3] ="MTreat"                                             #Rename Mortierella spp. dual-culture ave growth column                                       
colnames(Growthsummary)[4] ="PTreat"                                             #Rename Fusarium dual-culture ave growth column  
head(Growthsummary)                                                              #Show results in console
                                                                                 #Note: The Fusarium control (Pcon) is a single value, so not in the table!
#####Calculate relative inhibition (RI)----

RI <- Growthsummary %>% group_by(org) %>% mutate(MRatio = MTreat/MCon) %>%
  mutate(FRatio = PTreat/Pcon) %>% mutate(RI = FRatio/MRatio)                    #Calculate organism growth ratios, then relative inhibition
head(RI)                                                                         #Show results in console  
write.table(RI, file = "RI.csv", sep = ",", col.names = TRUE, row.names = TRUE)  #Export data frame as a CSV file

#####Calculate percent radial growth inhibition (PRGI) of Fusarium by Mortierella spp.----

PRGI_all <- tibbleIN %>% mutate(PRGI = ((Pcon-pgrow)/Pcon)*100)                  #Calculate PRGI
PRGI <- PRGI_all[!(PRGI_all$treat=="C" | PRGI_all$treat=="P5"),]                 #Remove controls
head(PRGI)                                                                       #Show results in console
write.table(PRGI, file = "PRGI.csv", sep = ",", 
            col.names = TRUE, row.names = TRUE)                                  #Export data frame as a CSV file

avePRGI <- PRGI %>% group_by(treat, org) %>% mutate(avePRGI = mean(PRGI)) %>%
  distinct(avePRGI, .keep_all = TRUE) %>%
  select(-c(rep, round, catorg, pgrow, PRGI, mgrow))                             #Calculate average PRGI by Mortierella spp. isolate
head(avePRGI)                                                                    #Show results in console
write.table(avePRGI, file = "avePRGI.csv", sep = ",",
            col.names = TRUE, row.names = TRUE)                                  #Export data frame as a CSV file                                                               #Show results in console       

#####Tukey tests & Scott-Knott for growth differences----

aovPRGI <- aov(formula=PRGI ~ org, data=PRGI)                                    #Run ANOVA to prep for Tukey
capture.output(summary(aovPRGI), file = "aovPRGI.txt")                           #Save summary of results in text file  
summary(aovPRGI)                                                                 #Show results in console

tky_PRGI <- TukeyHSD(aovPRGI)                                                    #Run Tukey to compare across treatments
capture.output(tky_PRGI, file = "tky_PRGI.txt")                                  #Save results in text file     

sk1 <- SK(PRGI ~org, data = PRGI, which ="org", sig.level=.05,)                  #Identify groups based on Scott-Knott mean clustering
print(sk1, digits = 2L, )                                                        #Show details on clusters in console
                                                                                 
#####Plot PRGI----

ggplot(PRGI, aes(x=org, y=PRGI, fill=org))+ 
  labs(x= ~italic("Mortierella")~ "Isolate",
  y= ~italic("F. solani")~ "Percent Radial Growth Inhibition", fill="Isolate", 
  title = ~italic("Fusarium")~ "Growth Inhibition") + 
  theme(legend.text = element_text(face = "italic")) +
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust=0.5, face = "bold",)) + geom_boxplot()   #Plot PRGI box plot, without Scott-Knott clustering

boxplot(sk1, args.legend=list(x='bottomleft'), adj=-0.5, 
  main= ~italic("F. solani")~ "PRGI by" ~italic("Mortierella")~ "Isolate",
  xlab= ~italic("Mortierella")~ "Isolate",
  ylab= ~italic("F. solani")~ "PRGI", margin=margin(30,0,0,0))                   #Scott-Knott clustered box plot
                                                                                 #Note: You will get an error message, but the graph is fine
###Part 2: Mortierella spp. temperature experiment----

#####Calculate percent diameter growth reduction (PDGR) by temperature----

aveTemps <- tibbletemp %>% group_by(temp, organism) %>%
 mutate(aveGrowth = mean(measurement)) %>% distinct(aveGrowth, .keep_all = TRUE) #Calculate average growth by organisms & temp, then remove some unnecessary columns
head(aveTemps)                                                                   #Show results in console

aveTempTreatments <- aveTemps[!(aveTemps$treat=="C" | aveTemps$temp=="25"),]     #Remove control temp (25) from the table; this is the baseline for comparison
head(aveTempTreatments)                                                          #Show results in console

aveTempControl <- aveTemps %>% subset(temp > 15) %>%
  rename(aveConGrowth = aveGrowth)                                               #Create table with control values, renaming column to prevent a duplicate when bind
head(aveTempControl)                                                             #Show results in console

repaveconsum <- aveTempControl[rep(seq_len(nrow(aveTempControl)), each = 2), ]   #Duplicate control values so column sizes match for binding/calculations
head(repaveconsum)                                                               #Show results in console

masterTemp <- cbind(aveTempTreatments, repaveconsum[8])                          #Bind the duplicated average control values to the treatment table currently without controls
head(masterTemp)                                                                 #Show results in console

TempPGIsummary <- masterTemp %>% mutate(PDGR = ((aveConGrowth-aveGrowth)/
  (aveConGrowth)*100)) %>% select(-c(treat, aveConGrowth, aveGrowth))            #Create table summarizing growth reduction by temp 
head(TempPGIsummary)                                                             #Show results in console

repavecontrols <- aveTempControl[rep(seq_len(nrow(aveTempControl)), each = 9), ] #Create 9 rows of each control value for the correct binding column size
head(repavecontrols)                                                             #Show results in console
masterTempfile <- cbind(tibbletemp, repavecontrols[8])                           #Bind the controls to the original data file (in tibble form)
head(masterTempfile)                                                             #Show results in console

TempPGIanalysis <- masterTempfile %>%
  mutate(PDGR = (((aveConGrowth-measurement)/
  (aveConGrowth)*100)*-1)) %>% select(-c(treat, aveConGrowth))                   #Calculate PDGR to be used in further analysis, making negative for plot
head(TempPGIanalysis)                                                            #Show results in console

#####Calculate growth rate (GR) per degree change (mm/1*C increase)----

ratedeg <- masterTempfile %>% mutate(rate = ((measurement)/10)) %>%
  select(-c(treat, aveConGrowth))                                                #Calculate the growth rate (10*C change between each temp)
head(ratedeg)                                                                    #Show results in console

averatedeg <- ratedeg  %>% group_by(temp, organism) %>%
  mutate(averate = mean(rate))                                                   #Calculate the average rate by org and temperature
head(averatedeg)                                                                 #Show results in console

addj <- averatedeg[!(averatedeg$replicate=="R2" | averatedeg$replicate=="R3"),]  #Remove duplicate rows
SumTempGR <- addj %>% select(-c(replicate, round, rate, measurement, catorg)) %>%            
mutate_if(is.numeric, ~round(., 3))                                              #Remove unnecessary columns & round values
head(SumTempGR)                                                                  #Show results in console

x <- SumTempGR[!(SumTempGR$temp=="15"| SumTempGR$temp=="25"),]                   #Create new column for growth rate change between 15*C & 25*C
y <- SumTempGR[!(SumTempGR$temp=="5" | SumTempGR$temp=="25"),]                   #Create new column for growth rate change between 5*C & 25*C
z <- SumTempGR[!(SumTempGR$temp=="5" | SumTempGR$temp=="15"),]                   #Create new column for growth rate change between 5*C & 15*C
v <- cbind(x, y[3], z[3])                                                        #Bind new columns into new table
SummaryRates <- v %>% group_by(organism) %>% select(-c(temp))                    #Organize table by organism
colnames(SummaryRates)[2] ="five"                                                #Rename column
colnames(SummaryRates)[3] ="fifteen"                                             #Rename column
colnames(SummaryRates)[4] ="twentyfive"                                          #Rename column
SummaryRates <- SummaryRates[!duplicated(SummaryRates), ]                        #Remove duplicate rows
head(SummaryRates)                                                               #Show results in console

Ratetable <- SummaryRates %>% group_by(organism) %>%
  mutate(GRchange25 = twentyfive-fifteen) %>% mutate(GRchange5 = fifteen-five)   #Calculate growth rate change between temperatures
head(Ratetable)                                                                  #Show results in console

#####Tukey tests for differences in growth----

aovtemp <- aov(formula=PDGR ~ catorg, data=TempPGIanalysis)                      #Run ANOVA on PDGR by treatment 
capture.output(summary(aovtemp), file = "aovTemp.txt")                           #Save summary of results in text file  
print(aovtemp)                                                                   #Show results in console

tky_temp <- TukeyHSD(aovtemp)                                                    #Run Tukey to compare across treatments
capture.output(tky_temp, file = "tky_pdgrtemp.txt")                              #Save results in text file  

#####Plot temp growth----

ggplot(TempPGIanalysis, aes(x=factor(catorg,level=c('25A', '25G', '25X', 
  '25E', '25F', '15A', '15G', '15X', '15E', '15F', '5A', '5G', '5X', '5E', 
  '5F')), y=PDGR, fill=organism)) + 
  labs(x= ~italic("Mortierella & Fusarium")~ "Isolate By Temperature", 
       y= "Percent Diameter Growth Inhibition", fill="Isolate", 
       title = "Temperature Growth Inhibition") + 
  theme(legend.text = element_text(face = "italic")) +  
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust=0.5, face = "bold",)) + geom_boxplot()   #Make a box plot of PDGR (standardized to 25*C control)

ggplot(data=aveTemps, aes(x=temp, y=aveGrowth, group=organism)) +
  geom_line(aes(linetype=organism,  color = organism,))+ 
  geom_point(aes(shape=organism)) +
  labs(x= "Temperature", y= "Average Isolate Growth (mm)", 
       title = ~italic("Mortierella & Fusarium")~ "Growth By Temperature") + 
  theme(legend.text = element_text(face = "italic")) + 
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust=0.5, face = "bold",))                    #Make a line graph of average growth at each temperature

###Part 3: Fungicide sensitivity experiment----

sig_stars <- function(p) {cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", ""), right = TRUE)}   #Save significance star values for results

perform_dunnett_from_file <- function(file_list) {
  data <- read.csv(file_list)                                                    #Read the data in each individual fungicide isolate file
  
  data$Catorg <- as.factor(data$Catorg)                                          #Set group as factor
  data$Catorg <- relevel(data$Catorg, ref = "Control")                           #Set control as the reference
  
  aov_result <- aov(GRate ~ Catorg, data = data)                                 #Run ANOVA
  aov_table <- as.data.frame(summary(aov_result)[[1]])                           #Make ANOVA data frame
  aov_table$file <- basename(file_list)                                          #Name data frame based on input file name
  aov_table$signif <- sig_stars(aov_table$`Pr(>F)`)                              #Add significance stars to ANOVA results in output file
  
  dunnett_result <- DunnettTest(x = data$GRate, g = data$Catorg)                 #Perform Dunnett's test
  
  result_df <- as.data.frame(dunnett_result$Control)                             #Create a data frame with the Dunnett's test results
  result_df$file <- basename(file_list)                                          #Create file named by isolate for results
  result_df$signif <- sig_stars(result_df$pval)                                  #Add significance stars to results 
  
  return(list(anova   = aov_table, dunnett = result_df))                         #Set output want to save
}

all_results <- purrr::map(file_list, perform_dunnett_from_file)                  #Apply the function to all files & combine results
anova_all <- dplyr::bind_rows(purrr::map(all_results, "anova"))                  #Combine all ANOVA results in the same file
dunnett_all <- dplyr::bind_rows(purrr::map(all_results, "dunnett"))              #Combine all Dunnett's test results in the same file

print(all_results)                                                               #View the combined results
capture.output(all_results, file = "allDunnettResults.txt")                      #Save results in text file  

plots <- list()                                                                  #Create a list for plots made in loop

for(file_path in file_list) {
  
  data <- read.csv(file_path)                                                    #Read the data in each individual fungicide isolate file
  
  data$Catorg <- factor(trimws(data$Catorg))                                     #Set group as factor
  data$Catorg <- relevel(data$Catorg, ref = "Control")                           #Set control as the reference
  
  p <- ggplot(data, aes(x = Catorg, y = GRate, fill = Treat)) + geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) + theme_minimal(base_size = 14) +
    theme(plot.background = element_rect(fill = "white", color = NA),  
          panel.background = element_rect(fill = "white", color = NA),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +                                            #Note: Remove if want plot legends
    labs(x = "Treatment", y = "Growth Rate", fill = "Treatment")                 #Generate box plot
                                                                                 #Note: If want plot titles, add this in labs() title = paste(tools::file_path_sans_ext(basename(file_path))), 
  out_file <- paste0(tools::file_path_sans_ext(basename(file_path)),
                     "_boxplot.png")                                             #Create PNG file named by isolate for box plot

  ggsave(out_file, plot = p, width = 6, height = 4, bg = "white")                #Save individual isolate box plot PNG files
  plots[[basename(file_path)]] <- p                                              #Add box plot to the list to combine
}

plots                                                                            #Show individual isolate box plots
combined_plot <- wrap_plots(plots, ncol = 2)                                     #Combine all plots in one figure with 2 columns
combined_plot                                                                    #Show plots combined
ggsave("all_isolates_boxplots.png", combined_plot, width = 8, 
       height = 12, bg = "white")                                                #Save combined box plot as a PNG file
