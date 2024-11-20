#Analysis script for wave 2 COVCOG data, written by Josefina Weinerova (2024)
#Bayesian analysis is done in JASP

renv::restore() #restores packages from the renv file. For help look here: https://rstudio.github.io/renv/articles/renv.html

# renv::install("naniar")
# renv::install("brms")
# renv::install("MASS")
# renv::install("GGally")
# renv::install("ggcorrplot")
# renv::install("corrtable")
# renv::install("gridExtra")
# renv::install("gtable")
# renv::install("cowplot")
#renv::install("zoo")
library(tidyverse)
library(lmerTest)
library(tidyselect)
library(pwr)
library(tidyverse)
library(easystats)
library(car)
library(rstatix)
library(data.table)
library(effectsize)
library(brms)
library(MASS)
library(GGally)
library(ggcorrplot)
library(corrtable)
library(Hmisc)
library(gridExtra)
library(gtable)
library(naniar)
library(cowplot)
library(naniar)

df<- read_csv("COVCOG_data_140524_2.csv") #loading data
#View(df) #checking the csv was corectly parsed
excl <- df %>% dplyr::select(starts_with("Q2.04")) #subsetting columns relevant to exclusion criteria

#data cleaning:
#excluding everyone with any mentioned psychiatric or neurodevelopmental condition affecting cognitive function or major medical procedure within the last 6 months
psych_dis<-unique(excl$Q2.04.23OtherMentalPsychiatric_Specify)
psych_dis<-psych_dis[c(2,3,5,6,9:11,14,17,19,20)]
neurod_dis <-unique(excl$Q2.04.25Neurodevelopmental_Specify)
neurod_dis <- neurod_dis[c(2:10,12:17)]
noncov_dis<-unique(excl$Q2.04.1NonCOVIDMedical_Specify)
noncov_dis <- noncov_dis[9:11]
df_cl <- df %>% subset(!(Q2.04.23OtherMentalPsychiatric_Specify %in% psych_dis)) %>% subset(!(Q2.04.25Neurodevelopmental_Specify %in% neurod_dis)) %>% subset(!(Q2.04.1NonCOVIDMedical_Specify %in% noncov_dis))
nrow(df_cl) #getting final number of rows = 296

table(df_cl$COVID_Grouping) #group numbers
#setting variable types
df_cl$COVID_Grouping<-as.factor(df_cl$COVID_Grouping) #Covid group as factor
df_cl$Q1.01Age_Text<-as.numeric(df_cl$Q1.01Age_Text) #Age as numeric
df_cl.cols <- c("Q1.02Sex_Quantised", "Q1.03Country_Quantised","Q1.08.Education_Quantised")
df_cl[df_cl.cols]<-lapply(df_cl[df_cl.cols],factor) #factors: sex, country, edu, country
#replacing all 9999 values with NA
df_cl<-df_cl %>% replace_with_na_all(condition = ~.x == 9999)
df_cl<-df_cl %>% replace_with_na_all(condition = ~.x == "9999")
df_cl<-df_cl %>% replace_with_na_all(condition = ~.x == "#N/A")

#demographics for sample section
table(df_cl$COVID_Grouping,df_cl$Q1.02Sex_Quantised)
table(df_cl$COVID_Grouping)
summary(df_cl$Q1.01Age_Text)
tapply(df_cl$Q1.01Age_Text,df_cl$COVID_Grouping,summary)
tapply(df_cl$Q1.01Age_Text,df_cl$COVID_Grouping,sd)
table(df_cl$Q1.03Country_Specify)
df_cl$Country_binary <- if_else(df_cl$Q1.03Country_Quantised==1,1,0) 
sd(df_cl$Q1.01Age_Text)

#colourblind friendly palette for plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#ANALYSIS 1 - replication of Guo et al. - looking at word list item and pictorial associative memory tasks
#ANCOVA replication
##verbal item memory task percent correct
options(scipen=100, digits=4) #just sets number of decimal points shown
aovreport <- aov(VMemoryRecog_CorrectPercent ~ COVID_Grouping + Q1.01Age_Text + Q1.02Sex_Quantised + Q1.03Country_Quantised + Q1.08.Education_Quantised, data = df_cl)
Anova(aovreport, type="III")
effectsize::eta_squared(aovreport,partial=TRUE)
#nonverbal associative memory task percent correct
aovreport <- aov(NVMemoryAsso_CorrectPercent ~ COVID_Grouping + Q1.01Age_Text + Q1.02Sex_Quantised + Q1.03Country_Quantised + Q1.08.Education_Quantised, data = df_cl)
Anova(aovreport, type="III")
effectsize::eta_squared(aovreport,partial=TRUE)
#verbal item memory RTs
aovreport <- aov(VMemoryRecog_RTOverall ~ COVID_Grouping + Q1.01Age_Text + Q1.02Sex_Quantised + Q1.03Country_Quantised + Q1.08.Education_Quantised, data = df_cl)
Anova(aovreport, type="III")
effectsize::eta_squared(aovreport,partial=TRUE)
#nonverbal associative memory RT overall
aovreport <- aov(NVMemoryAsso_RTOverall ~ COVID_Grouping + Q1.01Age_Text + Q1.02Sex_Quantised + Q1.03Country_Quantised + Q1.08.Education_Quantised, data = df_cl)
Anova(aovreport, type="III")
effectsize::eta_squared(aovreport,partial=TRUE)

#descriptives for the variables
df_cl%>%group_by(COVID_Grouping)%>% summarise(mean = mean(VMemoryRecog_CorrectPercent,na.rm=TRUE),sd=sd(VMemoryRecog_CorrectPercent,na.rm=TRUE))
df_cl%>%group_by(COVID_Grouping)%>% summarise(mean = mean(NVMemoryAsso_CorrectPercent,na.rm=TRUE),sd=sd(NVMemoryAsso_CorrectPercent,na.rm=TRUE))
df_cl%>%group_by(COVID_Grouping)%>% summarise(mean = mean(VMemoryRecog_RTOverall,na.rm=TRUE),sd=sd(VMemoryRecog_RTOverall,na.rm=TRUE))
df_cl%>%group_by(COVID_Grouping)%>% summarise(mean = mean(NVMemoryAsso_RTOverall,na.rm=TRUE),sd=sd(NVMemoryAsso_RTOverall,na.rm=TRUE))

#JASP does not allow for files with more than 1000 columns to be uploaded so creating a file with less columns
df_JASP<- df_cl%>%dplyr::select(Participant_Private_ID, Completion_Date, COVID_Grouping, Q1.01Age_Text,Q1.02Sex_Quantised,Q1.03Country_Quantised,Q1.08.Education_Quantised,
                         VMemoryAsso_CorrectPercent,VMemoryAsso_RTOverall,VMemoryRecog_CorrectPercent,VMemoryRecog_RTOverall,NVMemoryRecog_CorrectPercent,NVMemoryRecog_RTOverall,NVMemoryAsso_CorrectPercent,NVMemoryAsso_RTOverall)
#file for ANALYSIS 2
write.csv(df_JASP,"COVCOG_wave2_JASP.csv")


#pivoting data for Analysis 3 and 4
#scaling and pivoting accuracy
df_JASP$VMemoryRecog_CorrectPercent<-scale(df_JASP$VMemoryRecog_CorrectPercent, center = TRUE, scale=TRUE)
df_JASP$NVMemoryAsso_CorrectPercent<-scale(df_JASP$NVMemoryAsso_CorrectPercent,center=TRUE,scale=TRUE)
df_JASP1 <- pivot_longer(df_JASP, cols=c(VMemoryRecog_CorrectPercent,NVMemoryAsso_CorrectPercent),
                             names_to = "MemTest",
                             values_to = "MemAccuracy",
                             values_drop_na = TRUE)

#Plot for ANALYSIS 3 
ggplot(df_JASP1, aes(y=MemAccuracy, x=MemTest,color=COVID_Grouping)) + geom_boxplot(notch=TRUE, lwd=1.5,fatten=1) + scale_color_manual(values=cbPalette) + theme_cowplot(12) + labs(title="Scaled memory accuracy for verbal item memory \n and nonverbal associative memory task", y="Z scores (Accuracy)", x="Memory task", color="Covid group") +
  theme(plot.title = element_blank(),legend.position="none", axis.title=element_text(size=18, face="bold"), axis.text.x=element_text(size=17),axis.text.y=element_text(size=17)) + scale_x_discrete(labels= c("Nonverbal Associative", "Verbal Item")) + scale_y_continuous(limits=c(-3,3),breaks = c(-3,-2,-1,0,1,2,3))

#CSV for ANALYSIS 3 
write.csv(df_JASP1,"COVCOG_wave2_piv.csv")

#scaling and pivoting for RTs
df_JASP$VMemoryRecog_RTOverall <- scale(df_JASP$VMemoryRecog_RTOverall,center = TRUE, scale=TRUE)
df_JASP$NVMemoryAsso_RTOverall <- scale(df_JASP$NVMemoryAsso_RTOverall,center = TRUE, scale = TRUE)
df_JASP2 <- pivot_longer(df_JASP, cols=c(VMemoryRecog_RTOverall,NVMemoryAsso_RTOverall),
                         names_to = "MemTest",
                         values_to = "MemRT",
                         values_drop_na = TRUE)
#CSV for ANALYSIS 4
write.csv(df_JASP2, "COVCOG_wave2_pivRT.csv")

#pivoting for Accu and RTs for all memory  tasks for ANALYSIS 5
#accuracy
df_JASP$VMemoryAsso_CorrectPercent<-scale(df_JASP$VMemoryAsso_CorrectPercent, center = TRUE, scale=TRUE)
df_JASP$NVMemoryRecog_CorrectPercent<-scale(df_JASP$NVMemoryRecog_CorrectPercent,center=TRUE,scale=TRUE)
df_JASP3 <- pivot_longer(df_JASP, cols=c(VMemoryRecog_CorrectPercent,NVMemoryAsso_CorrectPercent,VMemoryAsso_CorrectPercent,NVMemoryRecog_CorrectPercent),
                         names_to = "MemTest",
                         values_to = "MemAccuracy",
                         values_drop_na = TRUE)
#creating columns with verbal vs nonverbal and associative vs item
df_JASP3<- df_JASP3 %>% mutate(MemType = if_else(df_JASP3$MemTest==c("VMemoryRecog_CorrectPercent","NVMemoryRecog_CorrectPercent"),"Recog","Asso"))
df_JASP3$MemStim <- if_else(startsWith(df_JASP3$MemTest,"N"),"Nonverbal","Verbal")
#CSV file for ANALYSIS 5 for accuracy 
write.csv(df_JASP3, "COVCOG_wave2_accuPivo.csv")

#RTs for ANALYSIS 5
df_JASP$VMemoryAsso_RTOverall<-scale(df_JASP$VMemoryAsso_RTOverall, center = TRUE, scale=TRUE)
df_JASP$NVMemoryRecog_RTOverall<-scale(df_JASP$NVMemoryRecog_RTOverall,center=TRUE,scale=TRUE)
df_JASP4 <- pivot_longer(df_JASP, cols=c(VMemoryRecog_RTOverall,NVMemoryAsso_RTOverall,VMemoryAsso_RTOverall,NVMemoryRecog_RTOverall),
                         names_to = "MemTest",
                         values_to = "MemRT",
                         values_drop_na = TRUE)

df_JASP4<- df_JASP4 %>% mutate(MemType = if_else(df_JASP4$MemTest==c("VMemoryRecog_RTOverall","NVMemoryRecog_RTOverall"),"Recog","Asso"))
df_JASP4$MemStim <- if_else(startsWith(df_JASP4$MemTest,"N"),"Nonverbal","Verbal")
#CSV for ANALYSIS 5
write.csv(df_JASP4, "COVCOG_wave2_RTPivo.csv")

#plot for accuracy on all 4 tasks for ANALYSIS 5 section
ggplot(df_JASP3,aes(y=MemAccuracy, x= MemType, colour=COVID_Grouping)) + geom_boxplot(notch=TRUE, lwd=1.5,fatten=1) + facet_grid(~MemStim) + theme_cowplot(12) +
  labs(title="The scaled accuracy scores for Memory tasks", y="Z scores (Accuracy)",x="Memory task", color="Covid group") + scale_color_manual(values=cbPalette,labels=c("No-Covid", "Covid")) +
  theme(plot.title = element_blank(), legend.position = "none", axis.title=element_text(size=18,face="bold"), axis.text.x=element_text(size=17),axis.text.y=element_text(size=17), strip.text = element_text(
    size = 18)) + scale_x_discrete(labels= c("Associative", "Item")) + scale_y_continuous(limits=c(-5,5),breaks=c(-5,-2.5,0,2.5,5))

#creating csv for Additional ANALYSES SECTION (registered)
df_cl$Syn_ErrorPerc <- as.numeric(df_cl$Syn_ErrorPerc)
df_cl$Syn_TimeOutsPerc<-as.numeric(df_cl$Syn_TimeOutsPerc)
df_cl$SynCorrPerc <- 100-(df_cl$Syn_ErrorPerc + df_cl$Syn_TimeOutsPerc)
df_JASP_all<- df_cl%>%dplyr::select(Participant_Private_ID, Completion_Date, COVID_Grouping, Q1.01Age_Text,Q1.02Sex_Quantised,Q1.03Country_Quantised,Q1.08.Education_Quantised,
                         VMemoryAsso_CorrectPercent,VMemoryAsso_RTOverall,VMemoryRecog_CorrectPercent,VMemoryRecog_RTOverall,NVMemoryRecog_CorrectPercent,NVMemoryRecog_RTOverall,NVMemoryAsso_CorrectPercent,NVMemoryAsso_RTOverall,
                         DS_MaxSpan,DS_Score, CF_CorrectExclRep, SynCorrPerc,Syn_RToverall,WCST_Correct,WCST_RT_correct)


cols_all <- c("VMemoryAsso_CorrectPercent","VMemoryAsso_RTOverall","VMemoryRecog_CorrectPercent","VMemoryRecog_RTOverall","NVMemoryRecog_CorrectPercent","NVMemoryRecog_RTOverall","NVMemoryAsso_CorrectPercent","NVMemoryAsso_RTOverall",
              "DS_MaxSpan","DS_Score", "CF_CorrectExclRep", "SynCorrPerc","Syn_RToverall","WCST_Correct","WCST_RT_correct")
df_JASP_all[cols_all]<-sapply(df_JASP_all[cols_all],as.numeric)
sapply(df_JASP_all,"class")
#CSV file for ADDITIONAL ANALYSES section (registered)
write.csv(df_JASP_all,"COVCOG_wave2_allCog.csv")

#preparing the csv for the ANALYSIS 6 and 7 - the effect of vaccination
df_vacc <- df_cl %>% filter(COVID_Grouping==1) %>% filter(Q3.06COVIDTimes_Quantised==1 & Q3.05.1COVICMoreThanOnce_Quantised %in% c(3,4)) #filtering in only people who have had covid once during two questions -how many times do you think you have had covid? and do you think you have had covid more than once? 
#assigning groups - vaccinated before or after first infection
df_vacc$VACC_Groupin <- NA
df_vacc$VACC_Groupin[df_vacc$Q3.01.4VaccineType1_Quantised == 4] <- 1 #if Johnson & Johnson vaccine assign to vacc group (I checked in the file beforehand that all of those with Johnson Johnson vacc had it more than 3 weeks before infection)
df_vacc$VACC_Groupin[df_vacc$Q3.05.2.1BeforeVaccine == 1]<-0 #had covid before vaccine

df_vacc$VACC_Groupin[df_vacc$Q3.01.1WhetherHadVaccine_Quantised == 2] <- 0
df_vacc$VACC_Groupin[df_vacc$Q3.05.2.5After2More3 == 1 | df_vacc$Q3.05.2.6After3Less3 == 1 | df_vacc$Q3.05.2.7After3More3 ==1] <- 1
df_vacc$VACC_Groupin[df_vacc$Q3.05.3.1BeforeVaccine ==1]<-0 #think they had covid before vaccine - second variable if they did not have test confirmed covid - but people who are truly unsure about whether they had covid have been excluded at this point for previous analysis
df_vacc$VACC_Groupin[df_vacc$Q3.05.2.5After2More3==1 & df_vacc$Q3.05.3.1BeforeVaccine ==1] <- NA #one person answered both that they had Covid before being vacc and that they had after 2 doses, so I excluded them

#selecting only necessary columns so that the file is not too big for JASP
df_vacc_all<- df_vacc%>%dplyr::select(VACC_Groupin, Participant_Private_ID, Completion_Date, COVID_Grouping, Q1.01Age_Text,Q1.02Sex_Quantised,Q1.03Country_Quantised,Q1.08.Education_Quantised,
                             VMemoryAsso_CorrectPercent,VMemoryAsso_RTOverall,VMemoryRecog_CorrectPercent,VMemoryRecog_RTOverall,NVMemoryRecog_CorrectPercent,NVMemoryRecog_RTOverall,NVMemoryAsso_CorrectPercent,NVMemoryAsso_RTOverall,
                             DS_MaxSpan,DS_Score, CF_CorrectExclRep, SynCorrPerc,Syn_RToverall,WCST_Correct,WCST_RT_correct)

#CSV for ANALYSIS 6: effect of vaccination
write.csv(df_vacc_all,"COVCOG_wave2_vacc.csv")



#UNREGISTERED ANALYSES
#RT analysis without the outliers (to double check conclusions for Analysis 2)
#for verbal item memory RT removing outliers
summary(df_cl$VMemoryRecog_RTOverall) # finding the IQR
Q <- quantile(df_cl$VMemoryRecog_RTOverall, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(df_cl$VMemoryRecog_RTOverall)
up <-  Q[2]+3*iqr # Upper Range  
low<- Q[1]-3*iqr # Lower Range
df_RT <-df_cl
df_RT$VMemoryRecog_RTOverall<- ifelse(df_RT$VMemoryRecog_RTOverall<low | df_RT$VMemoryRecog_RTOverall>up,NA,df_RT$VMemoryRecog_RTOverall)
sum(is.na(df_RT$VMemoryRecog_RTOverall))
aovreport <- aov(VMemoryRecog_RTOverall ~ COVID_Grouping + Q1.01Age_Text + Q1.02Sex_Quantised + Q1.03Country_Quantised + Q1.08.Education_Quantised, data = df_RT)
Anova(aovreport, type="III")
effectsize::eta_squared(aovreport,partial=TRUE)
#for nonverbal associative memory RTs removing outliers
boxplot(df_cl$NVMemoryAsso_RTOverall)
Q2 <- quantile(df_cl$NVMemoryAsso_RTOverall, probs=c(.25, .75), na.rm = FALSE)
iqr2 <- IQR(df_cl$NVMemoryAsso_RTOverall)
up2 <-  Q2[2]+3*iqr2 # Upper Range  
low2<- Q2[1]-3*iqr2 # Lower Range

df_RT$NVMemoryAsso_RTOverall<- ifelse(df_RT$NVMemoryAsso_RTOverall<low2 | df_RT$NVMemoryAsso_RTOverall>up2,NA,df_RT$NVMemoryAsso_RTOverall)
sum(is.na(df_RT$NVMemoryAsso_RTOverall))
aovreport <- aov(NVMemoryAsso_RTOverall ~ COVID_Grouping + Q1.01Age_Text + Q1.02Sex_Quantised + Q1.03Country_Quantised + Q1.08.Education_Quantised, data = df_RT)
Anova(aovreport, type="III")
effectsize::eta_squared(aovreport,partial=TRUE)
#creating a scaled csv files to do the same in JASP
df_RT$VMemoryRecog_RTOverall <- scale(df_RT$VMemoryRecog_RTOverall,center = TRUE, scale=TRUE)
df_RT$NVMemoryAsso_RTOverall <- scale(df_RT$NVMemoryAsso_RTOverall,center = TRUE, scale = TRUE)

df_RT<- df_RT %>% dplyr::select(Participant_Private_ID, Completion_Date, COVID_Grouping, Q1.01Age_Text,Q1.02Sex_Quantised,Q1.03Country_Quantised,Q1.08.Education_Quantised,
       VMemoryAsso_CorrectPercent,VMemoryAsso_RTOverall,VMemoryRecog_CorrectPercent,VMemoryRecog_RTOverall,NVMemoryRecog_CorrectPercent,NVMemoryRecog_RTOverall,NVMemoryAsso_CorrectPercent,NVMemoryAsso_RTOverall)


df_RT2 <- pivot_longer(df_RT, cols=c(VMemoryRecog_RTOverall,NVMemoryAsso_RTOverall),
                         names_to = "MemTest",
                         values_to = "MemRT",
                         values_drop_na = FALSE)
#CSV for analysing the effect of Covid status on Verbal Item Memory and Nonverbal Associative memory RTs and interaction with memory task (ANALYSIS 4) after outlier removal
write.csv(df_RT2,"COVCOG_wave2_RT_nooutliers.csv")

#plot with 2 memory tasks RTs with outlier removal applied (presented in ANALYSIS 4 section)
df_RT2$COVID_Grouping <- as.factor(df_RT2$COVID_Grouping)
ggplot(df_RT2, aes(y=MemRT, x=MemTest,color=COVID_Grouping)) + geom_boxplot(notch=TRUE,lwd=1.5,fatten=1) + scale_color_manual(values=cbPalette, labels=c("No-Covid", "Covid")) + theme_cowplot(12) + labs(title="Scaled memory accuracy for verbal item memory \n and nonverbal associative memory task", y="Z scores (RTs)", x="Memory task", color="Covid group") +
  theme(plot.title = element_blank(), legend.title = element_text(size =18, face="bold"),legend.text = element_text(size = 17), axis.title=element_text(size=18,face="bold"), axis.text.x=element_text(size=17),axis.text.y=element_text(size=17)) + scale_x_discrete(labels= c("Nonverbal Associative", "Verbal Item")) + scale_y_continuous(limits=c(-3,3),breaks = c(-3,-2,-1,0,1,2,3))


#Time Since Infection 
df_cl$Completion_Date <- as.Date(df_cl$Completion_Date, "%d/%m/%Y")
df_cl$DateFirstCovid <- paste(df_cl$Q3.06.1.1COVIDDate1Day,df_cl$Q3.06.1.1COVIDDate1Month,df_cl$Q3.06.1.1COVIDDate1Year) #putting together date of first covid into one column from 3
df_cl$DateFirstCovid <- as.Date(df_cl$DateFirstCovid, "%d %m %Y")
#
difftime(df_cl$Completion_Date,df_cl$DateFirstCovid, units = "days")

#number of days since last Covid infection
df_cl$Q3.06.1.2COVIDDate2 <- as.Date(df_cl$Q3.06.1.2COVIDDate2, "%d/%m/%Y")
df_cl$Q3.06.1.3COVIDDate3 <- as.Date(df_cl$Q3.06.1.3COVIDDate3, "%d/%m/%Y")
df_cl$Q3.06.1.4COVIDDate4 <- as.Date(df_cl$Q3.06.1.4COVIDDate4, "%d/%m/%Y")

#calculating the days since last infection
df_cl$TimeSinceLastInf <- ifelse(df_cl$Q3.06COVIDTimes_Quantised == 1, difftime(df_cl$Completion_Date,df_cl$DateFirstCovid, units = "days"),
                                 ifelse(df_cl$Q3.06COVIDTimes_Quantised == 2, difftime(df_cl$Completion_Date,df_cl$Q3.06.1.2COVIDDate2, units = "days"),
                                        ifelse(df_cl$Q3.06COVIDTimes_Quantised == 3, difftime(df_cl$Completion_Date,df_cl$Q3.06.1.3COVIDDate3, units = "days"),
                                               ifelse(df_cl$Q3.06COVIDTimes_Quantised == 4, difftime(df_cl$Completion_Date,df_cl$Q3.06.1.4COVIDDate4, units = "days"),NA))))

#excluding unrealistic time delay
is.na(df_cl$TimeSinceLastInf)<- df_cl$DateFirstCovid<"2020-01-01" #assigning NAs to dates that predate the spread of the pandemic - first reported patient outside of China was on 13 January 2020
is.na(df_cl$TimeSinceLastInf)<- df_cl$TimeSinceLastInf<0 # assigning NAs if the time since infection is negative - i.e. participant put date of most recent infection after date of data collection

#vaccination but with combined memory tasks z scores
df_vacc_all_expl<-df_vacc_all
df_vacc_all_expl$VMemoryAsso_CorrectPercent<-scale(df_vacc_all_expl$VMemoryAsso_CorrectPercent, center = TRUE, scale=TRUE)
df_vacc_all_expl$NVMemoryRecog_CorrectPercent<-scale(df_vacc_all_expl$NVMemoryRecog_CorrectPercent,center=TRUE,scale=TRUE)
df_vacc_all_expl$NVMemoryAsso_CorrectPercent<-scale(df_vacc_all_expl$VMemoryAsso_CorrectPercent, center = TRUE, scale=TRUE)
df_vacc_all_expl$VMemoryRecog_CorrectPercent<-scale(df_vacc_all_expl$NVMemoryRecog_CorrectPercent,center=TRUE,scale=TRUE)
df_vacc_all_expl <- pivot_longer(df_vacc_all_expl, cols=c(VMemoryRecog_CorrectPercent,NVMemoryAsso_CorrectPercent,VMemoryAsso_CorrectPercent,NVMemoryRecog_CorrectPercent),
                         names_to = "MemTest",
                         values_to = "MemAccuracy",
                         values_drop_na = TRUE)
#creating csv file for the analysis of the Effect of Vaccination status on accuracy across all memory tasks
write.csv(df_vacc_all_expl, "COVCOG_wave2_vacc_memtasks.csv")


# re-run the analyses with depression and anxiety as covariates and with category fluency number of repetitions - ADDING DEPRESSION AND ANXIETY AS COVARIATES
JASP_extended <- df_cl%>%dplyr::select(Participant_Private_ID, Completion_Date, COVID_Grouping, Q1.01Age_Text,Q1.02Sex_Quantised,Q1.03Country_Quantised,Q1.08.Education_Quantised,Q2.04.21Depression,Q2.04.22Anxiety,
               TimeSinceLastInf,VMemoryAsso_CorrectPercent,VMemoryAsso_RTOverall,VMemoryRecog_CorrectPercent,VMemoryRecog_RTOverall,NVMemoryRecog_CorrectPercent,NVMemoryRecog_RTOverall,NVMemoryAsso_CorrectPercent,NVMemoryAsso_RTOverall,
               DS_MaxSpan,DS_Score, CF_CorrectExclRep, CF_repititions, SynCorrPerc,Syn_RToverall,WCST_Correct,WCST_RT_correct)

#Outlier removal applied across all tasks
mod_data <- JASP_extended%>%dplyr::select(Participant_Private_ID, Completion_Date, COVID_Grouping, Q1.01Age_Text,Q1.02Sex_Quantised,Q1.03Country_Quantised,Q1.08.Education_Quantised,Q2.04.21Depression,Q2.04.22Anxiety,
                                   TimeSinceLastInf,VMemoryAsso_CorrectPercent,VMemoryAsso_RTOverall,VMemoryRecog_CorrectPercent,VMemoryRecog_RTOverall,NVMemoryRecog_CorrectPercent,NVMemoryRecog_RTOverall,NVMemoryAsso_CorrectPercent,NVMemoryAsso_RTOverall,
                                   DS_MaxSpan,DS_Score, CF_CorrectExclRep, CF_repititions, SynCorrPerc,Syn_RToverall,WCST_Correct,WCST_RT_correct)
cols_all <- c("VMemoryAsso_CorrectPercent","VMemoryAsso_RTOverall","VMemoryRecog_CorrectPercent","VMemoryRecog_RTOverall","NVMemoryRecog_CorrectPercent","NVMemoryRecog_RTOverall","NVMemoryAsso_CorrectPercent","NVMemoryAsso_RTOverall",
              "DS_MaxSpan","DS_Score", "CF_CorrectExclRep", "CF_repititions", "SynCorrPerc","Syn_RToverall","WCST_Correct","WCST_RT_correct")

mod_data[cols_all]<-sapply(mod_data[cols_all],as.numeric) #all task columns as numeric

#how many NAs were there prior to outlier removal (tracking how many data points were removed)
na_prior<-colSums(is.na(mod_data))

colSums(!is.na(mod_data)) #calculating numbers prior to outlier removal for each task

#removing outliers
outliers<-function(variable){
  iqr<-IQR(variable,na.rm=TRUE)
  Q <- quantile(variable, probs=c(.25, .75),na.rm=TRUE)
  mild_low<-Q[1]-(3*iqr)
  mild_high<-Q[2]+(3*iqr)
  new_variable<-ifelse(variable < mild_low | variable > mild_high,NA, variable)
  return(new_variable)
}

mod_data[cols_all] <- sapply(mod_data[cols_all], outliers) #applying outlier removal

#how many data points removed per variable
na_post <- colSums(is.na(mod_data))
#scaling after outlier removal
mod_data[cols_all] <- sapply(mod_data[cols_all], scale, center=TRUE,scale=TRUE)

#file for the analysis of depression and anxiety as covariates and the analysis of memory task RTs for Analysis 2 and for check of additional registered analyses without outliers
write.csv(mod_data, "COVCOG_wave2_extended.csv")

#categorising Time Since Infection so that it can be analysed within ANCOVA
quantile(mod_data$TimeSinceLastInf,na.rm=TRUE)

mod_data$TimeSinceLastInfCat <- ifelse(mod_data$TimeSinceLastInf<= 84, 1, 
                                       ifelse(mod_data$TimeSinceLastInf > 84 & mod_data$TimeSinceLastInf <=197,2,
                                              ifelse(mod_data$TimeSinceLastInf > 197 & mod_data$TimeSinceLastInf<=443,3,
                                                     ifelse(mod_data$TimeSinceLastInf >443 & mod_data$TimeSinceLastInf <= 941,4,NA))))

write.csv(mod_data, "COVCOG_wave2_tSinceInfCat.csv")

#calculate how many rows there are for each cognitive variable for the time since infection analysis
numbersTimSInf<-mod_data %>% drop_na(TimeSinceLastInfCat)
data.frame(colSums(!is.na(numbersTimSInf)))


#plotting the RTs for memory tasks after outlier removal
mod_data$COVID_Grouping <- as.factor(mod_data$COVID_Grouping)
df_allMemRT <- pivot_longer(mod_data, cols=c(VMemoryAsso_RTOverall,VMemoryRecog_RTOverall,NVMemoryRecog_RTOverall,NVMemoryAsso_RTOverall),
                            names_to="MemTest",
                            values_to="MemRT",
                            values_drop_na=FALSE)

df_allMemRT$MemType<- if_else(grepl("Recog",df_allMemRT$MemTest),"Recog", 
                              if_else(grepl("Asso",df_allMemRT$MemTest),"Asso",NA))
                                                              
df_allMemRT$MemStim <- if_else(startsWith(df_allMemRT$MemTest,"N"),"Nonverbal","Verbal")

#plot of RTs shown in ANALYSIS 5 section
ggplot(df_allMemRT,aes(y=MemRT, x= MemType, colour=COVID_Grouping)) + geom_boxplot(notch=TRUE,fatten =1,lwd=1.5) + facet_grid(~MemStim) + theme_cowplot(12) +
  labs(title="The scaled accuracy scores for Memory tasks", y="Z scores (RTs)",x="Memory task", color="Covid group") + scale_color_manual(values=cbPalette,labels=c("No-Covid", "Covid")) +
  theme(plot.title = element_blank(), legend.title = element_text(size =18,face="bold"),legend.text = element_text(size = 17), axis.title=element_text(size=18,face="bold"), axis.text.x=element_text(size=17),axis.text.y=element_text(size=17), strip.text = element_text(
    size = 18)) + scale_x_discrete(labels= c("Associative", "Item")) + scale_y_continuous(limits=c(-5,5),breaks=c(-5,-2.5,0,2.5,5))


#creating a df for analysing the interaction between memory type and Covid status with depression and anxiety as covariates
df_allMemAcc <- pivot_longer(mod_data, cols=c(VMemoryAsso_CorrectPercent,VMemoryRecog_CorrectPercent,NVMemoryRecog_CorrectPercent,NVMemoryAsso_CorrectPercent),
                             names_to="MemTest",
                             values_to="MemAcc",
                             values_drop_na=FALSE)
df_allMemAcc$MemType<- if_else(grepl("Recog",df_allMemAcc$MemTest),"Recog", 
                          if_else(grepl("Asso",df_allMemAcc$MemTest),"Asso",NA))
df_allMemAcc$MemStim <- if_else(startsWith(df_allMemAcc$MemTest,"N"),"Nonverbal","Verbal")

#joining the pivoted df together
df_allMemAcc <- df_allMemAcc%>%dplyr::select("Participant_Private_ID","MemType","MemStim","MemAcc")
df_allMem <- merge(df_allMemRT,df_allMemAcc, by=c("Participant_Private_ID","MemType","MemStim"))
head(df_allMem)

write.csv(df_allMem,"COVCOG_wave2_explorative_reanalysis5.csv")



#analysing the count data
#apply log transformation
# Add 1 to avoid log(0) issues if there are zeros in your data
mod_data_nosc$log_CF_rep <- log(mod_data_nosc$CF_repititions + 1)
mod_data_nosc$log_DS_Score <- log(mod_data_nosc$DS_Score + 1)
mod_data_nosc$log_DS_MaxSpan <- log(mod_data_nosc$DS_MaxSpan + 1)
# #add the TimeSinceInfection column for the Time since infection analysis
# TimeSI <- mod_data %>% dplyr::select(Participant_Private_ID,TimeSinceLastInfCat)
# mod_data_nosc <- merge(mod_data_nosc, TimeSI, by = "Participant_Private_ID")
#write the csv for analysis of the log transformed columns
write.csv(mod_data_nosc,"COVCOG2_wave2_logtrans.csv")


#self reported memory correlational analysis
col_mem<-c("Q4.01.105SpatialMemoryProblem_Quantised","Q4.01.106ShortTermEventMemory_Quantised","Q4.01.107MediumTermEventMemory_Quantised","Q4.01.108LongTermEventMemory_Quantised",
           "Q4.11.105SpatialMemoryProblem_Quantised","Q4.11.106ShortTermEventMemory_Quantised","Q4.11.107MediumTermEventMemory_Quantised","Q4.11.108LongTermEventMemory_Quantised",
           "Q4.13.105SpatialMemoryProblem_Quantised","Q4.13.106ShortTermEventMemory_Quantised","Q4.13.107MediumTermEventMemory_Quantised","Q4.13.108LongTermEventMemory_Quantised")
df_cl[col_mem]<-sapply(df_cl[col_mem],as.numeric)
mem_data <- pivot_longer(df_cl, cols=c("Q4.01.105SpatialMemoryProblem_Quantised","Q4.01.106ShortTermEventMemory_Quantised","Q4.01.107MediumTermEventMemory_Quantised","Q4.01.108LongTermEventMemory_Quantised",
                                       "Q4.11.105SpatialMemoryProblem_Quantised","Q4.11.106ShortTermEventMemory_Quantised","Q4.11.107MediumTermEventMemory_Quantised","Q4.11.108LongTermEventMemory_Quantised",
                                       "Q4.13.105SpatialMemoryProblem_Quantised","Q4.13.106ShortTermEventMemory_Quantised","Q4.13.107MediumTermEventMemory_Quantised","Q4.13.108LongTermEventMemory_Quantised"),
             names_to = "Question",
             values_to = "Reported_memory",
             values_drop_na = FALSE)
#dataframe with the self-report memory measures
selfMemory <- df_cl %>% dplyr::select(Participant_Private_ID, Q4.01.105SpatialMemoryProblem_Quantised,Q4.01.106ShortTermEventMemory_Quantised,Q4.01.107MediumTermEventMemory_Quantised,Q4.01.108LongTermEventMemory_Quantised,
                                      Q4.11.105SpatialMemoryProblem_Quantised,Q4.11.106ShortTermEventMemory_Quantised,Q4.11.107MediumTermEventMemory_Quantised,Q4.11.108LongTermEventMemory_Quantised,
                                      Q4.13.105SpatialMemoryProblem_Quantised,Q4.13.106ShortTermEventMemory_Quantised,Q4.13.107MediumTermEventMemory_Quantised,Q4.13.108LongTermEventMemory_Quantised)


#for accuracy:
MemTests <- mod_data_nosc %>%dplyr::select(Participant_Private_ID,COVID_Grouping,Q1.01Age_Text,VMemoryAsso_CorrectPercent,VMemoryRecog_CorrectPercent,NVMemoryRecog_CorrectPercent,NVMemoryAsso_CorrectPercent)
onlyMemAcc <- merge(MemTests, selfMemory, by = "Participant_Private_ID")

onlyMemAcc$COVID_Grouping <- as.factor(onlyMemAcc$COVID_Grouping)

result <- onlyMemAcc %>%
  group_by(COVID_Grouping) %>%
  dplyr::summarize(spatial_3w = mean(Q4.11.105SpatialMemoryProblem_Quantised, na.rm = TRUE),
                   spatial_1to2days = mean(Q4.13.105SpatialMemoryProblem_Quantised, na.rm = TRUE),
                   short_3w = mean(Q4.11.106ShortTermEventMemory_Quantised,na.rm = T),
                   short_1to2days = mean(Q4.13.106ShortTermEventMemory_Quantised,na.rm=T),
                   medium_3w = mean(Q4.11.107MediumTermEventMemory_Quantised,na.rm=T),
                   medium_1to2days = mean(Q4.13.107MediumTermEventMemory_Quantised,na.rm=T),
                   long_3w = mean(Q4.11.108LongTermEventMemory_Quantised,na.rm=T),
                   long_1to2days = mean(Q4.13.108LongTermEventMemory_Quantised,na.rm=T))

print(result)
result_piv <- pivot_longer(result, cols=c(spatial_3w,spatial_1to2days, short_3w,short_1to2days,medium_3w,medium_1to2days,long_3w,long_1to2days),         
                           names_to = "Self_rep",
                           values_to = "Mean_mem",
                           values_drop_na = TRUE)
ggplot(result_piv, aes(y=Mean_mem, x=Self_rep, colour=COVID_Grouping)) + geom_point(size=6,shape=17) + theme_cowplot(12) +
  labs(title="Mean self-report ratings of memory by Covid status", y="Mean memory rating (scale 1-7)",x="Self-report memory type/time period", color="Covid group") + scale_color_manual(values=cbPalette,labels=c("No-Covid", "Covid")) +
  theme(plot.title = element_blank(), legend.title = element_text(size =18,face="bold"),legend.text = element_text(size = 17), axis.title=element_text(size=18,face="bold"), axis.text.x=element_text(size=17),axis.text.y=element_text(size=17)) + 
  scale_x_discrete(labels = c(
    "long-term memory\n1-2d prior", 
    "long-term memory\n3w post", 
    "medium-term memory\n1-2d prior", 
    "medium-term memory\n3w post", 
    "short-term memory\n1-2d prior", 
    "short-term memory\n3w post", 
    "spatial memory\n1-2d prior", 
    "spatial memory\n3w post"
  ))





#grid
# column names without the "self-report" mentioned for shorter labels
colnames(onlyMemAcc) <- c("Participant_Private_ID", "COVID_Grouping", "Age", "Verbal Associative Memory", 
                          "Verbal Item Memory", "Nonverbal Item Memory", "Nonverbal Associative Memory", 
                          "Spatial Memory at Inf.", "Short-Term Memory at Inf.", "Medium-Term Memory at Inf.", 
                          "Long-Term Memory at Inf.", "Spatial Memory 3w Post", "Short-Term Memory 3w Post", 
                          "Medium-Term Memory 3w Post", "Long-Term Memory 3w Post", 
                          "Spatial Memory 1-2d Prior", "Short-Term Memory 1-2d Prior", 
                          "Medium-Term Memory 1-2d Prior", "Long-Term Memory 1-2d Prior")

# Define the colors for Self-Report and Non-Self-Report labels
label_faces_x <- c(rep("plain", 5), rep("bold.italic", 12))
label_faces_y <- c(rep("plain", 4), rep("bold.italic", 12))

#All Accuracy
# Run correlation analysis and p-value adjustments
dat1 <- onlyMemAcc[3:19]
dat1 <- Hmisc::rcorr(as.matrix(dat1), type = "pearson")
p_values_vector <- as.vector(dat1$P)
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "fdr")
adjusted_p_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(dat1$P), ncol = ncol(dat1$P))
rownames(adjusted_p_matrix) <- rownames(dat1$P)
colnames(adjusted_p_matrix) <- colnames(dat1$P)
dat1$P_corrected <- adjusted_p_matrix
dat1$P_binary <- ifelse(dat1$P_corrected < 0.05, 0.06, 0.01)

# Generate the correlation plot 
ggcorrplot(dat1$r, type="upper", p.mat=dat1$P_binary, sig.level = 0.05, outline.col = "white", pch = 8,
           ggtheme = theme_bw) + 
  ggtitle("Correlation Plot of Memory Task Accuracy") +  
  theme(plot.title = element_blank(),  
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = label_faces_x),  
        axis.text.y = element_text(size = 13, face = label_faces_y),
        legend.text = element_text(size = 11)) + 
  scale_fill_gradient2(low ="#0072B2", high ="#E69F00", mid = "white",
                       midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name="Correlation",
                       breaks = c( -1,-0.5,-0.3, 0 ,0.3, 0.5,1),  
                       labels = c("-1","-0.5","-0.3"," 0"," 0.3"," 0.5"," 1"),
                       guide = guide_colorbar(barwidth = 1, barheight = 10)) 

#Covid Accuracy
# Run correlation analysis and p-value adjustments
onlyMemAcc_cov <- onlyMemAcc %>% filter(COVID_Grouping==1)
dat2<-onlyMemAcc_cov[3:19]
dat2 <- Hmisc::rcorr(as.matrix(dat2),type = "pearson")
#correcting for multiple comparisons
# Convert the p-value matrix to a vector
p_values_vector <- as.vector(dat2$P)
# Apply FDR adjustment
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "fdr")
# Convert the adjusted p-values vector back to a matrix
adjusted_p_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(dat2$P), ncol = ncol(dat2$P))
# Assign row and column names to the adjusted matrix
rownames(adjusted_p_matrix) <- rownames(dat2$P)
colnames(adjusted_p_matrix) <- colnames(dat2$P)
# Store the adjusted p-values back in the original dat2 object
dat2$P_corrected <- adjusted_p_matrix
#creating a flipped p-value matrix - values that were below 0.05 after correction are now at 0.06 and values that were above 0.05 are now at 0.01. This so that on the plot its significant correlations that have a mark on them, default is for nonsignificant
dat2$P_binary <- ifelse(dat2$P_corrected<0.05,0.06,0.01) 


# Generate the correlation plot 
ggcorrplot(dat2$r, type="upper", p.mat=dat2$P_binary, sig.level = 0.05, outline.col = "white", pch = 8,
           ggtheme = theme_bw) + 
  ggtitle("Correlation Plot of Memory Task Accuracy") +  
  theme(plot.title = element_blank(),  
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = label_faces_x),  
        axis.text.y = element_text(size = 13, face = label_faces_y),
        legend.text = element_text(size = 11)) + 
  scale_fill_gradient2(low ="#0072B2", high ="#E69F00", mid = "white",
                       midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name="Correlation",
                       breaks = c( -1,-0.5,-0.3, 0 ,0.3, 0.5,1),  
                       labels = c("-1","-0.5","-0.3"," 0"," 0.3"," 0.5"," 1"),
                       guide = guide_colorbar(barwidth = 1, barheight = 10)) 



#No-Covid Accuracy
# Run correlation analysis and p-value adjustments
onlyMemAcc_nc <- onlyMemAcc %>% filter(COVID_Grouping==0)
dat3 <-onlyMemAcc_nc[3:19]

dat3<-Hmisc::rcorr(as.matrix(dat3), type="pearson")
#correcting for multiple comparisons
# Convert the p-value matrix to a vector
p_values_vector <- as.vector(dat3$P)
# Apply FDR adjustment
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "fdr")
# Convert the adjusted p-values vector back to a matrix
adjusted_p_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(dat3$P), ncol = ncol(dat3$P))
# Assign row and column names to the adjusted matrix
rownames(adjusted_p_matrix) <- rownames(dat3$P)
colnames(adjusted_p_matrix) <- colnames(dat3$P)
# Store the adjusted p-values back in the original dat1 object
dat3$P_corrected <- adjusted_p_matrix
#creating a flipped p-value matrix - values that were below 0.05 after correction are now at 0.06 and values that were above 0.05 are now at 0.01. This so that on the plot its significant correlations that have a mark on them, default is for nonsignificant
dat3$P_binary <- ifelse(dat3$P_corrected<0.05,0.06,0.01) 


# Generate the correlation plot 
ggcorrplot(dat3$r, type="upper", p.mat=dat3$P_binary, sig.level = 0.05, outline.col = "white", pch = 8,
           ggtheme = theme_bw) + 
  ggtitle("Correlation Plot of Memory Task Accuracy") +  
  theme(plot.title = element_blank(),  
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = label_faces_x),  
        axis.text.y = element_text(size = 13, face = label_faces_y),
        legend.text = element_text(size = 11)) + 
  scale_fill_gradient2(low ="#0072B2", high ="#E69F00", mid = "white",
                       midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name="Correlation",
                       breaks = c( -1,-0.5,-0.3, 0 ,0.3, 0.5,1),  
                       labels = c("-1","-0.5","-0.3"," 0"," 0.3"," 0.5"," 1"),
                       guide = guide_colorbar(barwidth = 1, barheight = 10)) 

#Fisher's Z transformation to compare the correlations between groups
dat2$z<-0.5*(log(1+dat2$r) - log(1-dat2$r))
dat3$z <-0.5*(log(1+dat3$r) - log(1-dat3$r))

z_diff_acc<- (dat2$z - dat3$z)/sqrt((1/(dat2$n-3))+(1/(dat3$n-3)))
p_diff_acc <- 2*pnorm(-abs(z_diff_acc))
p_diff_acc_ad <- p.adjust(p_diff_acc, method="fdr")
p_diff_acc_ad <-matrix(p_diff_acc_ad, nrow = nrow(p_diff_acc), ncol = ncol(p_diff_acc))
# Assign row and column names to the adjusted matrix
rownames(p_diff_acc_ad) <- rownames(p_diff_acc)
colnames(p_diff_acc_ad) <- colnames(p_diff_acc)

ggcorrplot(z_diff_acc, p.mat = p_diff_acc_ad, type="upper",insig="blank",lab=TRUE,ggtheme = ggplot2::theme_bw) + ggtitle("Z scores of differences between correlations in the Covid and No-Covid group") + 
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(size = 13),  # Adjust x-axis text size
        axis.text.y = element_text(size=13),
        legend.text=element_blank()) 


#RTS
MemTestsRT<- mod_data_nosc %>%dplyr::select(Participant_Private_ID,COVID_Grouping,Q1.01Age_Text,VMemoryAsso_RTOverall,VMemoryRecog_RTOverall,NVMemoryRecog_RTOverall,NVMemoryAsso_RTOverall)
onlyMemRT <- merge(MemTestsRT, selfMemory, by = "Participant_Private_ID")

colnames(onlyMemRT) <- c("Participant_Private_ID", "COVID_Grouping", "Age", "Verbal Associative Memory", 
                         "Verbal Item Memory", "Nonverbal Item Memory", "Nonverbal Associative Memory", 
                         "Spatial Memory at Inf.", "Short-Term Memory at Inf.", "Medium-Term Memory at Inf.", 
                         "Long-Term Memory at Inf.", "Spatial Memory 3w Post", "Short-Term Memory 3w Post", 
                         "Medium-Term Memory 3w Post", "Long-Term Memory 3w Post", 
                         "Spatial Memory 1-2d Prior", "Short-Term Memory 1-2d Prior", 
                         "Medium-Term Memory 1-2d Prior", "Long-Term Memory 1-2d Prior")
str(onlyMemRT)
dat1<-onlyMemRT[3:19]

dat1 <- Hmisc::rcorr(as.matrix(dat1),type = "pearson")
# Convert the p-value matrix to a vector
p_values_vector <- as.vector(dat1$P)
# Apply FDR adjustment
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "fdr")
# Convert the adjusted p-values vector back to a matrix
adjusted_p_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(dat1$P), ncol = ncol(dat1$P))
# Assign row and column names to the adjusted matrix
rownames(adjusted_p_matrix) <- rownames(dat1$P)
colnames(adjusted_p_matrix) <- colnames(dat1$P)
# Store the adjusted p-values back in the original dat1 object
dat1$P_corrected <- adjusted_p_matrix
#creating a flipped p-value matrix - values that were below 0.05 after correction are now at 0.06 and values that were above 0.05 are now at 0.01. This so that on the plot its significant correlations that have a mark on them, default is for nonsignificant
dat1$P_binary <- ifelse(dat1$P_corrected<0.05,0.06,0.01)

# Generate the correlation plot 
ggcorrplot(dat1$r, type="upper", p.mat=dat1$P_binary, sig.level = 0.05, outline.col = "white", pch = 8,
           ggtheme = theme_bw) + 
  ggtitle("Correlation Plot of Memory Task Accuracy") +  
  theme(plot.title = element_blank(),  
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = label_faces_x),  
        axis.text.y = element_text(size = 13, face = label_faces_y),
        legend.text = element_text(size = 11)) + 
  scale_fill_gradient2(low ="#0072B2", high ="#E69F00", mid = "white",
                       midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name="Correlation",
                       breaks = c( -1,-0.5,-0.3, 0 ,0.3, 0.5,1),  
                       labels = c("-1","-0.5","-0.3"," 0"," 0.3"," 0.5"," 1"),
                       guide = guide_colorbar(barwidth = 1, barheight = 10)) 

#for groups - COVID
onlyMemRT_cov <- onlyMemRT %>% filter(COVID_Grouping==1)
dat2<-onlyMemRT_cov[3:19]

dat2 <- Hmisc::rcorr(as.matrix(dat2),type = "pearson")
# Convert the p-value matrix to a vector
p_values_vector <- as.vector(dat2$P)
# Apply FDR adjustment
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "fdr")
# Convert the adjusted p-values vector back to a matrix
adjusted_p_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(dat2$P), ncol = ncol(dat2$P))
# Assign row and column names to the adjusted matrix
rownames(adjusted_p_matrix) <- rownames(dat2$P)
colnames(adjusted_p_matrix) <- colnames(dat2$P)
# Store the adjusted p-values back in the original dat1 object
dat2$P_corrected <- adjusted_p_matrix
#creating a flipped p-value matrix - values that were below 0.05 after correction are now at 0.06 and values that were above 0.05 are now at 0.01. This so that on the plot its significant correlations that have a mark on them, default is for nonsignificant
dat2$P_binary <- ifelse(dat2$P_corrected<0.05,0.06,0.01)

# Generate the correlation plot 
ggcorrplot(dat2$r, type="upper", p.mat=dat2$P_binary, sig.level = 0.05, outline.col = "white", pch = 8,
           ggtheme = theme_bw) + 
  ggtitle("Correlation Plot of Memory Task Accuracy") +  
  theme(plot.title = element_blank(),  
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = label_faces_x),  
        axis.text.y = element_text(size = 13, face = label_faces_y),
        legend.text = element_text(size = 11)) + 
  scale_fill_gradient2(low ="#0072B2", high ="#E69F00", mid = "white",
                       midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name="Correlation",
                       breaks = c( -1,-0.5,-0.3, 0 ,0.3, 0.5,1),  
                       labels = c("-1","-0.5","-0.3"," 0"," 0.3"," 0.5"," 1"),
                       guide = guide_colorbar(barwidth = 1, barheight = 10)) 


#for groups - NOCOVID
onlyMemRT_nc <- onlyMemRT %>% filter(COVID_Grouping==0)
dat3 <-onlyMemRT_nc[3:19]

dat3<-Hmisc::rcorr(as.matrix(dat3), type="pearson")
#correcting for multiple comparisons
# Convert the p-value matrix to a vector
p_values_vector <- as.vector(dat3$P)
# Apply FDR adjustment
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "fdr")
# Convert the adjusted p-values vector back to a matrix
adjusted_p_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(dat3$P), ncol = ncol(dat3$P))
# Assign row and column names to the adjusted matrix
rownames(adjusted_p_matrix) <- rownames(dat3$P)
colnames(adjusted_p_matrix) <- colnames(dat3$P)
# Store the adjusted p-values back in the original dat1 object
dat3$P_corrected <- adjusted_p_matrix
#creating a flipped p-value matrix - values that were below 0.05 after correction are now at 0.06 and values that were above 0.05 are now at 0.01. This so that on the plot its significant correlations that have a mark on them, default is for nonsignificant
dat3$P_binary <- ifelse(dat3$P_corrected<0.05,0.06,0.01)

# Generate the correlation plot 
ggcorrplot(dat3$r, type="upper", p.mat=dat3$P_binary, sig.level = 0.05, outline.col = "white", pch = 8,
           ggtheme = theme_bw) + 
  ggtitle("Correlation Plot of Memory Task Accuracy") +  
  theme(plot.title = element_blank(),  
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = label_faces_x),  
        axis.text.y = element_text(size = 13, face = label_faces_y),
        legend.text = element_text(size = 11)) + 
  scale_fill_gradient2(low ="#0072B2", high ="#E69F00", mid = "white",
                       midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name="Correlation",
                       breaks = c( -1,-0.5,-0.3, 0 ,0.3, 0.5,1),  
                       labels = c("-1","-0.5","-0.3"," 0"," 0.3"," 0.5"," 1"),
                       guide = guide_colorbar(barwidth = 1, barheight = 10)) 

#Fisher's z transformation for RTs
dat2$z<-0.5*(log(1+dat2$r) - log(1-dat2$r))
dat3$z <-0.5*(log(1+dat3$r) - log(1-dat3$r))

z_diff_RT<- (dat2$z - dat3$z)/sqrt((1/(dat2$n-3))+(1/(dat3$n-3)))
p_diff_RT <- 2*pnorm(-abs(z_diff_RT))
p_diff_RT_ad <- p.adjust(p_diff_RT, method="fdr")
p_diff_RT_ad <-matrix(p_diff_RT_ad, nrow = nrow(p_diff_RT), ncol = ncol(p_diff_RT))
# Assign row and column names to the adjusted matrix
rownames(p_diff_RT_ad) <- rownames(p_diff_RT)
colnames(p_diff_RT_ad) <- colnames(p_diff_RT)

ggcorrplot(z_diff_RT, p.mat = p_diff_RT_ad, type="upper",insig="blank",lab=TRUE,ggtheme = ggplot2::theme_bw) + ggtitle("Z scores of differences between correlations in the Covid and No-Covid group") + 
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(size = 13),  # Adjust x-axis text size
        axis.text.y = element_text(size=13),
        legend.text=element_blank())  # Adjust these values



#update library dependencies
#renv::snapshot()
  

