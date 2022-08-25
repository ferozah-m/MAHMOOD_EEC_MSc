#Clear workspace
rm(list=ls())

#set and check working directory where datafiles are saved
setwd('/Users/ferozah/Desktop/')

#load in required packages
require(ggplot2)
require(dplyr)
require(tidyverse)
library(tidyr)
library(lme4)
library(sjPlot)
library(MuMIn)
library(RColorBrewer)

#initial plotting of potential random effects of brood weight and number of bees
# but a started with same number of bees and brood weight was relatively 
# even in method, this seems adequately controlled
data1<-read.csv("pre post brood and worker data.csv", header=TRUE)
boxplot(pre_number_of_bees~treatment,data=data1)
boxplot(post_number_of_bees~treatment,data=data1)
boxplot(pre_brood_weight~treatment,data=data1)
boxplot(post_brood_weight~treatment,data=data1)



#BROOD TEMPERTAUTRE (H1)
#load thermal data 
temp<-read.csv("all temp data .csv", header=TRUE)
str(temp)

#mean, sd and se of brood temp of each treatment group
temp %>% group_by(treatmentgroup) %>%
  summarise(mean = mean(brood_temp, na.rm = TRUE),
            sd = sd(brood_temp, na.rm=TRUE),
            se = sd / (sqrt(length(brood_temp))))

#construct a simple linear model
#compare multiple models: null model, full model and different models with all biologically relevant variables
M0<-lm(brood_temp~1, data = temp)
M1<-lm(brood_temp~TREATMENT, data = temp) 
M2<-lm(brood_temp~TEMP, data = temp)
M3<-lm(brood_temp~TREATMENT+TEMP, data = temp)
M4<-lm(brood_temp~TREATMENT*TEMP+day,data = temp)
M5<-lmer(brood_temp~TREATMENT*TEMP+day+(1|colony_of_origin)+(1|microcolony),data = temp)

#MODEL SELECTION
#Preform model selection based on AIC criteria,
model.sel(M0, M1, M2, M3, M4, M5)
#this suggests M4 is the best model to use
summary(M4)
anova(M4, test = "Chisq")
#final temp model
tempmodel<-lm(hot_spot~TREATMENT*TEMP+day,data = temp)
summary(tempmodel)
anova(tempmodel, test = "Chisq")
plot(tempmodel) #check model diagnostics 
tab_model(tempmodel)#obtain R squared and 95 Confidence Interval values

#box plot of brood temperature of treatment groups
ggplot(data = temp, aes(x = as.factor(treatmentgroup), as.numeric(brood_temp))) +
  stat_boxplot(geom ='errorbar',width=0.5)+
  geom_boxplot(fill=colours) +
  theme_classic() +
  geom_signif(comparisons = list(c("controlhigh", "controllow"),c("controlhigh","imihigh"),
                                 c("controlhigh","imilow"),c("controllow","imihigh"),
                                 c("controllow","imilow"),c("imihigh","imilow")),
              map_signif_level = TRUE, textsize = 4, col = 'black', vjust = 0.1,y_position = c(36,36,36,34,34,32.5)) +
  ggtitle('') + labs(x = "Treatment group", y = "Brood temperature (°C)")+
  theme(plot.title = element_text(hjust =0.5))

#load dataset for brood temperature over time plot
temptime<-read.csv("temp over time new.csv", header=TRUE)
temptimeplot<-ggplot(data = temptime, aes(x = day, y = brood.temperature,shape=treatment_group,group=treatment_group))+geom_point()
#manually create a colour picker with preffered colour blind firendly colours for all grpahs from colour brewer package
colours=c("#1B9E77","#D95F02","#7570B3","#E7298A")
temptimeplot+geom_smooth(method = "lm",aes(colour=treatment_group))+theme_classic() +ylab('Brood Temperature (°C)')+xlab("Day")+scale_colour_manual(values=colours)+theme(legend.position="top")+labs(colour="Treatment Group",shape="Treatment Group")



#PROPOTION (H2)
#load propotion dataset
prop<-read.csv("draft prop data.csv", header=TRUE)
#mean, sd and se propotion of each treatment group
prop %>% group_by(treatmentgroup) %>%
  summarise(mean = mean(propotion, na.rm = TRUE),
            sd = sd(propotion, na.rm=TRUE),
            se = sd / (sqrt(length(propotion))))

#construct a generalized linear model 
#compare multiple models: null model, full model and different models with all biologically relevant variables
#binomial models as proportion of on vs off brood can be converted into 1 vs 0 data
M0<- glm(propotion~1, data = prop, family = "binomial")
M1<- glm(propotion~TREATMENT, data = prop,family = "binomial" )
M2<- glm(propotion~TEMP, data = prop,family = "binomial" )
M3<- glm(propotion~TREATMENT+TEMP, data = prop,family = "binomial" )
M4<- glm(propotion~TREATMENT*TEMP+day, data = prop,family = "binomial" )
M5<- glmer(propotion~TREATMENT*TEMP+day+(1|colony_of_origin)+(1|microcolony),data = prop,family = "binomial" )
#MODEL SELECTION
#Preform model selection based on AIC criteria,
model.sel(M0, M1, M2, M3, M4, M5)
#this suggests M4 is the best model to use
summary(M4)
anova(M4, test = "Chisq")
#final propotion model
propmodel<-glm(propotion~TREATMENT*TEMP,data=prop,family = binomial)
summary(propmodel)
anova(propmodel,test="Chisq")
plot(propmodel) #check model diagnostics 
#dispersion is 12.165/16.704 = 0.728 - quite good
tab_model(propmodel)#obtain R squared and 95 Confidence Interval values
#r sqaured is 27%

#boxlplot
ggplot(data = prop, aes(x = treatmentgroup, y=propotion)) +
  stat_boxplot(geom ='errorbar',width=0.5)+
  geom_boxplot(fill=colours) +
  theme_classic() +
  geom_signif(comparisons = list(c("controlhigh", "controllow"),c("controlhigh","imihigh"),
                                 c("controlhigh","imilow"),c("controllow","imihigh"),
                                 c("controllow","imilow"),c("imihigh","imilow")),
              map_signif_level = TRUE, textsize = 4, col = 'black', vjust = 0.1,y_position = c(1.5,1.5,1.5,1.25,1.25,1)) +
  ggtitle('') + labs(x = "Treatment group", y = "Propotion")+
  theme(plot.title = element_text(hjust =0.5))

#proportion change over time graph
propovertime<-read.csv("Prop over time new.csv", header=TRUE)
str(propovertime)

#building model up 
proptimegraph <- ggplot(data = propovertime, aes(x = day, y = propotion,shape=treatment_group,group=treatment_group)) +
  geom_point()
proptimegraph + geom_smooth(method = "lm")+theme_classic() +ylab('Propotion of worker bees on brood')+xlab("Day")+ scale_color_manual(values = colours)+theme(legend.position="top")
proptimegraph + geom_smooth(method = "lm",aes(colour=treatment_group))+theme_classic() +ylab('Propotion of worker bees on brood')+xlab("Day")
theme(legend.title = element_text(colour="blue", size=10, 
                                  face="bold"))
proptimegraph + geom_smooth(method = "lm",aes(colour=treatment_group))+theme_classic() +ylab('Propotion of worker bees on brood')+xlab("Day")+scale_colour_manual(values=colours)+theme(legend.position="top")+labs(colour="Treatment Group",shape="Treatment Group") 


#ACTIVITY (H3)
#load in acitivty data
activity<-read.csv("activity data.csv", header=TRUE)
#cbind glm for inactivity 
inactive<-glm(cbind(inactive,3-inactive)~treatmentgroup,data=activity,family=binomial)
summary(inactive)
anova(inactive)
tab_model(inactive)
#cbind glm for activity
active<-glm(cbind(active,3-active)~treatmentgroup,data=activity,family=binomial)
summary(active)
anova(active)
tab_model(active)
#error so remove total counts over 4 as this are errors (only 3 measures taken per video so cannot be over 3)
#then reran model
activity$active[which(activity$active > 1)] <- 3
#cbind glm for brood care and nursing activity
brood<-glm(cbind(brood,3-brood)~treatmentgroup,data=activity,family=binomial)
summary(brood)
tab_model(brood)

#dispersion parameter for all models is high

#activity graphs for day 1 and day 5
activityplotd1<-read.csv("activityplotd1.csv", header=TRUE)
activityplotd5<-read.csv("activityplotd5.csv", header=TRUE)

d1<-ggplot(data=activityplotd1, aes(x=treatment_group, y=count, fill=Behaviour)) +
  geom_bar(stat="identity")+xlab("Treatment Group")+ylab('Average count of behaviour state')
d1+scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+theme_classic()+theme(legend.position="top")
d5<-ggplot(data=activityplotd5, aes(x=treatment_group, y=count, fill=Behaviour)) +
  geom_bar(stat="identity")+xlab("Treatment Group")+ylab('Average count of behaviour state')
d5+scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+theme_classic()+theme(legend.position="top")
