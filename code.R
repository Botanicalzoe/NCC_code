
  # NCC_code
Warming causes contrasting behavioral responses of spiders by changing their prey size spectra 
 ############ data analysis  ########
### 2023.7.5 #########
library(readxl)
library(lme4) 
library(nlme)
library(car)
library(data.table)
library(vegan)
library(dplyr)
library(piecewiseSEM)
########  example _spider body size the same as others #########
df<- read_excel("E:\\2021-11-蜘蛛文章\\1-NCC返修-2023.3.11\\20230630\\data_for spider_behavior_and_abundance_2023.7.5-final.xlsx",sheet= "2_spider_abundance")
df$year<- as.factor(year(df$date))
df$chambers<- as.factor(df$chambers)
df_large<- subset(df,species %in% "large_spider")  # small_spider

library(lmerTest)		 
M1<- lmer(mass_mg~treatment*year+(1|chambers)+(1|year:chambers),data=df_large,
          control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
summary(M1)
aov <- anova(M1)


M2<- glmer.nb(number~treatment*year+(1|chambers)+(1|year:chambers),data=df_large,
              control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
car::Anova(M2)    # for spider and prey abundance


###### mesh size of small spider  ########
M2<- glmer(mesh_size_mm2~treatment*year+(1|chambers:year),data=df_large,family=Gamma(link=log),
     control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
summary(M2)
anova(M2)
car::Anova(M2)


################### piecewiseSEM #####################
# mean of each variables per chamber per year ####
df<- read_excel("E:\\2021-11-蜘蛛文章\\data-20211115.xlsx",sheet= "raw_data_20230620")
df_tr <- df[,c("temperature","lower_soil_moisture","biomass_per_m2","relative_biomass_graminoids","num_shared","num_specialized_large","num_specialized_small",
               "large_spider_abund","small_spider_abund","small_spider_mesh_size","small_spider_web_diamter","large_spider_mesh_size","large_spider_web_diamter")] %>%
  decostand("standardize", na.rm = TRUE) %>%   
  bind_cols(df[,c(1:4)]) 
df_tr<- na.omit(df_tr)

model_psem<- psem(  #full model
  lme(lower_soil_moisture~ temperature,random= list(year=~1,chambers=~1),data=df_tr),
  lme(relative_biomass_graminoids~  lower_soil_moisture+ temperature,random= list(year=~1,chambers=~1),data=df_tr),  
  lme(biomass_per_m2~ lower_soil_moisture+ temperature,random= list(year=~1,chambers=~1),data=df_tr),
  lme(num_shared~ lower_soil_moisture+ temperature+relative_biomass_graminoids+ biomass_per_m2,random= list(year=~1,chambers=~1),data=df_tr), 
  lme(num_specialized_large~ temperature+ lower_soil_moisture+ relative_biomass_graminoids+ biomass_per_m2,random= list(year=~1,chambers=~1),data=df_tr),
  lme(num_specialized_small~ temperature+ lower_soil_moisture+ relative_biomass_graminoids+ biomass_per_m2,random= list(year=~1,chambers=~1),data=df_tr),
  lme(large_spider_abund~ temperature+lower_soil_moisture+ num_specialized_large+ num_shared,random= list(year=~1,chambers=~1),data=df_tr),  
  lme(small_spider_abund~  temperature+lower_soil_moisture+ num_shared+  num_specialized_small,random= list(year=~1,chambers=~1),data=df_tr),
  lme(small_spider_mesh_size~ temperature+ lower_soil_moisture+ num_shared+  num_specialized_small,random= list(year=~1,chambers=~1),data=df_tr),
  lme(large_spider_mesh_size~ temperature+ lower_soil_moisture+ num_shared+ num_specialized_large,random= list(year=~1,chambers=~1),data=df_tr),  
  lme(small_spider_web_diamter~ temperature+lower_soil_moisture+  num_shared+  num_specialized_small,random= list(year=~1,chambers=~1),data=df_tr),
  lme(large_spider_web_diamter~ temperature+ lower_soil_moisture+ num_shared+ num_specialized_large,random= list(year=~1,chambers=~1),data=df_tr),
  (num_shared %~~% num_specialized_large),
  (num_shared %~~% num_specialized_small),
  (num_specialized_large %~~% num_specialized_small),
  (small_spider_abund %~~% large_spider_abund),
  (large_spider_abund %~~% large_spider_mesh_size),
  (large_spider_abund %~~% large_spider_web_diamter),
  (small_spider_abund %~~% small_spider_mesh_size),
  (small_spider_abund %~~% small_spider_web_diamter),
  (relative_biomass_graminoids %~~% biomass_per_m2),
   (large_spider_mesh_size %~~% large_spider_web_diamter),
  (small_spider_mesh_size%~~% small_spider_web_diamter)
)
summary(model_psem,.progressBar=F,standardize="scale")



###########  fitted model  ############
model_psem1<- psem(
  lme(lower_soil_moisture~ temperature,random= list(year=~1,chambers=~1),data=df_tr),
  lme(relative_biomass_graminoids~  lower_soil_moisture+ temperature,random= list(year=~1,chambers=~1),data=df_tr),  
  lme(biomass_per_m2~ temperature,random= list(year=~1,chambers=~1),data=df_tr),
  lme(num_shared~ relative_biomass_graminoids,random= list(year=~1,chambers=~1),data=df_tr), 
  lme(num_specialized_large~ lower_soil_moisture+ relative_biomass_graminoids,random= list(year=~1,chambers=~1),data=df_tr),
  lme(num_specialized_small~ temperature+ lower_soil_moisture+ biomass_per_m2,random= list(year=~1,chambers=~1),data=df_tr),
  lme(large_spider_abund~  num_specialized_large,random= list(year=~1,chambers=~1),data=df_tr),  
  lme(small_spider_abund~   num_specialized_small,random= list(year=~1,chambers=~1),data=df_tr),
  lme(small_spider_mesh_size~ temperature+  num_specialized_small,random= list(year=~1,chambers=~1),data=df_tr),
  lme(large_spider_mesh_size~  num_shared+ num_specialized_large,random= list(year=~1,chambers=~1),data=df_tr),  
  lme(small_spider_web_diamter~ num_specialized_small,random= list(year=~1,chambers=~1),data=df_tr),
  lme(large_spider_web_diamter~  num_shared+ num_specialized_large,random= list(year=~1,chambers=~1),data=df_tr),
  (num_shared %~~% num_specialized_small),  
  (small_spider_abund %~~% small_spider_mesh_size),
  (relative_biomass_graminoids %~~% biomass_per_m2)
)
summary(model_psem1,.progressBar=F,standardize="scale")

############  End    ##############
