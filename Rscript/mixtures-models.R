#--------------------------------------------------------#
# 				               ANALYSIS	               			   #
#--------------------------------------------------------#

rm(list=ls())

##########################################################
# Load packages:
library(ggplot2)
library(tidyr) # for gather
library(tidytext) # for reorder_within
library(dplyr) # for mutate
library(stringr)
library(scales)
library(mice)
library(forcats)
library(arm) # for bayesian glm
library(MASS) # for glm.nb

source("../00-functions.R")

##########################################################
# Load data

load("../data.RData")

# Prepare Intern and Extern scales 
nd <- nd %>%
  filter(test %in% "SDQ") %>%
  droplevels() %>%
  dplyr::select(nid, scale, value) %>%
  spread(scale, value) %>%
  mutate(Intern = Emotional + Peer,
         Extern = Conduct + Hyperactivity) %>%
  gather(scale, value, -nid)

# Select specific chemicals
pollutants <- pollutants %>%
  filter(compound %in% c("MMP", "MEP", "MBzP", 
                         "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP", 
                         "BPA", "BPF", 
                         "SumDINCH", "SumDEHTP")) %>%
  droplevels()

# Join data frames from chemicals, neurodevelopmental outcomes and demographics
d <- left_join(pollutants, nd) %>%
  left_join(demo)

# Remove data frames and keep final d data frame
rm(demo)
rm(nd)
rm(pollutants)

# Prepare binary variables
d <- d %>%
  mutate(EducationHigh = ifelse(educm == "2", 1, 0)) %>%
  mutate(HouseholdDivided = ifelse(household == "2", 1, 0)) %>% 
  mutate(HasSiblings = ifelse(siblings != "0", 1, 0)) %>%
  mutate(AgeSchoolOld = ifelse(age.school %in% c("2", "3"), 1, 0)) %>%
  mutate(Female = ifelse(sex == "female", 1, 0)) %>%
  mutate(TraumaticEvents = ifelse(traumatic == "1", 1, 0)) %>%
  mutate(ResidenceUrban = ifelse(resid == "1", 0, 1)) %>%
  mutate(SESlessAfluent = ifelse(ses == "2", 1, 0))

# Change of pro-social to anti-social
d %>%
  filter(scale == "Prosocial") %>%
  count(value)

d <- d %>%
  mutate(value = ifelse(scale == "Prosocial", 10 - value, value))


################################################################################## BAYES GWQS

#--------------------------------------------------------------------------- Load package
devtools::load_all("BayesGWQS")
library(ggmcmc)
library(coda)

# y vector of outcomes
# x matrix of component data
# z vector or matrix of covariates
# x.s vector of the number of components
# x.s number of components

############################################################# FINAL MODELS

#-------------------------------------------------------------- MODELS WITH 3 GROUPS
#-------------------------------------------------------------- PH/NON-PH/BPs

#-------------------------------------------------------- Linear GWQS with sex adjustment
#-------------------------------------------------------------- 3 families (ph/non-ph/bps)
#--------------------------------------------------------------------- 4 quantiles

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  # The reviewer says to focus on the following compounds
  # DINCH, DEHTP, 
  # DiDP, DiNP, DEHP, DBP and DiBP, and as individual chemicals, MMP, MEP, MBzP, 
  # BPA, BPS and BPF
  # But we avoid BPS because of low DF (20%)
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH"),
                          c("BPA", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=4, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "DINCH+DEHTP", "BPA+BPF"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="Final-BQWS-linear-3families-4quantiles-adjusted.RData")
load("Final-BQWS-linear-3families-4quantiles-adjusted.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "DINCH+DEHTP",
    Parameter == "beta3" ~ "BPA+BPF"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom")
ggsave("Final-BQWS-linear-3families-4quantiles-adjusted.pdf", height=6, width=8)

# Table of coefs
table <- SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "BPA",
    Parameter == "beta3" ~ "Substitutes"))) %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(median, 3), ci)) %>%
  mutate(pval.dic = ifelse(high < 0 | low > 0, "**",
                           ifelse(High <  0| Low > 0, "*", "-"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(Compound, Scale, beta.ci.pval) %>%
  spread(Scale, beta.ci.pval)
write.csv(table, "Final-BQWS-linear-3families-4quantiles-adjusted.csv")


#-------------------------------------------------------- Linear GWQS with sex stratification
#-------------------------------------------------------------- 3 families (ph/non-ph/bps)
#--------------------------------------------------------------------- 4 quantiles

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("SumDEHTP", "SumDINCH"),
                            c("BPA", "BPF"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=4, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "DINCH+DEHTP", "BPA+BPF"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="Final-BQWS-linear-3families-4quantiles-stratified.RData")
load("Final-BQWS-linear-3families-4quantiles-stratified.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("Final-BQWS-linear-3families-4quantiles-stratified.pdf", height=6, width=8)



# Table of coefs
table <- SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "BPA",
    Parameter == "beta3" ~ "Substitutes"))) %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(median, 3), ci)) %>%
  mutate(pval.dic = ifelse(high < 0 | low > 0, "**",
                           ifelse(High <  0| Low > 0, "*", "-"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(Compound, Scale, Sex, beta.ci.pval) %>%
  spread(Scale, beta.ci.pval)
write.csv(table, "Final-BQWS-linear-3families-4quantiles-stratified.csv")



#-------------------------------------------------------------- MODELS WITH 3 GROUPS
#-------------------------------------------------------------- PH/Substitutes/BPA

#-------------------------------------------------------- Linear GWQS with sex adjustment
#------------------------------------------------------- 3 families (ph/substitutes/BPA)
#--------------------------------------------------------------------- 4 quantiles

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear

  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("BPA"),
                          c("SumDEHTP", "SumDINCH", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=4, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "BPA", "Substitutes"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="Final-BQWS-linear-3families-4quantiles-adjusted-alternative.RData")
load("Final-BQWS-linear-3families-4quantiles-adjusted-alternative.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "BPA",
    Parameter == "beta3" ~ "Substitutes"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom")
ggsave("Final-BQWS-linear-3families-4quantiles-adjusted-alternative.pdf", height=6, width=8)

# Table of coefs
table <- SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "BPA",
    Parameter == "beta3" ~ "Substitutes"))) %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(median, 3), ci)) %>%
  mutate(pval.dic = ifelse(high < 0 | low > 0, "**",
                           ifelse(High <  0| Low > 0, "*", "-"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(Compound, Scale, beta.ci.pval) %>%
  spread(Scale, beta.ci.pval)
write.csv(table, "Final-BQWS-linear-3families-4quantiles-adjusted-alternative.csv")




#-------------------------------------------------------- Linear GWQS with sex stratification
#-------------------------------------------------------------- 3 families (ph/non-ph/bps)
#--------------------------------------------------------------------- 4 quantiles

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("BPA"),
                            c("SumDEHTP", "SumDINCH", "BPF"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=4, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "BPA", "Substitutes"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="Final-BQWS-linear-3families-4quantiles-stratified-alternative.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "BPA",
    Parameter == "beta3" ~ "Substitutes"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("Final-BQWS-linear-3families-4quantiles-stratified-alternative.pdf", height=6, width=8)



# Table of coefs
table <- SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "BPA",
    Parameter == "beta3" ~ "Substitutes"))) %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(median, 3), ci)) %>%
  mutate(pval.dic = ifelse(high < 0 | low > 0, "**",
                           ifelse(High <  0| Low > 0, "*", "-"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(Compound, Scale, Sex, beta.ci.pval) %>%
  spread(Scale, beta.ci.pval)
write.csv(table, "Final-BQWS-linear-3families-4quantiles-stratified-alternative.csv")

############################################################# END OF FINAL MODELS

#-------------------------------------------------------------- MODELS WITH 3 GROUPS
#-------------------------------------------------------------- PH/NON-PH/BPs

#-------------------------------------------------------- Linear GWQS with sex adjustment
#-------------------------------------------------------------- 3 families (ph/non-ph/bps)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH"),
                          c("BPA", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=3, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "Non-phthalates", "Bisphenols"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-3fam-adjusted.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-3fam-adjusted.pdf", height=6, width=8)


#-------------------------------------------------------- Linear GWQS with sex stratification
#-------------------------------------------------------------- 3 families (ph/non-ph/bps)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH"),
                          c("BPA", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=3, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "Non-phthalates", "Bisphenols"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-3fam-stratified.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-3fam-stratified.pdf", height=6, width=8)


#------------------------------------------------------------ SAME AS BEFORE BUT 2 QUANTILES

#--------------------------------------------------- Linear GWQS with sex adjustment 2 QUANT
#---------------------------------------------------------- 3 families (ph/non-ph/bps)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH"),
                          c("BPA", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=2, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "Non-phthalates", "Bisphenols"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-3fam-adjusted-2quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-3fam-adjusted-2quant.pdf", height=6, width=8)


#--------------------------------------------- Linear GWQS with sex stratification 2 QUANT
#--------------------------------------------------------- 3 families (ph/non-ph/bps)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("SumDEHTP", "SumDINCH"),
                            c("BPA", "BPF"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=2, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "Non-phthalates", "Bisphenols"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-3fam-stratified-2quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-3fam-stratified-2quant.pdf", height=6, width=8)



#------------------------------------------------------------ SAME AS BEFORE BUT 4 QUANTILES

#--------------------------------------------------- Linear GWQS with sex adjustment 4 QUANT
#---------------------------------------------------------- 3 families (ph/non-ph/bps)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH"),
                          c("BPA", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=4, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "Non-phthalates", "Bisphenols"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-3fam-adjusted-4quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-3fam-adjusted-4quant.pdf", height=6, width=8)


#--------------------------------------------- Linear GWQS with sex stratification 4 QUANT
#--------------------------------------------------------- 3 families (ph/non-ph/bps)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("SumDEHTP", "SumDINCH"),
                            c("BPA", "BPF"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=4, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "Non-phthalates", "Bisphenols"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-3fam-stratified-4quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-3fam-stratified-4quant.pdf", height=6, width=8)



#-------------------------------------------------------------- MODELS WITH 2 GROUPS
#-------------------------------------------------------------- ORIGINAL / SUBSTITUTES

#-------------------------------------------------------- Linear GWQS with sex adjustment
#------------------------------------------------------ 2 families (original / substitutes)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP", "BPA"),
                          c("SumDEHTP", "SumDINCH", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 2, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 2, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=3, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates+BPA", "Non-phthalates+BPF"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-2fam-adjusted.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates+BPA",
    Parameter == "beta2" ~ "Non-phthalates+BPF"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-2fam-adjusted.pdf", height=6, width=8)


#-------------------------------------------------------- Linear GWQS with sex stratification
#------------------------------------------------------ 2 families (original / substitutes)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP", "BPA"),
                            c("SumDEHTP", "SumDINCH", "BPF"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 2, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 2, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=3, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates+BPA", "Non-phthalates+BPF"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-2fam-stratified.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates+BPA",
    Parameter == "beta2" ~ "Non-phthalates+BPF"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-2fam-stratified.pdf", height=6, width=8)


#------------------------------------------------------------ SAME AS BEFORE BUT 2 QUANTILES

#--------------------------------------------------- Linear GWQS with sex adjustment 2 QUANT
#------------------------------------------------------ 2 families (original / substitutes)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP", "BPA"),
                          c("SumDEHTP", "SumDINCH", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 2, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 2, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=2, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates+BPA", "Non-phthalates+BPF"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-2fam-adjusted-2quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates",
    Parameter == "beta3" ~ "Bisphenols"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-2fam-adjusted-2quant.pdf", height=6, width=8)


#--------------------------------------------- Linear GWQS with sex stratification 2 QUANT
#------------------------------------------------------ 2 families (original / substitutes)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP", "BPA"),
                            c("SumDEHTP", "SumDINCH", "BPF"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 2, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 2, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=2, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates+BPA", "Non-phthalates+BPF"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-2fam-stratified-2quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates+BPA",
    Parameter == "beta2" ~ "Non-phthalates+BPF"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-2fam-stratified-2quant.pdf", height=6, width=8)



#------------------------------------------------------------ SAME AS BEFORE BUT 4 QUANTILES

#--------------------------------------------------- Linear GWQS with sex adjustment 4 QUANT
#------------------------------------------------------ 2 families (original / substitutes)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP", "BPA"),
                          c("SumDEHTP", "SumDINCH", "BPF"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 2, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 2, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=4, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates+BPA", "Non-phthalates+BPF"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-2fam-adjusted-4quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates+BPA",
    Parameter == "beta2" ~ "Non-phthalates+BPF"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-2fam-adjusted-4quant.pdf", height=6, width=8)


#--------------------------------------------- Linear GWQS with sex stratification 4 QUANT
#------------------------------------------------------ 2 families (original / substitutes)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP", "BPA"),
                            c("SumDEHTP", "SumDINCH", "BPF"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 2, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 2, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=4, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates+BPA", "Non-phthalates+BPF"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-2fam-stratified-4quant.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates+BPA",
    Parameter == "beta2" ~ "Non-phthalates+BPF"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-2fam-stratified-4quant.pdf", height=6, width=8)



#---------------------------------------------------------- MODELS WITH 3 GROUPS alternative2
#-------------------------------------------------------------- PH/NON-PH+BPF/BPA

#-------------------------------------------------------- Linear GWQS with sex adjustment
#-------------------------------------------------------------- 3 families (ph/non-ph+BPF/BPA)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH", "BPF"),
                          c("BPA"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=3, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "Non-phthalates", "Bisphenols"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-3fam-adjusted-alt2.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates+BPF",
    Parameter == "beta3" ~ "BPA"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-3fam-adjusted-alt2.pdf", height=6, width=8)


#-------------------------------------------------------- Linear GWQS with sex stratification
#-------------------------------------------------------------- 3 families (ph/non-ph+BPF/BPA)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("SumDEHTP", "SumDINCH", "BPF"),
                            c("BPA"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=3, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "Non-phthalates+BPF", "BPA"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-3fam-stratified-alt2.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates+BPF",
    Parameter == "beta3" ~ "BPA"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-3fam-stratified-alt2.pdf", height=6, width=8)


#------------------------------------------------------------ SAME AS BEFORE BUT 2 QUANTILES

#--------------------------------------------------- Linear GWQS with sex adjustment 2 QUANT
#-------------------------------------------------------------- 3 families (ph/non-ph+BPF/BPA)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH", "BPF"),
                          c("BPA"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=2, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "Non-phthalates+BPF", "BPA"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-3fam-adjusted-2quant-alt2.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates+BPF",
    Parameter == "beta3" ~ "BPA"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-3fam-adjusted-2quant-alt2.pdf", height=6, width=8)


#--------------------------------------------- Linear GWQS with sex stratification 2 QUANT
#-------------------------------------------------------------- 3 families (ph/non-ph+BPF/BPA)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("SumDEHTP", "SumDINCH", "BPF"),
                            c("BPA"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=2, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "Non-phthalates+BPF", "BPA"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-3fam-stratified-2quant-alt2.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates+BPF",
    Parameter == "beta3" ~ "BPA"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-3fam-stratified-2quant-alt2.pdf", height=6, width=8)


#------------------------------------------------------------ SAME AS BEFORE BUT 4 QUANTILES

#--------------------------------------------------- Linear GWQS with sex adjustment 4 QUANT
#-------------------------------------------------------------- 3 families (ph/non-ph+BPF/BPA)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF", "BPS")) %>%
    unique() %>%
    dplyr::select(nid, value,
                  compound, concentration, 
                  cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                  HasSiblings, age, Female, bmi) %>%
    group_by(compound) %>% 
    mutate(concentration = escalamenta(log(concentration))) %>%
    ungroup() %>% 
    mutate(edatm = escalamenta(edatm),
           age = escalamenta(age),
           cotinine.7y = escalamenta(log(cotinine.7y)),
           bmi = escalamenta(bmi)) %>%
    spread(compound, concentration) %>%
    filter(!is.na(BPF)) %>%
    filter(!is.na(BPA)) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(HasSiblings)) %>%
    filter(!is.na(EducationHigh)) %>%
    filter(!is.na(HouseholdDivided)) %>%
    filter(!is.na(cotinine.7y)) %>%
    as.data.frame()
  table(complete.cases(d2))
  
  #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
  #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
  outcome <- d2$value # CONTINUOUS outcome  for linear
  
  group_compounds <- list(c("MMP", "MEP", "MBzP", 
                            "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                          c("SumDEHTP", "SumDINCH", "BPF"),
                          c("BPA"))
  group_compounds
  
  group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                         "HasSiblings", "age", "Female", "bmi", "cotinine.7y"))
  group_covars
  
  component <- make.X(d2, 3, group_compounds)
  head(component)
  component.num <- make.x.s(d2, 3, group_compounds)
  component.num
  
  covariates <- make.X(d2, 1, group_covars)
  head(covariates)
  
  work_dir <- tempdir()
  
  message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
  fit <- bgwqs.fit(y=outcome, x=component, 
                   x.s=component.num, 
                   z = covariates,
                   n.quantiles=4, # 4 quantiles
                   working.dir = work_dir)
  
  par(mfrow=c(1,3))
  weight.plot(fit, group.names = c("Phthalates", "Non-phthalates+BPF", "BPA"),
              group.list = group_compounds, x.s = component.num)
  
  fit1 <- fit
  class(fit1) <- c("list", "mcmc.list")
  
  S <- ggs(as.mcmc.list(fit$Samples))
  SS <- bind_rows(SS, mutate(S, Scale = Scale))
}

save(SS, file="250505-bqws-linear-3fam-adjusted-4quant-alt2.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates+BPF",
    Parameter == "beta3" ~ "BPA"))) %>%
  ggplot(aes(x = median, y = Scale, color = Compound)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") + xlim(c(-3, 3))
ggsave("250505-BQWS-linear-3fam-adjusted-4quant-alt2.pdf", height=6, width=8)


#--------------------------------------------- Linear GWQS with sex stratification 4 QUANT
#-------------------------------------------------------------- 3 families (ph/non-ph+BPF/BPA)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("SumDEHTP", "SumDINCH", "BPF"),
                            c("BPA"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=4, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "Non-phthalates+BPF", "BPA"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-3fam-stratified-4quant-alt2.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates+BPF",
    Parameter == "beta3" ~ "BPA"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-3fam-stratified-4quant-alt2.pdf", height=6, width=8)


#--------------------------------------------- Linear GWQS with sex stratification 5 QUANT
#-------------------------------------------------------------- 3 families (ph/non-ph+BPF/BPA)

# Big loop for each of the scales
Scales <- unique(d$scale)
#Scales <- Scales[!str_detect(Scales, "Total")]
nS <- length(Scales)

SS <- NULL
for (G in c("female", "male")) {
  for (s in 1:nS) {
    Scale <- Scales[s]
    d2 <- d %>%
      filter(sex == G) %>%
      filter(scale == Scales[s]) %>%
      filter(compound %in% c("MMP", "MEP", "MBzP", 
                             "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                             "SumDEHTP", "SumDINCH",
                             "BPA", "BPF", "BPS")) %>%
      unique() %>%
      dplyr::select(nid, value,
                    compound, concentration, 
                    cotinine.7y, edatm, EducationHigh, HouseholdDivided, 
                    HasSiblings, age, bmi) %>%
      group_by(compound) %>% 
      mutate(concentration = escalamenta(log(concentration))) %>%
      ungroup() %>% 
      mutate(edatm = escalamenta(edatm),
             age = escalamenta(age),
             cotinine.7y = escalamenta(log(cotinine.7y)),
             bmi = escalamenta(bmi)) %>%
      spread(compound, concentration) %>%
      filter(!is.na(BPF)) %>%
      filter(!is.na(BPA)) %>%
      filter(!is.na(value)) %>%
      filter(!is.na(HasSiblings)) %>%
      filter(!is.na(EducationHigh)) %>%
      filter(!is.na(HouseholdDivided)) %>%
      filter(!is.na(cotinine.7y)) %>%
      as.data.frame()
    table(complete.cases(d2))
    
    #outcome <- d2$Emotional # BINARY outcome (0/1 based on clinical vs. normal scales)
    #outcome <- d2$value # CONTINUOUS / Counts outcome  for negative binomial
    outcome <- d2$value # CONTINUOUS outcome  for linear
    
    group_compounds <- list(c("MMP", "MEP", "MBzP", 
                              "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP"),
                            c("SumDEHTP", "SumDINCH", "BPF"),
                            c("BPA"))
    group_compounds
    
    group_covars <- list(c("edatm", "EducationHigh", "HouseholdDivided", 
                           "HasSiblings", "age", "bmi", "cotinine.7y"))
    group_covars
    
    component <- make.X(d2, 3, group_compounds)
    head(component)
    component.num <- make.x.s(d2, 3, group_compounds)
    component.num
    
    covariates <- make.X(d2, 1, group_covars)
    head(covariates)
    
    work_dir <- tempdir()
    
    message(paste0("Processing ", Scale, ", with N = ", length(outcome)))
    fit <- bgwqs.fit(y=outcome, x=component, 
                     x.s=component.num, 
                     z = covariates,
                     n.quantiles=5, # 4 quantiles
                     working.dir = work_dir)
    
    par(mfrow=c(1,3))
    weight.plot(fit, group.names = c("Phthalates", "Non-phthalates+BPF", "BPA"),
                group.list = group_compounds, x.s = component.num)
    
    fit1 <- fit
    class(fit1) <- c("list", "mcmc.list")
    
    S <- ggs(as.mcmc.list(fit$Samples))
    SS <- bind_rows(SS, mutate(S, Scale = Scale, Sex = factor(G)))
  }
}
save(SS, file="250505-bqws-linear-3fam-stratified-5quant-alt2.RData")

SS %>%
  filter(Parameter != "beta0") %>% 
  filter(str_detect(Parameter, "^beta")) %>% 
  droplevels() %>% 
  group_by(Parameter, Scale, Sex) %>%
  summarize(median = quantile(value, 0.5),
            low = quantile(value, 0.025),
            high = quantile(value, 0.975),
            Low = quantile(value, 0.05),
            High = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  mutate(Compound = factor(case_when(
    Parameter == "beta1" ~ "Phthalates",
    Parameter == "beta2" ~ "Non-phthalates+BPF",
    Parameter == "beta3" ~ "BPA"))) %>%
  mutate(Scale = factor(case_when(
    Scale == "Conduct" ~ "Conduct problems",
    Scale == "Emotional" ~ "Emotional symptoms",
    Scale == "Hyperactivity" ~ "Hyperactivity/Inattention",
    Scale == "Peer" ~ "Peer relationships problems",
    Scale == "Prosocial" ~ "Prosocial behavior",
    Scale == "Total" ~ "Total difficulties",
    Scale == "Intern" ~ "Internalizing score",
    Scale == "Extern" ~ "Externalizing score"))) %>%
  mutate(Scale = fct_relevel(Scale, c("Conduct problems", "Emotional symptoms",
                                      "Hyperactivity/Inattention", 
                                      "Peer relationships problems",
                                      "Prosocial behavior", "Total difficulties",
                                      "Internalizing score",
                                      "Externalizing score"))) %>%
  ggplot(aes(x = median, y = Compound, color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4), lwd = 1) +
  xlab(expression(paste("Standardized ", beta, "-coefficients"))) +
  ylab("") +
  theme_bw() +
  facet_wrap(~Scale, scales="free") +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Set2")
ggsave("250505-BQWS-linear-3fam-stratified-5quant-alt2.pdf", height=6, width=8)


