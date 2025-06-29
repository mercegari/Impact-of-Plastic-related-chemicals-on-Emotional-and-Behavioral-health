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

#---------------------------------- Comparison of Linear Models vs Negative Binomial Models

#------------ Linear model & Negative Binomial Model with sex adjustment ------------------

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4, 2),
           dimnames = list(Outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", 
                                       "Total difficulties",
                                       "Internalizing score",
                                       "Externalizing score"),
                           Compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval"),
                           Model = c("Linear", "Negative binomial")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
  for (c in 1:nC) {
    C <- compound.labels[c]
    message(paste0("  Compound: ", C))
      d.now <- d %>%
        filter(compound == C) %>%
        filter(scale == O) %>%
        dplyr::select(-nid, -compound, -scale)
    
    model <- bayesglm(value ~ escalamenta(edatm) + 
                      EducationHigh + HouseholdDivided + HasSiblings + Female +
                        escalamenta(age) + 
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(concentration)),
                      data = d.now, family = gaussian())
    model.nb <- glm.nb(value ~ escalamenta(edatm) + 
                      EducationHigh + HouseholdDivided + HasSiblings + Female +
                        escalamenta(age) + 
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(concentration)),
                      data = d.now)
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64, outerCI = 1.96, intercept = FALSE)
    sum <- summary(model)
    M[o,c,2,1] <- m$Value[1] # beta-coef
    M[o,c,1,1] <- m$LowOuter[1] # low CI (2 sigma)
    M[o,c,3,1] <- m$HighOuter[1] # high CI (2 sigma)
    M[o,c,4,1] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
    # Negative binomial
    m <- coefplot:::buildModelCI(model.nb, innerCI = 1.64, outerCI = 1.96, intercept = FALSE)
    sum <- summary(model.nb)
    M[o,c,2,2] <- m$Value[1] # beta-coef
    M[o,c,1,2] <- m$LowOuter[1] # low CI (2 sigma)
    M[o,c,3,2] <- m$HighOuter[1] # high CI (2 sigma)
    M[o,c,4,2] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
  }
}


Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$Compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(Outcome, beta.ci.pval)

# Plot
cbPalette <- c("#E69F00", "#56B4E9")
Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=reorder_within(Compound, beta, Outcome), 
                      y=beta, ymin=low, ymax=high, color = Model),
                  position = position_dodge(width = 0.6)) +
  facet_wrap(~Outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_x_reordered() +
  ggtitle("Model comparison: Linear vs. Negative binomial") +
  scale_color_brewer(palette="Set2") 

ggsave("Figure-SI-Comparison-linear-negbin.png", height=10, width=14)


# Correlation plot between linear and negbin models
Mbs %>%
  dplyr::select(Outcome, Compound, Model, beta) %>%
  pivot_wider(names_from = Model, values_from = beta) %>%
  ggplot(aes(x = Linear, y = `Negative binomial`)) +
  ggpubr::stat_cor(aes(label = after_stat(r.label))) +
  geom_point() +
  facet_wrap(~ Outcome) +
  theme_bw()

ggsave("Correlation-coefs-linear-negbin.png", height=6, width=10)

Mbs %>%
  dplyr::select(Outcome, Compound, Model, beta) %>%
  pivot_wider(names_from = Model, values_from = beta) %>%
  ggplot(aes(x = `Linear (before)`, y = `Linear`)) +
  ggpubr::stat_cor(aes(label = after_stat(r.label))) +
  geom_point() +
  facet_wrap(~ Outcome)

coefs.single <- Mbs %>%
  dplyr::select(Outcome, Compound, Model, beta) %>%
  pivot_wider(names_from = Model, values_from = beta) %>%
  mutate(model = "single-adjusted",
         sex = "both")

coefs.single

#############################################################################################

#------------------------------------------------------------------------------- MAIN MODELS

#------------------------------------------------------------------------------- ONE EXPOSURE

#--------------------------------------------------- Linear model with sex adjustment

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
  for (c in 1:nC) {
    C <- compound.labels[c]
    message(paste0("Compound: ", C))
    d.now <- d %>%  
      filter(compound == C) %>%
      filter(scale == O) %>%
      dplyr::select(-nid, -compound, -scale)
    
    model <- bayesglm(value ~ escalamenta(edatm) + EducationHigh + HouseholdDivided +
                        HasSiblings + Female + escalamenta(age) + 
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(concentration)),
                      data = d.now, family = gaussian())

    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    M[o,c,2] <- m$Value[1] # beta-coef
    M[o,c,1] <- m$LowOuter[1] # low CI (2 sigma)
    M[o,c,3] <- m$HighOuter[1] # high CI (2 sigma)
    M[o,c,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
  }
}
M

Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-1pollutant-sex-adjusted.csv")

# Plot
Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=reorder_within(compound, beta, outcome), 
                      y=beta, ymin=low, ymax=high)) +
  facet_wrap(~outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_x_reordered() +
  ggtitle("Linear Model (single exposure): Adjusted by Sex")
ggsave("250503-linear-1pollutant-sex-adjusted.pdf", width=12, height=10)

beta.lin.adj <- Mbs %>%
  dplyr::select(outcome, compound, beta) %>%
  mutate(model = "linear", 
         sex = "both")

#--------------------------------------------- Negative Binomial model with sex adjustment

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),                           compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
  for (c in 1:nC) {
    C <- compound.labels[c]
    message(paste0("Compound: ", C))
    d.now <- d %>%  
      filter(compound == C) %>%
      filter(scale == O) %>%
      dplyr::select(-nid, -compound, -scale)
    
    model <- glm.nb(value ~ escalamenta(edatm) + EducationHigh + HouseholdDivided +
                        HasSiblings + Female + escalamenta(age) + 
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(concentration)),
                      data = d.now)
    
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    M[o,c,2] <- m$Value[1] # beta-coef
    M[o,c,1] <- m$LowOuter[1] # low CI (2 sigma)
    M[o,c,3] <- m$HighOuter[1] # high CI (2 sigma)
    M[o,c,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
  }
}
M

Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)

write.csv(table, "250503-negbin-1pollutant-sex-adjusted.csv")

# Plot
Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=reorder_within(compound, beta, outcome), 
                      y=beta, ymin=low, ymax=high)) +
  facet_wrap(~outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_x_reordered() +
  ggtitle("Negative Binomial Model (single exposure): Adjusted by Sex")

ggsave("250503-negbin-1pollutant-sex-adjusted.pdf", width=12, height=10)


beta.nb.adj <- Mbs %>%
  dplyr::select(outcome, compound, beta) %>%
  mutate(model = "NB", 
         sex = "both")

#--------------------------------------------------------- Linear Model stratified by sex

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

# Sex
sexes <- unique(as.character(d$sex))
nS <- length(unique(d$sex))

M <- array(NA, dim=c(nO, nC, nS, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           sex = sexes, 
                           parameter = c("low", "beta", "high", "pval")))


for (s in 1:nS) {
  S <- sexes[s]
  for (o in 1:nO) {
    O <- outcomes[o]
    for (c in 1:nC) {
      C <- compound.labels[c]
      d.now <- d %>%
        filter(compound == C) %>%
        filter(scale == O) %>%
        filter(sex == S) %>%
        dplyr::select(-nid, -compound, -scale, -sex)
      model <- bayesglm(value ~ escalamenta(edatm) + EducationHigh + 
                          HouseholdDivided + HasSiblings + 
                          escalamenta(age) + 
                          escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                          escalamenta(log(concentration)),
                        data = d.now, family = gaussian())
      print(C)
      print(AIC(model))
      print(summary(model))
      m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                   outerCI = 1.96,
                                   intercept = FALSE)
      sum <- summary(model)
      
      M[o,c,s,2] <- m$Value[1] # beta-coef
      M[o,c,s,1] <- m$LowOuter[1] # low CI (2 sigma)
      M[o,c,s,3] <- m$HighOuter[1] # high CI (2 sigma)
      M[o,c,s,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
    }
  }
}
M

Mb <- tbl_df(as.data.frame.table(M)) 
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table 
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-1pollutant-sex-stratified.csv")

# Plot
cbPalette <- c("#E69F00","#56B4E9")
ggplot(Mbs, aes(x=compound, 
                y=beta, ymin=low, ymax=high, color=sex)) +
  geom_pointrange(size=0.3, position=position_dodge(width=0.5)) +
  facet_wrap(~outcome, ncol=3, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_color_brewer(palette="Set2") +
  ggtitle("Linear Model (single exposure): Sex Stratification")
ggsave("250503-linear-1pollutant-sex-stratified.pdf", width=12, height=10)

beta.lin.str <- Mbs %>%
  dplyr::select(outcome, compound, beta, sex) %>%
  mutate(model = "linear")

#------------------------------------------------- Negative Binomial Model stratified by sex

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

# Sex
sexes <- unique(as.character(d$sex))
nS <- length(unique(d$sex))

M <- array(NA, dim=c(nO, nC, nS, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                                      # "Fluid IQ", "Crystalized IQ",
                           compound = compound.labels,
                           sex = sexes, 
                           parameter = c("low", "beta", "high", "pval")))


for (s in 1:nS) {
  S <- sexes[s]
  for (o in 1:nO) {
    O <- outcomes[o]
    for (c in 1:nC) {
      C <- compound.labels[c]
      d.now <- d %>%
        filter(compound == C) %>%
        filter(scale == O) %>%
        filter(sex == S) %>%
        dplyr::select(-nid, -compound, -scale, -sex)
      model <- glm.nb(value ~ escalamenta(edatm) + EducationHigh + 
                          HouseholdDivided + HasSiblings + 
                          escalamenta(age) + 
                          escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                          escalamenta(log(concentration)),
                        data = d.now)
      print(C)
      print(AIC(model))
      print(summary(model))
      # print(summary(model.int))
      m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                   outerCI = 1.96,
                                   intercept = FALSE)
      sum <- summary(model)
      
      M[o,c,s,2] <- m$Value[1] # beta-coef
      M[o,c,s,1] <- m$LowOuter[1] # low CI (2 sigma)
      M[o,c,s,3] <- m$HighOuter[1] # high CI (2 sigma)
      M[o,c,s,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
    }
  }
}


Mb <- tbl_df(as.data.frame.table(M)) 
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-negbin-1pollutant-sex-stratified.csv")

# Plot
cbPalette <- c("#E69F00","#56B4E9")
ggplot(Mbs, aes(x=compound, 
                y=beta, ymin=low, ymax=high, color=sex)) +
  geom_pointrange(size=0.3, position=position_dodge(width=0.5)) +
  facet_wrap(~outcome, ncol=3, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_color_brewer(palette="Set2") +
  ggtitle("Negative Binomial Model (single exposure): Sex Stratification")
ggsave("250503-negbin-1pollutant-sex-stratified.pdf", width=12, height=10)


beta.nb.str <- Mbs %>%
  dplyr::select(outcome, compound, beta, sex) %>%
  mutate(model = "NB")


#############################################################################################

# Comparison of Linear vs. NB models in single-pollutant (adjusted and stratified)

beta.lin.adj
beta.lin.str
beta.nb.adj
beta.nb.str

single.model <- bind_rows(
  beta.lin.adj, beta.lin.str, beta.nb.adj, beta.nb.str) %>%
  spread(model, beta)

single.model %>%
  rename(Approach = sex) %>%
  mutate(Approach = ifelse(Approach == "both", "Sex-adjusted",
                        ifelse(Approach == "female", "Sex-stratified (F)", 
                               "Sex-stratified (M)"))) %>%
  ggplot(aes(x=linear, y=NB)) +
  geom_point(aes(color=Approach), alpha=0.6) +
  ggpubr::stat_cor(aes(label=after_stat(r.label))) +
  facet_wrap(~outcome, scales="free") +
  theme_bw() +
  theme(legend.position="bottom") +
  #theme_bw() +
  ylab("Negative Binomial Model") +
  xlab("Multivariable Linear Model") +
  scale_color_manual(values=c("#999888", "#E69F00", "#56B4E9")) +
  ggtitle("A. Single exposure models")
ggsave("250505-Comparison-NB-Linear-coefs.png", height = 6, width=8)


#############################################################################################

#----------------------------------------------------------------------------- ALL EXPOSURES
d.wide <- d %>%
  dplyr::select(-`creatinine in g/l`, -concentration.cr) %>%
  spread(compound, concentration) %>%
  unique()

d.wide

#--------------------------------------------------- Linear model with sex adjustment

# Compounds (all together)
nC <- 12
compound.labels <- c("MMP", "MEP", "MBzP", 
                     "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                     "BPA", "BPF", "SumDINCH", "SumDEHTP")
# Outcomes
outcomes <- unique(d.wide$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
    d.now <- d.wide %>%  
      filter(scale == O) %>%
      unique() %>%
      dplyr::select(-nid, -scale)
    
    model <- bayesglm(value ~ escalamenta(edatm) + EducationHigh + HouseholdDivided +
                        HasSiblings + Female + escalamenta(age) + 
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(SumDEHTP)) + escalamenta(log(SumDINCH)) +
                        escalamenta(log(BPF)) +
                        escalamenta(log(BPA)) + escalamenta(log(SumDiBP)) +
                        escalamenta(log(SumDnBP)) + escalamenta(log(SumDEHP)) +
                        escalamenta(log(SumDiNP)) + escalamenta(log(SumDiDP)) +
                        escalamenta(log(MBzP)) + escalamenta(log(MEP)) +
                        escalamenta(log(MMP)),
                      data = d.now, family = gaussian())
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    M[o,1:12,2] <- m$Value[1:12] # beta-coef
    M[o,1:12,1] <- m$LowOuter[1:12] # low CI (2 sigma)
    M[o,1:12,3] <- m$HighOuter[1:12] # high CI (2 sigma)
    M[o,1:12,4] <- rev(sum$coefficients[,4])[1:12] # p-value (last covariate from summary(model))
 # }
}
M

Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-joint-sex-adjusted.csv")

# Plot
Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=compound, 
                      y=beta, ymin=low, ymax=high)) +
  facet_wrap(~outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  ggtitle("Linear Model (joint exposure): Adjusted by Sex")
ggsave("250503-linear-joint-sex-adjusted.pdf", width=12, height=10)


beta.linear.adj <- Mbs %>%
  dplyr::select(outcome, compound, beta) %>%
  mutate(sex = "both",
         model = "linear")

#--------------------------------------------- Negative Binomial model with sex adjustment

# Compounds (all together)
nC <- 12
compound.labels <- c("MMP", "MEP", "MBzP", 
                     "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                     "BPA", "BPF", "SumDINCH", "SumDEHTP")
# Outcomes
outcomes <- unique(d.wide$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
    d.now <- d.wide %>%  
      filter(scale == O) %>%
      unique() %>%
      dplyr::select(-nid, -scale)
    
    model <- glm.nb(value ~ escalamenta(edatm) + EducationHigh + HouseholdDivided +
                        HasSiblings + Female + escalamenta(age) + 
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(SumDEHTP)) + escalamenta(log(SumDINCH)) +
                        escalamenta(log(BPF)) +
                        escalamenta(log(BPA)) + escalamenta(log(SumDiBP)) +
                        escalamenta(log(SumDnBP)) + escalamenta(log(SumDEHP)) +
                        escalamenta(log(SumDiNP)) + escalamenta(log(SumDiDP)) +
                        escalamenta(log(MBzP)) + escalamenta(log(MEP)) +
                        escalamenta(log(MMP)),
                      data = d.now)
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    M[o,1:12,2] <- m$Value[1:12] # beta-coef
    M[o,1:12,1] <- m$LowOuter[1:12] # low CI (2 sigma)
    M[o,1:12,3] <- m$HighOuter[1:12] # high CI (2 sigma)
    M[o,1:12,4] <- rev(sum$coefficients[,4])[1:12] # p-value (last covariate from summary(model))
  }

Mb <- tbl_df(as.data.frame.table(M)) 
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-negbin-joint-sex-adjusted.csv")

# Plot
Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=compound, 
                      y=beta, ymin=low, ymax=high)) +
  facet_wrap(~outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  ggtitle("Negative Binomial Model (joint exposure): Adjusted by Sex")
ggsave("250503-negbin-joint-sex-adjusted.pdf", width=12, height=10)


beta.nb.adj <- Mbs %>%
  dplyr::select(outcome, compound, beta) %>%
  mutate(sex = "both",
         model = "NB")

#--------------------------------------------------------- Linear Model stratified by sex

# Compounds
nC <- 12
compound.labels <- c("MMP", "MEP", "MBzP", 
                     "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                     "BPA", "BPF", "SumDINCH", "SumDEHTP")
# Outcomes
outcomes <- unique(d.wide$scale)
nO <- length(outcomes)

# Sex
sexes <- unique(as.character(d$sex))
nS <- length(unique(d$sex))

M <- array(NA, dim=c(nO, nC, nS, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           sex = sexes, 
                           parameter = c("low", "beta", "high", "pval")))


for (s in 1:nS) {
  S <- sexes[s]
  for (o in 1:nO) {
    O <- outcomes[o]
      d.now <- d.wide %>%
        filter(scale == O) %>%
        filter(sex == S) %>%
        dplyr::select(-nid, -scale, -sex)
      model <- bayesglm(value ~ escalamenta(edatm) + EducationHigh + 
                          HouseholdDivided + HasSiblings + 
                          escalamenta(age) + 
                          escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                          escalamenta(log(SumDEHTP)) + escalamenta(log(SumDINCH)) +
                          escalamenta(log(BPF)) +
                          escalamenta(log(BPA)) + escalamenta(log(SumDiBP)) +
                          escalamenta(log(SumDnBP)) + escalamenta(log(SumDEHP)) +
                          escalamenta(log(SumDiNP)) + escalamenta(log(SumDiDP)) +
                          escalamenta(log(MBzP)) + escalamenta(log(MEP)) +
                          escalamenta(log(MMP)),
                        data = d.now, family = gaussian())
      
      print(C)
      print(AIC(model))
      print(summary(model))
      m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                   outerCI = 1.96,
                                   intercept = FALSE)
      sum <- summary(model)
      
      M[o,1:12,s,2] <- m$Value[1:12] # beta-coef
      M[o,1:12,s,1] <- m$LowOuter[1:12] # low CI (2 sigma)
      M[o,1:12,s,3] <- m$HighOuter[1:12] # high CI (2 sigma)
      M[o,1:12,s,4] <- rev(sum$coefficients[,4])[1:12] # p-value (last covariate from summary(model))
  }
}

Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-joint-sex-stratified.csv")

# Plot
cbPalette <- c("#E69F00","#56B4E9")
ggplot(Mbs, aes(x=compound, 
                y=beta, ymin=low, ymax=high, color=sex)) +
  geom_pointrange(size=0.3, position=position_dodge(width=0.5)) +
  facet_wrap(~outcome, ncol=3, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  #geom_hline(yintercept = 1, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_color_brewer(palette="Set2") +
  ggtitle("Linear Model (single exposure): Sex Stratification")
ggsave("250503-linear-joint-sex-stratified.png", width=12, height=10)

beta.linear.str <- Mbs %>%
  dplyr::select(outcome, compound, beta, sex) %>%
  mutate(model = "linear")

#------------------------------------------------- Negative Binomial Model stratified by sex

# Compounds
nC <- 12
compound.labels <- c("MMP", "MEP", "MBzP", 
                     "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                     "BPA", "BPF", "SumDINCH", "SumDEHTP")
# Outcomes
outcomes <- unique(d.wide$scale)
nO <- length(outcomes)

# Sex
sexes <- unique(as.character(d$sex))
nS <- length(unique(d$sex))

M <- array(NA, dim=c(nO, nC, nS, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           sex = sexes, 
                           parameter = c("low", "beta", "high", "pval")))

for (s in 1:nS) {
  S <- sexes[s]
  for (o in 1:nO) {
    O <- outcomes[o]
      d.now <- d.wide %>%
      filter(scale == O) %>%
      filter(sex == S) %>%
      dplyr::select(-nid, -scale, -sex)
    model <- bayesglm(value ~ escalamenta(edatm) + EducationHigh + 
                        HouseholdDivided + HasSiblings + 
                        escalamenta(age) + 
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(SumDEHTP)) + escalamenta(log(SumDINCH)) +
                        escalamenta(log(BPF)) +
                        escalamenta(log(BPA)) + escalamenta(log(SumDiBP)) +
                        escalamenta(log(SumDnBP)) + escalamenta(log(SumDEHP)) +
                        escalamenta(log(SumDiNP)) + escalamenta(log(SumDiDP)) +
                        escalamenta(log(MBzP)) + escalamenta(log(MEP)) +
                        escalamenta(log(MMP)),
                      data = d.now, family = gaussian())
    
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    
    M[o,1:12,s,2] <- m$Value[1:12] # beta-coef
    M[o,1:12,s,1] <- m$LowOuter[1:12] # low CI (2 sigma)
    M[o,1:12,s,3] <- m$HighOuter[1:12] # high CI (2 sigma)
    M[o,1:12,s,4] <- rev(sum$coefficients[,4])[1:12] # p-value (last covariate from summary(model))
  }
}

Mb <- tbl_df(as.data.frame.table(M)) 
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-negbin-joint-sex-stratified.csv")

# Plot
cbPalette <- c("#E69F00","#56B4E9")
ggplot(Mbs, aes(x=compound, 
                y=beta, ymin=low, ymax=high, color=sex)) +
  geom_pointrange(size=0.3, position=position_dodge(width=0.5)) +
  facet_wrap(~outcome, ncol=3, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_color_brewer(palette="Set2") +
  ggtitle("Linear Model (single exposure): Sex Stratification")
ggsave("250503-negbin-joint-sex-stratified.pdf", width=12, height=10)

beta.nb.str <- Mbs %>%
  dplyr::select(outcome, compound, beta, sex) %>%
  mutate(model = "NB")


#############################################################################################

# Comparison of Linear vs. NB models in joint-pollutant (adjusted and stratified)

beta.linear.adj
beta.linear.str
beta.nb.adj
beta.nb.str

joint.model <- bind_rows(
  beta.linear.adj, beta.linear.str, beta.nb.adj, beta.nb.str) %>%
  spread(model, beta)

joint.model %>%
  rename(Approach = sex) %>%
  mutate(Approach = ifelse(Approach == "both", "Sex-adjusted",
                        ifelse(Approach == "female", "Sex-stratified (F)", 
                               "Sex-stratified (M)"))) %>%
  ggplot(aes(x=linear, y=NB)) +
  geom_point(aes(color=Approach), alpha=0.6) +
  ggpubr::stat_cor(aes(label=after_stat(r.label))) +
  facet_wrap(~outcome, scales="free") +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Negative Binomial Model") +
  xlab("Multivariable Linear Model") +
  scale_color_manual(values=c("#999888", "#E69F00", "#56B4E9")) +
  ggtitle("B. Joint exposures models")
ggsave("250505-Comparison-NB-Linear-coefs-Joint.png", height = 6, width=8)


####################################################################### SENSITIVITY ANALYSES

# Since we have demonstrated that NB and linear are the same, the sensitivity will be only
# performed on linear

#------------------------------------------------------------------------------- ONE EXPOSURE

#---------------------------------------------- Linear model with sex adjustment SENSITIVITY
#----------------------------------------------------------------------- ALL COVARS (8 + 4)

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
  for (c in 1:nC) {
    C <- compound.labels[c]
    message(paste0("Compound: ", C))
    d.now <- d %>%  
      filter(compound == C) %>%
      filter(scale == O) %>%
      dplyr::select(-nid, -compound, -scale)
    
    model <- bayesglm(value ~ escalamenta(edatm) + EducationHigh + 
                        SESlessAfluent +
                        HouseholdDivided +
                        TraumaticEvents +
                        HasSiblings + Female + escalamenta(age) + 
                        AgeSchoolOld + ResidenceUrban +
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(concentration)),
                      data = d.now, family = gaussian())
    
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    M[o,c,2] <- m$Value[1] # beta-coef
    M[o,c,1] <- m$LowOuter[1] # low CI (2 sigma)
    M[o,c,3] <- m$HighOuter[1] # high CI (2 sigma)
    M[o,c,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
  }
}
M

Mb <- tbl_df(as.data.frame.table(M)) 
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-1pollutant-sex-adjusted-sensitivity-12covars.csv")

Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=reorder_within(compound, beta, outcome), 
                      y=beta, ymin=low, ymax=high)) +
  facet_wrap(~outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_x_reordered() +
  ggtitle("Linear Model (single exposure): Adjusted by Sex")
ggsave("250503-linear-1pollutant-sex-adjusted-sensitivity-12covars.pdf", width=12, height=10)


#--------------------------------------------- Linear Model stratified by sex SENSITIVITY
#----------------------------------------------------------------------- ALL COVARS (8 + 4)

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

# Sex
sexes <- unique(as.character(d$sex))
nS <- length(unique(d$sex))

M <- array(NA, dim=c(nO, nC, nS, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           # "Fluid IQ", "Crystalized IQ",
                           # "Mathematical skills", 
                           # "Cognition", 
                           # "Psychomotor skills", 
                           # "Language skills"),
                           compound = compound.labels,
                           sex = sexes, 
                           parameter = c("low", "beta", "high", "pval")))

for (s in 1:nS) {
  S <- sexes[s]
  for (o in 1:nO) {
    O <- outcomes[o]
    for (c in 1:nC) {
      C <- compound.labels[c]
      d.now <- d %>%
        filter(compound == C) %>%
        filter(scale == O) %>%
        filter(sex == S) %>%
        dplyr::select(-nid, -compound, -scale, -sex)
      
      model <- bayesglm(value ~ escalamenta(edatm) + EducationHigh + 
                          SESlessAfluent +
                          HouseholdDivided +
                          TraumaticEvents +
                          HasSiblings + escalamenta(age) + 
                          AgeSchoolOld + ResidenceUrban +
                          escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                          escalamenta(log(concentration)),
                        data = d.now, family = gaussian())
      
      print(C)
      print(AIC(model))
      print(summary(model))
      m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                   outerCI = 1.96,
                                   intercept = FALSE)
      sum <- summary(model)
      
      M[o,c,s,2] <- m$Value[1] # beta-coef
      M[o,c,s,1] <- m$LowOuter[1] # low CI (2 sigma)
      M[o,c,s,3] <- m$HighOuter[1] # high CI (2 sigma)
      M[o,c,s,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
    }
  }
}

Mb <- tbl_df(as.data.frame.table(M)) 
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-1pollutant-sex-stratified-sensitivity-12covars.csv")

# Plot
cbPalette <- c("#E69F00","#56B4E9")
ggplot(Mbs, aes(x=compound, 
                y=beta, ymin=low, ymax=high, color=sex)) +
  geom_pointrange(size=0.3, position=position_dodge(width=0.5)) +
  facet_wrap(~outcome, ncol=3, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_color_brewer(palette="Set2") +
  ggtitle("Linear Model (single exposure): Sex Stratification")
ggsave("250503-linear-1pollutant-sex-stratified-sensitivity-12covars.pdf", width=12, height=10)


#---------------------------------------------- Linear model with sex adjustment SENSITIVITY
#----------------------------------------------------------------- REDUCED COVARS (5 + 4)

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
  for (c in 1:nC) {
    C <- compound.labels[c]
    message(paste0("Compound: ", C))
    d.now <- d %>%  
      filter(compound == C) %>%
      filter(scale == O) %>%
      dplyr::select(-nid, -compound, -scale)
    
    model <- bayesglm(value ~ escalamenta(edatm) + #EducationHigh + 
                        SESlessAfluent +
                        #HouseholdDivided +
                        TraumaticEvents +
                        HasSiblings + Female + #escalamenta(age) + 
                        AgeSchoolOld + ResidenceUrban +
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(concentration)),
                      data = d.now, family = gaussian())
    
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    M[o,c,2] <- m$Value[1] # beta-coef
    M[o,c,1] <- m$LowOuter[1] # low CI (2 sigma)
    M[o,c,3] <- m$HighOuter[1] # high CI (2 sigma)
    M[o,c,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
  }
}

Mb <- tbl_df(as.data.frame.table(M)) 
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-1pollutant-sex-adjusted-sensitivity-reduced-9covars.csv")

# Plot
Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=reorder_within(compound, beta, outcome), 
                      y=beta, ymin=low, ymax=high)) +
  facet_wrap(~outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_x_reordered() +
  ggtitle("Linear Model (single exposure): Adjusted by Sex")
ggsave("250503-linear-1pollutant-sex-adjusted-sensitivity-reduced-9covars.pdf", width=12, height=10)


#--------------------------------------------- Linear Model stratified by sex SENSITIVITY
#----------------------------------------------------------------- REDUCED COVARS (5 + 4)

# Compounds
nC <- length(unique(d$compound))
compound.labels <- unique(d$compound)

# Outcomes
outcomes <- unique(d$scale)
nO <- length(outcomes)

# Sex
sexes <- unique(as.character(d$sex))
nS <- length(unique(d$sex))

M <- array(NA, dim=c(nO, nC, nS, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           sex = sexes, 
                           parameter = c("low", "beta", "high", "pval")))


for (s in 1:nS) {
  S <- sexes[s]
  for (o in 1:nO) {
    O <- outcomes[o]
    for (c in 1:nC) {
      C <- compound.labels[c]
      d.now <- d %>%
        filter(compound == C) %>%
        filter(scale == O) %>%
        filter(sex == S) %>%
        dplyr::select(-nid, -compound, -scale, -sex)
      
      model <- bayesglm(value ~ escalamenta(edatm) + #EducationHigh + 
                          SESlessAfluent +
                          #HouseholdDivided +
                          TraumaticEvents +
                          HasSiblings + #escalamenta(age) + 
                          AgeSchoolOld + ResidenceUrban +
                          escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                          escalamenta(log(concentration)),
                        data = d.now, family = gaussian())
      
      print(C)
      print(AIC(model))
      print(summary(model))
      # print(summary(model.int))
      m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                   outerCI = 1.96,
                                   intercept = FALSE)
      sum <- summary(model)
      
      M[o,c,s,2] <- m$Value[1] # beta-coef
      M[o,c,s,1] <- m$LowOuter[1] # low CI (2 sigma)
      M[o,c,s,3] <- m$HighOuter[1] # high CI (2 sigma)
      M[o,c,s,4] <- rev(sum$coefficients[,4])[1] # p-value (last covariate from summary(model))
    }
  }
}


Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250503-linear-1pollutant-sex-stratified-sensitivity-reduced-9covars.csv")

# Plot
cbPalette <- c("#E69F00","#56B4E9")
ggplot(Mbs, aes(x=compound, 
                y=beta, ymin=low, ymax=high, color=sex)) +
  geom_pointrange(size=0.3, position=position_dodge(width=0.5)) +
  facet_wrap(~outcome, ncol=3, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_color_brewer(palette="Set2") +
  ggtitle("Linear Model (single exposure): Sex Stratification")
ggsave("250503-linear-1pollutant-sex-stratified-sensitivity-reduced-9covars.pdf", width=12, height=10)


#------------------------------------------------------------------------------- ALL EXPOSURES

#--------------------------------------------- Linear Model adjusted by sex SENSITIVITY
#----------------------------------------------------------------- REDUCED COVARS (5 + 4)

# Compounds (all together)
nC <- 12
compound.labels <- c("MMP", "MEP", "MBzP", 
                     "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                     "BPA", "BPF", "SumDINCH", "SumDEHTP")

# Outcomes
outcomes <- unique(d.wide$scale)
nO <- length(outcomes)

M <- array(NA, dim=c(nO, nC, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           parameter = c("low", "beta", "high", "pval")))

for (o in 1:nO) {
  O <- outcomes[o]
  message(paste0("Outcome: ", O))
  d.now <- d.wide %>%  
    filter(scale == O) %>%
    unique() %>%
    dplyr::select(-nid, -scale)
  
  model <- bayesglm(value ~ escalamenta(edatm) + #EducationHigh + 
                      SESlessAfluent +
                      #HouseholdDivided +
                      TraumaticEvents +
                      HasSiblings + Female + #escalamenta(age) + 
                      AgeSchoolOld + ResidenceUrban +
                      escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                      escalamenta(log(SumDEHTP)) + escalamenta(log(SumDINCH)) +
                      escalamenta(log(BPF)) +
                      escalamenta(log(BPA)) + escalamenta(log(SumDiBP)) +
                      escalamenta(log(SumDnBP)) + escalamenta(log(SumDEHP)) +
                      escalamenta(log(SumDiNP)) + escalamenta(log(SumDiDP)) +
                      escalamenta(log(MBzP)) + escalamenta(log(MEP)) +
                      escalamenta(log(MMP)),
                    data = d.now, family = gaussian())
  print(C)
  print(AIC(model))
  print(summary(model))
  m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                               outerCI = 1.96,
                               intercept = FALSE)
  sum <- summary(model)
  M[o,1:12,2] <- m$Value[1:12] # beta-coef
  M[o,1:12,1] <- m$LowOuter[1:12] # low CI (2 sigma)
  M[o,1:12,3] <- m$HighOuter[1:12] # high CI (2 sigma)
  M[o,1:12,4] <- rev(sum$coefficients[,4])[1:12] # p-value (last covariate from summary(model))
  # }
}
M

Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250508-linear-joint-sex-adjusted-sensitivity-reduced-covars.csv")

# Plot
Mbs %>%
  ggplot() +
  geom_pointrange(aes(x=compound, 
                      y=beta, ymin=low, ymax=high)) +
  facet_wrap(~outcome, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  ggtitle("Linear Model (joint exposure): Adjusted by Sex")
#ggsave("250503-linear-joint-sex-adjusted.pdf", width=12, height=10)

#--------------------------------------------------------- Linear Model stratified by sex

# Compounds
nC <- 12
compound.labels <- c("MMP", "MEP", "MBzP", 
                     "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                     "BPA", "BPF", "SumDINCH", "SumDEHTP")

# Outcomes
outcomes <- unique(d.wide$scale)
nO <- length(outcomes)

# Sex
sexes <- unique(as.character(d$sex))
nS <- length(unique(d$sex))

M <- array(NA, dim=c(nO, nC, nS, 4),
           dimnames = list(outcome = c("Conduct problems", "Emotional symptoms", 
                                       "Hyperactivity/Inattention",
                                       "Peer relationships problems",
                                       "Prosocial behavior", "Total difficulties",
                                       "Internalizing score", "Externalizing score"),
                           compound = compound.labels,
                           sex = sexes, 
                           parameter = c("low", "beta", "high", "pval")))

for (s in 1:nS) {
  S <- sexes[s]
  for (o in 1:nO) {
    O <- outcomes[o]
    d.now <- d.wide %>%
      filter(scale == O) %>%
      filter(sex == S) %>%
      dplyr::select(-nid, -scale, -sex)
    
    model <- bayesglm(value ~ escalamenta(edatm) + #EducationHigh + 
                        SESlessAfluent +
                        #HouseholdDivided + 
                        TraumaticEvents +
                        HasSiblings + #escalamenta(age) + 
                        AgeSchoolOld + ResidenceUrban +
                        escalamenta(bmi) + escalamenta(log(cotinine.7y)) +
                        escalamenta(log(SumDEHTP)) + escalamenta(log(SumDINCH)) +
                        escalamenta(log(BPF)) +
                        escalamenta(log(BPA)) + escalamenta(log(SumDiBP)) +
                        escalamenta(log(SumDnBP)) + escalamenta(log(SumDEHP)) +
                        escalamenta(log(SumDiNP)) + escalamenta(log(SumDiDP)) +
                        escalamenta(log(MBzP)) + escalamenta(log(MEP)) +
                        escalamenta(log(MMP)),
                      data = d.now, family = gaussian())
    
    print(C)
    print(AIC(model))
    print(summary(model))
    m <- coefplot:::buildModelCI(model, innerCI = 1.64,
                                 outerCI = 1.96,
                                 intercept = FALSE)
    sum <- summary(model)
    
    M[o,1:12,s,2] <- m$Value[1:12] # beta-coef
    M[o,1:12,s,1] <- m$LowOuter[1:12] # low CI (2 sigma)
    M[o,1:12,s,3] <- m$HighOuter[1:12] # high CI (2 sigma)
    M[o,1:12,s,4] <- rev(sum$coefficients[,4])[1:12] # p-value (last covariate from summary(model))
  }
}


Mb <- tbl_df(as.data.frame.table(M))
Mbs <- spread(Mb, parameter, Freq)
unique(Mbs$compound)

# Table
table <- Mbs %>%
  mutate(ci = paste("[", round(low, 3), ";", round(high, 3), "]")) %>%
  mutate(beta.ci = paste(round(beta, 3), ci)) %>%
  dplyr::select(-low, -beta, -high, -ci) %>%
  mutate(pval.dic = ifelse(pval < 0.05, "**",
                           ifelse(pval >=0.1, "-", "*"))) %>%
  mutate(beta.ci.pval = paste(beta.ci, " ", pval.dic)) %>%
  dplyr::select(-pval, -beta.ci, -pval.dic) %>%
  spread(outcome, beta.ci.pval)
write.csv(table, "250508-linear-joint-sex-stratified-sensitivity-reduced-covars.csv")

# Plot
cbPalette <- c("#E69F00","#56B4E9")
ggplot(Mbs, aes(x=compound, 
                y=beta, ymin=low, ymax=high, color=sex)) +
  geom_pointrange(size=0.3, position=position_dodge(width=0.5)) +
  facet_wrap(~outcome, ncol=3, scales="free") +
  geom_hline(yintercept = 0, lty=3) +
  xlab("") +
  ylab(expression(paste("Standardized ", beta, "-coefficients"))) +
  theme_classic() +
  coord_flip() +
  scale_color_brewer(palette="Set2") +
  ggtitle("Linear Model (single exposure): Sex Stratification")
#ggsave("250503-linear-joint-sex-stratified.png", width=12, height=10)

