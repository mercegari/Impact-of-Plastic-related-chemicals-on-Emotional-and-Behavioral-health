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

################################################################################## 
# Reviewer comments:
# Reviewer 2
# Thank you for addressing my previous comments. I only have a remaining major comment:
# In the environmental health field, the most common mixture models include WQS, BWQS, 
# Quantile G-computation and BKMR. The authors mentioned they implemented a new model. 
# However, they needed to introduce adaptations since the package only supported binary 
# outcomes and not linear ones.
# The question is: why use this new model that needs adaptation (and that should ideally #
# be validated before implementation), when you can easily analyze linear outcomes 
# with Quantile G-comp 
# (https://cran.r-project.org/web/packages/qgcomp/vignettes/qgcomp-vignette.html) 
# or with BWQS (https://rdrr.io/github/ElenaColicino/bwqs/man/bwqs.html)?
# I believe you need to validate your current mixture results using another mixture #
# model that has been previously and widely implemented.
# "We have applied a Bayesian Grouped Weighted¡ Quantile Sum Regression Model, using the 
# R package BayesGWQS. This package allowed us to include up to three groups of chemicals, 
# namely phthalate metabolites (MMP, MEP, MBzP and sums of DEHP, DiDP, DiNP, DiBP and DnBP), 
# non-phthalate metabolites (sums of DINCH and DEHTP) and bisphenols (BPA and BPF).
# However, the aforementioned package was only available for binary outcomes, and therefore, 
# we have adapted it to support linear and negative binomial models, using a modified 
# JAGS code. Although this is not the scope of this work, we will contact the authors 
# of this package to provide these adaptations."


################################################################################## BWQS

devtools::install_github("ElenaColicino/bwqs", build_vignettes = TRUE)
library(BWQS)
library(rstan)

browseVignettes("BWQS")

library(MASS)
library(knitr)
library(clusterGeneration)
library(kableExtra)

#----------------------------------------------------------------------- Sex adjusted
# Run BWQS model using a loop for each SDQ scale (sex-adjusted models)

Scales <- unique(d$scale)
nS <- length(Scales)

SS <- NULL
for (s in 1:nS) {
  Scale <- Scales[s]
  d2 <- d %>%
    filter(scale == Scales[s]) %>%
    filter(compound %in% c("MMP", "MEP", "MBzP", 
                           "SumDiDP", "SumDiNP", "SumDEHP", "SumDnBP", "SumDiBP",
                           "SumDEHTP", "SumDINCH",
                           "BPA", "BPF")) %>%
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
  
  message(paste0("Processing ", Scale, ", with N = ", length(d2$value)))
  
  fit <- bwqs(value ~ edatm + EducationHigh + HouseholdDivided +
                    HasSiblings + age + Female + bmi + cotinine.7y, 
                  mix_name = c('MMP','MEP','MBzP',
                               'SumDiDP','SumDiNP', 'SumDEHP','SumDnBP','SumDiBP',
                               'SumDEHTP','SumDINCH','BPA','BPF'),
                  data = d2, q = 4, family = "gaussian")
  
  SS <- bind_rows(SS, tibble(as.data.frame(fit$summary_fit)[2,], Scale = Scale))

}

SS
save(SS, file="BQWS-sex-adjusted-output.RData")

ggplot(SS, aes(x=mean, y=reorder(Scale, mean), xmin=`2.5%`, xmax=`97.5%`)) +
  geom_point() +
  geom_pointrange(position = position_dodge(width=0.4)) +
  geom_vline(xintercept=0, lty=3) +
  ylab("") + 
  xlab(expression(paste("Standardized ", beta, "-coefficients (95% CI)"))) +
  theme_bw()
ggsave("BWQS-plot-sex-adjusted.pdf", height=4, width=6)
ggsave("BWQS-plot-sex-adjusted.png", height=4, width=6)
  
#----------------------------------------------------------------------- Sex stratified
# Run BWQS model using a loop for each SDQ scale and sex (sex-stratified models)

Scales <- unique(d$scale)
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
                             "BPA", "BPF")) %>%
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
    
  message(paste0("Processing ", Scale, G, ", with N = ", length(d2$value)))
  
  fit <- bwqs(value ~ edatm + EducationHigh + HouseholdDivided +
                HasSiblings + age + bmi + cotinine.7y, 
              mix_name = c('MMP','MEP','MBzP',
                           'SumDiDP','SumDiNP', 'SumDEHP','SumDnBP','SumDiBP',
                           'SumDEHTP','SumDINCH','BPA','BPF'),
              data = d2, q = 4, family = "gaussian")
  
  SS <- bind_rows(SS, tibble(as.data.frame(fit$summary_fit)[2,], Scale = Scale, Sex = factor(G)))
  }
}

SS
save(SS, file="BQWS-sex-stratified-output.RData")

ggplot(SS, aes(x=mean, y=Scale, #y=reorder(Scale, mean), 
               xmin=`2.5%`, xmax=`97.5%`, color=Sex)) +
  geom_point(position=position_dodge(width=0.4)) +
  geom_pointrange(position = position_dodge(width=0.4)) +
  geom_vline(xintercept=0, lty=3) +
  ylab("") + 
  xlab(expression(paste("Standardized ", beta, "-coefficients (95% CI)"))) +
  theme_bw() +
  theme(legend.position="bottom")
  scale_color_brewer(palette="Set2")
ggsave("BWQS-plot-sex-stratified.pdf", height=4, width=6)
ggsave("BWQS-plot-sex-stratified.png", height=4, width=6)


