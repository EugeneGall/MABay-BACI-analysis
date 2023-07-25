# Bayes_BACI_MA_Bay_3.R
# Code originally written by Chris Haak, UMass Amherst
# brms BACI analyses by Eugene.Gallagher@umb.edu 6/2/2022, revised 8/15/22, 
# 5/1/23, 7/21,23
## Haak source('C:/Data/THESIS_STATS_R/thesis/data/biostats.R')

# setwd("O:/MABay alpha Haak")
library(bayesplot)
library(BHH2)
library(brms)
library(car)
library(cowplot)
library(dplyr)
library(emmeans)
library(forcats)
library(gganimate)
library(ggdist)
library(ggplot2)
library(ggrepel)
library(Hmisc)
library(lattice)
library(lme4)
library(lmerTest)
library(loo)
library(magrittr)
library(MCMCglmm)
library(modelr)  # part of the tidyverse package, for finding p values
library(pbkrtest)
library(posterior)
library(purrr)
library(RColorBrewer)
library(rethinking)
library(rstan)
library(rstanarm)
library(sciplot)
library(tidybayes)
library(tidyr)
library(tidyverse)

# These options make rstan run faster, Gallagher's PC has 8 cores
# http://mjskay.github.io/tidybayes/articles/tidybayes.html
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load data as a tibble and check structure -------------------------------
MWRA<-read_csv("../data/MWRA_MA_92-08new5.csv")  # read in data as a tibble
MWRA
str(MWRA)
summary(MWRA)

# Convert some numerics to factors ----------------------------------------

MWRA$fYear<-as.factor(MWRA$year)
MWRA$Case<-as.factor(MWRA$Case)
MWRA$rep<-as.factor(MWRA$rep)
MWRA$period<-as.factor(MWRA$period)

# relevel the order of period, so pre precedes post in plots ---------------

MWRA$period <- relevel(MWRA$period, ref = "pre")

# Plot Fisher's log series alpha --------------------------------------

ggplot(data=MWRA) + 
  geom_histogram(mapping = aes(x=alpha), binwidth = 1)
ggplot(data=MWRA, mapping = aes (x=alpha, color=period)) +
  geom_freqpoly(binwidth = 1)
ggplot(data=MWRA, mapping = aes(x=period, y = alpha)) + geom_boxplot()

boxplot(alpha ~ region*period, data=MWRA) #heteroscedastic
leveneTest(alpha ~ region*period, data=MWRA) #heteroscedastic p = 7e-08


# Add log-transformed alpha -----------------------------------------------

MWRA$lnalpha<-log(MWRA$alpha)

# Check log-transformed response ------------------------------------------
ggplot(data=MWRA) + 
  geom_histogram(mapping = aes(x=lnalpha), binwidth = 0.1)
ggplot(data=MWRA, mapping = aes (x=lnalpha, color=period)) +
  geom_freqpoly(binwidth = 0.1)
ggplot(data=MWRA, mapping = aes(x=period, y = lnalpha)) + geom_boxplot()
boxplot(lnalpha ~ region*period, data=MWRA) #better
leveneTest(lnalpha ~ region*period, data=MWRA) #still heteroscedastic p = 2e-04

# Note Andrews says homoscedasticity is not a problem for mixed effects models


# log(Fisher's alpha) over time by station --------------------------------

xyplot(alpha ~ fYear | factor(Station),type = "o", group=MWRA$region,col=c(1,4),
       data=MWRA)

ggplot(MWRA,
       aes(x = year+8/12, y = lnalpha)) + geom_point() +
  stat_smooth(method = 'lm', se = FALSE) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station)
#no clear trend - Nearfield seems slightly more variable than farfield
#overall trend seems to be slightly upwards with time

ggplot(MWRA,
       aes(x = year+8/12, y = lnalpha)) + geom_point() +
  geom_smooth(data = MWRA, mapping = aes(x = year, y = lnalpha), se = TRUE) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station)
# curvilinearities evident in many plots: FF06, FF07. FF12, FF14... indicating
# that year should be either treated as a factor or a curvilinear regression
# should be used: beta spline or GAM. Also, the following stations should be
# deleted: FF01, NF01, NF03, NF06, and NF11.

# Use the tidyverse filter operator to delete 5 stations, could be better ways
mwra <- filter(MWRA, Station != 'FF01' & Station != 'NF01' & Station != 'NF03' &
               Station != 'NF06' & Station != 'NF11')
str(MWRA) 
str(mwra) # Number of samples deleted from 864 to 854.

# Now redo the ggplot facet plot
ggplot(mwra,
       aes(x = year+8/12, y = lnalpha)) + geom_point() +
  geom_smooth(data = mwra, mapping = aes(x = year, y = lnalpha), se = TRUE) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station)

## Nested ANOVA (fixed factors)
M_1<-lm(lnalpha~region*period + Station%in%region + year%in%period, 
                data=mwra)
summary(M_1)
# BACI effect (regionnear:periodpre) -0.053 se=0.023 p=0.024

# Add the final level of nesting: replicates within stations:
M_2<-lm(lnalpha~region*period + Station%in%region + year%in%period +
                  rep%in%Station, data=mwra)
summary(M_2)
anova(M_2)
# BACI effect (regionnear:periodpre) -0.056 se=0.023 p=0.018, but this is
# not an appropriate model, as it is testing the F with 1 & 777 df.

# Mixed effects ANOVA using lme4 (lmer function) --------------------------

# FULL MODEL (REGION + PERIOD + REGION:PERIOD INTERACTION),
# station random nested in region, year random nested in period)\
# Need to use restricted maximum likelihood for model comparison with Wilks's 
# theorem, i.e., the likelihood ratio test

M_bad<- lmer(lnalpha ~ region*period +(1|region:Station) + 
                       (1|period:fYear) + (1|Station:rep), 
                       data=mwra, REML=TRUE)
# Can't use the above model, the result is singular, I can fit with brms
# but the fit adds little at all to the model

M_3<- lmer(lnalpha ~ region*period +(1|region:Station) + 
             (1|period:fYear), 
             data=mwra, REML=TRUE)

summary(M_3) # BACI effect 0.053, se= 0.021, p =0.01317

#Based on t-values, interaction should be significant
#strong correlation of intercept w/ region, but this is not a problem
## Type III anova table with p-values for F-tests based on Satterthwaite's
## method:

add_predictions(mwra, M_3) %>% 
  ggplot(aes(x = year+8/12, y = lnalpha, shape = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model')

# Don't do separate facet plots for stations

add_predictions(mwra, M_3) %>% 
  ggplot(aes(x = year+8/12, y = lnalpha, shape = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  ggtitle('lme4 BACI Model')

# This is a key figure, with the CI for the smoother diverging.


# Add back facets colour by region
add_predictions(mwra, M_3) %>% 
  ggplot(aes(x = year+8/12, y = lnalpha, colour = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model')

# Plot a line plot, not a lowess smooth plot, note somewhat bad fit to NF08
add_predictions(mwra, M_3) %>% 
  ggplot(aes(x = year+8/12, y = lnalpha, colour = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model')

# Find the p value using Kenward-Rogers & Wilks's Theorem -----------------

## Choose type II anova table with Kenward-Roger method for the F-test:
if(requireNamespace("pbkrtest", quietly = TRUE))
  anova(M_3, type=2, ddf="Kenward-Roger")
# p = 0.01317

# Test with Wilks's theorem, Bates's & Andrews' recommendation; dropping
# only the BACI term region:period

M_4 <- lmer(lnalpha ~  1 + region + period + (1|region:Station) +
              (1|period:fYear), data=mwra, REML=TRUE)
summary(M_4)
anova(M_3,M_4)
# BACI effect p = 0.0131 using Wilks's Theorem (likelihood ratio test)


# Plot the results using the frequentist BACI model -----------------------

###plot means and bootstrapped .95 CIs (using Hmisc 1000 bootstrapped samples) 
# ?smean.cl.boot)

#define CI bootstrap function
# define getCI function for bootstrapped data using the smean.cl.boot from
# Mmisc package (default is 1000 boostrapped samples)
getCI<-function(x)smean.cl.boot(x)[-1]

P_2 <- lineplot.CI(year,fitted(M_3), group = region, 
            data = mwra, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE, , xlab="Year",
            ylab = "log(Fisher's alpha)")
            abline(v=9.5,lty=3,col="red")

P_2       
# Plot AMO  --------------------------------------------
            # Read in AMO monthly mean data (1992 - 2008) from
# https://psl.noaa.gov/data/climateindices/
            
amo<-read_csv("../data/AMO_Monthly_mean_1992_2008.csv")
glimpse(amo)
            
P_3<-ggplot(amo, aes(x = Year, y = AMO)) + geom_point() +
              geom_line(data = amo, mapping = aes(x = Year, y = rm_8)) +
              geom_vline (aes(xintercept = 2000+9.9/12))
P_3


# plot Aug1948-2022 -------------------------------------------------------

amo_48<-read_csv("../data/AMO_sm8_48_22.csv")
glimpse(amo_48)
P_3a<-ggplot(amo_48, aes(x = Year, y = AMO)) + geom_point() +
  geom_line(data = amo_48, mapping = aes(x = Year, y = rm_8)) +
  geom_vline (aes(xintercept = 2000+9.5/12))
P_3a
amo <- amo[1:204,1:4]

glimpse(mwra)

# Drop the individual points

P_4<-ggplot(amo, aes(x = Year, y = AMO)) + 
  geom_line(data = amo, mapping = aes(x = Year, y = rm_8))
P_4

# composite contains the far and near predictions by year
composite <- tibble(Year=seq(1992,2008),far=P_2$vals[,1],near=P_2$vals[,2])
P_5<-ggplot(composite, aes(x = Year, y = far)) +
  geom_line(data = composite, mapping = aes(x = Year, y = far)) +
  aes(x = Year, y = near) + 
  geom_line(data = composite, mapping = aes(x = Year, y = near)) +
  geom_vline (aes(xintercept = 2000+9.5/12))
P_5

library(patchwork)
P_4/P_5

# create a ggplot graphic for lnalpha for composite plot.

# Calculate the correlation between AMO and lnalpha
# 8-month average AMO (Jan-Aug) for Aug 1992-2008 from psl.noaa.gov
AMO_Aug<- c(-0.2035, -0.224375, -0.260625, 0.126375, -0.06325, -0.0135, 0.3485,
            0.1195, 0.020625, 0.0285, 0.033875, 0.144375, 0.161, 0.274, 0.218,
            0.116125, 0.13)

# Use dplyr pipe to extract mean for lnAug across all stations
# Analyze correlation with 8-mo smoothed AMO

lnalphaAug <- mwra %>%
  group_by(year) %>%
  summarise(mean = mean(lnalpha), n = n())
lnalphaAug$mean
cor(AMO_Aug, lnalphaAug$mean,  method = "pearson")
# r = 0.702

cordata <- tibble(AMO_Aug=AMO_Aug,lnalphaAug=lnalphaAug$mean)
cordata

ggplot(cordata, aes(x = AMO_Aug, y = lnalphaAug)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "August AMO (8 mo smoother)", y = "log(Fisher's alpha)")

# Fit the model

fit_AMO_alpha <- lm(lnalphaAug ~ AMO_Aug , data = cordata)
summary(fit_AMO_alpha)
plot(fit_AMO_alpha)


# Bayesian analysis with brms, including nested reps within stations-----------

M_5 <- brm(lnalpha ~ 1 + region*period +(1|region:Station) +
             (1|period:fYear) + (1|Station:rep), data=mwra)
M_5  
plot(M_5)


# There appears to be little value in the reps within station variance, the
# the posterior appears to include 0, so let's drop that term
# with that term the BACI effect was 0.05, se=0.02
# and credibility interval of 0.01 to 0.09

M_6 <- brm(lnalpha ~ region*period + (1|region:Station) +
             (1|period:fYear), data=mwra)
M_6
plot(M_6)
# M_6 is superior to M_5 in that it involves one less variance estimator, reps
# within stations, which had a credible interval including 0
# The BACI effect was -0.05, se=0.02 and credibility interval of -0.01 to -0.10

# Diagnostic plots for Bayesian fitted model ----------------------------------

# install.packages("bayesplot") # installed 6/13/22
?? bayesplot

# See the following reference on diagnostic plots for Bayesian regression
# https://towardsdatascience.com/evaluating-bayesian-mixed-models-in-r-python-27d344a03016
# evaluate simulated data with observed data using Bayesplot

# Use M_6 as the model to check for diagnostics

# posterior predicted samples, ndraws = no. of datasets to generate
pp_data <- posterior_predict(M_6, ndraws = 200) 

color_scheme_set("gray" )
pdense2 <- ppc_dens_overlay(y = mwra$lnalpha,
                            yrep = pp_data) +
  labs(x = "lnalpha", title = "BACI model") +
  theme(title = element_text(size=10))
pdense2
# Caption, 200 simulated data drawn from the posterior (light lines) fit the
# the observed data well. 

# Observed and simulated mean SKEW metrics  (1000 simulated datasets)
ppc_stat(y = mwra$lnalpha, 
         yrep = posterior_predict(M_6, ndraws = 1000),
         stat = "mean")

# Check skewness
#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  (1000 simulated datasets)
ppc_stat(y = mwra$lnalpha, 
         yrep = posterior_predict(M_6, ndraws = 1000),
         stat = "skew")

# observed and simulated MEDIAN metric per region (1000 simulated datasets) 
ppc_stat_grouped(y = mwra$lnalpha,
                 yrep = posterior_predict(M_6, ndraws = 500),
                 group = mwra$region, 
                 stat = "median")
# Very close, but median seems underestimated for near region

# Use LOO, leave one out validation:

# It's important to set save_psis as TRUE
model_loo <- loo(M_6, save_psis = TRUE)
w <- weights(model_loo$psis_object)

color_scheme_set("gray")
ppc_loo_pit_overlay(mwra$lnalpha, 
                    posterior_predict(M_6), 
                    lw = w)
# Should approximate a uniform distribution, but the maximum at 0.5 indicates
# a slight tendency to overdispersion.


# Look for influential single observation using loo -----------------------

k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Since all observations have a k̂ <0.5 then none are considered influential in
# the posterior given our model is able to characterize them correctly. 
# In other words, our model’s LOO-PSIS estimates are reliable.

# Overall Model comparison, with and without BACI term -------------------------

# fit model with brms without BACI term
M_7 <- brm(lnalpha ~ 1 + region + period +(1|region:Station) +
             (1|period:fYear), data=mwra)
M_7

# Takes about 1 minute
plot(M_7)


# Kfold CV (not run)
# kf <- kfold(bayes_rintercept)

# Compute LOO CV for each model
# base model with BACI term
loo1 <- loo(M_6, save_psis = TRUE, cores = 8)

# Model BACI term and reps within stations
loo2 <- loo(M_5, save_psis = TRUE, cores = 8) # includes nested reps within stations
loo2

# Reduced model with no BACI term
loo3 <- loo(M_7, save_psis = TRUE, cores = 8)
loo3

# Comparison
loo_compare(loo1, loo2, loo3)
#     elpd_diff se_diff
# M_6  0.0       0.0   
# M_5 -1.0       0.5   
# M_7 -1.6       2.6  

# M_6 provides the best fit, slightly better than M_5(reps within stations) and 
# M_7 (no BACI term)

#-----------------------------------------

#####  Pointwise comparison ------------

# Obtain pointwise ELPD values
elpd1 <- loo1$pointwise[,"elpd_loo"]  # Base model with BACI
elpd2 <- loo3$pointwise[,"elpd_loo"]  # no BACI term

# Build differences dataframe, looking at Stations
elpd_df <- tibble(Region = mwra$region,
                  diff12 = elpd2 - elpd1) %>% 
  mutate(idx = 1:n())

# Plot comparing with and without BACI model
ggplot(elpd_df, aes(x = idx, y = diff12, color = mwra$region)) +
  geom_point(alpha=0.7) +
  geom_hline(aes(yintercept=0)) +
  theme_bw()
# It's tough to see one model fitting the data better.


# Final plot with tidybayes from
# http://mjskay.github.io/tidybayes/articles/tidybayes.html


# from the tidybayes vignette
# http://mjskay.github.io/tidybayes/

# Need to drop the call to data_grid for multilevel model

# From: http://mjskay.github.io/tidybayes/articles/tidy-brms.html
mwra %>%
  group_by(region) %>%
  add_epred_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")
# 95% credible interval is too light to be seen


mwra %>%
  group_by(region) %>%
  add_epred_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")
# The key 95% credible interval is not visible.


# Plot posterior predictions instead of means -----------------------------

mwra %>%
  group_by(region) %>%
  add_predicted_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region, fill = region)) +
  stat_lineribbon(aes(y = .prediction), .width = .95, alpha = 1/4) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")
# This is an effective plot, in showing the subtle differences between regions


# facet plot --------------------------------------------------------------

mwra %>%
  group_by(region) %>%
  add_predicted_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region, fill = region)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5), color = brewer.pal(5, "Blues")[[5]]) +
  geom_point(data = mwra) +
  scale_fill_brewer() +
  facet_grid(. ~ period, space = "free_x", scales = "free_x")
# This plot is kind of a mess

mwra %>%
  group_by(region) %>%
  add_epred_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region, fill = region)) +
  stat_lineribbon(aes(y = .epred), .width = .95, color = brewer.pal(5, "Blues")[[5]]) +
  geom_point(data = mwra) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer() +
  facet_grid(. ~ period, space = "free_x", scales = "free_x")
# Can't distinguish far and near lines.

# from http://mjskay.github.io/tidybayes/articles/tidybayes.html
mwra %>%
  group_by(region) %>%
#  data_grid(hp = seq_range(hp, n = 51)) %>%
  add_epred_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")
# Can't see the 95% credible levels

mwra %>%
  group_by(region) %>%
  #  data_grid(hp = seq_range(hp, n = 51)) %>%
  add_epred_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .epred), .width = .95, alpha = 1/4) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set2")  
# tells the story, but with horrible Irish colors

# from https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html
mwra %>%
#  data_grid(agegp) %>%
#  add_epred_draws(M_6, dpar = TRUE, category = "tobgp") %>%
add_epred_draws(M_6, dpar = TRUE, category = "region") %>%
  ggplot(aes(x = year, y = .epred, color = region)) +
  stat_pointinterval(position = position_dodge(width = .4)) +
  geom_line(aes(y = .epred)) +
#  stat_lineribbon(aes(y = .epred), .width = .95) +
  scale_size_continuous(guide = "none") +
  xlab("\nYear") +
  ylab("log (Fisher's alpha)\n") +  
  geom_vline (aes(xintercept = 2000.25)) +
  scale_color_manual(values = brewer.pal(6, "Blues")[-c(1,2)])
# won't plot

# From https://ourcodingclub.github.io/tutorials/brms/
(model_fit <- mwra %>%
    add_predicted_draws(M_6) %>%  # adding the posterior distribution
    ggplot(aes(x = year, y = lnalpha, color = region)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = mwra, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    xlab("\nYear") +
    ylab("log (Fisher's alpha)\n") +  
    geom_vline (aes(xintercept = 2000.25)) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))
# Lousy graph


# Box, Hunter, Hunter (2005) graphical ANOVA ----------------------------------

# Fit model without nested factors (years within period, stations within regions)

M_8<-lm(lnalpha~region*period, data=mwra)

anovaPlot(M_8, stacked = TRUE, base = TRUE, axes = TRUE,
          faclab = TRUE, labels = TRUE, cex = par("cex"),
          cex.lab = par("cex.lab"))
summary(M_8)
#very low r-squared (0.053) and very high SSR (51.96) suggest lots of unexplained variance
anova(M_8) 

#Interaction term is not significant F 1,860 = 2.31, P=.129
#insignificant interaction suggest there is not divergence b/w regions over time.  However the amount of unexplained variance
#and within group variability (SSR) are very high, interfering with our ability to detect an effect.  
#Also residual DF (850)is inflated because replicate samples within stations and within years are not independent of one another.

plot(M_8) #residuals look ok
hist(residuals(M_8)) #resids fairly normal


# Power Analysis: would 1 sample (no replicates) suffice? -----------------

# Use the tidyverse filter operator to replicates at stations, could be better ways
mwra_1 <- filter(mwra, rep != '2' & rep != 3)
str(mwra_1) # Number of samples deleted from 864 to 854 to 437

ggplot(mwra_1,
       aes(x = year+8/12, y = lnalpha)) + geom_point() +
  geom_smooth(data = mwra_1, mapping = aes(x = year, y = lnalpha), se = TRUE) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station)
# with just 1 replicate, 31 of 36 stations sampled, but BACI effect undetectable

M_9<- lmer(lnalpha ~ region*period +(1|region:Station) + 
             (1|period:fYear), data=mwra_1, REML=TRUE)
summary(M_9) 
# BACI effect for full replicates: 0.053, se= 0.021, p =0.01317
# BACI effect for no replicates   -0.046  se= 0.033  p =0.172

#Based on t-values, interaction should be significant
#strong correlation of intercept w/ region, but I don't think this is a problem
## Type III anova table with p-values for F-tests based on Satterthwaite's
## method:

add_predictions(mwra_1, M_9) %>% 
  ggplot(aes(x = year, y = lnalpha, shape = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model, no replicates')

# Don't do separate facet plots for stations

add_predictions(mwra_1, M_9) %>% 
  ggplot(aes(x = year, y = lnalpha, shape = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  ggtitle('lme4 BACI Model, no replicates')

# Add back facets colour by region
add_predictions(mwra_1, M_9) %>% 
  ggplot(aes(x = year, y = lnalpha, colour = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model, no replicates')

# Plot a line plot, not a lowess smooth plot, note bad fit to NF08
add_predictions(mwra_1, M_9) %>% 
  ggplot(aes(x = year, y = lnalpha, colour = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model, no replicates')

# Find the p value using Kenward-Rogers & Wilks's Theorem -----------------

## Choose type II anova table with Kenward-Roger method for the F-test:
if(requireNamespace("pbkrtest", quietly = TRUE))
  anova(M_9, type=2, ddf="Kenward-Roger")
# p = 0.01317 for full replicates
# p = 0.17219 for no replicates

# Test with Wilks's theorem, Bates's & Andrews' recommendation; dropping
# only the BACI term region:period

M_10 <- lmer(lnalpha ~  1 + region + period + (1|region:Station) +
              (1|period:fYear), data=mwra_1, REML=TRUE)
summary(M_10)
anova(M_9,M_10)
# BACI effect p = 0.0131 for full replicates
# BACI effect p = 0.1705 for no replicates


# Plot fitted frequentist model, no replicates ----------------------------

###plot means and bootstrapped .95 CIs (using Hmisc (1000 bootstrapped samples) 
# ?smean.cl.boot)

lineplot.CI(year,fitted(M_9), group = region, 
            data = mwra_1, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE, xlab = "Year", 
            ylab = "log(Fisher's alpha)")
            abline(v=9.25,lty=3,col="red")

# Bayesian analysis with brms, no reps within stations --------------------

M_11 <- brm(lnalpha ~ region*period + (1|region:Station) +
             (1|period:fYear), data=mwra_1)
M_11
plot(M_11)
# Little BACI effect evident, 95% credible interval contains 0
#                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# regionnear:periodpost   -0.05      0.03    -0.11     0.02  1.00     2980     2398


# Plot the results with no replicates, Bayesian model ---------------------

mwra_1 %>%
  group_by(region) %>%
  #  data_grid(hp = seq_range(hp, n = 51)) %>%
  add_epred_draws(M_11) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .epred), .width = .95, alpha = 1/4) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set2")  
# tells the story, but with horrible Irish colors

# Analysis with 3 replicates reduced to 2 ---------------------------------

# Use the tidyverse filter operator to replicates at stations, could be better ways
mwra_2 <- filter(mwra, rep != '2')
str(mwra_2) # Number of samples deleted from 864 to 854 to 645 (437 with no reps)

ggplot(mwra_2,
       aes(x = year+8/12, y = lnalpha)) + geom_point() +
  geom_smooth(data = mwra_2, mapping = aes(x = year, y = lnalpha), se = TRUE) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station)
# with just 2 replicate, 31 of 36 stations sampled, and BACI effect detectable

M_12<- lmer(lnalpha ~ region*period +(1|region:Station) + 
             (1|period:fYear), data=mwra_2, REML=TRUE)
summary(M_12)
# BACI effect for full replicates:-0.053  se= 0.021, p =0.01311
# BACI effect for 2 replicates    -0.058  se= 0.025  p =0.02004
# BACI effect for no replicates   -0.046  se= 0.033  p =0.172

#Based on t-values, interaction should be significant
#strong correlation of intercept w/ region, but I don't think this is a problem
## Type III anova table with p-values for F-tests based on Satterthwaite's
## method:

add_predictions(mwra_2, M_12) %>% 
  ggplot(aes(x = year, y = lnalpha, shape = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model, 2 replicates')

# Don't do separate facet plots for stations

add_predictions(mwra_2, M_12) %>% 
  ggplot(aes(x = year, y = lnalpha, shape = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  ggtitle('lme4 BACI Model, 2 replicates')

# Add back facets colour by region
add_predictions(mwra_2, M_12) %>% 
  ggplot(aes(x = year, y = lnalpha, colour = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model, 2 replicates')

# Plot a line plot, not a lowess smooth plot, note bad fit to NF08
add_predictions(mwra_2, M_12) %>% 
  ggplot(aes(x = year, y = lnalpha, colour = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model, 2 replicates')

# Find the p value using Kenward-Rogers & Wilks's Theorem -----------------

## Choose type II anova table with Kenward-Roger method for the F-test:
if(requireNamespace("pbkrtest", quietly = TRUE))
  anova(M_12, type=2, ddf="Kenward-Roger")
# p = 0.01317 for full replicates
# p = 0.02004 for 2 replicates
# p = 0.17219 for no replicates

# Test with Wilks's theorem, Bates's & Andrews' recommendation; dropping
# only the BACI term region:period

M_13 <- lmer(lnalpha ~  1 + region + period + (1|region:Station) +
               (1|period:fYear), data=mwra_2, REML=TRUE)
summary(M_13)
anova(M_12,M_13)
# BACI effect p = 0.0131 for full replicates
# BACI effect p = 0.01986 for 2 replicates
# BACI effect p = 0.1705 for no replicates


# Plot fitted frequentist model, no replicates ----------------------------

###plot means and bootstrapped .95 CIs (using Hmisc (1000 bootstrapped samples) 

lineplot.CI(year,fitted(M_12), group = region, 
            data = mwra_2, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE)
            abline(v=9.25,lty=3,col="red")

# Bayesian analysis with brms, 1 and 2 reps, not 3, within stations -----------

M_14 <- brm(lnalpha ~ region*period + (1|region:Station) +
              (1|period:fYear), data=mwra_2)
M_14
plot(M_14)
# Strong BACI effect evident, with 2, not with no replication
#                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# regionnear:periodpost    -0.05      0.03    -0.11     0.02 1.00     2980     2398 1  rep M_11
# regionnear:periodpost    -0.06      0.02    -0.10    -0.01 1.00     4046     2818 2 reps M_14
# regionnear:periodpost    -0.05      0.02    -0.09    -0.01 1.00     3953     2721 3 reps M_6

# Plot the results with 2 replicates, Bayesian model ---------------------

mwra_2 %>%
  group_by(region) %>%
  #  data_grid(hp = seq_range(hp, n = 51)) %>%
  add_epred_draws(M_14) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .epred), .width = .95, alpha = 1/4) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set2")  
# tells the story, but with horrible Irish colors

# DONE! 8/5/22 12:21 AM

# Fit the generalized linear mixed effects model for species --------

# Note that there are more animals in the nearfield, which isn't a problem
# for Fisher's alpha, but will create a near vs. far difference in the 
# pre outfall period.
str(mwra)
# fit species with lme4's model
# see Bruckner brms: An R Package for Bayesian Generalized
# Linear Mixed Models using Stan

M_15 <- glmer(spp ~ region*period + (1|region:Station) +
             (1|period:fYear), family = poisson("log"), data=mwra)
M_15
summary(M_15)

# Determine p-value using Wilks's theorem
M_16 <- glmer(spp ~ region + period + (1|region:Station) +
                (1|period:fYear), family = poisson("log"), data=mwra)
M_16
summary(M_16)

# regionnear:periodpost -0.10495    0.01749  -5.999 1.98e-09 ***

anova(M_15,M_16) # p = 2.03e-9

# Plot the results of the lme4 model, using 1000 bootsrapped samples.
lineplot.CI(mwra$year,fitted(M_15), group = mwra$region, 
            data = mwra, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE, xlab = "Year", 
            ylab = "Species per grab")
            abline(v=9.25,lty=3,col="red")


# Analyze species with Bayesian glmm; brms --------------------------------

M_17 <- brm(spp ~ region*period + (1|region:Station) +
                (1|period:fYear), family = poisson("log"), data=mwra)
M_17
plot(M_17)
# Since the log link function is being used, one must take antilogarithms
# to express results.
exp(c(4.08,3.95, 4.2))
exp(c(0.1,-.02, .22))
exp(c(.18,.05, .31))
exp(c(-.1,-.14, -.06))
# On average there were 59 species per grab with 95% credible interval 51 to 67
# species

# The nearfield had 11% more species than the farfield with a 95% credible
# interval of (-2% to 25% more)

# The post outfall period had 20% more species than the pre-outfall period with
# a 95% credible interval of (5% to 37% more).

# The nearfield/farfield species ratio in the post outfall period was 10% fewer
# nearfield species in the post-outfall period, relative to the farfield, with a
# 95% credible interval of (14% to 5% less)

# Use M_17 as the model to check for diagnostics for the Bayesian glmm --------


# posterior predicted samples, ndraws = no. of datasets to generate
pp_data <- posterior_predict(M_17, ndraws = 200) 

color_scheme_set("gray" )
pdense3 <- ppc_dens_overlay(y = mwra$spp,
                            yrep = pp_data) +
  labs(x = "Species per grab", title = "BACI model, glmm Poisson") +
  theme(title = element_text(size=10))
pdense3
# Reasonable fit

# Observed and simulated mean SKEW metrics  (1000 simulated datasets)
ppc_stat(y = mwra$spp, 
         yrep = posterior_predict(M_17, ndraws = 1000),
         stat = "mean")

# Check skewness
#  Fisher-Pearson Skew function based on NIST definition
color_scheme_set("gray")

# observed and simulated SKEW metric  (1000 simulated datasets)
ppc_stat(y = mwra$spp, 
         yrep = posterior_predict(M_17, ndraws = 1000),
         stat = "skew")
# Skewness is poorly modeled by the glmm, but the formula for skewness is for a
# normally distributed response and the Poisson is positively skewed!
# The skewness of a Poisson distribution is 1/sqrt(lambda), the normal is zero

# observed and simulated MEDIAN metric per region (1000 simulated datasets) 
ppc_stat_grouped(y = mwra$spp,
                 yrep = posterior_predict(M_17, ndraws = 500),
                 group = mwra$region, 
                 stat = "median")
# Excellent fit for both the near and far regions


# Use LOO for model validation --------------------------------------------


# It's important to set save_psis as TRUE
model_loo <- loo(M_17, save_psis = TRUE)
w <- weights(model_loo$psis_object)

color_scheme_set("gray")
ppc_loo_pit_overlay(mwra$spp, 
                    posterior_predict(M_17), 
                    lw = w)
# Should approximate a uniform distribution, but the minimum at 0.25 indicates
# a slight tendency to underdisperson?

# Look for influential single observation using loo -----------------------

k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)


ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Since all observations have a k̂ <0.5 then none are considered influential in
# the posterior given our model is able to characterize them correctly. 
# In other words, our model’s LOO-PSIS estimates are reliable.

# Overall Model comparison, with and without BACI term -------------------------

# fit model with brms without BACI term
M_18 <- brm(spp ~ 1 + region + period +(1|region:Station) +
              (1|period:fYear), family = poisson("log"), data=mwra)
M_18
# Takes about 1 minute
plot(M_18)

# Compute LOO CV for each model
# base model with BACI term
loo1 <- loo(M_17, save_psis = TRUE, cores = 8)

# Reduced model with no BACI term
loo2 <- loo(M_18, save_psis = TRUE, cores = 8)

# Comparison
loo_compare(loo1, loo2)
#     elpd_diff se_diff
# M_6  0.0       0.0   
# M_5 -1.0       0.5   
# M_7 -1.6       2.6  

# M_6 provides the best fit, slightly better than M_5(reps within stations) and 
# M_7 (no BACI term)

# Model validation
loo(M_17, M_18)

# Try a negative binomial model from Kurz statistical rethinking book --------
# Kurz Statistical Rethinking with brms, ggplot2, and the tidyverse
# https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/

M_19 <- brm(data = mwra, spp ~ region*period + (1|region:Station) +
                    (1|period:fYear), family = negbinomial,
    prior = c(prior(normal(0, 10), class = Intercept),
              prior(normal(0, 1), class = b),
              prior(gamma(0.01, 0.01), class = shape)),  # this is the brms default
    iter = 4000, warmup = 1000, cores = 2, chains = 2,
    seed = 11) %>% 
  print()
# plot(M_19)

# now produces an error, and I don't know why
# > plot(M_19)
# Error in xy.coords(x, y, xlabel, ylabel, log) : 
#   'x' is a list, but does not have components 'x' and 'y'
# Since the negative binomial uses the log link, one needs to exponentiate to
# get back into the count metric, i.e., number of species per grab


# Power analysis: Would 5 y pre and post been sufficient?  -----------

mwra_5y <- filter(mwra, year > 1995 & year < 2006)

str(mwra_5y) # Number of samples deleted from 864 to 854 to 536

ggplot(mwra_5y,
       aes(x = year+8/12, y = lnalpha)) + geom_point() +
  geom_smooth(data = mwra_5y, mapping = aes(x = year, y = lnalpha), se = TRUE) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station)

# Now redo the ggplot facet plot
ggplot(mwra_5y,
       aes(x = year, y = lnalpha)) + geom_point() +
  geom_smooth(data = mwra_5y, mapping = aes(x = year, y = lnalpha), se = TRUE) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station)

# Mixed effects ANOVA using lme4 (lmer function) --------------------------

# FULL MODEL (REGION + PERIOD + REGION:PERIOD INTERACTION),
# station random nested in region, year random nested in period)\
# Need to use restricted maximum likelihood for model comparison with Wilks's 
# theorem, i.e., the likelihood ratio test

M_20<- lmer(lnalpha ~ region*period +(1|region:Station) + 
             (1|period:fYear), data=mwra_5y, REML=TRUE)

summary(M_20) # BACI effect -0.040, se= 0.024, p =0.0943
#Fixed effects:
#  Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)             2.654719   0.066980  35.093610  39.634   <2e-16 ***
#  regionnear              0.009533   0.073081  30.049650   0.130   0.8971    
# periodpost              0.039817   0.038688  10.731710   1.029   0.3260    
# regionnear:periodpost  -0.040224   0.023997 495.914324  -1.676   0.0943 .

add_predictions(mwra_5y, M_20) %>% 
  ggplot(aes(x = year, y = lnalpha, shape = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model: 5-y pre & post')

# Don't do separate facet plots for stations

add_predictions(mwra_5y, M_20) %>% 
  ggplot(aes(x = year, y = lnalpha, shape = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  ggtitle('lme4 BACI Model: 5-y pre & post')
# This is a key figure, with the CI for the smoother diverging.

# Add back facets colour by region
add_predictions(mwra_5y, M_20) %>% 
  ggplot(aes(x = year, y = lnalpha, colour = region)) +
  geom_point() +
  geom_smooth(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model: 5-y pre & post')

# Plot a line plot, not a lowess smooth plot, note bad fit to NF08
add_predictions(mwra_5y, M_20) %>% 
  ggplot(aes(x = year, y = lnalpha, colour = region)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  geom_vline (aes(xintercept = 2000.25)) +
  facet_wrap(~Station) +
  ggtitle('lme4 BACI Model: 5-y pre & post')
# Note that years are scrunched together

# Find the p value using Kenward-Rogers & Wilks's Theorem -----------------

## Choose type II anova table with Kenward-Roger method for the F-test:
if(requireNamespace("pbkrtest", quietly = TRUE))
  anova(M_20, type=2, ddf="Kenward-Roger")
# p = 0.094

# Test with Wilks's theorem, Bates's & Andrews' recommendation; dropping
# only the BACI term region:period

M_21 <- lmer(lnalpha ~  1 + region + period + (1|region:Station) +
              (1|period:fYear), data=mwra_5y, REML=TRUE)
summary(M_21)
anova(M_20,M_21)
# Model AIC difference < 1 indicating little discernible difference
# BACI effect p = 0.09375 using Wilks's Theorem (likelihood ratio test)

# Plot the results using the frequentist BACI model -----------------------

lineplot.CI(year,fitted(M_20), group = region, 
            data = mwra_5y, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE)
            abline(v=5.25,lty=3,col="red")

# Bayesian analysis with brms, including nested reps within stations-----------

M_22 <- brm(lnalpha ~ region*period + (1|region:Station) +
             (1|period:fYear), data=mwra_5y)
M_22
# 1:03 to compile, 29 seconds to sample
plot(M_22)


# The BACI effect was -0.04, se=0.02 and credibility interval of -0.09 to 0.01

# Diagnostic plots for Bayesian fitted model ----------------------------------

# See the following reference on diagnostic plots for Bayesian regression
# https://towardsdatascience.com/evaluating-bayesian-mixed-models-in-r-python-27d344a03016
# evaluate simulated data with observed data using Bayesplot

# Use M_22 as the model to check for diagnostics

# posterior predicted samples, ndraws = no. of datasets to generate
pp_data <- posterior_predict(M_22, ndraws = 200) 

color_scheme_set("gray" )
pdense2 <- ppc_dens_overlay(y = mwra_5y$lnalpha,
                            yrep = pp_data) +
  labs(x = "lnalpha", title = "BACI model: 5 y pre & post") +
  theme(title = element_text(size=10))
pdense2
# Caption, 200 simulated data drawn from the posterior (light lines) fit the
# the observed data very well. 

# Observed and simulated mean SKEW metrics  (1000 simulated datasets)
ppc_stat(y = mwra_5y$lnalpha, 
         yrep = posterior_predict(M_22, ndraws = 1000),
         stat = "mean")
# Escellent fit

# Check skewness

color_scheme_set("gray")
# observed and simulated SKEW metric  (1000 simulated datasets)
ppc_stat(y = mwra_5y$lnalpha, 
         yrep = posterior_predict(M_22, ndraws = 1000),
         stat = "skew")
# The skewness is a bit off, and I don't know why

# observed and simulated MEDIAN metric per region (1000 simulated datasets) 
ppc_stat_grouped(y = mwra_5y$lnalpha,
                 yrep = posterior_predict(M_22, ndraws = 500),
                 group = mwra_5y$region, 
                 stat = "median")
# median seems overestimated for far region

# Use LOO, leave one out validation:

# It's important to set save_psis as TRUE
model_loo <- loo(M_22, save_psis = TRUE)
w <- weights(model_loo$psis_object)

color_scheme_set("gray")
ppc_loo_pit_overlay(mwra_5y$lnalpha, 
                    posterior_predict(M_22), 
                    lw = w)
# Should approximate a uniform distribution, but the maximum at 0.5 indicates
# a slight tendency to overdispersion.

# Look for influential single observation using loo -----------------------

k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)


ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Since only one observations has a k̂ >0.5, but is less than 0.7
# the posterior given our model is able to characterize them correctly. 
# In other words, our model’s LOO-PSIS estimates are reliable.

# Overall Model comparison, with and without BACI term -------------------------

# fit model with brms without BACI term
M_23 <- brm(lnalpha ~ 1 + region + period +(1|region:Station) +
             (1|period:fYear), data=mwra_5y)
M_23

# Takes about 1 minute
plot(M_23)

# Compute LOO CV for each model
# base model with BACI term
loo1 <- loo(M_23, save_psis = TRUE, cores = 8)

# Reduced model with no BACI term
loo2 <- loo(M_22, save_psis = TRUE, cores = 8)

# Comparison
loo_compare(loo1, loo2)
#     elpd_diff se_diff
# M_22  0.0       0.0   
# M_23 -0.2       1.7

# M_22 is the better model, but the difference is much less than the standard
# error


#####  Pointwise comparison ------------

# Obtain pointwise ELPD values
elpd1 <- loo1$pointwise[,"elpd_loo"]  # Base model with BACI
elpd2 <- loo2$pointwise[,"elpd_loo"]  # no BACI term

# Build differences dataframe, looking at Stations
elpd_df <- tibble(Region = mwra_5y$region,
                  diff12 = elpd2 - elpd1) %>% 
  mutate(idx = 1:n())

# Plot comparisons
ggplot(elpd_df, aes(x = idx, y = diff12, color = mwra_5y$region)) +
  geom_point(alpha=0.7) +
  geom_hline(aes(yintercept=0)) +
  theme_bw()
# It's tough to see one model fitting the data better than the other


# Final plot with tidybayes from
# http://mjskay.github.io/tidybayes/articles/tidybayes.html
# from the tidybayes vignette
# http://mjskay.github.io/tidybayes/
# From: http://mjskay.github.io/tidybayes/articles/tidy-brms.html

mwra_5y %>%
  group_by(region) %>%
  add_predicted_draws(M_22) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.66, .95), alpha = 1/4) +
  geom_point(data = mwra_5y) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")
# This conveys the key information, considerable overlap among regions, but
# a slight divergence

# Plot posterior predictions instead of means -----------------------------

mwra_5y %>%
  group_by(region) %>%
  add_predicted_draws(M_6) %>%
  ggplot(aes(x = year, y = lnalpha, color = region, fill = region)) +
  stat_lineribbon(aes(y = .prediction), .width = .95, alpha = 1/4) +
  geom_point(data = mwra_5y) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")
# This is an effective plot, in showing the subtle differences between regions


# from http://mjskay.github.io/tidybayes/articles/tidybayes.html

mwra_5y %>%
  group_by(region) %>%
  add_epred_draws(M_22) %>%
  ggplot(aes(x = year, y = lnalpha, color = region)) +
  stat_lineribbon(aes(y = .epred), .width = .95, alpha = 1/4) +
  geom_point(data = mwra_5y) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set2")  
# tells the story, but with horrible Irish colors
# The level legend is distracting

# The MWRA cut the sampling design after 2003 because there was no
# effect detected. Would the present design have detected an effect?

mwra_2003 <- filter(mwra, year < 2003)

str(mwra_2003) # Number of samples deleted from 854 to 633 grabs 

ggplot(mwra_2003,
       aes(x = year+8/12, y = lnalpha)) + geom_point() +
  geom_smooth(data = mwra_2003, mapping = aes(x = year, y = lnalpha), se = FALSE) +
  geom_vline (aes(xintercept = 2000+9.5/12)) +
  facet_wrap(~Station)


# Test MWRA conclusion that no effect was evident by 2003 -----------------

M_24 <- lmer(lnalpha ~  1 + region*period + (1|region:Station) +
               (1|period:fYear), data=mwra_2003, REML=TRUE)
summary(M_24)
# p=.687

M_25 <- lmer(lnalpha ~  1 + region + period + (1|region:Station) +
               (1|period:fYear), data=mwra_2003, REML=TRUE)
summary(M_25)
# Test with Wilks's theorem
anova(M_24,M_25)
# p = 0.687, difference in AIC 1.83

# Model AIC difference < 2 indicating little discernible difference
# BACI effect p = 0.687 using Wilks's Theorem (likelihood ratio test)

# Plot the results using the frequentist BACI model -----------------------

lineplot.CI(year,fitted(M_24), group = region, 
            data = mwra_2003, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE, xlab = "Year", 
            ylab = "log (Fisher's alpha)")
            abline(v=9+9/12,lty=3,col="red")

# Bayesian analysis with brms, including nested reps within stations-----------

M_26 <- brm(lnalpha ~ region*period + (1|region:Station) +
              (1|period:fYear), data=mwra_2003)
M_26
# The estimated BACI effect was -0.01 (se=0.03) with a 95% credible interval of
# -0.08 to 0.05 

plot(M_26)


# Test for Pielou's J' effect ---------------------------------------------

M_27 <- lmer(JP ~  1 + region*period + (1|region:Station) +
               (1|period:fYear), data=mwra, REML=TRUE)
summary(M_27)
# p= 0.023 for 

M_28 <- lmer(JP ~  1 + region + period + (1|region:Station) +
               (1|period:fYear), data=mwra, REML=TRUE)
summary(M_28)
# Test with Wilks's theorem
anova(M_27,M_28)
# p = 0.02275, difference in AIC = 3.2, in the potentially better range. 

# Model AIC difference < 2 indicating little discernible difference
# BACI effect p = 0.687 using Wilks's Theorem (likelihood ratio test)

# Plot the results using the frequentist BACI model -----------------------

lineplot.CI(year,fitted(M_27), group = region, 
            data = mwra, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE, xlab = "Year", 
            ylab = "Pielou's J'")
            abline(v=9+9/12,lty=3,col="red")

# Bayesian analysis with brms, including nested reps within stations-----------

M_29 <- brm(JP ~ region*period + (1|region:Station) +
              (1|period:fYear), data=mwra)
M_29

plot(M_29)

mwra %>%
  group_by(region) %>%
  add_predicted_draws(M_29) %>%
  ggplot(aes(x = year, y = JP, color = region, fill = region)) +
  stat_lineribbon(aes(y = .prediction), .width = .95, alpha = 1/4) +
  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")
# This is an effective plot, in showing the subtle differences between regions

# plot without individual samples:
mwra %>%
  group_by(region) %>%
  add_predicted_draws(M_29) %>%
  ggplot(aes(x = year, y = JP, color = region, fill = region)) +
  stat_lineribbon(aes(y = .prediction), .width = .95, alpha = 1/4) +
#  geom_point(data = mwra) +
  geom_vline (aes(xintercept = 2000.25)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")
# This is an effective plot, in showing the subtle differences between regions

# I'm concerned about the odd post-outfall even-odd sampling, run  depth -------

M_30 <- lmer(depth ~  1 + region*period + (1|region:Station) +
               (1|period:fYear), data=mwra, REML=TRUE)
summary(M_30)
# p= 1.0 model unidentifiable for BACI effect 

# Plot the results using the frequentist BACI model -----------------------

lineplot.CI(year,fitted(M_30), group = region, 
            data = mwra, ci.fun=getCI, lwd=1,pch=c(1,17),col=c("blue","black"),
            lty=c(2,1),cex=1,err.width=.075,err.lty=c(2,1),
            err.col=c("blue","black"),xtick=TRUE, xlab = "Year", 
            ylab = "depth")
abline(v=9+9/12,lty=3,col="red")

# Bayesian analysis with brms, including nested reps within stations-----------

M_31 <- brm(depth ~ region*period + (1|region:Station) +
              (1|period:fYear), data=mwra)
M_31
# Model failed to converge, posterior distributions look hideous


# Analyze whether the regional effect is merely a depth effect ------------

ggplot(data=mwra, mapping = aes (x=lnalpha, color=period)) +
  geom_freqpoly(binwidth = 1)
ggplot(data=mwra, mapping = aes(x=period, y = lnalpha)) + geom_boxplot()

boxplot(lnalpha ~ region*(depth>35), data=mwra) 

# It doesn't appear that can account for the near vs. far effect.
