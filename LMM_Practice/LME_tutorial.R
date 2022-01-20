## Tutorial -  LME ##

# Environment: 
#  - R version 4.0.2 (2020-06-22)
#  - Platform: x86_64-w64-mingw32/x64 (64-bit)
#  - Running under: Windows 10 x64 (build 18362)
#
# Reference 
# 1) Mixed Effect Model
# - Singmann, & Kellen (2017). An Introduction to Mixed Models ~
# - Baayen, Davison, & Bates (2008). Mixed-effects modeling with ~
# - Lo, & Andrew (2015). To transform or not to transform ~
# - Jaeger (2008). Categorical data analysis: Away from ANOVAs ~
# - Mixed Model Tutorial (http://www.bodowinter.com/resources.html -> part one & two)
# - Mixed Model Reanalysis of RT data 
#   (https://cran.r-project.org/web/packages/afex/vignettes/afex_mixed_example.html#overview)
#
# 2) ANOVA, t-test, etc...
# - A Language, not a Letter: Learning Statistics in R (https://ademos.people.uic.edu/index.html)
# - ANOVA and Post-Hoc Contrasts: Reanalysis of Singmann and Klauer (2011)
#   (https://cran.r-project.org/web/packages/afex/vignettes/afex_anova_example.html)
# - Just Enough R 
#   (https://benwhalley.github.io/just-enough-r/contrasts-lmer.html)
#
# 3) power analysis in Linear Mixed Effects Model
# - https://jakewestfall.shinyapps.io/crossedpower/

# ---------------------------------------------------------- #
# 1. Preparing ####
# ---------------------------------------------------------- #

# get ready
rm(list=ls())
getwd()
setwd("C:/Users/sorel/Desktop/LMM/LMM_Practice") # working directory path
set.seed(17) # for replication
options(scipen=4)

# install // load packages 
# Some packages need to be loaded. 
# We use `pacman` as a package manager, which takes care of the other packages. 
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
if (!require("papaja", quietly = TRUE)) devtools::install_github("crsh/papaja")
if (!require("patchwork", quietly = TRUE)) devtools::install_github("thomasp85/patchwork")
if (!require("klippy", quietly = TRUE)) devtools::install_github("RLesur/klippy")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
if (!require("Rmisc", quietly = TRUE)) install.packages("Rmisc") # Never load it directly.
pacman::p_load(tidyverse, papaja, knitr, dplyr, car, psych, afex, lme4, lmerTest, 
               emmeans, ggplot2, ggpubr, lattice, latticeExtra, parallel, effects)
library("patchwork")
library("klippy")

options(knitr.kable.NA = '') # hide NA with knitr function
klippy::klippy()

# ---------------------------------------------------------- #
# 2. Loading Data ####
# ---------------------------------------------------------- #

## load data
dat <- read.csv("dat_tutorial.csv", header = T)
dim(dat) # 84 rows 5 columns
glimpse(dat, width=70)

# SN        <int> 1, 1, 1, 2, 2, 2, ...; 1 ~ 6, subject number
# Trial     <int> 1, 1, 1, 2, 2, 2, ...; 1 ~ 14, trial number
# Cond      <int> 1, 1, 1, 2, 2, 2, ...; 1 ~ 2, w/n factor
# Item      <int> 1, 2, 3, 4, 5, 6, 7, ...; 1 ~ 7, stimuli
# RT        <dbl> 2045, 2597, 2869, 2768, ...; RT, ms


## check data
view(dat)

table(dat$SN)
#  1  2  3  4  5  6 
# 14 14 14 14 14 14

table(dat$Cond, dat$SN) 
table(dat$Item, dat$SN) 

headTail(dat)

## change class of main factors: double to factor
dat$SN = factor(dat$SN)
dat$Cond = factor(dat$Cond, levels=c(1,2), labels=c("c1","c2"))
dat$Item = factor(dat$Item)
# dat$RT <- dat$RT*1000
# dat$Corr <- as.numeric(dat$Corr==1)

# ---------------------------------------------------------- #
# 3. Preprocessing ####
# ---------------------------------------------------------- #

# trimming
tdat <- dat %>% filter(RT > 200 & RT < 10000) %>%
  group_by(SN) %>% # grouping by participants
  nest() %>%
  mutate(lbound = map(data, ~mean(.$RT)-3*sd(.$RT)),
         ubound = map(data, ~mean(.$RT)+3*sd(.$RT))) %>% # make new data (3sd cut)
  unnest(c(lbound, ubound))%>% 
  unnest(data) %>% 
  mutate(Outlier = (RT < lbound)|(RT > ubound)) %>% # set outlier
  filter(Outlier == FALSE) %>% # filtering outlier
  ungroup() %>% 
  select(SN, Cond, Item, RT) # select variables to analyze

# outlier trial ratio
100-100*(nrow(tdat)/nrow(dat))

# before trimming
den1 <- ggplot(dat, aes(x=RT)) + 
  geom_density() + 
  theme_bw(base_size = 18) +
  labs(x = "Raw RT") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 

# after trimming
den2 <- ggplot(tdat, aes(x=RT)) + 
  geom_density() + 
  theme_bw(base_size = 18) + 
  labs(x = "Trimmed RT") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
# den1 + den2

# ---------------------------------------------------------- #
# 4. Descriptive Statistics ####
# ---------------------------------------------------------- #

# subject-level, long format
p1rtL <- tdat %>% group_by(SN, Cond) %>%
  dplyr::summarise(RT = mean(RT)) %>%
  ungroup()
p1rtL %>% kable(digits=2)

# subject-level, wide format
p1rtW <- p1rtL %>% spread(key=Cond, value = RT)
p1rtW %>% kable(digits=2)

# summary table (grand mean)
p1rtG <- p1rtL %>% group_by(Cond) %>%
  summarise(RT.m = mean(RT), RT.sd = sd(RT)) %>%
  ungroup()
p1rtG$RT.se <- Rmisc::summarySEwithin(data = p1rtL, measurevar = "RT", 
                                      idvar = "SN", withinvars = "Cond")$se
p1rtG$RT.ci <- Rmisc::summarySEwithin(data = p1rtL, measurevar = "RT", 
                                      idvar = "SN", withinvars = "Cond")$ci
p1rtG <- p1rtG %>% 
  mutate(lower.ci = RT.m-RT.ci,
         upper.ci = RT.m+RT.ci)
p1rtG %>% kable(digits=2)

# simple plot
apa_barplot(data=as.data.frame(p1rtL), 
            id="SN", dv="RT", 
            ylab = "RT (ms)", xlab = "Condition", 
            main = c("Response Time (ms)"),  
            dispersion =  within_subjects_conf_int, # w/n confidence interval
            factors=c("Cond"), 
            ylim = c(0, 3000),
            las=1)

a = ggplot(data=tdat, aes(x = SN, y = RT)) + 
  ggtitle("참가자별") +
  xlab("참가자") + ylab("RT") + theme_bw() +
  theme(text=element_text(size=14)) +
  theme(legend.title = element_text(face = 1,size = 15)) +
  geom_boxplot()

b = ggplot(data=tdat, aes(x = Item, y = RT)) + 
  ggtitle("자극별") +
  xlab("자극") + ylab("RT") + theme_bw() +
  theme(text=element_text(size=14)) +
  theme(legend.title = element_text(face = 1,size = 15)) +
  geom_boxplot()
ggarrange(a, b)

# ---------------------------------------------------------- #
# 5. Inference Statistics ####
# ---------------------------------------------------------- #

## Linear Mixed Effect Modeling

# fitting method: 
#   1) Maximum Likelihood Estimation(ML)
#   2) Restricted Maximum Likelihood Estimation(REML)

# significance test method: 
#   1) Kenward-Roger approximation : method ="KR"
#   2) Satterthwaite approximation : method = "S"
#   3) Likelihood Ratio Test : method = "LRT"
#   4) Parametric bootstrapping : method = "PB"

# function: lmer(lme4), mixed(afex)

# Model

# 1) By-subject random intercept model
# model <- lmer(y ~ x1 + (1 | SN), 
#               data = dataframe)

# 2) By-subject random intercept & slope model
# model <- lmer(y ~ x1 + (1 + x1 | SN), 
#               data = dataframe)

# 3) By-item random intercept model
# model <- lmer(y ~ x1 + (1 | Item), 
#               data = dataframe)

# 4) By-item random intercept & slope model
# model <- lmer(y ~ x1 + (1 + x1 | Item), 
#               data = dataframe)

# Cross-random effect model
# model <- lmer(y ~ x1 + (1 + x1 | SN) + (1 + x1 | Item), 
#               data = dataframe)


# Random Intercept only
rt.lmer.m0 <- lmer(RT ~ (1|SN), tdat)
summary(rt.lmer.m0)

# By-subject Random Intercept 
rt.lmer.m1 <- lmer(RT ~ Cond + (1|SN), tdat)
summary(rt.lmer.m1)
anova(rt.lmer.m1)
anova(rt.lmer.m0, rt.lmer.m1)

# By-subject Random Intercept & Slope
rt.lmer.m2 <- lmer(RT ~ Cond + (1+Cond|SN), tdat, REML = FALSE)
summary(rt.lmer.m2)
anova(rt.lmer.m2)
anova(rt.lmer.m1, rt.lmer.m2)

# By-item Random Intercept 
rt.lmer.m3 <- lmer(RT ~ Cond + (1|Item), tdat, REML = FALSE)
summary(rt.lmer.m3)
anova(rt.lmer.m3)

# By-item Random Intercept & Slope
rt.lmer.m4 <- lmer(RT ~ Cond + (1+Cond|Item), tdat, REML = FALSE)
summary(rt.lmer.m4)
anova(rt.lmer.m4)

# Cross-Random Effect Model (Full Model)
rt.lmer.full <- lmer(RT ~ Cond + (1 + Cond|SN) + (1 + Cond|Item), tdat, REML = FALSE)
summary(rt.lmer.full)
anova(rt.lmer.full)

# random effect 
coef(rt.lmer.full) # random intercept & slope
ranef(rt.lmer.full) # random effect

# coef(rt.lmer.full)$SN$`(Intercept)` - ranef(rt.lmer.full)$SN$`(Intercept)` # F1 Intercept
# coef(rt.lmer.full)$SN$`Condc2` - ranef(rt.lmer.full)$SN$`Condc2` # F1 Slope
# coef(rt.lmer.full)$Item$`(Intercept)` - ranef(rt.lmer.full)$Item$`(Intercept)` # F2 Intercept
# coef(rt.lmer.full)$Item$`Condc2` - ranef(rt.lmer.full)$Item$`Condc2` # F2 slope

### post-hoc
rt.posthoc <- emmeans(rt.lmer.full, pairwise ~ Cond, type = "response", adjust="bon") # adjust="bon"
rt.posthoc$contrasts %>% kable(digits=2)
plot(rt.posthoc, horizontal = FALSE, comparisons = T)
model <- rt.lmer.full

model_predicted_effects <- as.data.frame(effects::effect("Cond", model))
ggplot(data = model_predicted_effects,
       aes(x = Cond, y = fit)) +
  geom_pointrange(aes(ymax = upper,
                      ymin = lower),
                  position = position_dodge(width = 1)) +
  # geom_line(aes(x = Block, y = fit, group = Btw),
  # position = position_dodge(width = 1)) +
  ylab("RT") +
  xlab("Condition") +
  scale_colour_grey() +
  theme_classic() +
  theme(legend.justification=c(1,1), legend.position=c(1,1))
# this is estimated/predicted plot
# you can check exact numbers with model_predicted_effects

# Assumtion check
# 1) Linearity & Homoskedasticity(equal variance)
plot(fitted(rt.lmer.full),residuals(rt.lmer.full)) + 
  abline(h=0, col="red", lwd=1, lty=2)
# plot(residuals(rt.lmer.full))

# 2) Normality of residuals
hist(residuals(rt.lmer.full))
qqnorm(residuals(rt.lmer.full)) + 
  qqline(residuals(rt.lmer.full), col='2')

# model comparison
anova(rt.lmer.full, rt.lmer.m1) # LRT method
 
#### Afex
(nc <- detectCores())
cl <- makeCluster(rep("localhost", nc))
rt.afex.full <- afex::mixed(RT ~ Cond + (1+Cond|SN) + (1+Cond|Item),
                        data = tdat,
                        method = "S",  # "KR", "S", "LRT", "PB"
                        REML = FALSE,
                        cl = cl,
                        control = lmerControl(optCtrl = list(maxfun = 1e10)))

stopCluster(cl)

summary(rt.afex.full)
anova(rt.afex.full)

model <- rt.afex.full$full_model

# structure of random effect 
coef(model) # random intercept & slope
ranef(model) # random effect

### post-hoc
rt.posthoc <- emmeans(model, pairwise ~ Cond, type = "response", adjust="bon") # adjust="bon"
rt.posthoc$contrasts %>% kable(digits=2)
plot(rt.posthoc, horizontal = FALSE, comparisons = T)


model_predicted_effects <- as.data.frame(effects::effect("Cond", rt.lmer.full))
ggplot(data = model_predicted_effects,
       aes(x = Cond, y = fit)) +
  geom_pointrange(aes(ymax = upper,
                      ymin = lower),
                  position = position_dodge(width = 1)) +
  # geom_line(aes(x = Block, y = fit, group = Btw),
  # position = position_dodge(width = 1)) +
  ylab("RT") +
  xlab("Condition") +
  scale_colour_grey() +
  theme_classic() +
  theme(legend.justification=c(1,1), legend.position=c(1,1))


# ---------------------------------------------------------- #
# Appendix. GLMM ####
# ---------------------------------------------------------- #

#### GLMM for Correctness Analysis

### lme4 package
p1.acc.glmer <- glmer(Corr ~ Block*Btw + (Block|SN) + (1|IMname), p1,
                      family = binomial(link="logit"),
                      glmerControl(optimizer = c("bobyqa"),
                                   optCtrl = list(maxfun = 1e7)))
anova(p1.acc.glmer)
summary(p1.acc.glmer)

anova(p1.acc.glmer, p1.acc.glmer.reduced)

p1.acc.m7 <- emmeans(p1.acc.glmer.reduced1, pairwise ~ Block | Btw, type = "response") # adjust="bon"
p1.acc.m7$contrasts %>% kable(digits=2)

### afex package
# (nc <- detectCores())
# cl <- makeCluster(rep("localhost", nc))
# p1.acc.glmixed <- afex::mixed(Corr ~ Btw*Block + (Block|SN) + (1|IMname),
#                              data = p1,
#                              family = binomial(link="logit"),
#                              method = "LRT", cl = cl,
#                              expand_re = F,
#                              control = lmerControl(optimizer = "bobyqa",
#                                                    optCtrl = list(maxfun = 1e7)))
# stopCluster(cl)

#### GLMM for RT data

### lme4 package
# p1.rt.glmer.m3 <- glmer(RT ~ Block*Btw + (1|SN) + (1|IMname), tp1,
#                         family = inverse.gaussian(link="identity"),
#                         glmerControl(optimizer = c("bobyqa"),
#                                      optCtrl = list(maxfun = 1e7))) 
# 
# anova(p1.rt.glmer)
# summary(p1.rt.glmer)
# anova(p1.rt.glmer, p1.rt.glmer.reduced)
# 
# p1.rt.m7 <- emmeans(p1.rt.glmer.reduced1, pairwise ~ Block | Btw, type = "response") # adjust="bon"
# p1.rt.m7$contrasts %>% kable(digits=2)

### afec package
# (nc <- detectCores())
# cl <- makeCluster(rep("localhost", nc))
# p1.rt.glmixed <- afex::mixed(RT ~ Btw*Block + (1|SN) + (1|IMname),
#                              data = tp1,
#                              family = inverse.gaussian(link="identity"),
#                              method = "LRT", cl = cl,
#                              expand_re = F,
#                              control = lmerControl(optimizer = "bobyqa",
#                                                    optCtrl = list(maxfun = 1e7)))
# 
# # nested model을 비교하는 method (Singmann, & Kellen, 2017)
# # 1) Kenward-Roger approximation : method ="KR"
# # 2) Satterthwaite approximation : method = "S"
# # 3) Likelihood Ratio Test : method = "LRT"
# # 4) Parametric bootstrapping : method = "PB"
# stopCluster(cl)

