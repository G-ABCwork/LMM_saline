library(tidyverse)
library(ggpubr) # for ggqqplot()
library(rstatix) # for convert_as_factor()
library(Cairo)
ggplot2::theme_set(theme_classic()) # change the theme of ggplot2

# 1. One-way RM ANOVA -----------------------------------------------------
## (1) Data preparation
### Wide format
data("selfesteem", package = "datarium")
selfesteem
# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
selfesteem2 <- selfesteem %>% 
    pivot_longer(t1:t3, names_to = "time") %>% 
    convert_as_factor(id, time) %>% 
    rename(score = value)
selfesteem2

## (2) Summary staistics
### Compute some summary statistics of the self-esteem score by groups (time): mean and sd (standard deviation)
selfesteem2 %>%
    group_by(time) %>%
    get_summary_stats(value, type = "mean_sd")

## (3) Visualization
### Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(selfesteem2, x = "time", y = "score", add = "point")
CairoWin()
bxp

## (4) Check assumptions
### Outliers
selfesteem2 %>%
    group_by(time) %>%
    identify_outliers(value) # There were no extreme outliers.

### Normality assumption
selfesteem2 %>%
    group_by(time) %>%
    shapiro_test(score)
CairoWin()
ggqqplot(selfesteem2, "score", facet.by = "time")


### (5) Computation
### Assumption of sphericity:
### As mentioned in previous sections, the assumption of sphericity will be automatically checked during 
### the computation of the ANOVA test using the R function anova_test() [rstatix package]. 
### The Mauchly’s test is internally used to assess the sphericity assumption.
### By using the function get_anova_table() [rstatix] to extract the ANOVA table, 
### the Greenhouse-Geisser sphericity correction is automatically applied to factors 
### violating the sphericity assumption.
res_aov <- anova_test(data = selfesteem2, dv = score, wid = id, within = time)
get_anova_table(res_aov)
## -> The self-esteem score was statistically significantly different at the different time points during the diet, 
## F(2, 18) = 55.5, p < 0.0001, eta2[g] = 0.83.

## (6) Post-hoc tests
### u can perform multiple pairwise paired t-tests between the levels of the within-subjects factor (here time). 
### P-values are adjusted using the Bonferroni multiple testing correction method.
pwc <- selfesteem2 %>%
    pairwise_t_test(
        score ~ time, paired = TRUE,
        p.adjust.method = "bonferroni"
    )
pwc
# -> All the pairwise differences are statistically significant.

## (5) Report
### We could report the result as follow:
### The self-esteem score was statistically significantly different at the different time points, 
### F(2, 18) = 55.5, p < 0.0001, generalized eta squared = 0.82.
### Post-hoc analyses with a Bonferroni adjustment revealed that all the pairwise differences, 
### between time points, were statistically significantly different (p <= 0.05).

### Visualization: box plots with p-values
pwc2 <- pwc %>% add_xy_position(x = "time")
CairoWin()
bxp + 
    stat_pvalue_manual(pwc2) +
    labs(
        subtitle = get_test_label(res_aov, detailed = TRUE),
        caption = get_pwc_label(pwc)
    )

# 2. Two-way RM ANOVA -----------------------------------------------------
## (1) Data preparation
# Wide format
set.seed(123)
data("selfesteem2", package = "datarium")
# Load and show one random row by treatment group:
selfesteem2 %>% 
    sample_n_by(treatment, size = 1)

# Gather the columns t1, t2 and t3 into long format.
# Convert id and time into factor variables
selfesteem2 <- selfesteem2 %>%
    pivot_longer(t1:t3, "time") %>%
    convert_as_factor(id, time) %>% 
    rename(score = value)

# Inspect some random rows of the data by groups
set.seed(123)
selfesteem2 %>% sample_n_by(treatment, time, size = 1)

## (2) Summary staistics
selfesteem2 %>%
    group_by(treatment, time) %>%
    get_summary_stats(score, type = "mean_sd")

## (3) Visualization
bxp <- ggboxplot(
    selfesteem2, x = "time", y = "score",
    color = "treatment", palette = "jco"
)
CairoWin()
bxp

## (4) Check assumptions
## Outliers
selfesteem2 %>%
    group_by(treatment, time) %>%
    identify_outliers(score) # There were no extreme outliers.

## Normality assumption
selfesteem2 %>%
    group_by(treatment, time) %>%
    shapiro_test(score)
# -> The self-esteem score was normally distributed at each time point (p > 0.05), 
# except for ctr treatment at t1, as assessed by Shapiro-Wilk’s test.
CairoWin()
ggqqplot(selfesteem2, "score", ggtheme = theme_bw()) +
    facet_grid(time ~ treatment, labeller = "label_both")
# -> From the plot above, as all the points fall approximately along the reference line, we can assume normality.

## (5) Computation
res_aov <- anova_test(
    data = selfesteem2, dv = score, wid = id,
    within = c(treatment, time)
)
get_anova_table(res_aov)
# -> There is a statistically significant two-way interactions 
# between treatment and time, F(2, 22) = 30.4, p < 0.0001.

## (6) Post-hoc tests
###  A significant two-way interaction indicates that the impact that one factor (e.g., treatment) has on 
### the outcome variable (e.g., self-esteem score) depends on the level of the other factor (e.g., time) 
### (and vice versa). So, you can decompose a significant two-way interaction into:
###    Simple main effect: run one-way model of the first variable (factor A) at each level of 
###                        the second variable (factor B),
###    Simple pairwise comparisons: if the simple main effect is significant, run multiple pairwise comparisons 
###                                 to determine which groups are different.
### For a non-significant two-way interaction, you need to determine
### whether you have any statistically significant main effects from the ANOVA output.

### i. Procedure for a significant two-way interaction

#### Effect of treatment: In our example, we’ll analyze the effect of treatment on self-esteem score 
####  at every time point.
#### Note that, the treatment factor variable has only two levels (“ctr” and “Diet”); 
#### thus, ANOVA test and paired t-test will give the same p-values.
# Effect of treatment at each time point
one.way <- selfesteem2 %>%
    group_by(time) %>%
    anova_test(dv = score, wid = id, within = treatment) %>%
    get_anova_table() %>%
    adjust_pvalue(method = "bonferroni")
one.way

# -> Considering the Bonferroni adjusted p-value (p.adj), it can be seen that the simple main effect of treatment 
# was not significant at the time point t1 (p = 1). It becomes significant at t2 (p = 0.036) and t3 (p = 0.00051).

# Pairwise comparisons between treatment groups
pwc <- selfesteem2 %>%
    group_by(time) %>%
    pairwise_t_test(
        score ~ treatment, paired = TRUE,
        p.adjust.method = "bonferroni"
    )
pwc
# -> Pairwise comparisons show that the mean self-esteem score was significantly different between 
# ctr and Diet group at t2 (p = 0.12) and t3 (p = 0.00017) but not at t1 (p = 0.55).

#### Effect of time: Note that, it’s also possible to perform the same analysis 
#### for the time variable at each level of treatment. You don’t necessarily need to do this analysis.
# Effect of time at each level of treatment
one.way2 <- selfesteem2 %>%
    group_by(treatment) %>%
    anova_test(dv = score, wid = id, within = time) %>%
    get_anova_table() %>%
    adjust_pvalue(method = "bonferroni")
# -> After executing the R code above, you can see that the effect of time is significant only 
# for the control trial, F(2, 22) = 39.7, p < 0.0001. 

# Pairwise comparisons between time points
pwc2 <- selfesteem2 %>%
    group_by(treatment) %>%
    pairwise_t_test(
        score ~ time, paired = TRUE,
        p.adjust.method = "bonferroni"
    )
pwc2
# -> Pairwise comparisons show that all comparisons between time points were statistically significant 
# for control trial.

### ii. Procedure for non-significant two-way interaction
#### If the interaction is not significant, you need to interpret the main effects for each of the two variables: 
#### treatment and time. A significant main effect can be followed up with pairwise comparisons.
#### -> In our example (see ANOVA table in res.aov), there was a statistically significant main effects 
#### of treatment (F(1, 11) = 15.5, p = 0.002) and time (F(2, 22) = 27.4, p < 0.0001) on the self-esteem score.

# Pairwise paired t-test comparisons
# comparisons for treatment variable
selfesteem2 %>%
    pairwise_t_test(
        score ~ treatment, paired = TRUE, 
        p.adjust.method = "bonferroni"
    )
# comparisons for time variable
selfesteem2 %>%
    pairwise_t_test(
        score ~ time, paired = TRUE, 
        p.adjust.method = "bonferroni"
    )
# All pairwise comparisons are significant.
