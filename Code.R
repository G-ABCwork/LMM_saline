library(tidyverse)
library(tidymodels)
library(haven)
library(foreign) # for read.spss
library(MASS)
library(brant) # for testing the parallel regression assumption wit the brant test by Brant (1990) (평행선 검정)
library(ggpubr) # for ggqqplot()
# devtools::install_github("kassambara/rstatix")
library(rstatix) # for rm anova and more
library(naniar) # for missing data visualizations
library(Cairo)
library(conflicted)
conflict_prefer("t_test", "rstatix")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
ggplot2::theme_set(theme_classic()) # change the theme of ggplot2

rct0 <- read_sav("../완료_202105~202106_수액치료효과비교/data/초기통계데이터-BE, HCO3관련변경-22번제외.sav", encoding = "latin1")
rct_names <- read.spss("../완료_202105~202106_수액치료효과비교/data/초기통계데이터-BE, HCO3관련변경-22번제외.sav", reencode = "utf-8") %>% 
    names
colnames(rct0) <- rct_names

# 1. loading and wrangling data -------------------------------------------
rct <- rct0 %>% 
    select(
        # 독립변수
        배정수액2,
        # 반응변수
        starts_with("pH_"), starts_with("BE_"), starts_with("HCO3_"), starts_with("Cl_"),
        SBP_lowest, PR_lowest, Cl상승여부0hrto24hr, AKI발생여부, MAKE여부, 생존퇴원, 생존여부_6개월,
        퇴원시goodCPC, goodCPC_6개월,
        # 공변량
        성별, 나이, initial_shockableRhythm, bystanderCPR, no_flow_time
    ) %>% 
    rowwise() %>% 
    mutate(
        diff_pH = pH_24hr-pH_0hr,
        diff_BE = BE_24hr-BE_0hr,
        diff_HCO3 = HCO3_24hr-HCO3_0hr,
        diff_Cl = Cl_24hr-Cl_0hr
    )
glimpse(rct)
summary(rct)

rct2 <- rct %>% 
    mutate(
        # 0: NS, 1: PS
        배정수액2 = recode(as.numeric(배정수액2), `1` = "0", `2` = "1"),
        # 0: 하락, 1: 변화없음, 2: 상승
        Cl상승여부0hrto24hr = recode(as.numeric(Cl상승여부0hrto24hr), 
                                 `1` = "0", `2` = "1", `3` = "2"),
        # 0: 발생안함, 1: 발생
        AKI발생여부 = as.numeric(AKI발생여부) %>% as.character,
        # 0: 발생안함, 1: 발생
        MAKE여부 = as.numeric(MAKE여부) %>% as.character,
        # 0: 생존, 1: 퇴원
        생존퇴원 = recode(as.numeric(생존퇴원), `1` = "0", `2` = "1"),
        # 0: 생존, 1:사망
        생존여부_6개월 = as.numeric(생존여부_6개월),
        # 0: shockable, 1: nonshockable
        initial_shockableRhythm = recode(as.numeric(initial_shockableRhythm), 
                                         `1` = "0", `2` = "1"),
        # 0: 유, 1: 무, 2: 모름
        bystanderCPR = recode(as.numeric(bystanderCPR), 
                              `1` = "0", `2` = "1", `3` = "2"))

# 2. fitting proportional odds models -------------------------------------
cut_labels <- c("little", "a little", "big")
divide_group <- function(.data, .var){
    
    diff <- .data %>% pull({{.var}}) %>% na.omit() %>% as.numeric
    .data %>% 
        select({{.var}}, 배정수액2,
               # 공변량
               성별, 나이, initial_shockableRhythm, bystanderCPR, no_flow_time) %>% 
        na.omit() %>% 
        mutate(
            group = cut({{.var}},
                           breaks = c(min(diff), quantile(diff)[["25%"]],
                                      quantile(diff)[["75%"]], max(diff)+0.1),
                           right = FALSE, labels = cut_labels, ordered_result = TRUE)
        ) %>% 
        select(group, everything())
    
}
rct_ph <- divide_group(rct2, diff_pH)
rct_BE <- divide_group(rct2, diff_BE)
rct_HCO3 <- divide_group(rct2, diff_HCO3)
rct_Cl <- divide_group(rct2 %>% mutate(diff_Cl = abs(diff_Cl)), diff_Cl)
rct_SBP <- divide_group(rct2, SBP_lowest) %>% 
    mutate(group = recode_factor(group, "little" = "low", "a little" = "medium", "big" = "high",
                                 .ordered = TRUE))
rct_PR <- divide_group(rct2, PR_lowest) %>% 
    mutate(group = recode_factor(group, "little" = "low", "a little" = "medium", "big" = "high",
                                 .ordered = TRUE))
get_stat <- function(.data, var){
    .data %>% 
        group_by(배정수액2) %>% 
        summarize(
            n = n(),
            median = median({{var}}),
            IQR = IQR({{var}}),
            mean = mean({{var}})
        )
}
get_stat(rct_ph, diff_pH)
get_stat(rct_BE, diff_BE)
get_stat(rct_HCO3, diff_HCO3)
get_stat(rct_Cl, diff_Cl)
get_stat(rct_SBP, SBP_lowest)
get_stat(rct_PR, PR_lowest)

## (1) pH
# get_model <- function(.data, .target){ # recipe()의 모형식에서 작동 안하는듯..
#     .data %>%
#         select({{.target}}, 배정수액2:low_flow_time) %>%
#         recipe({{.target}} ~., data = .) %>%
#         step_normalize(all_predictors(), -all_nominal()) %>%
#         step_dummy(all_nominal(), -{{.target}}) %>%
#         prep %>%
#         juice %>%
#         polr({{.target}} ~ ., data = ., Hess = TRUE)
# }
# get_model(rct_odds, group_pH)

ph <- rct_ph %>%
    select(group, 배정수액2:no_flow_time) %>% 
    recipe(group ~., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -group, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2))
model_ph <- polr(group ~ ., data = ph, Hess = TRUE)
brant(model_ph)
summary(model_ph)
# Check the Overall Model Fit
anova(polr(group ~ 1, data = ph), model_ph)

get_odds <- function(.model){
    summary_table <- coef(summary(.model))
    pval <- pnorm(abs(summary_table[, "t value"]), lower.tail = FALSE)* 2
    summary_table <- cbind(summary_table, "p value" = round(pval, 3))
    summary_table %>% 
        as_tibble %>% 
        mutate(variables = rownames(summary_table),
               odds = exp(Value),
               lower_odds = exp(Value - `Std. Error`*qnorm(1-0.05/2)),
               upper_odds = exp(Value + `Std. Error`*qnorm(1-0.05/2))) %>% 
        select(variables, everything())
}
get_odds(model_ph) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)



## (2) BE
BE <- rct_BE %>% 
    select(group, 배정수액2:no_flow_time) %>% 
    recipe(group ~., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -group, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2))
model_BE <- polr(group ~ ., data = BE, Hess = TRUE)
brant(model_BE)
summary(model_BE)
# Check the Overall Model Fit
anova(polr(group ~ 1, data = rct_BE), model_BE)
get_odds(model_BE) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)


## (3) HCO3
HCO3 <- rct_HCO3 %>% 
    select(group, 배정수액2:no_flow_time) %>% 
    recipe(group ~., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -group, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2))
model_HCO3 <-  polr(group ~ ., data = HCO3, Hess = TRUE)
brant(model_HCO3)
summary(model_HCO3)
# Check the Overall Model Fit
anova(polr(group ~ 1, data = rct_HCO3), model_HCO3)
get_odds(model_HCO3) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

## (4) Cl
Cl <- rct_Cl %>% 
    select(group, 배정수액2:no_flow_time) %>% 
    recipe(group ~., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -group, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2))
model_Cl <- polr(group ~ ., data = Cl, Hess = TRUE)
brant(model_Cl)
summary(model_Cl)
# Check the Overall Model Fit
anova(polr(group ~ 1, data = rct_Cl), model_Cl)
get_odds(model_Cl) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

## (5) SBP_lowest
SBP <- rct_SBP %>%
    select(group, 배정수액2:no_flow_time) %>% 
    recipe(group ~., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -group, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2))
model_SBP <- polr(group ~ ., data = SBP, Hess = TRUE)
brant(model_SBP)
summary(model_SBP)
# Check the Overall Model Fit
anova(polr(group ~ 1, data = rct_SBP), model_SBP)
get_odds(model_SBP) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

## (6) PR
PR <- rct_PR %>% 
    select(group, 배정수액2:no_flow_time) %>% 
    recipe(group ~., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -group, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2)) 
model_PR <- polr(group ~ ., data = PR, Hess = TRUE)
brant(model_PR)
summary(model_PR)
# Check the Overall Model Fit
anova(polr(group ~ 1, data = rct_PR), model_PR)
get_odds(model_PR) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

# 3. fitting logistic regression models -------------------------------------
## (1) AKI 발생여부
rct_AKI <- rct2 %>% 
    select(AKI발생여부, 배정수액2, 성별:no_flow_time)
model_AKI <- rct_AKI %>% 
    recipe(AKI발생여부 ~ ., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -AKI발생여부, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2)) %>% 
    glm(AKI발생여부 ~ ., family = binomial, data = .)
get_odds <- function(.model){
    summary_table <- coef(summary(.model))
    pval <- pnorm(abs(summary_table[, "z value"]), lower.tail = FALSE)* 2
    summary_table <- cbind(summary_table, "p value" = round(pval, 3))
    summary_table %>% 
        as_tibble %>% 
        mutate(variables = rownames(summary_table),
               odds = exp(Estimate),
               lower_odds = exp(Estimate - `Std. Error`*qnorm(1-0.05/2)),
               upper_odds = exp(Estimate + `Std. Error`*qnorm(1-0.05/2))) %>% 
        select(variables, everything())
}
get_odds(model_AKI) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

## (2) MAKE 여부
rct_MAKE <- rct2 %>% 
    select(MAKE여부, 배정수액2, 성별:no_flow_time)
model_MAKE <- rct_MAKE %>% 
    recipe(MAKE여부 ~ ., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -MAKE여부, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2)) %>% 
    glm(MAKE여부 ~ ., family = binomial, data = .)
get_odds(model_MAKE) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

## (3) 생존퇴원: MAKE 여부와 정확하게 동일한 분포.
rct_surv <- rct2 %>% 
    select(생존퇴원, 배정수액2, 성별:no_flow_time)
model_surv <- rct_surv %>% 
    recipe(생존퇴원 ~ ., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -생존퇴원, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(배정수액2 = as.numeric(배정수액2)) %>% 
    glm(생존퇴원 ~ ., family = binomial, data = .)
get_odds(model_surv) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)
identical(rct_surv %>% pull(생존퇴원), rct_MAKE %>% pull(MAKE여부))

## (4) 생존여부_6개월
rct_6surv <- rct2 %>% 
    select(생존여부_6개월, 배정수액2, 성별:no_flow_time)
model_6surv <- rct_6surv %>% 
    recipe(생존여부_6개월 ~ ., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -생존여부_6개월, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(생존여부_6개월 = 생존여부_6개월-1,
           배정수액2 = as.numeric(배정수액2)) %>%
    glm(생존여부_6개월 ~ ., family = binomial, data = .)
get_odds(model_6surv) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

## (5) 퇴원시goodCPC
rct_CPC <- rct2 %>% 
    select(퇴원시goodCPC, 배정수액2, 성별:no_flow_time)
model_CPC <- rct_CPC %>% 
    recipe(퇴원시goodCPC ~ ., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -퇴원시goodCPC, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(퇴원시goodCPC = 퇴원시goodCPC-1,
           배정수액2 = as.numeric(배정수액2)) %>% 
    glm(퇴원시goodCPC ~ ., family = binomial, data = .)
get_odds(model_CPC) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

## (6) goodCPC_6개월
rct_6CPC <- rct2 %>% 
    select(goodCPC_6개월, 배정수액2, 성별:no_flow_time)
model_6CPC <- rct_6CPC %>% 
    recipe(goodCPC_6개월 ~ ., data = .) %>% 
    step_normalize(all_predictors(), -all_nominal()) %>% 
    step_dummy(all_nominal(), -goodCPC_6개월, -배정수액2) %>% 
    prep %>% 
    juice %>% 
    mutate(goodCPC_6개월 = goodCPC_6개월-1,
           배정수액2 = as.numeric(배정수액2)) %>% 
    glm(goodCPC_6개월 ~ ., family = binomial, data = .)
get_odds(model_6CPC) %>% 
    filter(variables == "배정수액2") %>% 
    select(`variables`, `p value`:upper_odds)

# 4. repeated measures with linear mixed models ------------------------------------------
glimpse(rct)

## (1) pH
repeated_pH <- rct %>% 
    rename(treatment = 배정수액2) %>% 
    select(treatment:pH_24hr)
CairoWin()
gg_miss_var(repeated_pH)
plot_miss_case <- function(.data){
    miss_case_summary(.data) %>% 
        filter(n_miss > 0) %>% 
        ggplot(aes(x = fct_reorder(factor(case), n_miss), y = n_miss)) +
        geom_bar(stat = "identity") +
        labs(x = "Case number", y = "The number of missing values") +
        coord_flip()
}
CairoWin()
plot_miss_case(repeated_pH)

ph <- repeated_pH %>% 
    pivot_longer(pH_0hr:pH_24hr, "time") %>% 
    mutate(
        time = recode(time, pH_0hr = "t1", pH_30min = "t2", pH_1hr = "t3", pH_2hr = "t4", 
                      pH_4hr = "t5", pH_6hr = "t6", pH_12hr = "t7", pH_18hr = "t8", pH_24hr = "t9")
    ) %>% 
    rename(concetration = value) %>% 
    arrange(time, treatment) %>% 
    mutate(id = rep(c(1:27, 28:53), 9)) %>% 
    select(id, everything()) %>% 
    mutate_at(vars(id:time), factor)

### summary statistics
ph %>%
    drop_na() %>% 
    group_by(treatment, time) %>%
    get_summary_stats(concetration, type = "mean_se") %>% 
    arrange(time)

## lmm
fit <- lmerTest::lmer(concetration  ~ time * treatment + (1|id), data = ph %>% drop_na())
fit
anova(fit)

#### 시간에 대한 사후 분석
ph %>%
    pairwise_t_test(
        concetration ~ time, paired = TRUE, 
        p.adjust.method = "bonferroni"
    ) %>% 
    filter(p.adj <=0.05)

## (2) Base excess
repeated_BE <- rct %>% 
    rename(treatment = 배정수액2) %>% 
    select(treatment, starts_with("BE_"))
CairoWin()
gg_miss_var(repeated_BE)
CairoWin()
plot_miss_case(repeated_BE)
BE <- repeated_BE %>% 
    pivot_longer(BE_0hr:BE_24hr, "time") %>% 
    mutate(
        time = recode(time, BE_0hr = "t1", BE_30min = "t2", BE_1hr = "t3", BE_2hr = "t4", 
                      BE_4hr = "t5", BE_6hr = "t6", BE_12hr = "t7", BE_18hr = "t8", BE_24hr = "t9")
    ) %>% 
    rename(concetration = value) %>% 
    arrange(time, treatment) %>% 
    mutate(id = rep(c(1:27, 28:53), 9)) %>% 
    select(id, everything()) %>% 
    mutate_at(vars(id:time), factor)

### summary statistics
BE %>%
    drop_na() %>% 
    group_by(treatment, time) %>%
    get_summary_stats(concetration, type = "mean_se") %>% 
    arrange(time)

## lmm
fit <- lmerTest::lmer(concetration  ~ time * treatment + (1|id), data = BE %>% drop_na())
fit
anova(fit)
# 시간, 처리에 관한 효과 유의하게 존재

#### 사후 분석
BE %>%
    pairwise_t_test(
        concetration ~ time, paired = TRUE, 
        p.adjust.method = "bonferroni"
    ) %>% 
    filter(p.adj<=0.05)

BE %>% 
    pairwise_t_test(
        concetration ~ treatment,
        p.adjust.method = "bonferroni"
    )
# -> 어차피 교호효과없으므로, group_by(time) 할 필요없이 다이렉트하게 비교하면 됨. 모든 포인트에서 PS의 BE가 더큼

## (3) HCO3
repeated_HCO3 <- rct %>% 
    rename(treatment = 배정수액2) %>% 
    select(treatment, starts_with("HCO3_"))
CairoWin()
gg_miss_var(repeated_HCO3)
CairoWin()
plot_miss_case(repeated_HCO3)

HCO3 <- repeated_HCO3 %>% 
    pivot_longer(HCO3_0hr:HCO3_24hr, "time") %>% 
    mutate(
        time = recode(time, HCO3_0hr = "t1", HCO3_30min = "t2", HCO3_1hr = "t3", HCO3_2hr = "t4", 
                      HCO3_4hr = "t5", HCO3_6hr = "t6", HCO3_12hr = "t7", HCO3_18hr = "t8", HCO3_24hr = "t9")
    ) %>% 
    rename(concetration = value) %>% 
    arrange(time, treatment) %>% 
    mutate(id = rep(c(1:27, 28:53), 9)) %>% 
    select(id, everything()) %>% 
    mutate_at(vars(id:time), factor)

### summary statistics
HCO3 %>%
    drop_na() %>% 
    group_by(treatment, time) %>%
    get_summary_stats(concetration, type = "mean_se") %>% 
    arrange(time)

## lmm
fit <- lmerTest::lmer(concetration  ~ time * treatment + (1|id), data = HCO3 %>% drop_na())
fit
anova(fit)
# 시간, 처리에 관한 효과 유의하게 존재

#### 사후 분석
HCO3 %>%
    pairwise_t_test(
        concetration ~ time, paired = TRUE, 
        p.adjust.method = "bonferroni"
    ) %>% 
    filter(p.adj<=0.05)

HCO3 %>% 
    pairwise_t_test(
        concetration ~ treatment,
        p.adjust.method = "bonferroni"
    )

## (4) Cl0
repeated_Cl <- rct %>% 
    rename(treatment = 배정수액2) %>% 
    select(treatment, starts_with("Cl_"))
miss_case_summary(repeated_Cl) %>% 
    filter(n_miss > 0)

Cl <- repeated_Cl %>% 
    pivot_longer(Cl_0hr:Cl_24hr, "time") %>% 
    mutate(
        time = recode(time, Cl_0hr = "t1", Cl_6hr = "t2", Cl_12hr = "t3", Cl_18hr = "t4", Cl_24hr = "t5")
    ) %>% 
    rename(concetration = value) %>% 
    arrange(time, treatment) %>% 
    mutate(id = rep(c(1:27, 28:53), 5)) %>% 
    select(id, everything()) %>% 
    mutate_at(vars(id:time), factor)

### summary statistics
Cl %>%
    group_by(treatment, time) %>%
    get_summary_stats(concetration, type = "mean_se") %>% 
    arrange(time) %>% 
    pull(se) %>% 
    round(2)

## lmm
fit <- lmerTest::lmer(concetration  ~ time * treatment + (1|id), data = Cl %>% drop_na())
fit
anova(fit)
# 시간, 처리, 교호효과에 관한 효과 유의하게 존재

#### 사후 분석
##### Simple main effect of group variable
# Pairwise comparisons between group levels
pwc <- Cl %>%
    group_by(time) %>%
    t_test(concetration ~ treatment) %>% 
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
pwc # # t1에서는 차이가 없으나, 나머지 time point에서는 차이가 있음. NS군의 Cl이 큼

##### Simple main effects of time variable
# Pairwise comparisons between group levels
one.way <- Cl %>%
    group_by(treatment) %>%
    anova_test(dv = concetration, wid = id, within = time) %>%
    get_anova_table() %>%
    adjust_pvalue(method = "bonferroni")
one.way # NS 군에서는 시간에 따른 Cl의 차이가 유의하나, PS에서는 존재하지 않음.

?pairwise_t_test
