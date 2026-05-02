#LOGISTIC REGRESSION MODEL FOR LANDSLIDE IN CATCHMENT USING HYPSOMETRIC INTEGRAL

library(lme4) # for multilevel models
library(tidyverse) # for data manipulation and plots
library(haven) #for reading sav data
library(sjstats) #for calculating intra-class correlation (ICC)
library(effects) #for plotting parameter effects
library(jtools) #for transformaing model summaries
library(ROCR) #for calculating area under the curve (AUC) statistics

df=read.csv("G:/My Drive/INVESTIGACION/POSDOC/GLM/df2.csv")
df2 <- df %>% mutate_at(c('hypso_inte', 'rainfall_cat'), ~(scale(.) %>% as.vector))
df_new <- df2 %>% mutate(y = if_else(landslides == 0, 0, 1))
df_new %>% group_by(y) %>% summarise(y = sum(y))
group_size(group_by(df_new,y))
#ONLY H_i
m1 <- glm(formula = y ~ hypso_inte, family = binomial(link = "logit"),data = df_new)
summary(m1)
summ(m1, exp = T)
plot(allEffects(m1))
Pred <- predict(m1, type = "response")
Pred <- if_else(Pred > 0.5, 1, 0)
ConfusionMatrix <- table(Pred, pull(df_new, y))
sum(diag(ConfusionMatrix))/sum(ConfusionMatrix)
ConfusionMatrix
Prob <- predict(m1, type="response")
df_new$prob = predict(m1, type="response")
Pred <- prediction(Prob, as.vector(pull(df_new, y)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC
#H_i & Rainfall
m2 <- glm(formula = y ~ hypso_inte+rainfall_cat, family = binomial(link = "logit"),data = df_new)
summary(m2)
summ(m2, exp = T)
plot(allEffects(m2))
Pred2 <- predict(m2, type = "response")
Pred2 <- if_else(Pred2 > 0.5, 1, 0)
ConfusionMatrix2 <- table(Pred2, pull(df_new, y))
sum(diag(ConfusionMatrix2))/sum(ConfusionMatrix2)
ConfusionMatrix2
Prob2 <- predict(m2, type="response")
df_new$prob2 = predict(m2, type="response")
Pred2 <- prediction(Prob2, as.vector(pull(df_new, y)))
AUC <- performance(Pred2, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC
write.csv(df_new, "G:/My Drive/INVESTIGACION/POSDOC/GLM/df4.csv", row.names=FALSE)
