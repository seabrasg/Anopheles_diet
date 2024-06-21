##########
# Author: Sofia Seabra (sgseabra@ihmt.nl.pt)
# This script uses input files from Zenodo ID: 
# Aim: Statistical analysis to compare each variable between type of diet, sex and/or generation 

# Used in analyses of the paper:
# "Long-term blood-free rearing of Anopheles mosquitoes with no effect on fitness, 
# Plasmodium infectivity nor microbiota composition."
# Authors: Joana Marques1*, Sofia G. Seabra1, Joana Gomes1, Ana Catarina Alves1, Henrique Silveira1*
#  1 Global Health and Tropical Medicine, GHTM, Associate Laboratory in Translation and Innovation Towards Global Health, LA-REAL, 
# Instituto de Higiene e Medicina Tropical, Universidade Nova de Lisboa, IHMT-NOVA, Rua da Junqueira 100, 1349-008 Lisboa, Portugal.
# *Corresponding authors

########################################################
# required R packages
library(openxlsx)
library(nlme)
library(FSA)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(dplyr)
library(plyr)
library(forcats)
library(cowplot)

## Set the working directory
setwd('')

###########################################################
###########################################################
# Proportion of feeding

df<-read.xlsx("./dados_bloodmeal_appetite.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

colnames(df)

sink("./R_output_appetite.txt", type=c("output", "message"), append = TRUE)

print(df)

print("########################")
# Perform t-test
t_test_result <- t.test(df[df$Diet == "B", "PFeeding"], df[df$Diet == "D", "PFeeding"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(PFeeding ~ Diet, data = df)
# Print the result
print(mwu_test_result)

print("########################")
# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(PFeeding ~ Generation, data = df)
# Print the result
print(kruskal_test_result)

print("########################")
print("# post-hoc Dunn test")
# Perform post-hoc Dunn test
dunn_test_result <- dunnTest(df$PFeeding, g = as.factor(df$Generation), method = "bonferroni")
# Print the post-hoc result
print(dunn_test_result)

print("########################")
print("ANCOVA model - dependent: PFeeding; categorical independent: Diet - ; covariate: Generation")
# Fit ANCOVA model - Diet - categorical independent; Generation - covariate (control for)
ancova_model <- lm(PFeeding ~ Diet + Generation, data = df)
# Print the summary of the model
summary(ancova_model)

########################################################
print("########################")
print("Linearity assumption between PFeeding and Generation for each Diet group")
# Separate data into two groups based on 'Diet' (B and D)
group_B <- subset(df, Diet == "B")
group_D <- subset(df, Diet == "D")

# Perform linear regression for Group B
lm_group_B <- lm(PFeeding ~ Generation, data = group_B)
summary(lm_group_B)

# Perform linear regression for Group D
lm_group_D <- lm(PFeeding ~ Generation, data = group_D)
summary(lm_group_D)

ggscatter(
  df, x = "Generation", y = "PFeeding",
  color = "Diet", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Diet)
  )

print("########################")
print("Homogeneity of regression slopes between the covariate and the grouping variable - interaction term should be not significant")
interaction_term<-df %>%anova_test(PFeeding ~ Diet*Generation)
print(interaction_term)

print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(PFeeding ~ Generation + Diet, data = df)
summary(model)

# Inspect the model diagnostic metrics
model.metrics <- augment(model) 
head(model.metrics, 3)
shapiro_test(model.metrics$.resid)

print("########################")
print("### Homogeneity of variances - ANCOVA assumes that the variance of the residuals is equal for all groups")
model.metrics %>% levene_test(.resid ~ Diet)

print("########################")
print("### Outliers - should not exist")
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

print("########################")
print("### Compute ANCOVA - covariate first")
res.aov <- df %>% anova_test(PFeeding ~  as.factor(Generation) + Diet)
print(res.aov)
res.aov <- df %>% anova_test(PFeeding ~  as.factor(Generation) + Diet + as.factor(Generation) * Diet)
print(res.aov)

res.aov <- df %>% anova_test(PFeeding ~  as.factor(Generation))
print(res.aov)

mymod <- aov(PFeeding ~  as.factor(Generation) + Diet , data=df)
summary(mymod)
TukeyHSD(mymod)

res.aov <- aov(PFeeding ~  as.factor(Diet), data=df)
summary(res.aov)

mymod <-  glm(PFeeding ~  as.factor(Generation) + Diet , data= df, family= "binomial")

library(emmeans)
em <- emmeans(mymod, "Generation")
contrast(em, "pairwise", adjust = "Sidak")

sink()


###########################################################
###########################################################
# LONGEVITY

df<-read.xlsx("./dados_longevidade.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

sink("./R_output_longevity.txt", type=c("output", "message"), append = TRUE)
print(df)
print("########################")
####
# Differences between Sex
t_test_result <- t.test(df[df$Sex == "M", "average_life_span"], df[df$Sex == "F", "average_life_span"])
# Print the result
print(t_test_result)

print("########################")
df %>% anova_test(average_life_span ~  Sex + Diet)

print("########################")
# Perform t-test
t_test_result <- t.test(df[df$Diet == "B", "average_life_span"], df[df$Diet == "D", "average_life_span"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(average_life_span ~ Diet, data = df)
# Print the result
print(mwu_test_result)

print("########################")
# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(average_life_span ~ Generation, data = df)
# Print the result
print(kruskal_test_result)

print("########################")
print("# post-hoc Dunn test")
# Perform post-hoc Dunn test
dunn_test_result <- dunnTest(df$average_life_span, g = as.factor(df$Generation), method = "bonferroni")
# Print the post-hoc result
print(dunn_test_result)

print("########################")
print("ANCOVA model - dependent: average_life_span; categorical independent: Diet - ; covariate: Generation")
# Fit ANCOVA model - Diet - categorical independent; Generation - covariate (control for)
ancova_model <- lm(average_life_span ~ Diet + Sex + Generation, data = df)
# Print the summary of the model
summary(ancova_model)

########################################################
print("########################")
print("Linearity assumption between average_life_span and Generation for each Diet group")
# Separate data into two groups based on 'Diet' (B and D)
group_B <- subset(df, Diet == "B")
group_D <- subset(df, Diet == "D")

# Perform linear regression for Group B
lm_group_B <- lm(average_life_span ~ Generation, data = group_B)
summary(lm_group_B)

# Perform linear regression for Group D
lm_group_D <- lm(average_life_span ~ Generation, data = group_D)
summary(lm_group_D)

ggscatter(
  df, x = "Generation", y = "average_life_span",
  color = "Diet", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Diet)
  )

print("########################")
print("Homogeneity of regression slopes between the covariate and the grouping variable - interaction term should be not significant")
interaction_term<-df %>%anova_test(average_life_span ~ Diet*Generation)
print(interaction_term)

print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(average_life_span ~ Generation + Diet, data = df)
summary(model)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(model.metrics, 3)
shapiro_test(model.metrics$.resid)

print("########################")
print("### Homogeneity of variances - ANCOVA assumes that the variance of the residuals is equal for all groups")
model.metrics %>% levene_test(.resid ~ Diet)

print("########################")
print("### Outliers - should not exist")
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

print("########################")
print("### Compute ANCOVA - covariate first")
res.aov <- df %>% anova_test(average_life_span ~  Diet + Sex + Generation)
print(res.aov)

sink()


###########################################################
###########################################################
# WING LENGTH

df<-read.xlsx("./dados_wing_length.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

sink("./R_output_wing_length.txt", type=c("output", "message"), append = TRUE)
print(df)
####
# Differences between Sex
t_test_result <- t.test(df[df$Sex == "M", "AVERAGE"], df[df$Sex == "F", "AVERAGE"])
# Print the result
print(t_test_result)

print("########################")
df %>% anova_test(AVERAGE ~  Sex + Diet)

print("########################")
# Perform t-test
t_test_result <- t.test(df[df$Diet == "B", "AVERAGE"], df[df$Diet == "D", "AVERAGE"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(AVERAGE ~ Diet, data = df)
# Print the result
print(mwu_test_result)

print("########################")
# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(AVERAGE ~ Generation, data = df)
# Print the result
print(kruskal_test_result)

print("########################")
print("# post-hoc Dunn test")
# Perform post-hoc Dunn test
dunn_test_result <- dunnTest(df$AVERAGE, g = as.factor(df$Generation), method = "bonferroni")
# Print the post-hoc result
print(dunn_test_result)

print("########################")
print("ANCOVA model - dependent: AVERAGE; categorical independent: Diet - ; covariate: Generation")
# Fit ANCOVA model - Diet - categorical independent; Generation - covariate (control for)
ancova_model <- lm(AVERAGE ~ Diet + Sex + Generation, data = df)
# Print the summary of the model
summary(ancova_model)

########################################################
print("########################")
print("Linearity assumption between AVERAGE and Generation for each Diet group")
# Separate data into two groups based on 'Diet' (B and D)
group_B <- subset(df, Diet == "B")
group_D <- subset(df, Diet == "D")

# Perform linear regression for Group B
lm_group_B <- lm(AVERAGE ~ Generation, data = group_B)
summary(lm_group_B)

# Perform linear regression for Group D
lm_group_D <- lm(AVERAGE ~ Generation, data = group_D)
summary(lm_group_D)

ggscatter(
  df, x = "Generation", y = "AVERAGE",
  color = "Diet", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Diet)
  )

print("########################")
print("Homogeneity of regression slopes between the covariate and the grouping variable - interaction term should be not significant")
interaction_term<-df %>%anova_test(AVERAGE ~ Diet*Generation)
print(interaction_term)

print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(AVERAGE ~ Generation + Diet, data = df)
summary(model)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(model.metrics, 3)
shapiro_test(model.metrics$.resid)

print("########################")
print("### Homogeneity of variances - ANCOVA assumes that the variance of the residuals is equal for all groups")
model.metrics %>% levene_test(.resid ~ Diet)

print("########################")
print("### Outliers - should not exist")
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

print("########################")
print("### Compute ANCOVA - covariate first")
res.aov <- df %>% anova_test(AVERAGE ~  Generation + Sex + Diet)
print(res.aov)

sink()

###########################################################
###########################################################
# INFECTION RATE

df<-read.xlsx("./dados_infection.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

colnames(df)

sink("./R_output_infection_rate.txt", type=c("output", "message"), append = TRUE)

print(df)

print("########################")
# Perform t-test
t_test_result <- t.test(df[df$Diet == "B", "Infection_Rate"], df[df$Diet == "D", "Infection_Rate"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(Infection_Rate ~ Diet, data = df)
# Print the result
print(mwu_test_result)

print("########################")
# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Infection_Rate ~ Generation, data = df)
# Print the result
print(kruskal_test_result)

print("########################")
print("# post-hoc Dunn test")
# Perform post-hoc Dunn test
dunn_test_result <- dunnTest(df$Infection_Rate, g = as.factor(df$Generation), method = "bonferroni")
# Print the post-hoc result
print(dunn_test_result)

print("########################")
print("ANCOVA model - dependent: Infection_Rate; categorical independent: Diet - ; covariate: Generation")
# Fit ANCOVA model - Diet - categorical independent; Generation - covariate (control for)
ancova_model <- lm(Infection_Rate ~ Diet + Generation, data = df)
# Print the summary of the model
summary(ancova_model)

########################################################
print("########################")
print("Linearity assumption between Infection_Rate and Generation for each Diet group")
# Separate data into two groups based on 'Diet' (B and D)
group_B <- subset(df, Diet == "B")
group_D <- subset(df, Diet == "D")

# Perform linear regression for Group B
lm_group_B <- lm(Infection_Rate ~ Generation, data = group_B)
summary(lm_group_B)

# Perform linear regression for Group D
lm_group_D <- lm(Infection_Rate ~ Generation, data = group_D)
summary(lm_group_D)

ggscatter(
  df, x = "Generation", y = "Infection_Rate",
  color = "Diet", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Diet)
  )

print("########################")
print("Homogeneity of regression slopes between the covariate and the grouping variable - interaction term should be not significant")
interaction_term<-df %>%anova_test(Infection_Rate ~ Diet*Generation)
print(interaction_term)

print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(Infection_Rate ~ Generation + Diet, data = df)
summary(model)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(model.metrics, 3)
shapiro_test(model.metrics$.resid)

print("########################")
print("### Homogeneity of variances - ANCOVA assumes that the variance of the residuals is equal for all groups")
model.metrics %>% levene_test(.resid ~ Diet)

print("########################")
print("### Outliers - should not exist")
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

print("########################")
print("### Compute ANCOVA - covariate first")
res.aov <- df %>% anova_test(Infection_Rate ~  Generation + Diet)
print(res.aov)

sink()


###########################################################
###########################################################
# INFECTION INTENSITY

df<-read.xlsx("./dados_infection.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

sink("./R_output_infection_intensity.txt", type=c("output", "message"), append = TRUE)
print(df)

print("########################")
# Perform t-test
t_test_result <- t.test(df[df$Diet == "B", "Infection_Intensity"], df[df$Diet == "D", "Infection_Intensity"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(Infection_Intensity ~ Diet, data = df)
# Print the result
print(mwu_test_result)

print("########################")
# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Infection_Intensity ~ Generation, data = df)
# Print the result
print(kruskal_test_result)

print("########################")
print("# post-hoc Dunn test")
# Perform post-hoc Dunn test
dunn_test_result <- dunnTest(df$Infection_Intensity, g = as.factor(df$Generation), method = "bonferroni")
# Print the post-hoc result
print(dunn_test_result)

print("########################")
print("ANCOVA model - dependent: Infection_Intensity; categorical independent: Diet - ; covariate: Generation")
# Fit ANCOVA model - Diet - categorical independent; Generation - covariate (control for)
ancova_model <- lm(Infection_Intensity ~ Diet + Generation, data = df)
# Print the summary of the model
summary(ancova_model)

########################################################
print("########################")
print("Linearity assumption between Infection_Intensity and Generation for each Diet group")
# Separate data into two groups based on 'Diet' (B and D)
group_B <- subset(df, Diet == "B")
group_D <- subset(df, Diet == "D")

# Perform linear regression for Group B
lm_group_B <- lm(Infection_Intensity ~ Generation, data = group_B)
summary(lm_group_B)

# Perform linear regression for Group D
lm_group_D <- lm(Infection_Intensity ~ Generation, data = group_D)
summary(lm_group_D)

ggscatter(
  df, x = "Generation", y = "Infection_Intensity",
  color = "Diet", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Diet)
  )

print("########################")
print("Homogeneity of regression slopes between the covariate and the grouping variable - interaction term should be not significant")
interaction_term<-df %>%anova_test(Infection_Intensity ~ Diet*Generation)
print(interaction_term)

print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(Infection_Intensity ~ Generation + Diet, data = df)
summary(model)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(model.metrics, 3)
shapiro_test(model.metrics$.resid)

print("########################")
print("### Homogeneity of variances - ANCOVA assumes that the variance of the residuals is equal for all groups")
model.metrics %>% levene_test(.resid ~ Diet)

print("########################")
print("### Outliers - should not exist")
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

print("########################")
print("### Compute ANCOVA - covariate first")
res.aov <- df %>% anova_test(Infection_Intensity ~  Generation + Diet)
print(res.aov)

sink()

##############################################
##############################################
# Meal preference

df<-read.xlsx("./Dados_meal_preference.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

sink("./R_output_meal_preference.txt", type=c("output", "message"), append = TRUE)

print(df)
print("# New dataset with blue only")
df_new <- df[df$Outcome == "blue", ]
print(df_new)

# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(PFeeding ~ Diet, data = df_new)
# Print the result
print(mwu_test_result)

print("# New dataset with red only")
df_new <- df[df$Outcome == "red", ]
print(df_new)

# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(PFeeding ~ Diet, data = df_new)
# Print the result
print(mwu_test_result)

print("# New dataset with none only")
df_new <- df[df$Outcome == "none", ]
print(df_new)

# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(PFeeding ~ Diet, data = df_new)
# Print the result
print(mwu_test_result)

sink()


############################################################################
##########################################################################
############################################################################
##########################################################################
# PLOTS MANUSCRIPT
############################################################################
##########################################################################



###########################################################
# Plot Pfeeding + Meal preference

png(file = "./plot_feeding_mealpreference.png", bg="white",width = 8, height = 4,  units = "in", res = 300)

df<-read.xlsx("./dados_bloodmeal_appetite.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

colnames(df)

df$Generation <- factor(df$Generation, levels = unique(df$Generation))

df$Diet <- factor(df$Diet, levels = unique(df$Diet))

# Calculate the average PFeeding values and standard deviation
avg_cq <- aggregate(PFeeding ~ Generation + Diet, data = df, FUN = mean)
std_dev <- aggregate(PFeeding ~ Generation + Diet, data = df, FUN = sd)

# Merge average PFeeding and standard deviation data
avg_stdev <- merge(avg_cq, std_dev, by = c("Generation", "Diet"))


# Create the grouped bar chart

p1 <- ggplot(avg_stdev, aes(x = Generation, y = PFeeding.x, fill = Diet)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(y = PFeeding.x, ymin = PFeeding.x + PFeeding.y, ymax = PFeeding.x + PFeeding.y), 
                position = position_dodge(width = 0.9), width = 0.25, colour="#636363") +
  geom_linerange(aes(ymin = PFeeding.x, ymax = PFeeding.x + PFeeding.y),position = position_dodge(width = 0.9), colour="#636363") +
  scale_fill_manual(name = "Type of diet", labels = c("Blood", "BLOODless"), values =c("#d6604d", "#fee090")) +
  labs(title = "",
       x = "Generation", y = "Percentage of females feeding \n on mouse",tag="A") +
  theme_minimal()


##############################################
##############################################
# Meal preference

df<-read.xlsx("./Dados_meal_preference.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

colnames(df)


#######################
# Plot Meal preference

df$Diet <- factor(df$Diet, levels = unique(df$Diet))

df$Outcome <- factor(df$Outcome, levels = c("red", "blue", "none"))

# Create the grouped bar chart
p2 <- ggplot(df, aes(x = Replicate, y = PFeeding, fill=Outcome)) +
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual(name = "Preference", labels = c("Blood", "BLOODless","none"), values =c("#d6604d", "#fee090","#bdbdbd")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "",
       x = "Experiments", y = "Percentage of Feeding",tag="B") +
  facet_wrap( ~ Diet, ncol=2, strip.position = "top") +
  theme_minimal() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

p_all <- plot_grid(p1, p2, nrow=1, ncol=2)

p_all

dev.off()


###########################################################
###########################################################
# LONGEVITY + WING LENGTH
png(file = "./plot_longevity_winglength.png", bg="white",width = 7, height = 5,  units = "in", res = 300)

df<-read.xlsx("./dados_longevidade.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)
#######################
# Plot Longevity
df$Generation <- factor(df$Generation, levels = unique(df$Generation))

df$Diet <- factor(df$Diet, levels = unique(df$Diet))

df$Sex <- factor(df$Sex, levels = unique(df$Sex))
levels(df$Sex)
levels(df$Sex) <- c("Female", "Male")
levels(df$Sex)

# Create the grouped bar chart

p1 <- ggplot(df, aes(x = Generation, y = average_life_span, fill = Diet)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(y = average_life_span, ymin = average_life_span + SD, ymax = average_life_span + SD), 
                position = position_dodge(width = 0.9), width = 0.25, colour="#636363") +
  geom_linerange(aes(ymin = average_life_span, ymax = average_life_span + SD),position = position_dodge(width = 0.9), colour="#636363") +
  scale_fill_manual(name = "Type of diet", labels = c("Blood", "BLOODless"), values =c("#d6604d", "#fee090")) +
  labs(title = "",
       x = "Generation", y = "Longevity (number of days)",tag="A") +
  facet_wrap( ~ Sex, ncol=1, strip.position = "top") +
  theme_minimal() +  theme(legend.position="none") 

###########################################################
###########################################################
# WING LENGTH

df<-read.xlsx("./dados_wing_length.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

#######################
# Plot Wing length

df$Generation <- factor(df$Generation, levels = unique(df$Generation))

df$Diet <- factor(df$Diet, levels = unique(df$Diet))

df$Sex <- factor(df$Sex, levels = unique(df$Sex))
levels(df$Sex)
levels(df$Sex) <- c("Female", "Male")
levels(df$Sex)

# Create the grouped bar chart
p2 <- ggplot(df, aes(x = Generation, y = AVERAGE, fill = Diet)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(y = AVERAGE, ymin = AVERAGE + STDEV, ymax = AVERAGE + STDEV), 
                position = position_dodge(width = 0.9), width = 0.25, colour="#636363") +
  geom_linerange(aes(ymin = AVERAGE, ymax = AVERAGE + STDEV),position = position_dodge(width = 0.9), colour="#636363") +
  scale_fill_manual(name = "Type of diet", labels = c("Blood", "BLOODless"), values =c("#d6604d", "#fee090")) +
  labs(title = "",
       x = "Generation", y = "Wing length",tag="B") +
  facet_wrap( ~ Sex, ncol=1, strip.position = "top") +
  theme_minimal() +  theme(legend.position="none") 

p_all <- plot_grid(p1, p2, nrow=1, ncol=2,
                   label_size = 14,
                   label_x = 0.92, label_y = 0.75,
                   hjust = -0.5, vjust = -0.5,
                   align="v")

legend <- get_legend(p1 + 
                       guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "bottom"))


plot_grid(p_all, legend, ncol = 1, rel_heights = c(1, .1))

dev.off()


############################
############################
# PLOT Infection rate + Infection intensity

df<-read.xlsx("./dados_infection.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

png(file = "./plot_infection_rate_and_intensity.png", bg="white",width = 5, height = 8,  units = "in", res = 300)

p1 <- ggplot(df, aes(x = as.factor(Generation), y = Infection_Rate, fill = Diet)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Type of meal", labels = c("Blood", "BLOODless"), values =c("#d6604d", "#fee090")) +
  labs(title = "",
       x = "Generation", y = "Infection rate",tag="A") +
  theme_minimal() +  theme(legend.position="none") 


p2 <- ggplot(df, aes(x = as.factor(Generation), y = Infection_Intensity, fill = Diet)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Type of meal", labels = c("Blood", "BLOODless"), values =c("#d6604d", "#fee090")) +
  labs(title = "",
       x = "Generation", y = "Infection intensity",tag="B") +
  theme_minimal() +  theme(legend.position="none") 

library(cowplot)
p_all <- plot_grid(p1, p2, nrow=2, ncol=1,
                   label_size = 14,
                   label_x = 0.92, label_y = 0.75,
                   hjust = -0.5, vjust = -0.5,
                   align="v")

legend <- get_legend(p1 + 
                       guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "bottom"))

plot_grid(p_all, legend, ncol = 1, rel_heights = c(1, .1))


get_legend_custom <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract the legend from p1
legend <- get_legend_custom(p1)

# Display the legend
plot_grid(NULL, legend, ncol = 1, rel_heights = c(1, 0.1))

dev.off()



#################
# NEW PLOT Infection rate and intensity (individual dots)
#################
png(file = "./plot_infection_rate_and_intensity_jitter.png", bg="white",width = 5, height = 6,  units = "in", res = 300)

# Infection rate
df<-read.xlsx("./dados_infection.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

p1 <- ggplot(df, aes(x = as.factor(Generation), y = Infection_Rate, fill = Diet)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Type of meal", labels = c("Blood", "BLOODless"), values =c("#d6604d", "#fee090")) +
  labs(title = "",
       x = "Generation", y = "Infection rate",tag="A") +
  theme_minimal() +  theme(legend.position="none") 


df<-read.xlsx("./dados_oocyst_number.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

# Calculate the average OocystNumber values and standard deviation
avg_cq <- aggregate(OocystNumber ~ Generation + Diet, data = df, FUN = mean)
std_dev <- aggregate(OocystNumber ~ Generation + Diet, data = df, FUN = sd)

# Merge average PFeeding and standard deviation data
avg_stdev <- merge(avg_cq, std_dev, by = c("Generation", "Diet"))

colnames(avg_stdev)=c("Generation","Diet","Oocyst_mean", "Ooscyt_sd")

avg_stdev$Diet <- factor(avg_stdev$Diet, levels = unique(avg_stdev$Diet))
levels(avg_stdev$Diet)
avg_stdev$Generation <- factor(avg_stdev$Generation, levels = unique(avg_stdev$Generation))
levels(avg_stdev$Generation)
avg_stdev$Oocyst_mean <- as.numeric(avg_stdev$Oocyst_mean)
avg_stdev$Ooscyt_sd <- as.numeric(avg_stdev$Ooscyt_sd)


df$Group <- factor(df$Group, levels = c("B_F23","D_F23","B_F27","D_F27","B_F38","D_F38"))

p2 <- ggplot(df, aes(x = as.factor(Generation), y = OocystNumber, fill = as.factor(Diet))) +
  geom_jitter(position= position_jitterdodge(jitter.width = 0.5, dodge.width = 0.6),pch=21, stroke=0, size=2) +
  scale_fill_manual(name = "Type of diet", labels = c("Blood", "BLOODless"), values =c("#d6604d", "#fee090"))+
  stat_summary(fun=median, geom ="crossbar", position = position_dodge(width = 0.6), 
               width = 0.4, size=0.2, show.legend = FALSE) +
  labs(title = "",
       x = "Generation", y = "Oocyst Number",tag="B") +
  theme_minimal() +  theme(legend.position="bottom") 


library(cowplot)
p_all <- plot_grid(p1, p2, nrow=2, ncol=1,
                   label_size = 14,
                   label_x = 0.92, label_y = 0.75,
                   hjust = -0.5, vjust = -0.5,
                   align="v", rel_heights = c(1, 1))
p_all

dev.off()



