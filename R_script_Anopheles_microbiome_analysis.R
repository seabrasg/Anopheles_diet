##########
# Author: Sofia Seabra (sgseabra@ihmt.nl.pt)
# This script uses input files from Zenodo ID: 
# Aim: Study the taxonomic composition and diversity of the microbial communities associated with 
# mid-guts and salivary glands of Anopheles stephensi mosquitoes fed either on blood or on artificial blood-free diet.

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
library(forcats)
library(ape)
library("metacoder")
library(vegan)
library(plyr)
library(phyloseq)
library(cowplot)


#################
# Set working directory
setwd('')

#################
#### PLOT RELATIVE ABUNDANCES OF TAXA - column Taxonomy_ID with classification until the ORDER level

# Input: abundance matrix (OTU table) from Qiime2
df_tax<-read.xlsx("./table_OTU.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df_tax)
nrow(df_tax)

# Reshape data for ggplot
df_long <- tidyr::gather(df_tax, key = "Sample", value = "Abundance", 10:Ng, factor_key=TRUE)

#### Add columns with: 
#diet (artificial or blood control) 

df_long <- df_long %>%
  mutate(Type_sample = NA, Type_diet = NA,Type_tissue = NA, Sex = NA, Replicates = NA)

df_long$Type_sample[df_long$Sample %in% c("1")]<-"Blood_Females_SG_1"

df_long$Type_sample[df_long$Sample %in% c("2")]<-"Blood_Females_MG_1"

df_long$Type_sample[df_long$Sample %in% c("3")]<-"Blood_Females_SG_2"

df_long$Type_sample[df_long$Sample %in% c("4")]<-"Blood_Females_SG_3"

df_long$Type_sample[df_long$Sample %in% c("5")]<-"Blood_Females_MG_2"

df_long$Type_sample[df_long$Sample %in% c("6")]<-"Blood_Females_MG_3"

df_long$Type_sample[df_long$Sample %in% c("7")]<-"Blood_Males_MG_1"
df_long$Type_sample[df_long$Sample %in% c("8")]<-"Blood_Males_MG_2"
df_long$Type_sample[df_long$Sample %in% c("9")]<-"Blood_Males_MG_3"

df_long$Type_sample[df_long$Sample %in% c("10")]<-"Diet_Females_SG_1"
df_long$Type_sample[df_long$Sample %in% c("11")]<-"Diet_Females_SG_2"
df_long$Type_sample[df_long$Sample %in% c("12")]<-"Diet_Females_SG_3"

df_long$Type_sample[df_long$Sample %in% c("13")]<-"Diet_Females_MG_1"
df_long$Type_sample[df_long$Sample %in% c("14")]<-"Diet_Females_MG_2"
df_long$Type_sample[df_long$Sample %in% c("15")]<-"Diet_Females_MG_3"

df_long$Type_sample[df_long$Sample %in% c("16")]<-"Diet_Males_MG_1"
df_long$Type_sample[df_long$Sample %in% c("17")]<-"Diet_Males_MG_2"
df_long$Type_sample[df_long$Sample %in% c("18")]<-"Diet_Males_MG_3"

df_long$Type_sample[df_long$Sample %in% c("Ng")]<-"Negative_control"
df_long$Type_sample[df_long$Sample %in% c("Ps")]<-"Positive_control"


df_long$Type_diet[df_long$Sample %in% c("10","11","12",
                                        "13","14","15",
                                        "16","17","18")]<-"Diet"

df_long$Type_diet[df_long$Sample %in% c("1","2","3",
                                        "4","5","6",
                                        "7","8","9")]<-"Blood"

#type (salivary glands or midguts), 
df_long$Type_tissue[df_long$Sample %in% c("1","3","4",
                                          "10","11","12")]<-"SG"

df_long$Type_tissue[df_long$Sample %in% c("2","5","6",
                                          "7","8","9",
                                          "13","14","15","16",
                                          "17","18")]<-"MG"

#sex (female or male)
df_long$Sex[df_long$Sample %in% c("1","2","3",
                                  "4","5","6",
                                  "10","11","12",
                                  "13","14","15")]<-"Females"

df_long$Sex[df_long$Sample %in% c("7","8","9",
                                  "16","17","18")]<-"Males"

#replicates (1,2 or 3)
df_long$Replicates[df_long$Sample %in% c("1","2","7",
                                         "10","13","16")]<-"1"

df_long$Replicates[df_long$Sample %in% c("3","5","8",
                                         "11","14","17")]<-"2"

df_long$Replicates[df_long$Sample %in% c("4","6","9",
                                         "12","15","18")]<-"3"

table(df_long$Type_sample, df_long$Sample, useNA = "always")
table(df_long$Type_diet, df_long$Sample, useNA = "always")
table(df_long$Type_tissue, df_long$Sample, useNA = "always")
table(df_long$Sex, df_long$Sample, useNA = "always")
table(df_long$Replicates, df_long$Sample, useNA = "always")



#############
# Group by and sum abundances

df_grouped <- df_long %>% 
  group_by(Taxonomy_ID,Type_diet, Type_tissue, Sex, Replicates) %>% 
  summarise(Abundance = sum(Abundance))

# Filter rows with abundance greater than or equal to 1000
df_filtered <- df_grouped %>%
  filter(`Abundance` >= 1000) %>%
  filter( !is.na(`Type_diet`)) 

df_filtered$Type_diet <- recode(df_filtered$Type_diet, Diet = 'BLOODless', 
                                Blood = 'Blood')

n <- length(unique(df_filtered$Taxonomy_ID))


library(randomcoloR)
color_palette <- distinctColorPalette(n)
color_palette_total <- c( "#E39551", "#DAE29E","#8085C7",  "#81E658", "#C4AA91", "#DA70BC",
                          "#D5D852",  "#C743E3", "#8C66DD", "#7DDC91", "#7FBBD7", "#6CE0CE","#DCB7D5","#D2606E","#CDE4DA")


# Create a percent stacked barplot
png(file = "./plot_OTU_abundance.png", bg="white",width = 10, height = 8,  units = "in", res = 300)

ggplot(df_filtered, aes(x = Replicates, y = Abundance, fill = fct_reorder(Taxonomy_ID, Abundance, .desc=TRUE))) +
  geom_bar(position="fill", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = color_palette_total) +
  labs(title = "",
       x = "Samples",
       y = "Relative frequency") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  theme(legend.title=element_blank()) +
  facet_grid(Sex+Type_tissue~Type_diet, scales="free", space="free_x")

dev.off()



###############################################
##############################################
## PLOT DIVERSITY MICROBIOME - OTU + Simpson + Shannon

df<-read.xlsx("./table_diversity.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)
nrow(df)

png(file = "./plot_diversity.png", bg="white",width = 4, height = 8,  units = "in", res = 300)

# To change labels in facet
diet_labs <- c("BLODDless", "Blood")
names(diet_labs) <- c( "artificial", "control")

body_labs <- c("Midguts", "Salivary Glands")
names(body_labs) <- c("MidGuts","Salivary_Gland")

sex_labs <- c("Females", "Males")
names(sex_labs) <- c("females","males")

p1 <- ggplot(df, aes(x = Diet, y = observed_features), fill= Diet) +
  geom_point(aes(col = Diet), size=3) +
  geom_line(aes(col = Diet), size=1.5) +
  scale_colour_manual(name="Type of meal", labels=c("BLOODless","Blood"),values = c("#fee090","#d6604d")) +
  theme_bw() +
  labs(title = "",
       x = "",
       y = "Number of OTUs",tag="A") +
  facet_grid(~ Sex + `body-site`, scales="free", space="free_x",labeller = labeller(Sex = sex_labs, `body-site` = body_labs)) +
  scale_y_continuous(limits = c(100, 300)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_blank()) +  
  theme(legend.position="none") 

p2 <- ggplot(df, aes(x = Diet, y = shannon_entropy), fill= Diet) +
  geom_point(aes(col = Diet), size=3) +
  geom_line(aes(col = Diet), size=1.5) +
  scale_colour_manual(name="Type of meal", labels=c("BLOODless","Blood"),values = c("#fee090","#d6604d")) +
  theme_bw() +
  labs(title = "",
       x = "",
       y = "Shannon Entropy",tag="B") +
  facet_grid(~ Sex + `body-site`, scales="free", space="free_x",labeller = labeller(Sex = sex_labs, `body-site` = body_labs)) +
  scale_y_continuous(limits = c(0.5, 2.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_blank()) +  
  theme(legend.position="none") 


p3 <- ggplot(df, aes(x = Diet, y = simpson), fill= Diet) +
  geom_point(aes(col = Diet), size=3) +
  geom_line(aes(col = Diet), size=1.5) +
  scale_colour_manual(name="Type of meal", labels=c("BLOODless","Blood"),values = c("#fee090","#d6604d")) +
  theme_bw() +
  labs(title = "",
       x = "",
       y = "Simpson index",tag="C") +
  facet_grid(~ Sex + `body-site`, scales="free", space="free_x",labeller = labeller(Sex = sex_labs, `body-site` = body_labs)) +
  scale_y_continuous(limits = c(0.1, 0.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_blank()) +  
  theme(legend.position="none") 


p_all <- plot_grid(p1, p2, p3, nrow=3, ncol=1,
                   label_size = 14,
                   label_x = 0.92, label_y = 0.75,
                   hjust = -0.5, vjust = -0.5,
                   align="v")

legend <- get_legend(p1 + 
                       guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "bottom"))

plot_grid(p_all, legend, ncol = 1, rel_heights = c(1, .1))

dev.off()

#################
# Statistical analysis

df<-read.xlsx("./table_diversity.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df)

sink("./R_output_diversity.txt", type=c("output", "message"), append = TRUE)

print(df)

### ### ### ### ### ### 
### Observed features
print("########################################################################")
print("# Observed features")
# Perform t-test
t_test_result <- t.test(df[df$Diet == "control", "observed_features"], df[df$Diet == "artificial", "observed_features"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(observed_features ~ Diet, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(observed_features ~ Sex, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(observed_features ~ `body-site`, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)


print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(observed_features ~ Diet + Sex + `body-site`, data = df)
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

print("# Factorial ANOVA")
lmodel <- aov(observed_features ~ Diet + Sex + `body-site`, data = df)
# Print the summary of the model
summary(lmodel)


### ### ### ### ### ### 
### Shannon index
print("########################################################################")
print("# Shannon index")

# Perform t-test
t_test_result <- t.test(df[df$Diet == "control", "shannon_entropy"], df[df$Diet == "artificial", "shannon_entropy"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(shannon_entropy ~ Diet, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(shannon_entropy ~ Sex, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(shannon_entropy ~ `body-site`, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)

print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(shannon_entropy ~ Diet + Sex + `body-site`, data = df)
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

print("# Factorial ANOVA")
lmodel <- aov(shannon_entropy ~ Diet + Sex + `body-site`, data = df)
# Print the summary of the model
summary(lmodel)

### ### ### ### ### ### 
### Simpson index
print("########################################################################")
print("# Simpson index")

# Perform t-test
t_test_result <- t.test(df[df$Diet == "control", "simpson"], df[df$Diet == "artificial", "simpson"])
# Print the result
print(t_test_result)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(simpson ~ Diet, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(simpson ~ Sex, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)

print("########################")
# Perform Mann-Whitney U test
mwu_test_result <- wilcox.test(simpson ~ `body-site`, data = df)
# Print the result
print(mwu_test_result)
print("bonferroni adjusted p-value")
p.adjust(mwu_test_result$p.value, method = "bonferroni", n= 3)


print("########################")
print("Normality of residuals - Shapiro Wilk test should be not significant")
# Fit the model, the covariate goes first
model <- lm(simpson ~ Diet + Sex + `body-site`, data = df)
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

print("# Factorial ANOVA")
lmodel <- aov(simpson ~ Diet + Sex + `body-site` , data = df)
# Print the summary of the model
summary(lmodel)

sink()


########################
# Principal Coordinate Analysis and PERMANOVA

# Adapted from: https://rstudio-pubs-static.s3.amazonaws.com/343284_cbadd2f3b7cd42f3aced2d3f42dc6d33.html#introduction

# Input: abundance matrix, with samples in columns and Operational Taxonomic Units (OTUs) in rows - output OTU from qiime2

df_otu<-read.xlsx("./table_OTU.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(df_otu)
nrow(df_otu)
head(df_otu)

# Remove rows with only "Bacteria" assignment
df_otu <- filter(df_otu,
                 !Family == " ")
colnames(df_otu)
nrow(df_otu)
unique(df_otu$Taxonomy_ID)

### Get the OTU TABLE
# Remove columns
df_otu_filtered <- select (df_otu,
                           !c(Taxonomy_ID:Species,Ps,Ng))
colnames(df_otu_filtered)
row.names(df_otu_filtered)

# transpose (OTU_ID in columns and samples in rows)
library(data.table)
table_otu<- data.table::transpose(df_otu_filtered, keep.names = "OTU_ID", make.names = "OTU_ID")
colnames(table_otu)
row.names(table_otu) 
table_otu[,"OTU_ID"]
row.names(table_otu) <- table_otu$OTU_ID  # now the OTU_ID column has the sample ids and we will put them as row names
row.names(table_otu) 


### Get the TAXONOMY TABLE (columns: OTU_ID and Taxonomy)
table_taxonomy <- select (df_otu,
                          c(OTU_ID:Species))

colnames(table_taxonomy)
nrow(table_taxonomy)
row.names(table_taxonomy)
row.names(table_taxonomy) <- table_taxonomy$OTU_ID
row.names(table_taxonomy)

# Remove all the OTUs that don’t occur in our table_otu data set
table_taxonomy = table_taxonomy[row.names(table_taxonomy) %in% colnames(table_otu),]
nrow(table_taxonomy)

# remove unecessary columns (OTU_ID is now the row.names)
table_taxonomy <- select (table_taxonomy,
                          !c(OTU_ID))


##### GET the Metadata

meta <-read.xlsx("./table_sample_id.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
row.names(meta) <- meta$sample_id
colnames(meta)

# Order OTU and meta files by sample id
otu.clean <- table_otu[order(row.names(table_otu)),]
row.names(otu.clean) 

meta.clean = meta[order(row.names(meta)),]
row.names(meta.clean) 

# convert to abundance values to numeric
otu.clean <- otu.clean %>% mutate_if(is.character, as.numeric)

# OTU-based metrics
BC.nmds = metaMDS(otu.clean, distance="bray", k=2, trymax=1000)

BC.nmds$points

df_BC_nmds<-data.frame(BC.nmds$points)

df_BC_nmds$sample<-meta.clean$Sex_Diet_Site
levels(as.factor(df_BC_nmds$sample))


######################################################################################
# Principal Coordinates Analysis
library(ape)

pcoa_BC <- pcoa(BC.dist, correction="none", rn=NULL)

pcoa_BC$values["Relative_eig"]
pcoa_BC_df <- data.frame(pcoa_BC$vectors[,1:3])
colnames(pcoa_BC_df) <- c("PCo1", "PCo2","PCo3")

meta.clean$Sex_Diet_Site <- recode_factor(meta.clean$Sex_Diet_Site, F_artificial_MG="Females_Midguts_BLOODless",
                                          F_control_SG="Females_SalivaryGlands_Blood",
                                          F_artificial_SG="Females_SalivaryGlands_BLOODless",
                                          M_artificial_MG="Males_Midguts_BLOODless",
                                          F_control_MG="Females_Midguts_Blood",
                                          M_control_MG="Males_Midguts_Blood")
meta.clean$Sex_Diet_Site <- factor(meta.clean$Sex_Diet_Site, levels= c("Females_Midguts_BLOODless","Females_Midguts_Blood",
                                                                       "Females_SalivaryGlands_BLOODless","Females_SalivaryGlands_Blood",
                                                                       "Males_Midguts_BLOODless","Males_Midguts_Blood"))

pcoa_BC_df$sample<-meta.clean$Sex_Diet_Site # give names of samples - meta.clean has the same order of rows
levels(as.factor(pcoa_BC_df$sample))

write.csv(pcoa_BC_df, file ="./pcoa_BC_coordinates.csv")

# plot coordinates 

png(file = "./plot_pcoa_braycurtis.png", bg="white",width = 6, height = 4,  units = "in", res = 300)

ggplot(pcoa_BC_df, aes(x = PCo1, y = PCo2)) + 
  geom_point(aes(shape =sample,color= sample),size = 4) +
  
  scale_shape_manual(values= c(16,16,
                               10,10,
                               17,17)) +
  scale_colour_manual(values = c("#e3b94b","#d6604d",
                                 "#e3b94b","#d6604d",
                                 "#e3b94b","#d6604d")) +
  labs(title = "",
       x = "PCo1 (40.8%)", y = "PCo2 (23.9%)",tag="") +
  theme(legend.position = "top",legend.title=element_blank) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  theme_minimal()

dev.off()


######################################################################################
# PERMANOVA 
BC.dist=vegdist(otu.clean, distance="bray")

sink(file="./PERMANOVA_Bray_curtis_distances.txt")
BC.dist
adonis2(BC.dist ~ Diet + Sex + bodysite, data = meta.clean, permutations = 1000)

sink()




