---
title: "SPAH eoPRED score and placental pathologies"
author: "Hannah Illing"
date: "2026 February 02"
#Date created: 2025,08,08
output: 
  html_document:
    keep_md: yes
    code_folding: show
    toc: true  
    toc_depth: 4
    toc_float: 
      collapsed: true 
      smooth_scroll: true
---

## Set-up

### Load Packages 



### Style Guide



### Load Data


``` r
metadata <- readRDS(here::here("D. Pathology Project/02_Outputs/AA. Cleaned/SPAH_metadata_ss_filtered.rds"))
dim(metadata) #493 x 162
```

```
## [1] 493 152
```

``` r
#fix some variable classes
metadata <- metadata %>%
  mutate(Sex_Predic = as.factor(Sex_Predic),
         primi_gravida = as.factor(pcr_gravida == "1"),
         primi_para = as.factor(parity == "0"),
         Placenta_SGAvAGAvLGA = as.factor(Placenta_SGAvAGAvLGA),
         Feto_placental_wt = pcr_birth_wt/Placental_wt,
         pcr_diabetes_any = as.factor(pcr_diabetes_any),
         Race_AmIndALNative = as.factor(Race_AmIndALNative),
         Race_Asian = as.factor(Race_Asian),
         Race_Black = as.factor(Race_Black),
         Race_HINativePacIsl = as.factor(Race_HINativePacIsl), 
         Race_white = as.factor(Race_white),
         Race_Other = as.factor(Race_Other),
         Ethnicity_HispYN = as.factor(Ethnicity_HispYN),
         acute_inflammation_3cat = as.factor(acute_inflammation_3cat),
         chronic_inflammation_3cat = as.factor(chronic_inflammation_3cat)) %>% 
#note that one sample is mistakenly labelled as FVM in the FETAL_VASC_PATH column
#Use the FETAL_VASC_PATH column to create an updated FETAL_VASC_PATH variable
  mutate(FETAL_VASC_PATH = as.factor(ifelse(Fetal_vasc_path_3cat==0,0,1)))
```

## Figure 1 


``` r
#Does the cutoff accurately assign EOPE class to cases of EOPE?
table(metadata$PE_Status, metadata$PE_GA_state) #5/7 EOPE, 9/42 LOPE, 5/40 nPTB, 14/402 nTB
```

```
##               
##                EOPE LOPE nPTB nTB
##   EOPE            5    9    5  14
##   Normotensive    2   33   35 388
```

``` r
stat.test.pegastate <- dunn_test(EOPE~PE_GA_state, data= metadata)
stat.test.pegastate <- stat.test.pegastate %>%
  mutate(p.adj=round(p.adj, digits=3))%>%
  filter(p.adj<0.05)%>%
  mutate(y.position = c(0.89,0.95,0.81),
         PE_GA_state = group1,
         p.adj = case_when(p.adj<0.001 ~ "p < 0.001", 
                           p.adj>=0.05 ~ "",
                           .default = paste("p = ",p.adj)))

fig_1_a <- metadata %>%
  filter(!is.na(PE_GA_state))%>%
  ggplot(aes(x=PE_GA_state, y = EOPE, color=PE_GA_state, fill=PE_GA_state))+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2), size=2.5)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  #geom_hline(yintercept = 0.55, color="grey", linetype="dashed", linewidth=1)+
  scale_color_manual(values = categories)+
  scale_fill_manual(values = categories)+
  labs(y = "eoPRED score", x = "Clinical Diagnosis")+
  stat_pvalue_manual(stat.test.pegastate, label="{p.adj}", tip.length=0.01)+
  stat_compare_means(aes(group=PE_GA_state), label.x = 0.95)+
  scale_x_discrete(labels=c("EOPE"="EOPE","LOPE"="LOPE","nPTB"="PTB","nTB"="TB"))+
  ylim(0,1)
fig_1_a
```

![](03_SPAH_eoPRED_score_2026_files/figure-html/fig-1-part-1-1.png)<!-- -->

``` r
fig_1_b <- metadata %>%
  filter(!is.na(new_pcr_hypertension_gestation)) %>%
  filter(PE_GA_state %in% c("nPTB","nTB"))%>%
  ggplot(aes(x=PE_GA_state, y = EOPE, color=new_pcr_hypertension_gestation, fill=new_pcr_hypertension_gestation))+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2), size=2.5)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  labs(y = "eoPRED score", x = "Gestational Age")+
  scale_x_discrete(labels=c("nPTB"="Preterm","nTB"="Term"))+
  scale_fill_discrete(labels=c("1"="Gestational Hypertension","0"="Normotensive"),
                      type =c("grey55","grey22"))+
  scale_color_discrete(labels=c("1"="Gestational Hypertension","0"="Normotensive"),
                       type = c("grey55","grey22"))+
  stat_compare_means(label.x = 0.95)+
  theme(legend.position.inside = c(0.17,0.11),
        legend.position="inside", legend.background= element_rect(color = "grey"),
        legend.title = element_blank())+
  ylim(0,1)
fig_1_b
```

![](03_SPAH_eoPRED_score_2026_files/figure-html/fig-1-part-1-2.png)<!-- -->


``` r
## Is eoPRED score higher in all samples with MVM?
wilcox_test(EOPE ~ mvm_di, data = metadata)
```

```
## # A tibble: 1 x 9
##   .y.   group1 group2    n1    n2 statistic           p       p.adj p.adj.signif
## * <chr> <chr>  <chr>  <int> <int>     <dbl>       <dbl>       <dbl> <chr>       
## 1 EOPE  0      1        340   144     17330 0.000000373 0.000000373 ****
```

``` r
metadata %>%
  group_by(mvm_di) %>%
  summarise_at(vars(EOPE), list(MEAN = mean))
```

```
## # A tibble: 3 x 2
##   mvm_di  MEAN
##   <fct>  <dbl>
## 1 0      0.347
## 2 1      0.418
## 3 <NA>   0.353
```

``` r
## Is eoPRED score higher in all samples with MVM? (without EOPE)
wilcox_test(EOPE ~ mvm_di, data = (metadata %>%
              filter(!PE_GA_state == "EOPE")))
```

```
## # A tibble: 1 x 9
##   .y.   group1 group2    n1    n2 statistic          p      p.adj p.adj.signif
## * <chr> <chr>  <chr>  <int> <int>     <dbl>      <dbl>      <dbl> <chr>       
## 1 EOPE  0      1        340   136     17000 0.00000637 0.00000637 ****
```

``` r
#MVM grade
stat.test.mvmgrade <- dunn_test(EOPE ~mvm_score_3cat, data=metadata)
stat.test.mvmgrade <- stat.test.mvmgrade %>%
  mutate(p.adj=round(p.adj, digits = 3))%>%
  mutate(y.position = c(0.82,0.88,0.97), mvm_score_3cat = group1,
         p.adj = ifelse(p.adj<0.001, "p < 0.001",paste("p = ",p.adj)))

fig_1_c <- metadata %>%
  filter(!is.na(mvm_score_3cat))%>%
  ggplot(aes(x=mvm_score_3cat, y=EOPE, color = mvm_score_3cat, fill= mvm_score_3cat))+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2), size=2.5)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_fill_manual(values = c(categories2[1],mvm_pal[2], mvm_pal[4]))+
  scale_color_manual(values = c(categories2[1],mvm_pal[2],mvm_pal[4]))+
  labs(x="MVM Grade", y="eoPRED score")+
  ylim(0,1)+
  scale_x_discrete(labels = c("0" = "No MVM","1" = "Low Grade","2" = "High Grade"))+
  stat_pvalue_manual(stat.test.mvmgrade, label="{p.adj}", tip.length=0.01)+
  stat_compare_means(aes(group=mvm_score_3cat), label.x = 0.9)
fig_1_c
```

![](03_SPAH_eoPRED_score_2026_files/figure-html/fig-1-part2-1.png)<!-- -->

``` r
fig_1_d <- metadata %>%
  filter(!is.na(mvm_score_3cat))%>%
  filter(!is.na(PE_GA_state)) %>%
  ggplot(aes(x=PE_GA_state, y=EOPE, color = mvm_score_3cat, fill= mvm_score_3cat))+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2), size=2.5)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, position =  position_dodge2(width = 0.75, preserve = "single"))+
  scale_fill_manual(values = c(categories2[1],mvm_pal[2], mvm_pal[4]),
                     labels=c("0" = "No MVM","1" = "Low Grade","2" = "High Grade"))+
  scale_color_manual(values = c(categories2[1],mvm_pal[2],mvm_pal[4]),
                     labels=c("0" = "No MVM","1" = "Low Grade","2" = "High Grade"))+
  labs(x="Clinical Diagnosis", y="eoPRED score", color="", fill="")+
  ylim(0,1)+
  #stat_compare_means(aes(label = paste("p = ",round(as.numeric(..p.format..),3))))+
  annotate(geom="text", x=2, y=1, col="black", size = 4, label="p=0.02")+
  annotate(geom="text", x=3, y=1, col="black", size = 4, label="p=0.02")+
  annotate(geom="text", x=4, y=1, col="black", size = 4, label="p=0.002")+
  theme(legend.position.inside = c(0.1,0.15),
        legend.position="inside", legend.background= element_rect(color = "grey"),
        legend.title = element_blank())+
scale_x_discrete(labels=c("EOPE"="EOPE","LOPE"="LOPE","nPTB"="PTB","nTB"="TB"))
fig_1_d
```

![](03_SPAH_eoPRED_score_2026_files/figure-html/fig-1-part2-2.png)<!-- -->

``` r
fig_1 <- plot_grid(fig_1_a, NULL, fig_1_b, 
                   NULL, NULL, NULL,
                   fig_1_c, NULL, fig_1_d,
                   labels = c("A.","","B.","","","","C.","","D."),
                   nrow = 3, ncol=3,
                   rel_widths = c(1,0.1,1), rel_heights = c(1,0.1,1))
fig_1
```

![](03_SPAH_eoPRED_score_2026_files/figure-html/fig-1-part2-3.png)<!-- -->

``` r
ggsave("Figure 1 - 2026.png", plot = fig_1, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/eoPRED Paper/"), width = 14, height = 8, device = "png")
```

## Table 1

Since gestational age and cell-type proportions are associated, here I tried to look at the association between eoPRED score and cell-type proportions, adjusting for gestational age.


``` r
## Individual Cell Types
#EOPE ~ GA + STB
fit_cells_stb <- lm(EOPE ~ ga_continuous + Syncytiotrophoblast, data = metadata)
fit_cells_stb_coef <- cbind(coef(summary(fit_cells_stb)))
signif(fit_cells_stb_coef, digits = 3)
```

```
##                     Estimate Std. Error t value Pr(>|t|)
## (Intercept)           0.9620     0.0966    9.96 2.12e-21
## ga_continuous        -0.0208     0.0024   -8.64 8.25e-17
## Syncytiotrophoblast   0.2520     0.0511    4.93 1.10e-06
```

``` r
#EOPE ~ GA + CTB
fit_cells_ctb <- lm(EOPE ~ ga_continuous + Trophoblasts, data = metadata)
fit_cells_ctb_coef <- cbind(coef(summary(fit_cells_ctb)))
signif(fit_cells_ctb_coef, digits = 3)
```

```
##               Estimate Std. Error t value Pr(>|t|)
## (Intercept)     1.1900    0.09980   12.00 4.01e-29
## ga_continuous  -0.0211    0.00254   -8.29 1.08e-15
## Trophoblasts   -0.2320    0.09140   -2.54 1.13e-02
```

``` r
#EOPE ~ GA + Endothelial
fit_cells_end <- lm(EOPE ~ ga_continuous + Endothelial, data = metadata)
fit_cells_end_coef <- cbind(coef(summary(fit_cells_end)))
signif(fit_cells_end_coef, digits = 3)
```

```
##               Estimate Std. Error t value Pr(>|t|)
## (Intercept)     1.1800    0.09190   12.80 1.05e-32
## ga_continuous  -0.0192    0.00236   -8.15 3.11e-15
## Endothelial    -1.1500    0.19400   -5.91 6.31e-09
```

``` r
#EOPE ~ GA + Stromal
fit_cells_str <- lm(EOPE ~ ga_continuous + Stromal, data = metadata)
fit_cells_str_coef <- cbind(coef(summary(fit_cells_str)))
signif(fit_cells_str_coef, digits = 3)
```

```
##               Estimate Std. Error t value Pr(>|t|)
## (Intercept)     1.1000    0.08980   12.20 4.41e-30
## ga_continuous  -0.0178    0.00233   -7.65 1.09e-13
## Stromal        -1.1400    0.16200   -7.06 5.57e-12
```

``` r
#EOPE ~ GA + Hofbauer
fit_cells_hbc <- lm(EOPE ~ ga_continuous + Hofbauer, data = metadata)
fit_cells_hbc_coef <- cbind(coef(summary(fit_cells_hbc)))
signif(fit_cells_hbc_coef, digits = 3)
```

```
##               Estimate Std. Error t value Pr(>|t|)
## (Intercept)     1.0600    0.09300   11.40 5.60e-27
## ga_continuous  -0.0174    0.00242   -7.20 2.27e-12
## Hofbauer       -2.3300    0.52600   -4.44 1.13e-05
```

``` r
#EOPE ~ GA + nRBC
fit_cells_rbc <- lm(EOPE ~ ga_continuous + nRBC, data = metadata)
fit_cells_rbc_coef <- cbind(coef(summary(fit_cells_rbc)))
signif(fit_cells_rbc_coef, digits = 3)
```

```
##               Estimate Std. Error t value Pr(>|t|)
## (Intercept)     1.0800    0.09210   11.70 4.49e-28
## ga_continuous  -0.0192    0.00238   -8.08 4.99e-15
## nRBC            2.1000    0.41400    5.07 5.58e-07
```

## Supplementary Figure 1

MVM Lesions:

-   `FN_AA` = Fibrinoid necrosis/Acute atherosis involving basal or parietal decidual vessels
-   `Musc_BP_arterioles` = Muscularized basal plate arterioles
-   `MHMA` = Mural hypertrophy of membrane arterioles
-   `Basal_Dec_vasc_thromb` = Basal decidual vascular thrombosis
-   `Infarct` = Any infarction of the villous parenchyma characterized by collapse of the intervillous space and coagulative necrosis. 
-   `Infarct_multiple` = >1 Infarct
-   `Increased_SK` = Increased syncytial knots for gestational age
-   `Villous_AG` = Villous aglutination
-   `Increased_PVF` = Increased perivillous fibrin deposition for gestational age
-   `DVH_STV` = Distal villous hypoplasia/small terminal villi for gestational age

-   `BP_acute_intradecid_hemmorhage` = Basal plate with acute intradecidual hemorrhage
-   `RP_blood_hematoma` = Basal plate with Retroplacental blood/hematoma
-   `RP_blood_hematoma_hemosid` = Basal plate with Retroplacental blood/hematoma and hemosiderin
-   `RP_blood_hematoma_infarct` = Basal plate with Retroplacental blood/hematoma and associated infarction


``` r
# Read in the extra histology data
histodata <- read.csv(here::here("B. Data/A. Cohort Data/histology_cleaned_11.15.23.csv"))
dim(histodata) #575 x 98
```

```
## [1] 575  98
```

``` r
histodata$ï..ID <- sub("^","SPAH_", histodata$ï..ID)
head(histodata)
```

```
##       ï..ID Placental_wt Placenta_SGA_AGA_LGA Wt_Pctile Bilobed Accessory_lobe
## 1 SPAH_2001          570                  AGA        10       0              0
## 2 SPAH_2002          420                  AGA         3       0              0
## 3 SPAH_2003          496                  AGA         5       0              0
## 4 SPAH_2004          432                  AGA         3       0              0
## 5 SPAH_2005          842                  LGA        11       0              0
## 6 SPAH_2007          570                  AGA         9       0              0
##   CORD_ABNORMALITY Decreased_Coiling Hyper_Coiling Velamentous_Insertion
## 1                0                 0             0                     0
## 2                0                 0             0                     0
## 3                0                 0             0                     0
## 4                0                 0             0                     0
## 5                0                 0             0                     0
## 6                1                 0             1                     0
##   Marginal_Cord_Insertion Single_Umbilical_artery Cord_knot Abnormal_length
## 1                       0                       0         0               0
## 2                       0                       0         0               0
## 3                       0                       0         0               0
## 4                       0                       0         0               0
## 5                       0                       0         0               0
## 6                       0                       0         0               0
##   Furcate_Insertion ACUTE_INFLAMMATION Maternal_inflamm Maternal_stage
## 1                 0                  0                0              0
## 2                 0                  1                1              2
## 3                 0                  0                0              0
## 4                 0                  0                0              0
## 5                 0                  1                1              1
## 6                 0                  1                1              1
##   Maternal_stage_3cat Maternal_high_stage Fetal_inflamm Fetal_stage
## 1                   0                   0             0           0
## 2                   2                   1             0           0
## 3                   0                   0             0           0
## 4                   0                   0             0           0
## 5                   1                   0             1           1
## 6                   1                   0             1           2
##   Fetal_stage_3cat Fetal_high_stage Funisitis Granular_fetal Periph_funisitis
## 1                0                0         0              0                0
## 2                0                0         0              0                0
## 3                0                0         0              0                0
## 4                0                0         0              0                0
## 5                1                0         0              1                0
## 6                2                1         1              4                0
##   acute_inflammation_3cat Villous_edema CHRONIC_INFLAMMATION
## 1                       0             0                    1
## 2                       2             0                    1
## 3                       0             0                    1
## 4                       0             0                    1
## 5                       1             0                    0
## 6                       2             0                    1
##   Chronic_inflammation_sum Chronic_villitis Chronic_villitis_di
## 1                        1                0                   0
## 2                        1                0                   0
## 3                        2                0                   0
## 4                        2                0                   0
## 5                        0                0                   0
## 6                        1                0                   0
##   Chronic_basal_villitis CDPC Chronic_chorioamnio Chronic_chorionitis
## 1                      1    0                   0                   0
## 2                      0    1                   0                   0
## 3                      0    1                   0                   0
## 4                      1    1                   0                   0
## 5                      0    0                   0                   0
## 6                      0    1                   0                   0
##   Chronic_marg_decid Chronic_decid_perivasc Chronic_intervillositis
## 1                  0                      0                       0
## 2                  0                      0                       0
## 3                  1                      0                       0
## 4                  0                      0                       0
## 5                  0                      0                       0
## 6                  0                      0                       0
##   Eosin_Tcell_vasculitis chronic_inflammation_comp chronic_inflammation_3cat
## 1                      0                         1                         1
## 2                      0                         1                         1
## 3                      0                         2                         2
## 4                      0                         2                         2
## 5                      0                         0                         0
## 6                      0                         1                         1
##   FETAL_VASC_PATH Fetal_vasc_thrombi FVT_Multifocal Thrombi_chorionic_vessel
## 1               1                  0              0                        0
## 2               0                  0              0                        0
## 3               0                  0              0                        0
## 4               1                  1              1                        0
## 5               1                  0              0                        0
## 6               1                  0              0                        0
##   TCV_Multifocal Thrombi_stem_villous_vessel TSV_Multifocal
## 1              0                           0              0
## 2              0                           0              0
## 3              0                           0              0
## 4              0                           1              1
## 5              0                           0              0
## 6              0                           0              0
##   Thrombi_umbilical_vessel TUV_Multifocal Avascular_villi AV_Multifocal
## 1                        0              0               1             1
## 2                        0              0               0             0
## 3                        0              0               0             0
## 4                        0              0               1             1
## 5                        0              0               1             0
## 6                        0              0               1             0
##   Fetal_thrombotic_vasculopathy Myonecrosis Fetal_vasc_path_3cat Vessel_path
## 1                             1           0                    2           1
## 2                             0           0                    0           0
## 3                             0           0                    0           0
## 4                             0           0                    1           0
## 5                             0           0                    1           0
## 6                             0           0                    1           0
##   FN_AA Musc_BP_arterioles MHMA Basal_Dec_vasc_thromb Villous_changes Infarct
## 1     0                  1    0                     0               0       0
## 2     0                  0    0                     0               0       0
## 3     0                  0    0                     0               0       0
## 4     0                  0    0                     0               0       0
## 5     0                  0    0                     0               0       0
## 6     0                  0    0                     0               1       1
##   Infarct_multiple Increased_SK Villous_AG Increased_PVF DVH_STV MVM_lesion_sum
## 1                0            0          0             0       0              1
## 2                0            0          0             0       0              0
## 3                0            0          0             0       0              0
## 4                0            0          0             0       0              0
## 5                0            0          0             0       0              0
## 6                0            0          0             0       0              1
##   MVM_score mvm_score_3cat mvm_di EV_OF_ABRUPTION
## 1         1              0      0               0
## 2         0              0      0               0
## 3         0              0      0               0
## 4         0              0      0               0
## 5         0              0      0               0
## 6         1              0      0               0
##   BP_acute_intradecid_hemmorhage RP_blood_hematoma RP_blood_hematoma_hemosid
## 1                              0                 0                         0
## 2                              0                 0                         0
## 3                              0                 0                         0
## 4                              0                 0                         0
## 5                              0                 0                         0
## 6                              0                 0                         0
##   RP_blood_hematoma_infarct Focal_acute_IV_hem OLD_MBP_BLEED Memhem
## 1                         0                  0             0      0
## 2                         0                  0             0      0
## 3                         0                  0             0      0
## 4                         0                  0             0      0
## 5                         0                  0             0      0
## 6                         0                  0             0      0
##   BP_hemosid_depo Chorioam_hemosiderosis SC_IV_hem_thrombi Amnion_Nodosum
## 1               0                      0                 1              0
## 2               0                      0                 1              0
## 3               0                      0                 0              0
## 4               0                      0                 0              0
## 5               0                      0                 0              0
## 6               0                      0                 0              0
##   Meconium_histiocytosis Villous_dysmaturity Villous_chorangiosis
## 1                      0                   0                    0
## 2                      0                   0                    0
## 3                      0                   0                    0
## 4                      0                   0                    0
## 5                      0                   0                    1
## 6                      1                   0                    0
##   Massive_PVF_deposition BPmyo other_chronic_path Phenotype Phenotype_Group
## 1                      0     1                  0        cF              3b
## 2                      0     1                  0        Ac              6b
## 3                      0     0                  0         C              4a
## 4                      0     0                  0        Cf              4a
## 5                      0     1                  1        af              6a
## 6                      0     0                  0       Acf              7c
##                                                                                                            phenotype_label
## 1     High-grade fetal vascular pathology with low-grade chronic inflammation and/or low-grade maternal vascular pathology
## 2                             High-grade acute inflammation and low-grade fetal vascular pathology or chronic inflammation
## 3                                     High-grade chronic inflammation (with or without low-grade fetal vascular pathology)
## 4                                     High-grade chronic inflammation (with or without low-grade fetal vascular pathology)
## 5                Low-grade fetal vascular pathology or chronic inflammation (with or without low-grade acute inflammation)
## 6 High-grade acute inflammation and low-grade maternal vascular pathology or any two or more low-grade chronic pathologies
##          Comments Covid_supplement SPAH_supplement Placenta_SGAvAGAvLGA
## 1              cF                1               1                    2
## 2                                1               0                    2
## 3                                1               1                    2
## 4                                1               0                    2
## 5              af                1               1                    3
## 6 infarct, 0.6 cm                1               0                    2
```

``` r
histodata <- metadata %>%
  dplyr::select(Sample_ID, EOPE)%>%
  left_join(histodata, by=c("Sample_ID"="ï..ID"))

cov_path <- histodata %>%
  dplyr::select(
    EOPE,
    mvm_di,
    mvm_score_3cat,
    MVM_score,
    MVM_lesion_sum,
    `FN_AA`,
    `Musc_BP_arterioles`,
    `MHMA`,
    `Basal_Dec_vasc_thromb`,
    `Infarct`,
    `Infarct_multiple`,
    `Increased_SK`,
    `Villous_AG`,
    `Increased_PVF`,
    `DVH_STV`,
    ACUTE_INFLAMMATION,
    CHRONIC_INFLAMMATION,
    FETAL_VASC_PATH) %>%
  dplyr::rename(
    eoPRED = EOPE,
    `MVM (Y/N)` = mvm_di,
    `MVM grade` = mvm_score_3cat,
    `MVM score` = MVM_score,
    `MVM lesions (n)` = MVM_lesion_sum,
    `Fibrinoid necrosis/Acute atherosis` = `FN_AA`,
    `Muscularized basal plate arterioles` = `Musc_BP_arterioles`,
    `Mural hypertrophy of membrane arterioles` = `MHMA`,
    `Basal decidual vascular thrombosis` = `Basal_Dec_vasc_thromb`,
    `Infarct >1` = `Infarct_multiple`,
    `Increased syncytial knots` = `Increased_SK`,
    `Villous aglutination` = `Villous_AG`,
    `Increased perivillous fibrin deposition` = `Increased_PVF`,
    `Distal villous hypoplasia/small terminal villi` = `DVH_STV`,
    AI = ACUTE_INFLAMMATION,
    CI = CHRONIC_INFLAMMATION,
    FVM = FETAL_VASC_PATH
  )%>%
  mutate_if(is.character, as.factor)%>%
  mutate_if(is.factor, as.numeric)


res <- rcorr(as.matrix(cov_path))

# Initialize file path
png(height=750, width=800, file=here::here("D. Pathology Project/02_Outputs/AA. Cleaned/eoPRED Paper/Figure S1 - 2026.png"))

corrplot(res$r, p.mat=res$P, method='color', type='lower', diag =F,
         insig='blank', #p-val<0.5 not shown
         tl.col='black', tl.srt = 55, #axis text colour and angle
         addgrid.col = 'grey', col = rev(COL2('RdBu', 10)), #grid colours
         addCoef.col = 'black', number.cex = 0.6) #add R values

dev.off()
```

```
## png 
##   2
```
