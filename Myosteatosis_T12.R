library(mice)
library(rms)
library(foreign)
library(sandwich)
library(lmtest)
library(mitools)


#name of the dataset = thoracic_nsclc_included_stageI_II
dat_s <- split(thoracic_nsclc_included_stageI_II, thoracic_nsclc_included_stageI_II$hospital)

# Add Nelson-Aalen by center
dat_s <- lapply(dat_s, function(ZZ) {
  ZZ$nelsonaalen2 <- nelsonaalen(ZZ, time_delta_surgery_years, final_vital_status_1_dead) 
  ZZ
})

thoracic_nsclc_included_stageI_II_2 <- do.call(rbind, dat_s)
rownames(thoracic_nsclc_included_stageI_II_2) <- seq(nrow(thoracic_nsclc_included_stageI_II_2))


dat_impute <- thoracic_nsclc_included_stageI_II_2 %>% dplyr::select(ageatsurgery, gender, height_cm, weight_kg, race_3_cat, cigarette_smoking, asa_3_cat, zubrod_score_3_cat, fev1_predicted, dlco_predicted, open_vs_closed, neoadjuvant_therapy_general_cat, stage_tnm,
                                                                    diabetes, laa_950, steroids, copd, peripheral_vascular_disease, coronary_artery_disease, hypertension, cerebrvascular_history_cat,
                                                                    area_t5m, area_t5f,mean_t5m_hu, mean_t5f_hu, area_t8m, area_t8f,mean_t8m_hu, mean_t8f_hu, area_t10m, area_t10f, mean_t10m_hu, mean_t10f_hu, area_t12m, area_t12sat, area_t12imat, mean_t12m_hu, mean_t12sat_hu, mean_t12imat_hu,
                                                                    congestive_heart_failure, 
                                                                    minor_events, major_events, pulmonary_postoperative_events, surgery_until_discharge,final_vital_status_1_dead, time_delta_surgery_years, nelsonaalen2, cancer_histology, contrast_iv_thoracic, hospital)


sum(complete.cases(dat_impute))


dat_impute <- dat_impute %>% mutate(
  
  stage_tnm_cat = stage_tnm,
  race_2_cat = race_3_cat
)

dat_impute <- dat_impute %>% mutate(

  stage_tnm_cat = recode(stage_tnm_cat, '0' = '<=I', 'IA' = '<=I', 'IB' = '<=I', 'IIA' ='II', 'IIB' = 'II', 'IIIA' = '>=III', 'IIIB' = '>=III', 'IV' = '>=III'),
  race_2_cat = recode(race_2_cat, 'Black' = 'Other')
)

dat_impute$cigarette_smoking <- factor(dat_impute$cigarette_smoking, levels = c('Never smoked', 'Current smoker', 'Past smoker'))
dat_impute$cerebrvascular_history_cat <- factor(dat_impute$cerebrvascular_history_cat, levels = c('No CVD history', 'Yes'))
dat_impute$neoadjuvant_therapy_general_cat <- factor(dat_impute$neoadjuvant_therapy_general_cat, levels = c('None', 'Neoadjuvant Therapy'))
dat_impute <- dat_impute %>% mutate(minor_major_cat = case_when((major_events == 0 & minor_events<= 1) ~ 0,
                                                                (major_events >= 1 | minor_events >= 2) ~ 1, 
                                                                
))
dat_impute <- dat_impute %>% mutate(any_postoperative_events_3cat = case_when(minor_events +major_events == 0 ~ '0',
                                                                              minor_events +major_events ==1 ~ '1',
                                                                              minor_events +major_events >= 2 ~ '>=2'))                                    
dat_impute <- dat_impute %>% mutate(pulmonary_postoperative_events_cat = case_when(pulmonary_postoperative_events == 0 ~ '0',
                                                                                   pulmonary_postoperative_events ==1 ~ '1',
                                                                                   pulmonary_postoperative_events >= 2 ~ '>=2'))  


meth <- make.method(dat_impute)

n_vars  <- ncol(dat_impute)

predMat <- matrix(1, nrow=n_vars, ncol=n_vars) - diag(n_vars)
predMat[46,] <- 0 # time delta surgery
predMat[,46] <- 0

predMat[c(51:55),] <- 0
predMat[,c(51:55)] <- 0

dat_impute<-droplevels(dat_impute)

#####for hospital-specific imputation
dat_impute_mgh <- dat_impute[dat_impute$hospital == 'MGH',]
dat_impute_uc <- dat_impute[dat_impute$hospital == 'UChicago',]
dat_impute_uc$cancer_histology<-droplevels(dat_impute_uc$cancer_histology)
dat_impute_wash <- dat_impute[dat_impute$hospital == 'WashU',]

idat_mgh_stageI_II <- mice(data=dat_impute_mgh, meth=meth, m=100, seed=3535, predictorMatrix=predMat, print=FALSE)
idat_uc_stageI_II<- mice(data=dat_impute_uc, meth=meth, m=100, seed=3535, predictorMatrix=predMat, print=FALSE)
idat_wash_stageI_II <- mice(data=dat_impute_wash, meth=meth, m=100, seed=3535, predictorMatrix=predMat, print=FALSE)
idat_hos_stageI_II_corrected <- rbind(idat_mgh_stageI_II, idat_uc_stageI_II)
idat_hos_stageI_II_corrected <- rbind(idat_hos_stageI_II_corrected, idat_wash_stageI_II)

##hospital-specific imp
idat_imp <- complete(idat_hos_stageI_II_corrected, action='long')
idat_imp$any_postoperative_events_3cat <- as.factor(idat_imp$any_postoperative_events_3cat)
idat_imp$pulmonary_postoperative_events_cat <- as.factor(idat_imp$pulmonary_postoperative_events_cat)
idat_imp$any_postoperative_events_3cat <- factor(idat_imp$any_postoperative_events_3cat, levels = c('0', '1', '>=2')) 
idat_imp$pulmonary_postoperative_events_cat <- factor(idat_imp$pulmonary_postoperative_events_cat, levels = c('0', '1', '>=2')) 
idat_imp <- idat_imp %>% mutate(pulmonary_postoperative_events_logistic = case_when(pulmonary_postoperative_events == 0 ~ '0',
                                                                                    pulmonary_postoperative_events >=1 ~ '1'))
idat_imp$pulmonary_postoperative_events_logistic <- factor(idat_imp$pulmonary_postoperative_events_logistic, levels=c('0','1'))

idat_imp <- idat_imp %>% mutate(any_postoperative_events_logistic = case_when(minor_events +major_events == 0 ~ '0',
                                                                              minor_events +major_events >=1 ~ '1'))
idat_imp$any_postoperative_events_logistic <- factor(idat_imp$any_postoperative_events_logistic, levels=c('0','1'))

idat_imp$gender <- factor(idat_imp$gender, levels = c('Female', 'Male')) 
idat_imp_f <- idat_imp %>% filter(gender == 'Female')
idat_imp_m <- idat_imp %>% filter(gender == 'Male')

idat_imp$index_t12m <- idat_imp$area_t12m*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$index_t12s <- (idat_imp$area_t12sat/(idat_imp$height_cm*idat_imp$height_cm))*10000
idat_imp$index_t12i <- (idat_imp$area_t12imat/(idat_imp$height_cm*idat_imp$height_cm))*10000

idat_imp$index_t5m <- idat_imp$area_t5m*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$index_t8m <- idat_imp$area_t8m*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$index_t10m <- idat_imp$area_t10m*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$index_t12m <- idat_imp$area_t12m*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$sum_thoracic_m <- idat_imp$area_t5m + idat_imp$area_t8m + idat_imp$area_t10m + idat_imp$area_t12m
idat_imp$sum_thoracic_m_mean_hu <- (idat_imp$mean_t5m_hu + idat_imp$mean_t8m_hu + idat_imp$mean_t10m_hu + idat_imp$mean_t12m_hu)/4
idat_imp$index_sum_thoracic_m <- idat_imp$sum_thoracic_m*10000/(idat_imp$height_cm*idat_imp$height_cm)

idat_imp$adipose_sum_thoracic <- idat_imp$area_t12sat + idat_imp$area_t12imat + idat_imp$area_t10f + idat_imp$area_t8f + idat_imp$area_t5f
idat_imp$index_adipose_sum_thoracic <- ((idat_imp$area_t12sat + idat_imp$area_t12imat + idat_imp$area_t10f + idat_imp$area_t8f + idat_imp$area_t5f)/(idat_imp$height_cm*idat_imp$height_cm))*10000
idat_imp$area_t12f <- idat_imp$area_t12sat + idat_imp$area_t12imat
idat_imp$index_t12f <- ((idat_imp$area_t12sat + idat_imp$area_t12imat )/(idat_imp$height_cm*idat_imp$height_cm))*10000
idat_imp$mean_t12f_hu <- (idat_imp$mean_t12sat_hu + idat_imp$mean_t12imat_hu)/2
idat_imp$index_t5f <- idat_imp$area_t5f*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$index_t8f <- idat_imp$area_t8f*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$index_t10f <- idat_imp$area_t10f*10000/(idat_imp$height_cm*idat_imp$height_cm)
idat_imp$sum_thoracic_f_mean_hu <- (idat_imp$mean_t5f_hu + idat_imp$mean_t8f_hu + idat_imp$mean_t10f_hu + idat_imp$mean_t12f_hu)/4
idat_imp$index_t12s <- (idat_imp$area_t12sat/(idat_imp$height_cm*idat_imp$height_cm))*10000
idat_imp$index_t12i <- (idat_imp$area_t12imat/(idat_imp$height_cm*idat_imp$height_cm))*10000


idat_imp$gauge_index_t12m <- idat_imp$mean_t12m_hu*idat_imp$index_t12m
idat_imp$gauge_index_t12s <- idat_imp$mean_t12sat_hu*idat_imp$index_t12s
idat_imp$gauge_index_t12i <-idat_imp$mean_t12imat_hu*idat_imp$index_t12i
idat_imp$sum_thoracic_m_mean_hu_weighted <- (idat_imp$mean_t5m_hu*idat_imp$area_t5m + idat_imp$mean_t8m_hu*idat_imp$area_t8m + idat_imp$mean_t10m_hu*idat_imp$area_t10m + idat_imp$mean_t12m_hu*idat_imp$area_t12m)/idat_imp$sum_thoracic_m
idat_imp$sum_thoracic_f_mean_hu_weighted <- (idat_imp$mean_t5f_hu*idat_imp$area_t5f + idat_imp$mean_t8f_hu*idat_imp$area_t8f + idat_imp$mean_t10f_hu*idat_imp$area_t10f + idat_imp$mean_t12f_hu*idat_imp$area_t12f)/idat_imp$sum_thoracic_f

idat_imp$gauge_index_t12m <- idat_imp$mean_t12m_hu*idat_imp$index_t12m
idat_imp$gauge_index_t12s <- idat_imp$mean_t12sat_hu*idat_imp$index_t12s
idat_imp$gauge_index_t12i <-idat_imp$mean_t12imat_hu*idat_imp$index_t12i
idat_imp$gauge_index_t12f <-idat_imp$mean_t12f_hu*idat_imp$index_t12f
idat_imp$gauge_index_t10m <- idat_imp$mean_t10m_hu*idat_imp$index_t10m
idat_imp$gauge_index_t10f <- idat_imp$mean_t10f_hu*idat_imp$index_t10f
idat_imp$gauge_index_t8m <- idat_imp$mean_t8m_hu*idat_imp$index_t8m
idat_imp$gauge_index_t8f <- idat_imp$mean_t8f_hu*idat_imp$index_t8f
idat_imp$gauge_index_t5m <- idat_imp$mean_t5m_hu*idat_imp$index_t5m
idat_imp$gauge_index_t5f <- idat_imp$mean_t5f_hu*idat_imp$index_t5f


idat_imp$index_t12m_5 <- idat_imp$index_t12m/5
idat_imp$mean_t12m_hu_5 <- idat_imp$mean_t12m_hu/5
idat_imp$index_t12i_5 <- idat_imp$index_t12i/5
idat_imp$mean_t12imat_hu_5 <- idat_imp$mean_t12imat_hu/5
idat_imp$index_t12s_5 <- idat_imp$index_t12s/5
idat_imp$mean_t12sat_hu_5 <- idat_imp$mean_t12sat_hu/5

idat_imp$gauge_index_t12m_50 <- idat_imp$gauge_index_t12m/50
idat_imp$gauge_index_t12s_50 <- idat_imp$gauge_index_t12s/50
idat_imp$gauge_index_t12i_50 <-idat_imp$gauge_index_t12i/50

idat_imp$gauge_index_t12m_200 <- idat_imp$gauge_index_t12m/200
idat_imp$gauge_index_t12s_200 <- idat_imp$gauge_index_t12s/200
idat_imp$gauge_index_t12i_200 <-idat_imp$gauge_index_t12i/200

idat_imp <- idat_imp %>% mutate(
  race_2_cat = race_3_cat,
  race_2_cat = recode(race_2_cat, 'Black' = 'Other'),
)


idat_imp_s <- split(idat_imp_f,idat_imp_f$.imp)

idat_imp_s <- split(idat_imp,idat_imp$.imp)
mod1 <- vector('list',length(idat_imp_s))
mod2 <- vector('list',length(idat_imp_s))


# Models 

# Survival
lhs <- 'Surv(time_delta_surgery_years, final_vital_status_1_dead)'


var_list <- list(var1 =c("index_t5m", "mean_t5m_hu",  "index_t8m", "mean_t8m_hu", "index_t10m", "mean_t10m_hu",
                         "index_t12m", "mean_t12m_hu","index_sum_thoracic_m", "sum_thoracic_m_mean_hu_weighted"),
                 var2 = c("index_t5f", "mean_t5f_hu" , "index_t8f", "mean_t8f_hu", "index_t10f", "mean_t10f_hu",
                          "index_t12f", "mean_t12f_hu", "index_adipose_sum_thoracic", "sum_thoracic_f_mean_hu_weighted")
)

var_list <- list(var1 =c("gauge_index_t5m", "gauge_index_t8m",  "gauge_index_t10m", "gauge_index_t12m"),
                 var2 = c("gauge_index_t5f", "gauge_index_t8f",  "gauge_index_t10f", "gauge_index_t12f")
)

var_list <- list(var1 =c("index_t5m", "index_t8m", "index_t10m", "index_t12m", "index_sum_thoracic_m"),
                 var2 = c( "mean_t5m_hu", "mean_t8m_hu", "mean_t10m_hu", "mean_t12m_hu", "sum_thoracic_m_mean_hu"),
                 var3 = c("index_t5f", "index_t8f", "index_t10f", "index_t12f", "index_adipose_sum_thoracic"), 
                 var4 = c("mean_t5f_hu", "mean_t8f_hu", "mean_t10f_hu", "mean_t12f_hu", "sum_thoracic_f_mean_hu")
)
var_list <- list(var1 =c("gauge_index_t12m_200"),
                 var2 = c( "gauge_index_t12i_200"),
                 var3 = c("gauge_index_t12s_200")
)


var_list <- list(var1 =c("index_t12m_5"),
                 var2 = c("mean_t12m_hu_5"),
                 var3 = c("index_t12i_5"), 
                 var4 = c("mean_t12imat_hu_5"),
                 var5 = c("index_t12s_5"),
                 var6 = c("mean_t12sat_hu_5")
)


var_list <- list(var1 = c("index_sum_thoracic_m"),
                 var2 = c( "sum_thoracic_m_mean_hu_weighted"),
                 var3 = c("index_adipose_sum_thoracic"), 
                 var4 = c("sum_thoracic_f_mean_hu_weighted")
)


var_list <- list("area_t5m", "mean_t5m_hu", "index_t5m", "area_t8m", "mean_t8m_hu", "index_t8m", "area_t10m", "mean_t10m_hu", "index_t10m",
                 "area_t12m", "mean_t12m_hu", "index_t12m",  
                 "sum_thoracic_m", "sum_thoracic_m_mean_hu","index_sum_thoracic_m"
)
var_list <- list("area_t5f", "mean_t5f_hu", "index_t5f", "area_t8f", "mean_t8f_hu", "index_t8f", "area_t10f", "mean_t10f_hu", "index_t10f",
                 "area_t12f", "mean_t12f_hu", "index_t12f",  
                 "adipose_sum_thoracic", "sum_thoracic_f_mean_hu","index_adipose_sum_thoracic")



var_list <- list("area_t12sat", "mean_t12sat_hu", "index_t12s", "sat_area_cm2", "l3_sat_mean", "index_l3s", 
                 "sat_sum", "mean_sat_sum_hu","index_s_sum")
var_list <- list("area_t12imat", "mean_t12imat_hu", "index_t12i", "imat_area_cm2", "l3_imat_mean", "index_l3i", 
                 "imat_sum", "mean_imat_sum_hu","index_i_sum")

idat_imp_f <- idat_imp %>% filter(gender == 'Female')
idat_imp_m <- idat_imp %>% filter(gender == 'Male')
idat
idat_imp_c <- idat_imp %>% filter(contrast_iv_thoracic == 1)
idat_imp_nc <- idat_imp %>% filter(contrast_iv_thoracic == 0)

idat_imp_stageI <- idat_imp %>% filter(stage_tnm_cat == '<=I')
idat_imp_stageII <- idat_imp %>% filter(stage_tnm_cat == 'II')

idat_imp_s <- split(idat_imp_f,idat_imp_f$.imp)
idat_imp_s <- split(idat_imp_m,idat_imp_m$.imp)

idat_imp_s <- split(idat_imp_c,idat_imp_c$.imp)
idat_imp_s <- split(idat_imp_nc,idat_imp_nc$.imp)

idat_imp_s <- split(idat_imp_stageI,idat_imp_stageI$.imp)
idat_imp_s <- split(idat_imp_stageII,idat_imp_stageII$.imp)

idat_imp_s <- split(idat_imp,idat_imp$.imp)
mod1 <- vector('list',length(idat_imp_s))
mod2 <- vector('list',length(idat_imp_s))
mod1_out_bind <- data.frame()



for (i in 1:length(var_list$var1)){
  for(ix in seq_along(idat_imp_s)) {
    
    sfit <- coxph(as.formula(paste(lhs, paste(c(var_list$var1[i],var_list$var2[i], var_list$var5[i],var_list$var6[i],
                                                #"gender",
                                                "contrast_iv_thoracic", 
                                                "ageatsurgery", 
                                                #"race_3_cat","cigarette_smoking","dlco_predicted",
                                                # "fev1_predicted","diabetes","copd","hypertension","cerebrvascular_history_cat",
                                                "strata(hospital)","cluster(hospital)"),
                                              collapse=" + "), sep=" ~ ")), data=idat_imp_s[[ix]])
    
    ph_p <- cox.zph(sfit)$table
    ph_p <- ph_p['GLOBAL','p']
    
    sfit2 <- coxph(as.formula(paste(lhs, paste(c(paste(c(var_list$var1[i],   "contrast_iv_thoracic"), collapse=" * "), 
                                                 paste(c(var_list$var2[i],"contrast_iv_thoracic"), collapse=" * "), 
                                                 # paste(c(var_list$var3[i], "gender", "contrast_iv_thoracic"), collapse=" * "), 
                                                 paste(c(var_list$var5[i], "contrast_iv_thoracic"), collapse=" * "), 
                                                 paste(c(var_list$var6[i],  "contrast_iv_thoracic"), collapse=" * "), 
                                                 #"gender", 
                                                 "contrast_iv_thoracic",
                                                 #"ageatsurgery", 
                                                 #"race_3_cat","cigarette_smoking","dlco_predicted","fev1_predicted","diabetes","copd",
                                                 #"hypertension","cerebrvascular_history_cat", 
                                                 "strata(hospital)","cluster(hospital)"), collapse=" +"), sep=" ~ ")),
                   data=idat_imp_s[[ix]])
    
    
    
    mod1[[ix]] <- list('beta'=coef(sfit),'vcov'=vcov(sfit),'zph'=ph_p,'lrt_p'=lrtest(sfit, sfit2)[2,5])
    
  }
  mod1_zph <- sum(sapply(mod1, function(ZZ) ZZ$zph)<0.05)
  mod1_int <- sum(sapply(mod1, function(ZZ) ZZ$lrt_p)<0.05)
  mod1_est <- MIcombine(results=lapply(mod1, function(ZZ) ZZ$beta), variances=lapply(mod1, function(ZZ) ZZ$vcov))
  
  mod1_out     <- data.frame('est'=mod1_est$coefficients)
  mod1_out$CIl <- exp( mod1_est$coefficients - 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$CIu <- exp( mod1_est$coefficients + 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$HR  <- exp( mod1_est$coefficients )
  mod1_out$p   <- pnorm( abs( mod1_est$coefficients/ sqrt(diag(mod1_est$variance)) ), lower.tail=FALSE)*2
  mod1_out$p   <- format(mod1_out$p, scientific=FALSE)
  mod1_out$zph <- mod1_zph
  mod1_out$lrt_p <- mod1_int
  mod1_out$HR <- format(round(mod1_out$HR, 2), nsmall = 2)
  mod1_out$CIl <- format(round(mod1_out$CIl, 2), nsmall = 2)
  mod1_out$CIu <- format(round(mod1_out$CIu, 2), nsmall = 2)
  mod1_out$CI_combined <- paste(mod1_out$CIl,mod1_out$CIu,sep=",")
  mod1_out     <-  mod1_out[,c('HR','CI_combined','p', 'zph', 'lrt_p')]
  
  
  mod1_out_bind <- rbind(mod1_out_bind, mod1_out)
}



mod1_out_bind$p <- as.numeric(mod1_out_bind$p)
mod1_out_bind$p <- format(round(mod1_out_bind$p, 3), nsmall = 3)
mod1_out_bind_new <- mod1_out_bind[seq(17, nrow(mod1_out_bind), 30), ]
mod1_out_bind_table <- htmlTable::htmlTable(mod1_out_bind)
mod1_out_bind_table

######################test rcs

var2 = c("mean_t12m_hu_5"),
var3 = c("index_t12i_5"), 
var4 = c("mean_t12imat_hu_5"),
var5 = c("index_t12s_5"),
var6 = c("mean_t12sat_hu_5")
var_list <- list(var1 =c("index_t5m", "index_t8m", "index_t10m", "index_t12m", "index_sum_thoracic_m"),
                 var2 = c( "mean_t5m_hu", "mean_t8m_hu", "mean_t10m_hu", "mean_t12m_hu", "sum_thoracic_m_mean_hu_weighted"),
                 var3 = c("index_t5f", "index_t8f", "index_t10f", "index_t12f", "index_adipose_sum_thoracic"), 
                 var4 = c("mean_t5f_hu", "mean_t8f_hu", "mean_t10f_hu", "mean_t12f_hu", "sum_thoracic_f_mean_hu_weighted")
)

for(ix in seq_along(idat_imp_s)) {
  
  sfit <- coxph(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ index_t12m_5+mean_t12m_hu_5+index_t12s_5+mean_t12sat_hu_5+contrast_iv_thoracic+ageatsurgery+gender+race_3_cat+cigarette_smoking+dlco_predicted+fev1_predicted+diabetes+copd+hypertension+cerebrvascular_history_cat , data=idat_imp_s[[ix]])
  sfit_nonlinear <- coxph(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ rcs(index_t12m_5,5)+rcs(index_t12i_5,5)+rcs(mean_t12imat_hu_5,5)+rcs(index_t12s_5,5)+rcs(mean_t12sat_hu_5,5)+contrast_iv_thoracic+ageatsurgery+gender+race_3_cat+cigarette_smoking+dlco_predicted+fev1_predicted+diabetes+copd+hypertension+cerebrvascular_history_cat , data=idat_imp_s[[ix]])
  anova_lin <- anova(sfit, sfit_nonlinear)[2,4]
  mod1[[ix]] <- list('anova_lin'=anova_lin)

  
}
mod1_lin<- sum(sapply(mod1, function(ZZ) ZZ$anova_lin)<0.05)

idat_imp_s <- split(idat_imp,idat_imp$.imp)
mod1 <- vector('list',length(idat_imp_s))
mod2 <- vector('list',length(idat_imp_s))
mod1_out_bind <- data.frame()


# Binary - marginal model
bfit <- glm(as.formula(paste(lhs,rhs,sep='~')), data=tmp[tmp$gender=='Male',], family=binomial)
beta <- coef(bfit)
vcov <- vcovCL(bfit,cluster = ~hospital)
bfit <- lmtest:::coeftest(bfit, vcov = vcovCL, cluster = ~hospital)
mod_marg_m[[ix]] <- list('beta'=coef(bfit),'vcov'=vcov)

lhs <- 'pulmonary_postoperative_events_logistic'
lhs <- 'any_postoperative_events_logistic'
for (i in 1:length(var_list$var1)){
  for(ix in seq_along(idat_imp_s)) {
    bfit <- glm(as.formula(paste(lhs, paste(c(var_list$var1[i],var_list$var3[i],var_list$var4[i], var_list$var5[i],var_list$var6[i], 
                                              "contrast_iv_thoracic" ,
                                              "gender", 
                                              "ageatsurgery" 
                                              #,"race_3_cat","cigarette_smoking","dlco_predicted","fev1_predicted","diabetes","copd","hypertension","cerebrvascular_history_cat"
    ),collapse=" + "), sep=" ~ ")),
    data=idat_imp_s[[ix]], family=binomial)
    
    bfit2 <- glm(as.formula(paste(lhs, paste(c(paste(c(var_list$var1[i], "contrast_iv_thoracic", "gender"), collapse=" * "),
                                               #paste(c(var_list$var2[i], "contrast_iv_thoracic", "gender"), collapse=" * "),
                                               paste(c(var_list$var3[i], "contrast_iv_thoracic", "gender"), collapse=" * "),
                                               paste(c(var_list$var4[i], "contrast_iv_thoracic", "gender"), collapse=" * "),
                                               paste(c(var_list$var5[i], "contrast_iv_thoracic", "gender"), collapse=" * "), 
                                               paste(c(var_list$var6[i],"contrast_iv_thoracic", "gender"), collapse=" * "), 
                                               "contrast_iv_thoracic", "gender", "ageatsurgery", "race_3_cat","cigarette_smoking","dlco_predicted","fev1_predicted","diabetes","copd","hypertension","cerebrvascular_history_cat"), collapse=" +"), sep=" ~ ")), data=idat_imp_s[[ix]], family=binomial)
    
    anova_p <- anova(bfit, bfit2)
    anova_p <- as.data.frame(anova_p)
    anova_p <- anova_p$P[nrow(anova_p)]
    bfit_beta <- coef(bfit)
    
    bfit_vcov <- vcovCL(bfit,cluster = ~hospital)

    bfit <- lmtest:::coeftest(bfit, vcov = vcovCL, cluster = ~hospital)
    
    
    mod1[[ix]] <- list('beta'=bfit_beta,'vcov'=bfit_vcov)
     }
  
 
  mod1_est <- MIcombine(results=lapply(mod1, function(ZZ) ZZ$beta), variances=lapply(mod1, function(ZZ) ZZ$vcov))
  
  mod1_out     <- data.frame('est'=mod1_est$coefficients)
  mod1_out$CIl <- exp( mod1_est$coefficients - 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$CIu <- exp( mod1_est$coefficients + 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$OR  <- exp( mod1_est$coefficients )
  mod1_out$p   <- pnorm( abs( mod1_est$coefficients/ sqrt(diag(mod1_est$variance)) ), lower.tail=FALSE)*2
  mod1_out$p   <- format(mod1_out$p, scientific=FALSE)
  
  
  mod1_out$OR <- format(round(mod1_out$OR, 2), nsmall = 2)
  mod1_out$CIl <- format(round(mod1_out$CIl, 2), nsmall = 2)
  mod1_out$CIu <- format(round(mod1_out$CIu, 2), nsmall = 2)
  mod1_out$CI_combined <- paste(mod1_out$CIl,mod1_out$CIu,sep=",")
  mod1_out     <-  mod1_out[,c('OR','CI_combined','p')]
  
  
  mod1_out_bind <- rbind(mod1_out_bind, mod1_out)
}

mod1_out_bind$p <- as.numeric(mod1_out_bind$p)
mod1_out_bind$p <- format(round(mod1_out_bind$p, 3), nsmall = 3)
mod1_out_bind_new <- mod1_out_bind[seq(17, nrow(mod1_out_bind), 30), ]
mod1_out_bind_table <- htmlTable::htmlTable(mod1_out_bind)
mod1_out_bind_table


# Ordinal
idat_imp_s <- split(idat_imp,idat_imp$.imp)
mod1 <- vector('list',length(idat_imp_s))
mod2 <- vector('list',length(idat_imp_s))
mod1_out_bind <- data.frame()

for (i in 1:length(var_list$var1)){
  for(ix in seq_along(idat_imp_s)) {
    
    ofit <- MASS:::polr(as.formula(paste(lhs, paste(c(var_list$var1[i],var_list$var2[i], var_list$var5[i],var_list$var6[i], "contrast_iv_thoracic", "gender", "ageatsurgery", "race_3_cat","cigarette_smoking","dlco_predicted","fev1_predicted","diabetes","copd","hypertension","cerebrvascular_history_cat"), collapse=" + "), sep=" ~ ")),
                        data=idat_imp_s[[ix]],Hess=TRUE)
    
    
    ofit2 <- MASS:::polr(as.formula(paste(lhs, paste(c(paste(c(var_list$var1[i], "contrast_iv_thoracic", "gender"), collapse=" * "),
                                                       paste(c(var_list$var2[i], "contrast_iv_thoracic", "gender"), collapse=" * "),
                                                       #paste(c(var_list$var3[i], "contrast_iv_thoracic", "gender"), collapse=" * "),paste(c(var_list$var4[i], "contrast_iv_thoracic", "gender"), collapse=" * "),
                                                       paste(c(var_list$var5[i], "contrast_iv_thoracic", "gender"), collapse=" * "), paste(c(var_list$var6[i],"contrast_iv_thoracic", "gender"), collapse=" * "), 
                                                       "contrast_iv_thoracic", "gender", "ageatsurgery", "race_3_cat","cigarette_smoking","dlco_predicted","fev1_predicted","diabetes","copd","hypertension","cerebrvascular_history_cat"), collapse=" +"), sep=" ~ ")), data=idat_imp_s[[ix]],Hess=TRUE)
   
    anova_p <- anova(ofit, ofit2)
    anova_p <- as.data.frame(anova_p)
    anova_p <- anova_p$P[nrow(anova_p)]
    ofit_beta <- coef(ofit)
   
    ofit_vcov <- vcovCL(ofit,cluster = ~hospital)
    cnms <- colnames(ofit_vcov) 
    cnms_out <- grep('\\|', cnms) 
    ofit_vcov <- ofit_vcov[-cnms_out, -cnms_out]
    
    
    mod1[[ix]] <- list('beta'=ofit_beta,'vcov'=ofit_vcov,'anova_p'=anova_p)
    
  }
  
  mod1_int <- sum(sapply(mod1, function(ZZ) ZZ$anova_p)<0.05)
  mod1_est <- MIcombine(results=lapply(mod1, function(ZZ) ZZ$beta), variances=lapply(mod1, function(ZZ) ZZ$vcov))
  
  mod1_out     <- data.frame('est'=mod1_est$coefficients)
  mod1_out$CIl <- exp( mod1_est$coefficients - 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$CIu <- exp( mod1_est$coefficients + 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$OR  <- exp( mod1_est$coefficients )
  mod1_out$p   <- pnorm( abs( mod1_est$coefficients/ sqrt(diag(mod1_est$variance)) ), lower.tail=FALSE)*2
  mod1_out$p   <- format(mod1_out$p, scientific=FALSE)
  
  mod1_out$anova_p <- mod1_int
  mod1_out$OR <- format(round(mod1_out$OR, 2), nsmall = 2)
  mod1_out$CIl <- format(round(mod1_out$CIl, 2), nsmall = 2)
  mod1_out$CIu <- format(round(mod1_out$CIu, 2), nsmall = 2)
  mod1_out$CI_combined <- paste(mod1_out$CIl,mod1_out$CIu,sep=",")
  mod1_out     <-  mod1_out[,c('OR','CI_combined','p', 'anova_p')]
  
  mod1_out_bind <- rbind(mod1_out_bind, mod1_out)
}


mod1_out_bind$p <- as.numeric(mod1_out_bind$p)
mod1_out_bind$p <- format(round(mod1_out_bind$p, 3), nsmall = 3)
mod1_out_bind_new <- mod1_out_bind[seq(17, nrow(mod1_out_bind), 30), ]
mod1_out_bind_table <- htmlTable::htmlTable(mod1_out_bind)
mod1_out_bind_table


lhs <- 'pulmonary_postoperative_events_cat'

for (i in 1:length(var_list$var1)){
  for(ix in seq_along(idat_imp_s)) {
    
    ofit <- MASS:::polr(as.formula(paste(lhs, paste(c(var_list$var1[i],var_list$var2[i],"gender", "contrast_iv_thoracic", "ageatsurgery", "race_3_cat","cigarette_smoking","dlco_predicted","fev1_predicted","diabetes","copd","hypertension","cerebrvascular_history_cat"), collapse=" + "), sep=" ~ ")), data=idat_imp_s[[ix]],Hess=TRUE)
    
    
    beta <- coef(ofit)
    cl_vcov <- vcovCL(ofit,cluster = ~hospital)
    
    ofit <- lmtest:::coeftest(ofit, vcov = cl_vcov, cluster = ~hospital)
    
    
    mod1[[ix]] <- list('beta'=beta,'vcov'=cl_vcov)
   
  }
  
  
  mod1_est <- MIcombine(results=lapply(mod1, function(ZZ) ZZ$beta), variances=lapply(mod1, function(ZZ) ZZ$vcov))
  
  mod1_out     <- data.frame('est'=mod1_est$coefficients)
  mod1_out$CIl <- exp( mod1_est$coefficients - 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$CIu <- exp( mod1_est$coefficients + 1.96*sqrt(diag(mod1_est$variance)) )
  mod1_out$HR  <- exp( mod1_est$coefficients )
  mod1_out$p   <- pnorm( abs( mod1_est$coefficients/ sqrt(diag(mod1_est$variance)) ), lower.tail=FALSE)*2
  mod1_out$p   <- format(mod1_out$p, scientific=FALSE)
  
  mod1_out$lrtest_p <- mod1_int
  mod1_out$HR <- format(round(mod1_out$HR, 3), nsmall = 3)
  mod1_out$CIl <- format(round(mod1_out$CIl, 3), nsmall = 3)
  mod1_out$CIu <- format(round(mod1_out$CIu, 3), nsmall = 3)
  mod1_out$CI_combined <- paste(mod1_out$CIl,mod1_out$CIu,sep=",")
  mod1_out     <-  mod1_out[,c('HR','CI_combined','p', 'lrtest_p')]
  
  mod1_out_bind <- rbind(mod1_out_bind, mod1_out)
}


mod1_out_bind$p <- as.numeric(mod1_out_bind$p)
mod1_out_bind$p <- format(round(mod1_out_bind$p, 3), nsmall = 3)

mod1_out_bind_new <- mod1_out_bind[seq(17, nrow(mod1_out_bind), 30), ]
mod1_out_bind_table <- htmlTable::htmlTable(mod1_out_bind_new)



############################Correlation
install.packages("RVAideMemoire")
library("RVAideMemoire")
install.packages("PerformanceAnalytics")
install.packages("remotes")
remotes::install_remote("PerformanceAnalytics")
library("PerformanceAnalytics")
spearman.ci(thoracic_nsclc_included_stage_I_II$index_t12m, thoracic_nsclc_included_stage_I_II$mean_t12m_hu, nrep = 1000, conf.level = 0.95)

thoracic_nsclc_included_stage_I_II$index_t12s <- (thoracic_nsclc_included_stage_I_II$area_t12sat/(thoracic_nsclc_included_stage_I_II$height_cm*thoracic_nsclc_included_stage_I_II$height_cm))*10000
thoracic_nsclc_included_stage_I_II$index_t12i <- (thoracic_nsclc_included_stage_I_II$area_t12imat/(thoracic_nsclc_included_stage_I_II$height_cm*thoracic_nsclc_included_stage_I_II$height_cm))*10000

thoracic_nsclc_included_stage_I_II_female <- thoracic_nsclc_included_stage_I_II %>% filter(gender == "Female")
thoracic_nsclc_included_stage_I_II_male <- thoracic_nsclc_included_stage_I_II %>% filter(gender == "Male")
thoracic_nsclc_included_stage_I_II_nc <- thoracic_nsclc_included_stage_I_II %>% filter(contrast_iv_thoracic == 0)
thoracic_nsclc_included_stage_I_II_c <- thoracic_nsclc_included_stage_I_II %>% filter(contrast_iv_thoracic == 1)
fat_muscle_data <- thoracic_nsclc_included_stage_I_II %>% select(hospital, area_t5m, mean_t5m_hu, area_t8m, mean_t8m_hu, area_t10m, mean_t10m_hu, 
                                                                 area_t12m,mean_t12m_hu,area_t5f, mean_t5f_hu, area_t8f, mean_t8f_hu, area_t10f, mean_t10f_hu, 
                                                                 area_t12f,mean_t12f_hu)
T12fat_imat_data <- thoracic_nsclc_included_stage_I_II %>% select(index_t12m, mean_t12m_hu, index_t12i, mean_t12imat_hu, index_t12s, mean_t12sat_hu)

chart.Correlation(T12fat_imat_data, histogram=TRUE, method = "spearman", pch=19)

thoracic_nsclc_included_stage_I_II$gauge_index_t12s <- thoracic_nsclc_included_stage_I_II$mean_t12sat_hu*thoracic_nsclc_included_stage_I_II$index_t12s
thoracic_nsclc_included_stage_I_II$gauge_index_t12i <-thoracic_nsclc_included_stage_I_II$mean_t12imat_hu*thoracic_nsclc_included_stage_I_II$index_t12i

gauge_data <- thoracic_nsclc_included_stage_I_II %>% select(gauge_index_t12m,gauge_index_t12s, gauge_index_t12i )
thoracic_nsclc_included_stage_I_II$contrast_iv_thoracic
chart.Correlation(gauge_data, histogram=TRUE, method = "spearman", pch=19)



library(corrplot)

# do not edit
corrplot2 <- function(data,
                      method = "pearson",
                      sig.level = 0.05,
                      order = "original",
                      diag = FALSE,
                      type = "upper",
                      tl.srt = 90,
                      number.font = 1,
                      number.cex = 1,
                      tl.cex = 1,
                      mar = c(0, 0, 0, 0)) {
  library(corrplot)
  data_incomplete <- data
  data <- data[complete.cases(data), ]
  mat <- cor(data, method = method)
  cor.mtest <- function(mat, method) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = method)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(data, method = method)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(mat,
           method = "color", col = col(200), number.font = number.font,
           mar = mar, number.cex = number.cex,
           type = type, order = order,
           addCoef.col = "black", # add correlation coefficient
           tl.col = "black", tl.srt = tl.srt, # rotation of text labels
           # combine with significance level
           p.mat = p.mat, sig.level = sig.level, insig = "blank",
           # hide correlation coefficients on the diagonal
           diag = diag
  )
}
colnames(T12fat_imat_data) <- c('SMI', 'SMD', 'IMAT Index', 'IMAT Density', 'SAT Index', 'SAT Density')
T12fat_imat_data_f <- thoracic_nsclc_included_stage_I_II_female  %>% select(index_t12m, mean_t12m_hu, index_t12i, mean_t12imat_hu, index_t12s, mean_t12sat_hu)
T12fat_imat_data_m <- thoracic_nsclc_included_stage_I_II_male  %>% select(index_t12m, mean_t12m_hu, index_t12i, mean_t12imat_hu, index_t12s, mean_t12sat_hu)
T12fat_imat_data_nc <- thoracic_nsclc_included_stage_I_II_nc  %>% select(index_t12m, mean_t12m_hu, index_t12i, mean_t12imat_hu, index_t12s, mean_t12sat_hu)
T12fat_imat_data_c <- thoracic_nsclc_included_stage_I_II_c  %>% select(index_t12m, mean_t12m_hu, index_t12i, mean_t12imat_hu, index_t12s, mean_t12sat_hu)

colnames(T12fat_imat_data_f) <- c('SMI', 'SMD', 'IMAT Index', 'IMAT Density', 'SAT Index', 'SAT Density')
colnames(T12fat_imat_data_m) <- c('SMI', 'SMD', 'IMAT Index', 'IMAT Density', 'SAT Index', 'SAT Density')
corr_all <- corrplot2(
  data = T12fat_imat_data,
  method = "spearman",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75, 
  number.cex = 1.3
)
corr_female <- corrplot2(
  data = T12fat_imat_data_f,
  method = "spearman",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75
)
corr_male <- corrplot2(
  data = T12fat_imat_data_m,
  method = "spearman",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75
)


# One figure in row 1 and two figures in row 2
attach(mtcars)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
layout(matrix(c(1), 1, 1, byrow = TRUE))
corrplot2(
  data = T12fat_imat_data,
  method = "spearman",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75, 
  number.cex = 1.3
)
corrplot2(
  data = T12fat_imat_data_f,
  method = "spearman",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75, 
  number.cex = 1.3
)
corrplot2(
  data = T12fat_imat_data_m,
  method = "spearman",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75, 
  number.cex = 1.3
)

####wilcoxon test 
wilcox.test(thoracic_nsclc_included_stage_I_II$mean_t12sat_hu, thoracic_nsclc_included_stage_I_II$mean_t12imat_hu)
wilcox.test(thoracic_nsclc_included_stage_I_II_c$mean_t12sat_hu, thoracic_nsclc_included_stage_I_II_c$mean_t12imat_hu)
wilcox.test(thoracic_nsclc_included_stage_I_II_nc$mean_t12sat_hu, thoracic_nsclc_included_stage_I_II_nc$mean_t12imat_hu)
t.test(thoracic_nsclc_included_stage_I_II$mean_t12sat_hu, thoracic_nsclc_included_stage_I_II$mean_t12imat_hu)

###########################missing
library(naniar)
T12fat_imat_data_hos <- thoracic_nsclc_included %>% select(hospital, index_t12m, mean_t12m_hu, index_t12i, mean_t12imat_hu, index_t12s, mean_t12sat_hu)
colnames(T12fat_imat_data_hos) <- c('hospital', 'Muscle Index', 'Muscle Density', 'IMAT Index', 'IMAT Density', 'SAT Index', 'SAT Density')
vis_miss(fat_muscle_data)

gg_miss_var(T12fat_imat_data_hos, facet = hospital, show_pct = TRUE)




##################histograms sat imat density 
thoracic_nsclc_included_contrast <- thoracic_nsclc_included_stage_I_II %>% filter(contrast_iv_thoracic == 1)
thoracic_nsclc_included_noncontrast <- thoracic_nsclc_included_stage_I_II %>% filter(contrast_iv_thoracic == 0)


# 4 figures arranged in 2 rows and 2 columns
attach(mtcars)
par(mfrow=c(1,2))


hist(thoracic_nsclc_included_contrast$mean_t12imat_hu,ylim= c(0,70),xlim=c(-150,-20),col = "black", density = 30, angle = 135, breaks =25,  main = paste("Histogram of IMAT & SAT density on contrast-enhanced CT scans"), xlab = "Density [HU]", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(thoracic_nsclc_included_contrast$mean_t12sat_hu, add=T, col = "black", density = 10, angle = 45, breaks = 35 )
hist(thoracic_nsclc_included_noncontrast$mean_t12imat_hu,ylim= c(0,70),xlim=c(-150,-20),col = "black", density = 30, angle = 135, breaks =30,  main = paste("Histogram of IMAT & SAT density on non-contrast CT scans"), xlab = "Density [HU]", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(thoracic_nsclc_included_noncontrast$mean_t12sat_hu, add=T, col = "black", density = 10, angle = 45, breaks = 30 )
legend("topright", c("IMAT", "SAT"), col=c("black", "black"),density=c(30, 10), angle = c(135, 45), cex = 1.5)


hist(thoracic_nsclc_included_contrast$mean_t12imat_hu,ylim= c(0,70), xlim=c(-150,-20), col="red", breaks =sqrt(nrow(thoracic_nsclc_included_contrast)), main = paste("Histogram of IMAT & SAT density on contrast-enhanced CT scans"), xlab = "Density [HU]")
hist(thoracic_nsclc_included_contrast$mean_t12sat_hu, add=T, col=rgb(0, 1, 0, 0.5), breaks = 40 )
hist(thoracic_nsclc_included_noncontrast$mean_t12imat_hu,ylim= c(0,70),xlim=c(-150,-20), col="red", breaks =30,  main = paste("Histogram of IMAT & SAT density on non-contrast CT scans"), xlab = "Mean radioattenuation [HU]")
hist(thoracic_nsclc_included_noncontrast$mean_t12sat_hu, add=T, col=rgb(0, 1, 0, 0.5), breaks = 30 )
summary(thoracic_nsclc_included_contrast$mean_t12sat_hu)
summary(thoracic_nsclc_included_contrast$mean_t12imat_hu)

summary(thoracic_nsclc_included_noncontrast$mean_t12imat_hu)
summary(thoracic_nsclc_included_noncontrast_male$mean_t12sat_hu)
summary(thoracic_nsclc_included_contrast$mean_t12imat_hu)
summary(thoracic_nsclc_included_contrast_male$mean_t12sat_hu)

hist(thoracic_nsclc_included$mean_t12imat_hu,xlim=c(-150,-20), col="red", breaks =30,  main = paste("Histogram of IMAT & SAT density"), xlab = "Mean radioattenuation [HU]")
hist(thoracic_nsclc_included$mean_t12sat_hu, add=T, col=rgb(0, 1, 0, 0.5), breaks = 30 )
hist(thoracic_nsclc_included_stage_I_II$mean_t12imat_hu,ylim= c(0,250),xlim=c(-150,-20),col = "black", density = 20, angle = 135, breaks =30,  main = paste("Histogram of IMAT & SAT density"), xlab = "Mean radioattenuation [HU]")
hist(thoracic_nsclc_included_stage_I_II$mean_t12sat_hu, add=T, col = "black", density = 10, angle = 45, breaks = 30 )

###########################kaplan meier by derstine cut offs - 2018 - noncontrast 
install.packages('ggeasy')
library(ggeasy)

thoracic_nsclc_included_contrast <- thoracic_nsclc_included_stageI_II %>% filter(contrast_iv_thoracic == 1)
thoracic_nsclc_included_noncontrast <- thoracic_nsclc_included_stageI_II %>% filter(contrast_iv_thoracic == 0)
thoracic_nsclc_included_noncontrast_female <- thoracic_nsclc_included_noncontrast[thoracic_nsclc_included_noncontrast$gender == "Female",]
thoracic_nsclc_included_noncontrast_male <- thoracic_nsclc_included_noncontrast[thoracic_nsclc_included_noncontrast$gender == "Male",]
thoracic_nsclc_included_contrast_female <- thoracic_nsclc_included_contrast[thoracic_nsclc_included_contrast$gender == "Female",]
thoracic_nsclc_included_contrast_male <- thoracic_nsclc_included_contrast[thoracic_nsclc_included_contrast$gender == "Male",]
median(thoracic_nsclc_included_contrast$mean_t12imat_hu, na.rm = TRUE)

thoracic_nsclc_included_noncontrast_female$cat_smd <- cut(thoracic_nsclc_included_noncontrast_female$mean_t12m_hu,
                                                          breaks=c(-30, 31.3, 150),
                                                          labels=c('Low SMD', 'Normal SMD'))
thoracic_nsclc_included_noncontrast_male$cat_smd <- cut(thoracic_nsclc_included_noncontrast_male$mean_t12m_hu,
                                                        breaks=c(-30, 37.5, 150),
                                                        labels=c('Low SMD', 'Normal SMD'))
thoracic_nsclc_included_noncontrast_cat_smd <- rbind(thoracic_nsclc_included_noncontrast_female,thoracic_nsclc_included_noncontrast_male)

km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ cat_smd, data = thoracic_nsclc_included_noncontrast_cat_smd) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  )+ 
  geom_step(size = 2)+ 
  add_confidence_interval() +
  add_risktable( risktable_height = 0.33,
                 size = 7, # increase font size of risk table statistics
                 theme =   # increase font size of risk table title and y-axis label
                   list(
                     theme_risktable_default(axis.text.y.size = 22,
                                             
                                             plot.title.size = 22),
                     theme(plot.title = element_text(face = "bold"),
                           legend.text = element_text(size = 20)))) 
km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ cat_smd, data = thoracic_nsclc_included_noncontrast_cat_smd) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) +
  geom_step(size = 2)+ 
  add_confidence_interval()


km_derstine_smd <- km90_default + coord_cartesian(xlim = c(0, 8)) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent, 
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(breaks = 0:9, expand = c(0.02, 0)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 1)) +
  ggeasy::easy_y_axis_labels_size(size = 30) +
  ggeasy::easy_x_axis_labels_size(size = 30) +
  ggeasy::easy_y_axis_title_size(size = 30) +
  ggeasy::easy_x_axis_title_size(size = 30) +
  theme(legend.text=element_text(size=30))

km_derstine_smd
survdiff(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ cat_smd, data = thoracic_nsclc_included_noncontrast_cat_smd)

####################km for smi

thoracic_nsclc_included_noncontrast_female <- thoracic_nsclc_included_noncontrast[thoracic_nsclc_included_noncontrast$gender == "Female",]
thoracic_nsclc_included_noncontrast_male <- thoracic_nsclc_included_noncontrast[thoracic_nsclc_included_noncontrast$gender == "Male",]
thoracic_nsclc_included_noncontrast_female$cat_smi <- cut(thoracic_nsclc_included_noncontrast_female$index_t12m,
                                                          breaks=c(0, 20.8, 150),
                                                          labels=c('Low SMI', 'Normal SMI'))
thoracic_nsclc_included_noncontrast_male$cat_smi <- cut(thoracic_nsclc_included_noncontrast_male$index_t12m,
                                                        breaks=c(0, 28.8, 200),
                                                        labels=c('Low SMI', 'Normal SMI'))
thoracic_nsclc_included_noncontrast_cat_smi <- rbind(thoracic_nsclc_included_noncontrast_female,thoracic_nsclc_included_noncontrast_male)
summary(thoracic_nsclc_included_noncontrast_female$cat_smi)
km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ cat_smi, data = thoracic_nsclc_included_noncontrast_cat_smi) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  )+ 
  geom_step(size = 2)+ 
  add_confidence_interval() +
  add_risktable( risktable_height = 0.33,
                 size = 7, # increase font size of risk table statistics
                 theme =   # increase font size of risk table title and y-axis label
                   list(
                     theme_risktable_default(axis.text.y.size = 22,
                                             
                                             plot.title.size = 22),
                     theme(plot.title = element_text(face = "bold"),
                           legend.text = element_text(size = 20))))
km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ cat_smi, data = thoracic_nsclc_included_noncontrast_cat_smi) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) +
  geom_step(size = 2)+ 
  add_confidence_interval()

install.packages('ggeasy')
library(ggeasy)
km_derstine_smi <- km90_default + coord_cartesian(xlim = c(0, 8)) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent, 
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(breaks = 0:9, expand = c(0.02, 0)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 1)) +
  ggeasy::easy_y_axis_labels_size(size = 30) +
  ggeasy::easy_x_axis_labels_size(size = 30) +
  ggeasy::easy_y_axis_title_size(size = 30) +
  ggeasy::easy_x_axis_title_size(size = 30) +
  theme(legend.text=element_text(size=30))

km_derstine_smi
survdiff(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ cat_smi, data = thoracic_nsclc_included_noncontrast_cat_smi)

#####km of lowest quartile stratified by sex and contrast

thoracic_nsclc_included_noncontrast_female$quatile_smd <- cut(thoracic_nsclc_included_noncontrast_female$mean_t12m_hu,
                                                              breaks=c(-30, quantile(thoracic_nsclc_included_noncontrast_female$mean_t12m_hu, probs = 0.25, na.rm=TRUE), 150),
                                                              labels=c('Lowest Quartile - SMD', 'Other Quartiles  - SMD'))
thoracic_nsclc_included_noncontrast_male$quatile_smd <- cut(thoracic_nsclc_included_noncontrast_male$mean_t12m_hu,
                                                            breaks=c(-30, quantile(thoracic_nsclc_included_noncontrast_male$mean_t12m_hu, probs = 0.25, na.rm=TRUE), 150),
                                                            labels=c('Lowest Quartile - SMD', 'Other Quartiles  - SMD'))
thoracic_nsclc_included_contrast_female$quatile_smd <- cut(thoracic_nsclc_included_contrast_female$mean_t12m_hu,
                                                           breaks=c(-30, quantile(thoracic_nsclc_included_contrast_female$mean_t12m_hu, probs = 0.25, na.rm=TRUE), 150),
                                                           labels=c('Lowest Quartile - SMD', 'Other Quartiles  - SMD'))
thoracic_nsclc_included_contrast_male$quatile_smd <- cut(thoracic_nsclc_included_contrast_male$mean_t12m_hu,
                                                         breaks=c(-30, quantile(thoracic_nsclc_included_contrast_male$mean_t12m_hu, probs = 0.25, na.rm=TRUE), 150),
                                                         labels=c('Lowest Quartile - SMD', 'Other Quartiles  - SMD'))

thoracic_nsclc_included_quartile_smd <- rbind(thoracic_nsclc_included_noncontrast_female,thoracic_nsclc_included_noncontrast_male, thoracic_nsclc_included_contrast_female, thoracic_nsclc_included_contrast_male)



km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ quatile_smd, data = thoracic_nsclc_included_quartile_smd) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  )+ 
  geom_step(size = 2)+
  add_confidence_interval() +
  add_risktable( risktable_height = 0.33,
                 size = 7, # increase font size of risk table statistics
                 theme =   # increase font size of risk table title and y-axis label
                   list(
                     theme_risktable_default(axis.text.y.size = 22,
                                             
                                             plot.title.size = 22),
                     theme(plot.title = element_text(face = "bold"),
                           legend.text = element_text(size = 20)))) 
km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ quatile_smd, data = thoracic_nsclc_included_quartile_smd) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  )+
  geom_step(size = 2)+  
  add_confidence_interval() 

km_smd_quartile <- km90_default + coord_cartesian(xlim = c(0, 8)) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent, 
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(breaks = 0:9, expand = c(0.02, 0)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 1)) +
  ggeasy::easy_y_axis_labels_size(size = 30) +
  ggeasy::easy_x_axis_labels_size(size = 30) +
  ggeasy::easy_y_axis_title_size(size = 30) +
  ggeasy::easy_x_axis_title_size(size = 30) +
  theme(legend.text=element_text(size=30))

km_smd_quartile
survdiff(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ quatile_smd, data = thoracic_nsclc_included_quartile_smd)


###########km for imat index 


thoracic_nsclc_included_noncontrast_female$quatile_imat <- cut(thoracic_nsclc_included_noncontrast_female$index_t12imat,
                                                               breaks=c(-30, quantile(thoracic_nsclc_included_noncontrast_female$index_t12imat, probs = 0.75, na.rm=TRUE), 150),
                                                               labels=c('Other Quartiles  - IMAT index', 'Highest Quartile - IMAT index'))
thoracic_nsclc_included_noncontrast_male$quatile_imat <- cut(thoracic_nsclc_included_noncontrast_male$index_t12imat,
                                                             breaks=c(-30, quantile(thoracic_nsclc_included_noncontrast_male$index_t12imat, probs = 0.75, na.rm=TRUE), 150),
                                                             labels=c('Other Quartiles  - IMAT index', 'Highest Quartile - IMAT index'))
thoracic_nsclc_included_contrast_female$quatile_imat <- cut(thoracic_nsclc_included_contrast_female$index_t12imat,
                                                            breaks=c(-30, quantile(thoracic_nsclc_included_contrast_female$index_t12imat, probs = 0.75, na.rm=TRUE), 150),
                                                            labels=c('Other Quartiles  - IMAT index', 'Highest Quartile - IMAT index'))
thoracic_nsclc_included_contrast_male$quatile_imat <- cut(thoracic_nsclc_included_contrast_male$index_t12imat,
                                                          breaks=c(-30, quantile(thoracic_nsclc_included_contrast_male$index_t12imat, probs = 0.75, na.rm=TRUE), 150),
                                                          labels=c('Other Quartiles  - IMAT index', 'Highest Quartile - IMAT index'))

thoracic_nsclc_included_quartile_imat <- rbind(thoracic_nsclc_included_noncontrast_female,thoracic_nsclc_included_noncontrast_male, thoracic_nsclc_included_contrast_female, thoracic_nsclc_included_contrast_male)

thoracic_nsclc_included_quartile_imat$quatile_imat <- factor(thoracic_nsclc_included_quartile_imat$quatile_imat, levels = c('Highest Quartile - IMAT index', 'Other Quartiles  - IMAT index'))

km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ quatile_imat, data = thoracic_nsclc_included_quartile_imat) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  )+ 
  geom_step(size = 2)+
  add_confidence_interval() +
  add_risktable( risktable_height = 0.33,
                 size = 7, # increase font size of risk table statistics
                 theme =   # increase font size of risk table title and y-axis label
                   list(
                     theme_risktable_default(axis.text.y.size = 22,
                                             
                                             plot.title.size = 22),
                     theme(plot.title = element_text(face = "bold"),
                           legend.text = element_text(size = 20)))) 
km90_default <- survfit2(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ quatile_imat, data = thoracic_nsclc_included_quartile_imat) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  )+ 
  geom_step(size = 2)+  
  add_confidence_interval() 

km_imat_quartile <- km90_default + coord_cartesian(xlim = c(0, 8)) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent, 
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(breaks = 0:9, expand = c(0.02, 0)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 1)) +
  ggeasy::easy_y_axis_labels_size(size = 30) +
  ggeasy::easy_x_axis_labels_size(size = 30) +
  ggeasy::easy_y_axis_title_size(size = 30) +
  ggeasy::easy_x_axis_title_size(size = 30) +
  theme(legend.text=element_text(size=30))

km_imat_quartile
survdiff(Surv(time_delta_surgery_years, final_vital_status_1_dead) ~ quatile_imat, data = thoracic_nsclc_included_quartile_imat)


# 4 figures arranged in 2 rows and 2 columns
attach(mtcars)
par(mfrow=c(1,2))
km_derstine_smd
km_derstine_smi

attach(mtcars)
par(mfrow=c(1,2))
km_smd_quartile
km_imat_quartile

####################descriptives
library(ggpubr)
thoracic_nsclc_included_stage_I_II$contrast_iv_thoracic_word <- thoracic_nsclc_included_stage_I_II$contrast_iv_thoracic
thoracic_nsclc_included_stage_I_II$contrast_iv_thoracic_word <- recode(thoracic_nsclc_included_stage_I_II$contrast_iv_thoracic_word, '0' = 'Non-Contrast', '1' = 'Contrast')
# grouped boxplot pericardial

set.seed(1) # ensure jitter pts are the same position
smi_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=index_t12m)) + 
  labs(y = expression(paste("Skeletal Muscle Index [",cm^2,'/',m^2,']')), x = NULL) + geom_boxplot(outlier.alpha=0) + 
  geom_jitter(alpha=0.1,width=0.1) +  theme_pander(base_size =15) +
  facet_wrap(~gender)

smi_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=index_t12m)) + 
  labs(y = expression(paste("Skeletal Muscle Index [",cm^2,'/',m^2,']')), x = NULL) + geom_boxplot(outlier.alpha=0) + 
  geom_jitter(alpha=0.1,width=0.1) +  theme_pander(base_size =15) +
  facet_wrap(~gender)




attach(mtcars)
par(mfrow=c(3,2))




ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=index_t12m, fill=gender)) + labs(y = "Skeletal Muscle Index [cm^2/m^2]", x = NULL)  + theme(text = element_text(size=30)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())
ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=mean_t12m_hu, fill=gender)) + labs(y = "Skeletal Muscle Density [HU]", x = NULL)  + theme(text = element_text(size=30)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())
ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=index_t12i, fill=gender)) + labs(y = "IMAT Index [HU]", x = NULL)  + theme(text = element_text(size=30)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())
ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=mean_t12imat_hu, fill=gender)) + labs(y = "IMAT Density [HU]", x = NULL)  + theme(text = element_text(size=30)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())
ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=index_t12s, fill=gender)) + labs(y = "SAT Index [HU]", x = NULL)  + theme(text = element_text(size=30)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())
ggplot(thoracic_nsclc_included_stage_I_II, aes(x=contrast_iv_thoracic_word, y=mean_t12sat_hu, fill=gender)) + labs(y = "SAT Density [HU]", x = NULL)  + theme(text = element_text(size=30)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())

smi_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=gender, y=index_t12m, fill=contrast_iv_thoracic_word)) + labs(y = "Skeletal Muscle Index [cm^2/m^2]", x = NULL)  + theme(text = element_text(size=15)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())

smd_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=gender, y=mean_t12m_hu, fill=contrast_iv_thoracic_word)) + labs(y = "Skeletal Muscle Density [HU]", x = NULL)  + theme(text = element_text(size=15)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())

imat_i_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=gender, y=index_t12i, fill=contrast_iv_thoracic_word)) + labs(y = "IMAT Index [HU]", x = NULL)  + theme(text = element_text(size=15)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())

imat_hu_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=gender, y=mean_t12imat_hu, fill=contrast_iv_thoracic_word)) + labs(y = "IMAT Density [HU]", x = NULL)  + theme(text = element_text(size=15)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())

sat_i_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=gender, y=index_t12s, fill=contrast_iv_thoracic_word)) + labs(y = "SAT Index [HU]", x = NULL)  + theme(text = element_text(size=15)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())

sat_hu_boxplot <- ggplot(thoracic_nsclc_included_stage_I_II, aes(x=gender, y=mean_t12sat_hu, fill=contrast_iv_thoracic_word)) + labs(y = "SAT Density [HU]", x = NULL)  + theme(text = element_text(size=15)) + 
  geom_boxplot()  +  
  theme(legend.title = element_blank())

figure <- ggarrange(smi_boxplot, smd_boxplot, imat_i_boxplot, imat_hu_boxplot, sat_i_boxplot, sat_hu_boxplot,
                    labels = c("Muscle Index", "Muscle Density", "IMAT Index", "IMAT Density", "SAT Index", "SAT Density"),
                    ncol = 2, nrow = 3,
                    
                    vjust = 1.5)


thoracic_nsclc_included %>%
  tbl_summary(
    by = hospital,
    include = c( index_t12m, mean_t12m_hu, index_t12i, mean_t12imat_hu, index_t12s, mean_t12sat_hu),
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{N_nonmiss}/{N_obs} [{p_nonmiss}%]",
                                     "{median} ({p25}, {p75})", 
                                     "{mean} ({sd})"),
    
    label = c(index_t12m ~ "Muscle Index T12[cm^2]",
              mean_t12m_hu ~ "Muscle Density T12[HU]", 
              index_t12i ~ "IMAT Index T12[cm^2]", 
              mean_t12imat_hu ~ "IMAT Density T12[HU]", 
              index_t12s ~ "SAT Index T12[cm^2]", 
              mean_t12sat_hu ~ "SAT Density T12[HU]"
              
    ),
    missing = "no"
    
    
    
  ) %>%
  
  add_overall() %>%
  
  modify_table_styling(
    columns = label) %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Hospital**") %>%
  bold_labels()


l3mghpat_packyears <- thoracic_nsclc_included_stageI_II 
l3mghpat_packyears$packyears[l3mghpat_packyears$packyears == 0] <- NA
l3mghpat_packyears$hospital <- factor(l3mghpat_packyears$hospital, levels = c('MGH', 'WashU', 'UChicago'))
l3mghpat_packyears <- l3mghpat_packyears %>% mutate(
  cerebrvascular_history_cat = cerebrvascular_history
)

l3mghpat_packyears <- l3mghpat_packyears %>% mutate(
  cerebrvascular_history_cat = recode(cerebrvascular_history_cat, 'Any reversible event' = 'Yes', 'Cerebrovascular Accident - CVA' = 'Yes', 'Transient Ischemic Attack - TIA' = 'Yes', 'Any irreversible event' ='Yes')
)

droplevels(thoracic_nsclc_included_stage_I_II$stage_tnm_cat)
l3mghpat_packyears %>%
  tbl_summary(
    by = hospital,
    include = c( gender,
                 ageatsurgery,
                 race_3_cat,
                 bmi_calculated,
                 cigarette_smoking, 
                 packyears, 
                 zubrod_score_3_cat, 
                 asa_3_cat, 
                 fev1_predicted, 
                 dlco_predicted, 
                 stage_tnm,
                 open_vs_closed,
                 neoadjuvant_therapy_general, 
                 neoadjuvant_therapy_general_cat,
                 cancer_histology, 
                 diabetes, 
                 copd, 
                 
                 hypertension,
                 coronary_artery_disease,
                 
                 peripheral_vascular_disease, 
                 cerebrvascular_history,
                 cerebrvascular_history_cat,
                 laa_950),
    statistic = list(all_continuous() ~ "{median} [{p25},{p75}]",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    label = c(gender ~ "Gender",
              ageatsurgery ~ "Age at Surgery", 
              race_3_cat ~ "Race", 
              bmi_calculated ~ "BMI", 
              cigarette_smoking ~ "Smoking Status", 
              packyears ~ "Packyears", 
              zubrod_score_3_cat ~ "Zubrod Score",
              asa_3_cat ~ "ASA Classification", 
              fev1_predicted ~ "FEV1% predicted", 
              dlco_predicted ~ "DLCO% predicted",
              stage_tnm ~ "Tumor Stage",
              open_vs_closed ~ "Surgial Approach", 
              neoadjuvant_therapy_general ~ "Neoadjuvant Therapy",
              cancer_histology ~ "Lung cancer histology",
              diabetes ~ "Diabetes", 
              copd ~ "COPD",
              
              hypertension ~ "Hypertension",
              coronary_artery_disease ~ "Coronary Artery Disease", 
              peripheral_vascular_disease ~ "Peripheral Vascular Disease", 
              cerebrvascular_history ~ "Cerebrovascular History", 
              
              
              laa_950 ~ " %LAA-950 (Emphysema)"
    ),
    sort = c(cancer_histology ~ "frequency"),
    
    missing_text = "unknown"
    
    
  ) %>%
  add_n() %>%
  add_overall() %>%
  modify_table_styling(
    columns = label,
    rows = label == "Packyears",
    footnote = "Conditional statement: Only Current smokers and Past smokers considered n = 790") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Hospital**") %>%
  bold_labels()

thoracic_nsclc_included_stageI_II <- thoracic_nsclc_included_stageI_II %>% mutate(any_postoperative_events_3cat = case_when(minor_events +major_events == 0 ~ '0',
                                                                                                                            minor_events +major_events ==1 ~ '1',
                                                                                                                            minor_events +major_events >= 2 ~ '>=2'))                                    
thoracic_nsclc_included_stageI_II <- thoracic_nsclc_included_stageI_II %>% mutate(pulmonary_postoperative_events_cat = case_when(pulmonary_postoperative_events == 0 ~ '0',
                                                                                                                                 pulmonary_postoperative_events ==1 ~ '1',
                                                                                                                                 pulmonary_postoperative_events >= 2 ~ '>=2'))
thoracic_nsclc_included_stageI_II <- thoracic_nsclc_included_stageI_II %>% mutate(minor_events_cat = case_when(minor_events == 0 ~ '0',
                                                                                                               minor_events ==1 ~ '1',
                                                                                                               minor_events ==2 ~ '2',
                                                                                                               minor_events ==3 ~ '3',
                                                                                                               minor_events >= 4 ~ '>=4'))
thoracic_nsclc_included_stageI_II <- thoracic_nsclc_included_stageI_II %>% mutate(major_events_cat = case_when(major_events == 0 ~ '0',
                                                                                                               major_events ==1 ~ '1',
                                                                                                               major_events >= 2 ~ '>=2'))

thoracic_nsclc_included_outcome_des <- thoracic_nsclc_included_stageI_II
thoracic_nsclc_included_outcome_des <- thoracic_nsclc_included_outcome_des %>% mutate(
  final_vital_status_1_dead = recode(final_vital_status_1_dead, '1' = "No", '2' = 'Yes'))
thoracic_nsclc_included_outcome_des$any_postoperative_events_3cat <- factor(thoracic_nsclc_included_outcome_des$any_postoperative_events_3cat, levels = c('0', '1', '>=2'))
thoracic_nsclc_included_outcome_des$hospital <- factor(thoracic_nsclc_included_outcome_des$hospital, levels = c('MGH', 'WashU', 'UChicago'))

thoracic_nsclc_included_outcome_des %>%
  tbl_summary(
    by = stage_tnm_cat,
    include = c( final_vital_status_1_dead, minor_events_cat, major_events_cat, any_postoperative_events_3cat, pulmonary_postoperative_events_cat, status_30daysaftersurgery),
    
    statistic = all_categorical() ~ "{n} ({p}%)",
    
    label = c(final_vital_status_1_dead ~ "Death",
              any_postoperative_events_3cat ~ "Any Postoperative Events", 
              pulmonary_postoperative_events_cat ~ "Pulmonary Postoperative Events"),
    sort = c( pulmonary_postoperative_events_cat~ "frequency",  minor_events_cat~ "frequency", major_events_cat ~ "frequency"),
    missing = "no"
    
    
    
  ) %>%
  
  add_overall() %>%
  
  modify_table_styling(
    columns = label) %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Hospital**") %>%
  bold_labels()

#########################image aqiasiton

# dataname : image_df 
colnames(image_df)[1] <- "id_cmw"
image_df <- left_join(thoracic_nsclc_included_stage_I_II, image_df, by = "id_cmw")
image_df$SliceThickness
colnames(image_df)

tabel <- image_df %>%
  tbl_summary(
    
    include = c( Manufacturer, ManufacturerModelName, KVP, SliceThickness),
    
    statistic = all_categorical() ~ "{n} ({p}%)",
    #sort = c( pulmonary_postoperative_events_cat ~ "frequency")
    
  ) %>%
  
  modify_table_styling(
    columns = label) %>%
  
  bold_labels()


#####percentiles 

thoracic_nsclc_included_contrast <- thoracic_nsclc_included_stage_I_II %>% filter(contrast_iv_thoracic == 1)
thoracic_nsclc_included_noncontrast <- thoracic_nsclc_included_stage_I_II %>% filter(contrast_iv_thoracic == 0)
thoracic_nsclc_included_noncontrast_female <- thoracic_nsclc_included_noncontrast[thoracic_nsclc_included_noncontrast$gender == "Female",]
thoracic_nsclc_included_noncontrast_male <- thoracic_nsclc_included_noncontrast[thoracic_nsclc_included_noncontrast$gender == "Male",]
thoracic_nsclc_included_contrast_female <- thoracic_nsclc_included_contrast[thoracic_nsclc_included_contrast$gender == "Female",]
thoracic_nsclc_included_contrast_male <- thoracic_nsclc_included_contrast[thoracic_nsclc_included_contrast$gender == "Male",]

mean(thoracic_nsclc_included_noncontrast_female$mean_t12m_hu, na.rm = TRUE)
sd(thoracic_nsclc_included_noncontrast_female$mean_t12m_hu, na.rm = TRUE)
quantile(thoracic_nsclc_included_noncontrast_female$mean_t12m_hu, c(.05, .10, .25, .50, .75, .90, .95), na.rm = TRUE)




