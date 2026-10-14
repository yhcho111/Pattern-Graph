rm(list = ls(all = TRUE))
set.seed(2026, kind = "L'Ecuyer-CMRG")
all(sapply(c("lme4","sampling","data.table","plyr","dplyr","doParallel","fastGHQuad","merDeriv","MASS","numDeriv","Matrix","tidyr"), require, character.only = TRUE))
source("fns_DataAnalysis.R")

#################################################
####### Step 1: Data Exploration:
#################################################

data = read.csv("Myocardial infarction complications Database.csv") %>% 
  mutate(S_AD_KBRIG = ifelse(S_AD_KBRIG <= 90, 1, 0))%>%
  mutate(DEAD = ifelse (LET_IS == 0, 0, 1)) %>%

samp_ID = sample(1:nrow(data), 1700)

data.edit = data %>%
  filter(ID %in% samp_ID) %>%
  dplyr::select(no.missing.vars[k],   NA_KB  ,  S_AD_KBRIG) %>%
  mutate(L1 = eval(parse(text = no.missing.vars[k])),  L2 =   NA_KB  , L3 = S_AD_KBRIG) %>%
  filter(!is.na(L1)) %>%
  dplyr::select(L1:L3) %>%
  group_by(across(L1:L3)) %>% # add argument: make NA as a number for grouping
  summarize(count = n()) %>%
  ungroup(.) %>%
  mutate(S1 = ifelse(!is.na(L2) & !is.na(L3), 1, 0),
         S2 = ifelse(!is.na(L2) & is.na(L3), 1, 0),
         S3 = ifelse(is.na(L2) & !is.na(L3), 1, 0),
         S4 = ifelse(is.na(L2) | is.na(L3), 1, 0)) %>%
  mutate(pt = case_when(
    (S1 == 1) ~ "S1",
    (S2 == 1) ~ "S2",
    (S3 == 1) ~ "S3",
    (S4 == 1) ~ "S4"
  )) %>%
  dplyr::select(pt, count, L1:L3)


n = sum(data.edit$count)


all.pts = expand.grid(L1= c(1,0), L2= c(1,0), L3= c(1,0)) %>% arrange(desc(across(L1:L3))) %>%
  mutate(pt = "S1") %>% dplyr::select(pt, L1:L3) # all patterns pt, L1:L3

#################################################
####### Step 2: Calculate pattern probability:
#################################################
Orlr_PG1 = OrLr(data.edit, PG="PG1")
Orlr_PG2 = OrLr(data.edit, PG="PG2")
Orlr_PG3 = OrLr(data.edit, PG="PG3")
Orlr_PG4 = OrLr(data.edit, PG="PG4")
Orlr_PG5 = OrLr(data.edit, PG="PG5")


IPW.res_PG1 = IPW(data.edit, Orlr_PG1, PG="PG1", all.pts)
IPW.res_PG2 = IPW(data.edit, Orlr_PG2, PG="PG2", all.pts)
IPW.res_PG3 = IPW(data.edit, Orlr_PG3, PG="PG3", all.pts)
IPW.res_PG4 = IPW(data.edit, Orlr_PG4, PG="PG4", all.pts)
IPW.res_PG5 = IPW(data.edit, Orlr_PG5, PG="PG5", all.pts)

##########################################
####### Step 3: Each Method:
##########################################

########################
## Complete Case ##
########################
com_case.ns = all.pts %>% left_join(data.edit %>% filter(pt=="S1")) %>% mutate(count = replace_na(count, 0)) %>%
  mutate(t.pijk = count/sum(count))

cc_est = com_case.ns$count/sum(com_case.ns$count);cc_est

OR_CC = fit.odds(com_case.ns)

########################
## FI: Saturated Model##
########################
IPW.res_FI1 = IPW.res_cal(IPW.res_PG1) #L1:L3, P_S1:P_S4
IPW.res_FI2 = IPW.res_cal(IPW.res_PG2)
IPW.res_FI3 = IPW.res_cal(IPW.res_PG3)
IPW.res_FI4 = IPW.res_cal(IPW.res_PG4)
IPW.res_FI5 = IPW.res_cal(IPW.res_PG5)

ini.norm_pijk = com_case.ns$t.pijk ; names(ini.norm_pijk) = NULL


for(i in 3){ # i for interaction term
  for(j in  1:5){ # j for PG model
    S1 = c(T,T,T); S2 = c(T,T,F); S3 = c(T,F,T); S4 = c(T,F,F);
    
    eval(parse(text = paste0("PG", j, "_I",i," = EM_pijk(ini.norm_pijk, com_case.ns, all.pts, data.edit,
                               S1, S2, S3, S4, IPW.res_FI = IPW.res_FI",j,", n, interaction_terms =", i,")")))
    #est_category_prob:
    
    #est.pijk_PG
    eval(parse(text = paste0("est.pijk_PG", j, "_I", i," = PG", j, "_I", i,"$new.pijk")))
    eval(parse(text = paste0("ML_PG", j, "_I",i," = PG", j, "_I", i,"$new.Q")))
    
    
  }
}


# Make Graph
res.df.save = NULL
for(i in 3){
  for(j in 1:5){
    eval(parse(text = paste0("est.pijk = est.pijk_PG",j,"_I", i)))  
    eval(parse(text = paste0("ML =ML_PG", j, "_I",i)))
    res = c(est.pijk, ML, j, i)
    res.df.save = rbind(res.df.save, res)
  }
}
res.df.save 

baseline_PG = 1
res.df.long = res.df.save %>% as.data.frame(.) %>% 
  setNames(c(paste0("c",1:8),"lik", "PG", "I")) %>% 
  gather(., cate, cate.prob, c1:c8)
res.df.long2 = res.df.long %>% filter(PG !=baseline_PG) %>% dplyr::select(-lik) %>%  # RB of each category compared to PG3
  left_join(., res.df.long %>% filter(PG == baseline_PG) %>% 
              rename(PG_base.prob = cate.prob) %>% dplyr::select(-PG) ) %>% 
  mutate(RB = 100*(cate.prob-PG_base.prob)/PG_base.prob) %>% 
  mutate(PG = case_when( (PG == 2) ~ "PG2",
                         (PG == 3) ~ "PG3",
                         (PG == 4) ~ "PG4",
                         (PG == 5) ~ "PG5")) %>% dplyr::select(-lik)

res.df.long3 = res.df.long %>% filter(PG !=baseline_PG) %>% dplyr::select(-cate.prob) %>% 
  left_join(., res.df.long %>% filter(PG == baseline_PG) %>% rename(PG_base.Lik= lik) %>% dplyr::select(-cate.prob, -PG) ) %>% 
  mutate(Lik.Diff = PG_base.Lik-lik) %>% mutate(PG = case_when( (PG == 2) ~ "PG2",
                                                                (PG == 3) ~ "PG3",
                                                                (PG == 4) ~ "PG4",
                                                                (PG == 5) ~ "PG5")) %>% dplyr::select(-lik, -PG_base.Lik)
res.df.final = res.df.long2 %>% left_join(., res.df.long3)
library(ggplot2)
res.df.final %>% ggplot2::ggplot(., aes(x = PG_base.prob, y = RB, color = PG)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~PG) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),                          
    strip.background = element_rect(fill = "lightgray", color = "black", size = 1), 
    axis.title.x = element_text(size = 12, face = "bold"),                      
    axis.title.y = element_text(size = 12, face = "bold"),                      
    axis.text.x = element_text(size = 10),                         
    axis.text.y = element_text(size = 10),                         
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20)),  # 그래프 제목 크기 및 가운데 정렬
    plot.margin = margin(t = 30, b = 20, r = 10, l = 10),
    panel.spacing = unit(1.5, "lines")  
  ) +
  labs(
    title = "Relative Difference of Estimated Category Probabilities", 
    x = "Baseline PG",                                                
    y = "RD (%)"                                                      
  )


for(i in c(3)){ # i for interaction term
  for(j in  c(1,2,3,4,5)){ # j for PG model
    
    ## Estimation:
    eval(parse(text = paste0("plot(est.pijk_PG2_I", i, ",(est.pijk_PG", j,"_I", i, "-est.pijk_PG3_I", i,")/est.pijk_PG3_I",i,", main = 'PG",j, "_I",i,"'); abline(0,0,col='red')")))
    print(eval(parse(text = paste0("ML_PG", j, "_I",i ))))
  }
}


for(i in 3){ # i for interaction term
  for(j in 1:5){ # j for PG model
    print(paste0("PG",j, "_I",i))
    ML = eval(parse(text = paste0("ML_PG", j, "_I",i)))
    print(ML)
  }
}

# }


OR_IPW = fit.odds(IPW.res_PG1 %>% mutate(count = count/pi_L))
FI_EM.df1 = FI_EM.df(est.pijk_PG1_I3, data.edit, S1, S2, S3,S4, IPW.res_FI1, n) %>%
  mutate(count = count *final_prob )
OR_FI = fit.odds(FI_EM.df1)


Var.OR_FI.BT = mclapply(1:1000, function(.){
  est.FI.EM.b(est.pijk = est.pijk_PG1_I3,  com_case.ns= com_case.ns, 
              all.pts = all.pts, Orlr_PG = Orlr_PG1, 
              PG = "PG1", n, interaction_terms = 3)}, mc.cores = 5) %>% Reduce("rbind", .) #;Var.OR_FI.BT
SD.OR_FI.BT= apply(Var.OR_FI.BT,2,sd, na.rm = T);SD.OR_FI.BT
est = c(OR_CC, OR_IPW,  OR_FI)
CI_U =  est + 1.96*SD.OR_FI.BT
CI_L = est - 1.96*SD.OR_FI.BT
CI = cbind(est,CI_L, CI_U)

# CI_Long
CI = cbind(est, CI_L, CI_U) %>% 
  as.data.frame() %>%              
  setNames(c("Est.", "LL", "UL"))

CI

###############################################################
####### Linear Combination for Research Question #######
###############################################################
## Est for CC ase:OR_CC
vcov_CC = cov(Var.OR_FI.BT[,1:4])

## Est for IPW ase:OR_IPW
vcov_IPW = cov(Var.OR_FI.BT[,5:8])

## Est for IPW ase:OR_FI
vcov_FI = cov(Var.OR_FI.BT[,9:12])
stopifnot(all(c("(Intercept)","L2","L3","L2:L3") %in% names(OR_FI)))
stopifnot(all(colnames(vcov_FI) %in% names(OR_FI)))
nm <- names(OR_FI)

alpha <- 0.05
zcrit <- qnorm(1 - alpha/2)

L_norm <- setNames(rep(0, length(nm)), nm); L_norm["L2"] <- 1                 # Opioid effect at L3=0 (normotension)
L_hypo <- setNames(rep(0, length(nm)), nm); L_hypo[c("L2","L2:L3")] <- 1      # Opioid effect at L3=1 (hypotension)
L_L3_0 <- setNames(rep(0, length(nm)), nm); L_L3_0["L3"] <- 1                 # Hypotension effect at L2=0 (no opioid)
L_L3_1 <- setNames(rep(0, length(nm)), nm); L_L3_1[c("L3","L2:L3")] <- 1      # Hypotension effect at L2=1 (opioid)
L_int  <- setNames(rep(0, length(nm)), nm); L_int["L2:L3"] <- 1               # Interaction (L2×L3)

tab <- rbind(
  lin_summary(L_norm, OR_FI, vcov_FI, "Opioid effect at L3=0 (normotension)"),
  lin_summary(L_hypo, OR_FI, vcov_FI, "Opioid effect at L3=1 (hypotension)"),
  lin_summary(L_L3_0, OR_FI, vcov_FI, "Hypotension effect at L2=0 (no opioid)"),
  lin_summary(L_L3_1, OR_FI, vcov_FI, "Hypotension effect at L2=1 (opioid)"),
  lin_summary(L_int , OR_FI, vcov_FI, "Interaction (L2×L3)")
) |> as.data.frame()

# Clearer Version
tab$sig <- ifelse(as.numeric(tab$p) < alpha, "✴︎", "")
print(tab[, c("contrast","OR","LCL","UCL","p","sig")], row.names = FALSE)

#######################################
######## IPW Estimator ################
#######################################

OR_IPW = IPW.res_PG4  %>% mutate(p_11. = L1*L2*count/pi_L, 
                                 p_10. = L1*(1-L2)*count/pi_L, 
                                 p_01. = (1-L1)*count/pi_L, 
                                 p_00. = (1-L1)*(1-L2)*count/pi_L) %>%
  dplyr::select(contains("p_")) %>% apply(.,2,sum)
OR_IPW = (OR_IPW["p_11."]*OR_IPW["p_00."])/(OR_IPW["p_01."]*OR_IPW["p_10."]); OR_IPW


############################
## Missing Mechanism Check##
#############################
#----------------------------------------------------------------------
# Missingness: logistic regression model with L1 & L2
#----------------------------------------------------------------------
R_NA_KB <- ifelse(is.na(data$NA_KB),1,0) 
fit_missing_NA_KB <- glm(R_NA_KB ~ data$DEAD, family = "binomial")
summary(fit_missing_NA_KB)
exp(cbind(OR = coef(fit_missing_NA_KB), confint(fit_missing_NA_KB)))





