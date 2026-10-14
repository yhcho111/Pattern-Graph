
# Data Generation:
m1 = Scenario[s, "m1"]; m2 = Scenario[s, "m2"]; 
s11 = Scenario[s, "s11"]; rho = Scenario[s, "rho"]; s22 = Scenario[s, "s22"];

o_s2_0 = Scenario[s, "o_s2_0"]; o_s2_L1 = Scenario[s,"o_s2_L1"];
o_s3_0 = Scenario[s, "o_s3_0"]; o_s3_L2 = Scenario[s,"o_s3_L2"];
o_s4_0 = Scenario[s, "o_s4_0"]; 

n = Scenario[s, "n"];
s12 = rho*sqrt(s11*s22)

T.mean = c(m1, m2)
T.Sigma = matrix(c(s11,s12,s12,s22), ncol = 2)
t.mean = c(m1, m2, m1+m2) 
names(t.mean ) = paste0("t_",c("L1","L2","L1L2"))

#Create Sample n= 1000 ->nc
samp = mvrnorm(n = n, mu = T.mean, Sigma = T.Sigma) %>% as.data.frame(.) %>% setNames(paste0("L",1:2)) %>% mutate(ID = row_number()) 
samp.mean = samp %>% mutate(samp_L1 = L1, samp_L2 = L2, samp_L1L2 = L1+L2) %>% dplyr::select(contains("samp")) %>% apply(.,2,mean)


# Create Missing:
pic_by_mis = samp %>%
  mutate(O_S2 = exp(o_s2_0+o_s2_L1*L1), 
         O_S3 =exp(o_s3_0+o_s3_L2*L2), 
         O_S4 =exp(o_s4_0) ) %>% # Add intercept: more realistic model assumption
  mutate(Q_S1 = 1, Q_S2 = O_S2, Q_S3 = O_S3, Q_S4 = O_S4*(Q_S2+Q_S3)) %>%
  mutate(pi_L = 1/(Q_S1+Q_S2+Q_S3+Q_S4)) %>%
  mutate(P_S1 = pi_L, P_S2 = Q_S2*pi_L, P_S3 = Q_S3*pi_L, P_S4 = Q_S4*pi_L) ## create P(R=r|L) This is real probability for creating missingness

nc_by_mis = pic_by_mis %>% 
  apply(.,1,function(.){
    T.probc = .[c("P_S1","P_S2","P_S3","P_S4")] %>% unlist(.)
    rmultinom(1,size=1, T.probc)}) %>% as.data.frame(.) %>% t(.) %>% as.data.frame(.) %>%
  setNames(paste0("S",1:4))# %>% mutate(nc = S1+S2+S3+S4)

nc_by_mis = cbind(pic_by_mis, nc_by_mis) %>% as.data.frame(.) %>% dplyr::select(ID, L1:L2,S1:S4) 
apply(nc_by_mis,2,mean)

## Dataframe with missingness:
C_R = nc_by_mis%>%
  mutate(L2 = ifelse(S2==1 |S4 ==1, NA,L2),
         L1 = ifelse(S3==1 |S4 ==1, NA,L1)) %>% 
  mutate(pts = case_when((S1==1) ~ "S1",
                         (S2==1) ~ "S2",
                         (S3==1) ~ "S3",
                         (S4==1)  ~ "S4")) #%>% 
#mutate(ID = row_number())

S1 = c(T,T); S2 = c(T,F); S3 = c(F,T);S4 = c(F,F);


#Inverse probability weighting:
Orlr = OrLr(C_R)
for(i in 1:length(Orlr)){assign(names(Orlr)[[i]], Orlr[[i]])} #"coef.S2"    "coef.S3"    "coef.S4"  
IPW.res = IPW(C_R, Orlr)

IPW.sum = IPW.res %>% mutate(IPW_L1 = (L1)/pi_L,
                             IPW_L2 = (L2)/pi_L,
                             IPW_L1L2 = (L1+L2)/pi_L)%>%
  dplyr::select(contains("IPW_")) %>% apply(.,2,sum)
IPW.denom = IPW.res %>% mutate(denom = 1/pi_L) %>% dplyr::select(denom) %>% sum(.) %>% unlist(.)
IPW.mean = IPW.sum/IPW.denom


M = 100 # number of imputation

#Fractional Imputation:
# #EM algorithm:
ini.norm_pijk = EM_initials(C_R)
imp.Y.vals.kim = impY_hm_piL_kim(ini.norm_pijk, C_R, coef.S2, coef.S3, coef.S4, S1, S2, S3, S4,M)
est.pijk.kim = EM_pars.kim(ini.norm_pijk, imp.Y.vals.kim, C_R, M);
FI.EM.kim = FI_EM.df.kim(est.pijk.kim, imp.Y.vals.kim, coef.S2, coef.S3, coef.S4, S1, S2, S3, S4,M);
FI.EM.sum.kim = FI.EM.kim%>% mutate(FI.EM.K_L1 = (L1)*final_wt.t,
                                    FI.EM.K_L2 =L2*final_wt.t,
                                    FI.EM.K_L1L2 =(L1+L2)*final_wt.t)%>%
  dplyr::select(contains("FI.EM.K")) %>% apply(.,2,sum)
FI.EM.mean.kim = FI.EM.sum.kim/n

## RA
RAs = RA_imp.df(C_R, S1, S2, S3,S4, M)
RA.mean = RAs$RA_mean
RA_dat = RAs$RA_dat

## AIPW of Chen:
SEMI.mean = Semi_imp(C_R, RA_dat,coef.S2, coef.S3, coef.S4, M, n, IPW.res, ini.norm_pijk)


### Variance Estimation:
FI.EM.K.sd =lapply(1:100, function(.){return(est.FI.EM.b(est.pijk.kim, coef.S2, coef.S3, coef.S4, C_R, n, MODEL="K"))}) %>% Reduce("rbind",.) %>% apply(.,2,sd)


mean_res = c(t.mean, samp.mean, IPW.mean, RA.mean, FI.EM.mean.yang, FI.EM.mean.T, FI.EM.mean.kim, SEMI.mean, SEMI2.mean.yang, SEMI2.mean.kim, FI.EM.Y.sd, FI.EM.K.sd)

