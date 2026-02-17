all(sapply(c("lme4","sampling","data.table","dplyr","doParallel","fastGHQuad","merDeriv","MASS","numDeriv","Matrix","tidyr","reldist","plyr"), require, character.only = TRUE))
source("fns_Algorithm.R")
source("SEMI.R")
source("SEMI_FI.R")
# source("MEAN_IMPUTATION.R")
s = 8; T.seed = 31;
# Parameters:

lambda = 0.75;
a_i = 0.125; b_j = 0.12; c_k = 0.15;
ab_ij = -0.1; ac_ik=-0.1; bc_jk = 0.12;
o_s2_0 = -0.7;o_s2_L1 = 0.5;o_s2_L2 = 0.5;
o_s3_0 = -0.7; o_s3_L1 = 0.55;  o_s3_L3 = 0.45;
o_s4_0 = -0.5; o_s4_L1 = 0.5;
n = 1000;

res_list = NULL
methods = paste0(c("t.","pop.","samp","IPW.","FI.","RA."),"mean")
parameters = c("Y3","Y2Y3","Y1Y2Y3")
variable_names = expand_grid(methods,parameters) %>% apply(.,1,function(.){paste0(.,collapse = "_")})
n_iter = 1000

mean_res_list <- vector("list", n_iter)

for(m in 1:1000){
  m = m + 1;print(m);
  L1= c(-1,1);L2 = c(-1,1);L3 = c(-1,1); # log-linear: sum-to-zero constraints
  
  # Create Population:
  C = expand.grid(L1,L2,L3) %>% setNames(paste0("L",1:3)) %>% 
    mutate(lambda = lambda, a_i = L1*a_i, b_j = L2*b_j, c_k = L3*c_k, 
           ab_ij = L1*L2*ab_ij, ac_ik = L1*L3*ac_ik, bc_jk = L2*L3*bc_jk,
           pijk = exp(lambda + a_i+b_j+c_k+ab_ij+ac_ik+bc_jk)) %>%
    mutate(pijk = pijk/sum(pijk)) %>% mutate(L1 = ifelse(L1==-1,0,L1),
                                             L2 = ifelse(L2==-1,0,L2),
                                             L3 = ifelse(L3==-1,0,L3)) %>% 
    arrange(desc(L1), desc(L2), desc(L3)) 
  C$nc =rmultinom(1,size=n, C$pijk)
  
  t.mean <-  C  %>% mutate(t_Y2 = (L2)*pijk,
                           t_Y3 = (L3)*pijk,
                           t_Y2Y3= (L2*L3)*pijk,
                           t_Y1Y2Y3= (L1*L2*L3)*pijk,
                           t_Y10Y2Y3= ((1-L1)*L2*L3)*pijk,
                           t_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*pijk) %>% 
    dplyr::select(contains("t")) %>% apply(.,2,sum)
  
  # Create Missing:
  pic_by_mis = C %>%
    mutate(O_S2 = exp(o_s2_0+o_s2_L1*L1+o_s2_L2*L2), 
           O_S3 = exp(o_s3_0+o_s3_L1*L1+o_s3_L3*L3),
           O_S4 = exp(o_s4_0+o_s4_L1*L1) ) %>% # Add intercept: more realistic model assumption
    mutate(Q_S1 = 1, Q_S2 = O_S2, Q_S3 = O_S3, Q_S4 = O_S4*(Q_S2+Q_S3)) %>%
    mutate(pi_L = 1/(Q_S1+Q_S2+Q_S3+Q_S4)) %>%
    mutate(P_S1 = pi_L, P_S2 = Q_S2*pi_L, P_S3 = Q_S3*pi_L, P_S4 = Q_S4*pi_L) ## create P(R=r|L) This is real probability for creating missingness
  nc_by_mis = pic_by_mis %>% 
    apply(.,1,function(.){
      T.nc = .["nc"] %>% unlist(.)
      T.probc = .[c("P_S1","P_S2","P_S3","P_S4")] %>% unlist(.)
      rmultinom(1,size=T.nc, T.probc)}) %>% as.data.frame %>% t(.) %>% as.data.frame(.) %>%
    setNames(paste0("S",1:4))# %>% mutate(nc = S1+S2+S3+S4)
  nc_by_mis = cbind(pic_by_mis, nc_by_mis) %>% as.data.frame(.) %>% dplyr::select(L1:L3,nc, S1:S4)
  
  sum(nc_by_mis$S1)/n
  
  
  # nc_by_mis %>% apply(.,1,sum)
  C_R  = C %>% left_join(., nc_by_mis) %>% 
    mutate(cate = paste0("c",1:nrow(.))) %>% dplyr::select(L1:L3,nc,S1:S4,cate)
  
  S1 = c(T,T,T); S2 = c(T,T,F); S3 = c(T,F,T); S4 = c(T,F,F);
  
  
  samp.m.sum = C_R  %>% mutate(samp_Y2 = (L2)*nc,
                               samp_Y3 = (L3)*nc,    
                               samp_Y2Y3 = (L2*L3)*nc,
                               samp_Y1Y2Y3 = (L1*L2*L3)*nc,
                               samp_Y10Y2Y3= ((1-L1)*L2*L3)*nc,
                               samp_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*nc) %>% 
    dplyr::select(contains("samp_")) %>% apply(.,2,sum)
  samp.mean = samp.m.sum/n
  
  
  
  ## long_form for EM algorithm:
  C_R.pt = C_R %>% dplyr::select(cate, L1, L2, L3)
  C_R.cnt = C_R %>% dplyr::select(cate, S1:S4) %>% gather(., pt, count, S1:S4)
  C_R.cnt.pt = C_R.cnt %>% left_join(., C_R.pt) %>% 
    mutate(L3 = ifelse(pt == "S2"|pt == "S4",NA,L3),
           L2 = ifelse(pt == "S3"|pt == "S4",NA,L2))
  # C_R %>% .[,-9] %>% apply(.,2,sum)
  # sum(C_R$S1)/n
  
  #Inverse probability weighting:
  Orlr = OrLr(C_R); for(i in 1:length(Orlr)){assign(names(Orlr)[[i]], Orlr[[i]])};
  IPW.res = IPW(C_R, Orlr)
  
  IPW.sum = IPW.res %>% mutate(IPW_Y2 = (L2)*S1/pi_L,
                               IPW_Y3 = (L3)*S1/pi_L,
                               IPW_Y2Y3 = (L2*L3)*S1/pi_L,
                               IPW_Y1Y2Y3 = (L1*L2*L3)*S1/pi_L,
                               IPW_Y10Y2Y3= ((1-L1)*L2*L3)*S1/pi_L,
                               IPW_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*S1/pi_L) %>%
    dplyr::select(contains("IPW_")) %>% apply(.,2,sum)
  IPW.denom = IPW.res %>% mutate(ipw.denom = S1/pi_L) %>% dplyr::select(ipw.denom) %>% sum(.) %>% unlist
  IPW.mean = IPW.sum/IPW.denom
  # IPW.mean = IPW.sum/n
  # cbind(t.mean, IPW.mean)
  # # Model Parameter Estimation
  IPW.res_FI = IPW.res %>% 
    mutate(P_S2 = Q_S2.L*pi_L, P_S3 = Q_S3.L*pi_L, P_S4 = Q_S4.L*pi_L, P_S1 = pi_L, lat = P_S1 + P_S2 + P_S3 + P_S4) %>%
    dplyr::select(L1:L3,P_S1, P_S2, P_S3, P_S4)
  
  long.pop = C_R.cnt.pt[rep(1:nrow(C_R.cnt.pt), times = C_R.cnt.pt$count),] %>% dplyr::select(pt, L1:L3)
  short.pop = long.pop %>% group_by(pt,L1,L2,L3) %>% dplyr::summarize(count = n())
  
  
  # M = 100 # number of imputation
  SEMI_RA = RA_SEMI.imp(C_R, short.pop, S1, S2, S3, S4, Orlr, n, IPW.res, C_R.cnt.pt, hlist)
  SEMI.mean = SEMI_RA$SEMI.mean
  RA.mean = SEMI_RA$RA.mean %>% unlist(.)
  m_LR_R = SEMI_RA$m_LR_R
  # print(2)
  #Fractional Imputation:
  # #EM algorithm:
  ini.norm_pijk = rep(1/8,8)#C_R$S1/sum(C_R$S1) 
  est.pijk = EM_pijk(ini.norm_pijk ,C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI,n)
  
  FI.EM = FI_EM.df(est.pijk, long.pop, S1, S2, S3,S4, IPW.res_FI, n)
  FI.EM.sum = FI.EM %>% mutate(FI.EM_Y2 = (L2)*final_prob,
                               FI.EM_Y3 =L3*final_prob,
                               FI.EM_Y2Y3 =(L2*L3)*final_prob,
                               FI.EM_Y1Y2Y3 =(L1*L2*L3)*final_prob,
                               FI.EM_Y10Y2Y3= ((1-L1)*L2*L3)*final_prob,
                               FI.EM_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*final_prob
  )%>%
    dplyr::select(contains("FI.EM")) %>% apply(.,2,sum)
  FI.EM.mean = FI.EM.sum/n
  
  ##FI.EM.ST
  est.pijk.ST = EM_pijk.ST(ini.norm_pijk ,C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI,n)
  FI.EM.ST = FI_EM.df(est.pijk.ST, long.pop, S1, S2, S3,S4, IPW.res_FI, n)
  FI.EM.sum.ST = FI.EM.ST %>% mutate(FI.EM.ST_Y2 = (L2)*final_prob,
                                     FI.EM.ST_Y3 =L3*final_prob,
                                     FI.EM.ST_Y2Y3 =(L2*L3)*final_prob,
                                     FI.EM.ST_Y1Y2Y3 =(L1*L2*L3)*final_prob,
                                     FI.EM.ST_Y10Y2Y3= ((1-L1)*L2*L3)*final_prob,
                                     FI.EM.ST_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*final_prob)%>%
    dplyr::select(contains("FI.EM.ST")) %>% apply(.,2,sum)
  FI.EM.ST.mean = FI.EM.sum.ST/n
  
  
  
  
  SEMI2.mean = SEMI_Kim.imp(C_R, short.pop, S1, S2, S4, Orlr, n, IPW.res, C_R.cnt.pt, hlist, est.pijk)  # SEMI.YANG없다!! 
  SEMI.mean2 = SEMI2.mean$SEMI.mean
  SEMI.mean3 = SEMI2.mean$SEMI.mean2
  
  # # Mean Imputation:
  # FI.imp.mean = MEAN_IMF_FI.df(est.pijk, short.pop, IPW.res_FI, hlist, n, S1, S2, S3, S4)
  # print(4)
  # RA.imp.mean = MEAN_IMF_RA.df(C_R, short.pop, n, hlist, S1, S2, S3, S4)
  # print(5)
  
  #######Variance Estimation:
  FI.EM.sd =lapply(1:100, function(.){return(est.FI.EM.b(est.pijk, coef.S2, coef.S3, coef.S4, C_R, n, MODEL="TR"))}) %>% Reduce("rbind",.) %>% apply(.,2,sd)
  FI.EM.ST.sd =lapply(1:100, function(.){return(est.FI.EM.b(est.pijk.ST, coef.S2, coef.S3, coef.S4, C_R, n, MODEL="TR"))}) %>% Reduce("rbind",.) %>% apply(.,2,sd)
  
  mean_res = c(t.mean,samp.mean, IPW.mean, RA.mean, FI.EM.mean,FI.EM.ST.mean, SEMI.mean, SEMI.mean2, SEMI.mean3,FI.EM.sd,FI.EM.ST.sd) #, SEMI.FI.mean,FI.EM.real.mean, RA.FI.mean,RA.mean)#, FI.mean, FI.PA.mean, RA.mean)
  mean_res_list[[m]] <- mean_res
}
