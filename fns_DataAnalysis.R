OrLr= function(data.edit, PG){
  
  uniq.pt_S2 = data.edit %>% filter(pt %in% c("S1", "S2")) %>% dplyr::select(pt, count, L1, L2) %>% 
               group_by(pt,L1, L2) %>% summarize_at(vars(count), sum) %>% spread(., pt, count) %>%
               mutate(L1 = factor(L1), L2 = factor(L2)) %>% 
               mutate(across(S1:S2, ~ replace_na(., 0))) %>% as.data.frame(.)
  logis.S2 = glm(cbind(S2,S1)~L1+L2, family="binomial", data=uniq.pt_S2)
  coef.S2 = coef(logis.S2)
  Sel.Odds.S2 <- apply((uniq.pt_S2 %>% mutate(intercept = 1)  %>% 
                       mutate_if(is.factor, as.character) %>% 
                       mutate_if(is.character, as.numeric) %>% 
                       dplyr::select(intercept, L1,L2) %>% t(.))*coef.S2 ,2,sum) %>% exp(.)
  uniq.pt.S2 = uniq.pt_S2 %>% mutate(Sel.Odds = Sel.Odds.S2)
  coef.S2 = coef.S2
  
  
  uniq.pt_S3 = data.edit %>% filter(pt %in% c("S1", "S3")) %>% dplyr::select(pt, count, L1, L3) %>% 
               group_by(pt,L1, L3) %>% summarize_at(vars(count), sum) %>% spread(., pt, count) %>%
               mutate(L1 = factor(L1), L3 = factor(L3)) %>% 
               mutate(across(S1:S3, ~replace_na(., 0))) %>% as.data.frame(.)
  logis.S3 = glm(cbind(S3,S1)~L1+L3, family="binomial", data=uniq.pt_S3)
  coef.S3 = coef(logis.S3)
  Sel.Odds.S3 <- apply((uniq.pt_S3 %>% mutate(intercept = 1)  %>% 
                          mutate_if(is.factor, as.character) %>% 
                          mutate_if(is.character, as.numeric) %>% 
                          dplyr::select(intercept, L1,L3) %>% t(.))*coef.S3 ,2,sum) %>% exp(.)
  uniq.pt.S3 = uniq.pt_S3 %>% mutate(Sel.Odds = Sel.Odds.S3)
  coef.S3 = coef.S3
  
  
  if(PG=="PG1"){ #Fig 1.b (p14)
    
    uniq.pt_S4 = data.edit %>% filter(pt %in% c("S1", "S4")) %>% dplyr::select(pt, count, L1) %>% 
                 group_by(pt,L1) %>% summarize_at(vars(count), sum) %>% spread(., pt, count) %>%
                 mutate(L1 = factor(L1)) %>% 
                 mutate(across(S1:S4, ~replace_na(., 0))) %>% as.data.frame(.) %>% mutate(Sel.Odds = S4/S1)
    
    logis.S4 = glm(cbind(S4,S1)~L1, family="binomial", data=uniq.pt_S4)
    coef.S4 = coef(logis.S4)
    # Sel.Odds.S4 <- apply((uniq.pt_S4 %>% mutate(intercept = 1)  %>% 
    #                         mutate_if(is.factor, as.character) %>% 
    #                         mutate_if(is.character, as.numeric) %>% 
    #                         dplyr::select(intercept, L1) %>% t(.))*coef.S4 ,2,sum) %>% exp(.)
    # uniq.pt.S4 = uniq.pt_S4 %>% mutate(Sel.Odds = Sel.Odds.S4)
      
  } else if(PG=="PG2"){ #Fig 1.c (p14)
    
    uniq.pt_S4 = data.edit %>% filter(pt != "S1") %>% dplyr::select(pt, count, L1) %>% 
                 group_by(pt,L1) %>% summarize_at(vars(count), sum) %>% spread(., pt, count) %>%
                 mutate(L1 = factor(L1)) %>%  mutate(S2= replace_na(S3,0), S3= replace_na(S3,0))  %>%
                 mutate(S23 = S2+S3) %>% as.data.frame(.)  %>% mutate(Sel.Odds = S4/(S23)) %>% dplyr::select(-S2,-S3) %>%
                 mutate(Sel.Odds = ifelse(is.infinite(Sel.Odds), 1e6, Sel.Odds)) 
    
    logis.S4 = glm(cbind(S4,S23)~L1, family="binomial", data=uniq.pt_S4)
    coef.S4 = coef(logis.S4)
     
  } else if (PG== "PG3"){ #Fig 1.d (p14)
    
    uniq.pt_S4 = data.edit %>%  dplyr::select(pt, count, L1) %>% 
                 group_by(pt,L1) %>% summarize_at(vars(count), sum) %>% spread(., pt, count) %>%
                 mutate(L1 = factor(L1)) %>% mutate(S2= replace_na(S2,0))  %>% 
                 mutate(S123 = S1+S2+S3) %>% as.data.frame(.) %>% mutate(Sel.Odds = S4/(S123)) %>% dplyr::select(-S1,-S2,-S3) %>%
                 mutate(Sel.Odds = ifelse(is.infinite(Sel.Odds), 1e6, Sel.Odds)) 
    
    logis.S4 = glm(cbind(S4,S123)~L1, family="binomial", data=uniq.pt_S4)
    coef.S4 = coef(logis.S4)
  }else if (PG == "PG4") {
    
    uniq.pt_S4 = data.edit %>%  dplyr::select(pt, count, L1) %>% 
                 group_by(pt,L1) %>% summarize_at(vars(count), sum) %>% spread(., pt, count) %>%
                 mutate(L1 = factor(L1)) %>%  as.data.frame(.) %>% mutate(S2= replace_na(S2,0))  %>% mutate(Sel.Odds = S4/(S2)) %>% dplyr::select(-S1,-S3) %>%
                 mutate(Sel.Odds = ifelse(is.infinite(Sel.Odds), 1e6, Sel.Odds)) 
    
    logis.S4 = glm(cbind(S4,S2)~L1, family="binomial", data=uniq.pt_S4)
    coef.S4 = coef(logis.S4)
  }else if(PG == "PG5") {
    
    uniq.pt_S4 = data.edit %>%  dplyr::select(pt, count, L1) %>% 
      group_by(pt,L1) %>% summarize_at(vars(count), sum) %>% spread(., pt, count) %>%
      mutate(L1 = factor(L1)) %>%  as.data.frame(.) %>% mutate(S3= replace_na(S3,0))  %>% mutate(Sel.Odds = S4/(S2)) %>% dplyr::select(-S1,-S2) %>%
      mutate(Sel.Odds = ifelse(is.infinite(Sel.Odds), 1e6, Sel.Odds)) 
    
    logis.S4 = glm(cbind(S4,S3)~L1, family="binomial", data=uniq.pt_S4)
    coef.S4 = coef(logis.S4)
  }
  
  
  
  return(list(uniq.pt.S2 = uniq.pt.S2
              ,uniq.pt.S3 = uniq.pt.S3
              , uniq.pt.S4 = uniq.pt_S4
              , coef.S2 = coef.S4
              , coef.S3 = coef.S3
              , coef.S4 = coef.S4
              )
  )
}

IPW = function(data.edit, Orlr, PG, all.pts){
  core.samp = all.pts %>% left_join(data.edit %>%  filter(pt=="S1") ) %>% mutate(count = replace_na(count, 0))
  AL.1 = lapply(1:nrow(core.samp), function(i){
    T.L = core.samp[i,]
    T.L <- T.L %>% dplyr::mutate(Q_S1.L = 1,
                                 Q_S2.L = Orlr$uniq.pt.S2 %>% filter(L1==T.L$L1 & L2==T.L$L2) %>% dplyr::select(Sel.Odds) %>% unlist(.), #O_{110} = Q_{110}
                                 Q_S3.L = Orlr$uniq.pt.S3 %>% filter(L1==T.L$L1 & L3==T.L$L3) %>% dplyr::select(Sel.Odds) %>% unlist(.)) #O_{101} = Q_{101}
                                 # O_S4.L = Orlr$uniq.pt.S4 %>% filter(L1==T.L$L1) %>% dplyr::select(Sel.Odds) %>% unlist(.), 
                                 # Q_S4.L = O_S4.L*(Q_S2.L+Q_S3.L),
                                 # pi_L = 1/sum(Q_S1.L+Q_S2.L+Q_S3.L+Q_S4.L))})%>% Reduce("rbind",.)
    if(PG== "PG1"){
      
      T.L <- T.L %>% mutate(Q_S4.L = Orlr$uniq.pt.S4 %>% filter(L1==T.L$L1) %>% dplyr::select(Sel.Odds) %>% unlist(.)) 
      
    }else if (PG == "PG2"){
      
      T.L <- T.L %>% mutate(O_S4.L = Orlr$uniq.pt.S4 %>% filter(L1==T.L$L1) %>% dplyr::select(Sel.Odds) %>% unlist(.), 
                            Q_S4.L = O_S4.L*(Q_S2.L+Q_S3.L)) 
      
    }else if (PG== "PG3") {
      
      T.L <- T.L %>% mutate(O_S4.L = Orlr$uniq.pt.S4 %>% filter(L1==T.L$L1) %>% dplyr::select(Sel.Odds) %>% unlist(.), 
                            Q_S4.L = O_S4.L*(Q_S1.L + Q_S2.L + Q_S3.L)) 
    } else if (PG == "PG4"){
      
      T.L <- T.L %>% mutate(O_S4.L = Orlr$uniq.pt.S4 %>% filter(L1==T.L$L1) %>% dplyr::select(Sel.Odds) %>% unlist(.), 
                            Q_S4.L = O_S4.L*( Q_S2.L)) 
    }else if (PG == "PG5") {
      
      T.L <- T.L %>% mutate(O_S4.L = Orlr$uniq.pt.S4 %>% filter(L1==T.L$L1) %>% dplyr::select(Sel.Odds) %>% unlist(.), 
                            Q_S4.L = O_S4.L*( Q_S3.L)) 
    }
    
    T.L <- T.L %>% mutate(pi_L = 1/sum(Q_S1.L+Q_S2.L+Q_S3.L+Q_S4.L))
    
    
    # %>% mutate(P_S1 = 1/sum(Q_S1+Q_S2+Q_S3+Q_S4),
    #            P_S2 = Q_S2/sum(Q_S1+Q_S2+Q_S3+Q_S4),
    #            P_S3 = Q_S3/sum(Q_S1+Q_S2+Q_S3+Q_S4),
    #            P_S4 = Q_S4/sum(Q_S1+Q_S2+Q_S3+Q_S4)
    # )
  }) %>% Reduce("rbind",.)
  
  return(AL.1)
}

## Pattern Prob Calculation:

IPW.res_cal = function(IPW.res){
  
              IPW.res %>% 
              mutate(P_S2 = Q_S2.L*pi_L, P_S3 = Q_S3.L*pi_L, P_S4 = Q_S4.L*pi_L, P_S1 = pi_L, lat = P_S1 + P_S2 + P_S3 + P_S4) %>%
              dplyr::select(L1:L3,P_S1, P_S2, P_S3, P_S4) # L1:L3, P_S1:P_S4
            
          }





## EM algorithm for FI

Q_fn.ST= function(ini.norm_pijk, com_case.ns, all.pts, data.edit, S1, S2, S3, S4, IPW.res_FI, n, interaction_terms){
  
  T.cate.pijk = com_case.ns %>% dplyr::select(-count, pt) 
  
  mis_samp = sapply(paste0("S",2:4), function(pts){ 
    
    T.IPW.res_FI = IPW.res_FI %>% dplyr::select(L1:L3,paste0("P_",pts)) %>% setNames(c(paste0("L",1:3), "P_pt")) # L1:L3, \hat{P}(R=pts|L)
    T.C_R = all.pts %>% dplyr::select (L1:L3)  # Real information we only know marginal information
    T.mrg_cnt = data.edit %>% filter(pt == pts) %>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% 
                group_by(across(contains("L"))) %>% summarise_at(vars(count),sum) # Marginal count by observed variables (Real Situation)
    T.cnt = T.C_R %>% left_join(., T.mrg_cnt) %>% 
            left_join(., T.cate.pijk) %>% #f(L1,L2,L3|\hat{\theta}(t))
            left_join(., T.IPW.res_FI) %>% #P(R=pts|L1,L2, L3)
            mutate(P_pt = replace_na(P_pt,0)) %>%
            mutate(count = replace_na(count,0)) %>%
            mutate(pmis_obs = t.pijk*P_pt) # Calculate pmis_obs = P(L1=l1,L2=l2,L3=l3)*P(R=pts|L) numerator
    T.mrg.cnt = T.cnt %>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% # Remove unobserved variables
                group_by(across(contains("L"))) %>% summarise_at(vars(pmis_obs),sum) %>% # Summation wrt unobserved variables
                dplyr::rename(mrg = pmis_obs)  # calculate denominator: \sum_{unobs}P(L1=l1,L2=l2,L3=l3)*P(R=pts|L)
    T.Final.pt =  T.cnt %>% left_join(.,T.mrg.cnt) %>% dplyr::mutate(pmis_obs = pmis_obs/mrg) %>% # p(l_{unobs}|l_{obs},R=r)
                  mutate_at(vars(pmis_obs), ~ replace(., is.nan(.), 0)) %>%
                  mutate(new.cnt =count*pmis_obs)
    T.Final.pt$new.cnt})
  
  oneiter_C_R = com_case.ns %>% dplyr::select(-t.pijk,) %>% mutate(S1 = count) %>% dplyr::select(-count)
  oneiter_C_R = cbind(oneiter_C_R, mis_samp)
  oneiter_C_R = oneiter_C_R %>% mutate(new.nc = S1+S2+S3+S4)
  
  
  #Maximize:
  Max_1iter_C_R = oneiter_C_R %>% dplyr::select(L1:L3,new.nc) %>% dplyr::mutate(mle.nc = 1)
  Max_1iter_C_R = Max_1iter_C_R %>% dplyr::mutate(across(c(L1,L2,L3), factor, levels = c(1,0)))
  if(interaction_terms ==1){
    fit = (loglm(new.nc ~ L1+L2+L3, data = Max_1iter_C_R, param = T, fitted =T))
  }else if(interaction_terms == 2){
    fit = (loglm(new.nc ~ L1+L2+L3+L1:L2+L1:L3+L2:L3, data = Max_1iter_C_R, param = T, fitted =T))
  }else if(interaction_terms == 3){
    fit = (loglm(new.nc ~ L1+L2+L3+L1:L2+L1:L3+L2:L3 + L1:L2:L3, data = Max_1iter_C_R, param = T, fitted =T))
  }
  
  new.p.ijk = c(aperm(fit$fitted , aperm = c("L1","L2","L3")))/n
  new.Q.t = sum(oneiter_C_R$new.nc*log(new.p.ijk))
  prev.Q.t = sum(oneiter_C_R$new.nc*log(ini.norm_pijk))
  
  
  return(list( new.p.ijk = new.p.ijk,
               new.Q.t = new.Q.t,
               prev.Q.t = prev.Q.t))
}



EM_pijk = function(ini.norm_pijk, com_case.ns, all.pts, data.edit, S1, S2, S3, S4, IPW.res_FI, n, interaction_terms){
  
  old.pijk = ini.norm_pijk
  i = 0
  repeat{
    i = i+1; print(i)
    Step1_EM = Q_fn.ST(ini.norm_pijk = old.pijk, com_case.ns, all.pts, data.edit, 
                       S1, S2, S3, S4, IPW.res_FI, n, interaction_terms)
    new.Q = Step1_EM$new.Q.t
    new.pijk = Step1_EM$new.p.ijk
    # print(new.Q); print(Step1_EM$prev.Q.t); print(new.Q > Step1_EM$prev.Q.t);
    diff = sum(abs(new.pijk-old.pijk))
    if(diff>1e-3){ #change from 1e-8 to 1e-4 #12252022
      # if(!(all((abs(new.pijkl-old.pijk))<1e-8))){
      if (i<200){
      # if((abs((new.Q-old.Q)/old.Q))>10e-8){
      # if(all(abs((new.pijk-old.pijk)/new.pijk)>1e-8)){
      old.pijk = new.pijk
      old.Q = new.Q
    } else{break}
    }else{break}
    }
    
  return(list(new.pijk = new.pijk, new.Q = new.Q))
}




FI_EM.df = function(est.pijk, data.edit, S1, S2, S3,S4, IPW.res_FI, n){
  new.IPW.res_FI = IPW.res_FI; new.IPW.res_FI$est.pijk = est.pijk;
  # data.edit$id = 1:nrow(data.edit)
  imp.S1 = data.edit %>% filter(pt == "S1") %>% dplyr::select(pt,count,L1:L3) %>%mutate(final_prob = 1) %>% dplyr::select(pt, L1:L3, count,final_prob)
  
  ref_prob_pts = lapply(paste0("S",2:4), function(pts){
    T.IPW.res_FI = new.IPW.res_FI %>% dplyr::select(L1:L3,paste0("P_",pts),est.pijk) %>% setNames(c(paste0("L",1:3), "P_pt","est.pijk")) %>% mutate(t.new.prob = P_pt*est.pijk)
    T.mrg.prob = T.IPW.res_FI%>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% 
      group_by(across(contains("L"))) %>% summarise_at(vars(t.new.prob),sum) %>% dplyr::rename(t.mrg.prob = t.new.prob) 
    T.pt.mrg.prob = T.IPW.res_FI %>% left_join(., T.mrg.prob)%>% mutate(final_prob = t.new.prob/t.mrg.prob) %>% mutate(pt = pts)
    T.pt.mrg.prob%>% dplyr::select(pt,L1:L3,final_prob)}) %>% Reduce("rbind",.)
  # names( ref_prob_pts ) = paste0("S",2:4)
  imp_df = lapply(paste0("S", 2:4), function(pts){
    # T.pt.df = data.edit %>% filter(pt == pts) %>% group_by_all() %>% summarise(COUNT = n())
    T.pt.df =  data.edit %>% filter(pt == pts) %>% 
      dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))])%>% # L 선택 로직 (예: keep_vec="c(TRUE,TRUE,FALSE)")
      group_by(across(matches("^L\\d+$")), pt) %>%                  # 남아있는 L열들만 묶기
      summarise(count = sum(count), .groups = "drop") %>%
      left_join(., ref_prob_pts %>% filter(pt == pts)) %>% 
      dplyr::select(pt, L1,L2,L3,count, final_prob)
    T.pt.df}) %>% Reduce("rbind",.) %>% as.data.frame(.)
  return(as.data.frame(rbind(imp.S1, imp_df)))
}

######## Variance Estimation
est.FI.EM.b = function(est.pijk, com_case.ns, all.pts, Orlr_PG, PG,  n, interaction_terms){
  S1 = c(T,T,T); S2 = c(T,T,F); S3 = c(T,F,T); S4 = c(T,F,F);
  
  coef.S2= Orlr_PG$coef.S2; coef.S2[3] = ifelse(length(coef.S2)!=3,0, coef.S2[3])
  coef.S3= Orlr_PG$coef.S3;  coef.S4= Orlr_PG$coef.S4;
  
  # Data generate
  C_R.nc = rmultinom(1,size=n, est.pijk)
  C_R.b = all.pts %>% dplyr::select(L1:L3) %>% mutate(count = C_R.nc)
  
  # Missing Mechanism
  pic_by_mis.b = C_R.b %>%
    mutate(O_S2 = exp(coef.S2[1] + coef.S2[2]*L1 + coef.S2[3]*L2), 
           O_S3 = exp(coef.S3[1] + coef.S3[2]*L1 + coef.S3[3]*L3),
           O_S4 = exp(coef.S4[1] + coef.S4[2]*L1) ) %>% # Add intercept: more realistic model assumption
    mutate(Q_S1 = 1, Q_S2 = O_S2, Q_S3 = O_S3, Q_S4 = O_S4*(Q_S2+Q_S3)) %>%
    mutate(pi_L = 1/(Q_S1+Q_S2+Q_S3+Q_S4)) %>%
    mutate(P_S1 = pi_L, P_S2 = Q_S2*pi_L, P_S3 = Q_S3*pi_L, P_S4 = Q_S4*pi_L) ## create P(R=r|L) This is real probability for creating missingness
  
  nc_by_mis.b = pic_by_mis.b %>% 
                apply(.,1,function(.){
                  T.nc = .["count"] %>% unlist(.)
                  T.probc = .[c("P_S1","P_S2","P_S3","P_S4")] %>% unlist(.)
                  rmultinom(1,size=T.nc, T.probc)}) %>% as.data.frame %>% t(.) %>% as.data.frame(.) %>%
                setNames(paste0("S",1:4))# %>% mutate(nc = S1+S2+S3+S4)
  nc_by_mis.b = cbind(pic_by_mis.b, nc_by_mis.b) %>% as.data.frame(.) %>% dplyr::select(L1:L3, count, S1:S4)
  C_R.b  = nc_by_mis.b %>%  gather(pt,count, S1:S4) %>% dplyr::select(pt, count, L1:L3) #pt, count L1:L3
  
  
  # Create Estimates:
  Orlr.b = OrLr(C_R.b, PG)
  IPW.res.b = IPW(C_R.b, Orlr.b, PG, all.pts)
  IPW.res_FI.b = IPW.res_cal(IPW.res.b) 
  
  
  C_R.cnt.pt.b = C_R.b %>% 
                   mutate(L3 = ifelse(pt == "S2"|pt == "S4",NA,L3),
                          L2 = ifelse(pt == "S3"|pt == "S4",NA,L2)) %>% group_by(pt, L1, L2, L3)  %>%
                   summarise(count = sum(count, na.rm = T), .groups = 'drop')
  ## CC.b 
  
  cc_est = com_case.ns$count/sum(com_case.ns$count);cc_est
  
  OR_CC.df.b =  C_R.b %>% filter(pt == "S1") %>% mutate(t.pijk.b = count/sum(count)) %>% group_by(L2,L3) %>% summarise_at(vars(t.pijk.b), sum) 
  OR_CC.b = (OR_CC.df.b$t.pijk.b[OR_CC.df.b$L2 ==1 & OR_CC.df.b$L3 == 1] *OR_CC.df.b$t.pijk.b[OR_CC.df.b$L2 ==0 & OR_CC.df.b$L3 == 0]) / (OR_CC.df.b$t.pijk.b[OR_CC.df.b$L2 ==1 & OR_CC.df.b$L3 == 0] *OR_CC.df.b$t.pijk.b[OR_CC.df$L2 ==0 & OR_CC.df$L3 == 1])
  
  ## IPW
  
  OR_IPW.b = IPW.res.b  %>% mutate(p_11. = L1*L2*count/pi_L, 
                                      p_10. = L1*(1-L2)*count/pi_L, 
                                      p_01. = (1-L1)*count/pi_L, 
                                      p_00. = (1-L1)*(1-L2)*count/pi_L) %>%
              dplyr::select(contains("p_")) %>% apply(.,2,sum)
  OR_IPW.b = (OR_IPW.b["p_11."]*OR_IPW.b["p_00."])/(OR_IPW.b["p_01."]*OR_IPW.b["p_10."]); OR_IPW.b
  
  
  
  ## FI:
  
  ini.norm_pijk = C_R.cnt.pt.b %>% filter(pt == "S1") %>% mutate(count = count/sum(count))#C_R$S1/sum(C_R$S1) 
  
  ini.norm_pijk = rep(1/8,8)
  
  est.pijk.b = EM_pijk(ini.norm_pijk ,com_case.ns, all.pts,  C_R.cnt.pt.b, S1, S2, S3, S4, IPW.res_FI.b, n, interaction_terms)
  
  FI_EM.df2 = FI_EM.df( est.pijk.b$new.pijk, C_R.b , S1, S2, S3,S4, IPW.res_FI.b, n)
  
  OR_FI = FI_EM.df2 %>% mutate(p_11. = L1*L2*final_prob*count, 
                               p_10. = L1*(1-L2)*final_prob*count, 
                               p_01. = (1-L1)*L2*final_prob*count, 
                               p_00. = (1-L1)*(1-L2)*final_prob*count) %>%
          dplyr::select(contains("p_")) %>% apply(.,2,sum)
  
  OR_FI.b = (OR_FI["p_11."]*OR_FI["p_00."])/(OR_FI["p_01."]*OR_FI["p_10."]); OR_FI.b
  
  
  
  
  # FI.EM.b = FI_EM.df(est.pijk.b, long.pop.b, S1, S2, S3,S4, IPW.res_FI.b, n)
  # FI.EM.sum.b = FI.EM.b %>% mutate(FI.EM.b_Y2 = (L2)*final_prob,
  #                                  FI.EM.b_Y3 =L3*final_prob,
  #                                  FI.EM.b_Y2Y3 =(L2*L3)*final_prob,
  #                                  FI.EM.b_Y1Y2Y3 =(L1*L2*L3)*final_prob,
  #                                  FI.EM.b_Y10Y2Y3= ((1-L1)*L2*L3)*final_prob,
  #                                  FI.EM.b_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*final_prob)%>%
  #   dplyr::select(contains("FI.EM.b")) %>% apply(.,2,sum)
  # FI.EM.mean.b = FI.EM.sum.b/n
  # 
  # FI.EM.means = FI.EM.mean.b
  
  return(c(OR_CC.b, OR_IPW.b, OR_FI.b))
  
  
}



######## Variance Estimation for Pattern graph 1:
est.FI.EM.b = function(est.pijk, com_case.ns, all.pts, Orlr_PG, PG,  n, interaction_terms){
  S1 = c(T,T,T); S2 = c(T,T,F); S3 = c(T,F,T); S4 = c(T,F,F);
  
  coef.S2= Orlr_PG$coef.S2; coef.S2[3] = ifelse(length(coef.S2)!=3,0, coef.S2[3])
  coef.S3= Orlr_PG$coef.S3;  coef.S4= Orlr_PG$coef.S4;
  
  # Data generate
  C_R.nc = rmultinom(1,size=n, est.pijk)
  C_R.b = all.pts %>% dplyr::select(L1:L3) %>% mutate(count = C_R.nc)
  
  # Missing Mechanism
  pic_by_mis.b = C_R.b %>%
    mutate(O_S2 = exp(coef.S2[1] + coef.S2[2]*L1 + coef.S2[3]*L2), 
           O_S3 = exp(coef.S3[1] + coef.S3[2]*L1 + coef.S3[3]*L3),
           O_S4 = exp(coef.S4[1] + coef.S4[2]*L1) ) %>% 
    mutate(Q_S1 = 1, Q_S2 = O_S2, Q_S3 = O_S3, Q_S4 = O_S4) %>% ## Q_S4 has been changed 09112024
    mutate(pi_L = 1/(Q_S1+Q_S2+Q_S3+Q_S4)) %>%
    mutate(P_S1 = pi_L, P_S2 = Q_S2*pi_L, P_S3 = Q_S3*pi_L, P_S4 = Q_S4*pi_L) ## create P(R=r|L) This is real probability for creating missingness
  
  nc_by_mis.b = pic_by_mis.b %>% 
    apply(.,1,function(.){
      T.nc = .["count"] %>% unlist(.)
      T.probc = .[c("P_S1","P_S2","P_S3","P_S4")] %>% unlist(.)
      rmultinom(1,size=T.nc, T.probc)}) %>% as.data.frame %>% t(.) %>% as.data.frame(.) %>%
    setNames(paste0("S",1:4))# %>% mutate(nc = S1+S2+S3+S4)
  nc_by_mis.b = cbind(pic_by_mis.b, nc_by_mis.b) %>% as.data.frame(.) %>% dplyr::select(L1:L3, count, S1:S4)
  C_R.b  = nc_by_mis.b %>%  gather(pt,count, S1:S4) %>% dplyr::select(pt, count, L1:L3) #pt, count L1:L3
  
  
  # Create Estimates:
  Orlr.b = OrLr(C_R.b, PG)
  IPW.res.b = IPW(C_R.b, Orlr.b, PG, all.pts)
  IPW.res_FI.b = IPW.res_cal(IPW.res.b) 
  
  
  C_R.cnt.pt.b = C_R.b %>% 
    mutate(L3 = ifelse(pt == "S2"|pt == "S4",NA,L3),
           L2 = ifelse(pt == "S3"|pt == "S4",NA,L2)) %>% group_by(pt, L1, L2, L3)  %>%
    summarise(count = sum(count, na.rm = T), .groups = 'drop')
  
  ## CC.b 
  com_case.ns.b = C_R.b %>% filter(pt == "S1") %>% dplyr::select(pt, L1:L3, count) %>% mutate(t.pijk = count/sum(count))
  OR_CC.b = fit.odds(com_case.ns.b) 
  

  ## IPW
  OR_IPW.b = IPW.res.b  %>% mutate(count = count/pi_L)
  OR_IPW.b = fit.odds(OR_IPW.b) 
  
  
  ## FI:
  # ini.norm_pijk = C_R.cnt.pt.b %>% filter(pt == "S1") %>% mutate(count = count/sum(count))#C_R$S1/sum(C_R$S1) 
  ini.norm_pijk = rep(1/8,8)
  
  est.pijk.b = EM_pijk(ini.norm_pijk ,com_case.ns = com_case.ns.b, all.pts = all.pts,
                       data.edit = C_R.cnt.pt.b, S1, S2, S3, S4, 
                       IPW.res_FI = IPW.res_FI.b, n, interaction_terms)
  
  FI_EM.df.b = FI_EM.df(est.pijk = est.pijk.b$new.pijk, data.edit = C_R.b, S1, S2, S3,S4, IPW.res_FI = IPW.res_FI.b, n) %>%
              mutate(count = count * final_prob)
  
  OR_FI.b  = fit.odds(FI_EM.df.b)
  # OR_FI = FI_EM.df2 %>% mutate(p_11. = L1*L2*final_prob*count, 
  #                              p_10. = L1*(1-L2)*final_prob*count, 
  #                              p_01. = (1-L1)*L2*final_prob*count, 
  #                              p_00. = (1-L1)*(1-L2)*final_prob*count) %>%
  #   dplyr::select(contains("p_")) %>% apply(.,2,sum)
  
  # OR_FI.b = (OR_FI["p_11."]*OR_FI["p_00."])/(OR_FI["p_01."]*OR_FI["p_10."]); OR_FI.b
  
  
  
  
  # FI.EM.b = FI_EM.df(est.pijk.b, long.pop.b, S1, S2, S3,S4, IPW.res_FI.b, n)
  # FI.EM.sum.b = FI.EM.b %>% mutate(FI.EM.b_Y2 = (L2)*final_prob,
  #                                  FI.EM.b_Y3 =L3*final_prob,
  #                                  FI.EM.b_Y2Y3 =(L2*L3)*final_prob,
  #                                  FI.EM.b_Y1Y2Y3 =(L1*L2*L3)*final_prob,
  #                                  FI.EM.b_Y10Y2Y3= ((1-L1)*L2*L3)*final_prob,
  #                                  FI.EM.b_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*final_prob)%>%
  #   dplyr::select(contains("FI.EM.b")) %>% apply(.,2,sum)
  # FI.EM.mean.b = FI.EM.sum.b/n
  # 
  # FI.EM.means = FI.EM.mean.b
  
  return(c(OR_CC.b, OR_IPW.b, OR_FI.b))
  
  
}
fit.odds= function(df){
   model <- glm(L1 ~ L2 * L3, data = df, 
                family = binomial, weights = count)
   
   summary(model)
   
   # 오즈비 및 CI
   # exp(cbind(OR = coef(model), confint(model)))
   (coef(model))
 }

lin_summary <- function(L, beta, V, label){
  est <- sum(L * beta)
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  z   <- est / se
  p   <- 2 * pnorm(-abs(z))
  c(contrast = label,
    logOR = est,
    SE = se,
    OR = exp(est),
    LCL = exp(est - zcrit*se),
    UCL = exp(est + zcrit*se),
    p = p)
}
