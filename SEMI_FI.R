# semi-imputation using FI 
h1 = function(L){L["L2"]}
h2 = function(L){L["L3"]}
h3 = function(L){L["L2"]*L["L3"]}
h4 = function(L){L["L1"]*L["L2"]*L["L3"]}
h5 = function(L){(1-L["L1"])*L["L2"]*L["L3"]}
h6 = function(L){(1-L["L1"])*(1-L["L2"])*(1-L["L3"])}
hlist = list(h1 = h1, h2=h2, h3=h3, h4=h4, h5 = h5, h6 = h6)

EIF = function(C_R.cnt.pt.EIF,m_LR_R, h_par){
  lapply(1:nrow(C_R.cnt.pt.EIF), function(i){
    row.i = C_R.cnt.pt.EIF[i,]
    L_i = row.i %>% dplyr::select(L1:L3) %>% select_if(~ !any(is.na(.)))
    R_i = row.i$pt
    
    mLS2 = tryCatch(expr = {m_LR_R$S2 %>% filter((L1 == L_i$L1) & (L2 == L_i$L2))%>% ungroup(.)%>% dplyr::select(h_par)%>% unlist(.)},error = function(e){0}) 
    OLS2 = tryCatch(expr = {Orlr$uniq.pt.S2 %>% filter(L1 == L_i$L1, L2 == L_i$L2) %>% dplyr::select(Sel.Odds)%>% unlist(.)},error = function(e){0}) 
    
    mLS3 = tryCatch(expr = {m_LR_R$S3 %>% filter(L1 == L_i$L1, L3 == L_i$L3) %>% ungroup(.)%>% dplyr::select(h_par)%>% unlist(.)},error = function(e){0}) 
    OLS3 = tryCatch(expr = {Orlr$uniq.pt.S3 %>% filter(L1 == L_i$L1 & L3 == L_i$L3) %>% dplyr::select(Sel.Odds)%>% unlist(.)},error = function(e){0}) 
    
    mLS4 = tryCatch(expr = {m_LR_R$S4 %>% filter(L1 == L_i$L1) %>% ungroup(.)%>% dplyr::select(h_par)%>% unlist(.)},error = function(e){0}) 
    OLS4 = tryCatch(expr = {Orlr$uniq.pt.S4 %>% filter(L1 == L_i$L1 ) %>% dplyr::select(Sel.Odds)%>% unlist(.)},error = function(e){0}) 
    
    EIF_110_110 = mLS2*((R_i=="S2") - (OLS2)*(R_i=="S1")) 
    EIF_101_101 = mLS3*((R_i=="S3") - (OLS3)*(R_i=="S1")) 
    EIF_110_100 = mLS2*((R_i=="S2") - (OLS2)*(R_i=="S1"))*OLS4 
    EIF_101_100 = mLS3*((R_i=="S3") - (OLS3)*(R_i=="S1"))*OLS4 
    EIF_100_100 = mLS4*((R_i == "S4")-OLS4*(R_i == "S2"|R_i == "S3"))
    
    EIF_L_R = EIF_110_100+EIF_101_101+EIF_110_100 +EIF_101_100 +EIF_100_100
    
    return(unlist(EIF_L_R))
  }) %>% Reduce("c",.)
  
}

Semi_imp.FI= function(C_R, long.pop, S1, S2, S3, S4, M, Orlr, n, IPW.res_FI, C_R.cnt.pt, est.pijk, IPW.res){
  
  # Step1: do FI to calculate 
  est.C = C_R
  est.C$pijk = est.pijk
  est.C = est.C %>% left_join(., IPW.res_FI) %>% mutate(pijk_PS2 = pijk*P_S2
                                                        , pijk_PS3 = pijk*P_S3
                                                        , pijk_PS4 = pijk*P_S4)
  
  m_LR_R_FI <- lapply(paste0("S",1:4), function(pts){
    if(pts == "S1"){
      cnt.pt = est.C %>% dplyr::select(L1:L3,pts) 
      thetas = sapply(hlist, function(h){
        apply(cnt.pt%>% dplyr::select(-S1),1,h)
      }) 
      cbind(cnt.pt, thetas)
      
    }else if (pts == "S2"){
        cnt.pt = est.C %>% dplyr::select(L1:L2, pts) %>% group_by(L1,L2) %>% summarise_at(vars(contains("S")), funs(sum(.)))
        cond.prob = est.C %>%
        left_join(., est.C %>% group_by(L1,L2)  %>% summarise_at(vars(contains("pijk_PS2")), funs(sum(.)))%>% setNames(c("L1","L2","marg"))) %>%
          mutate(cond = pijk_PS2/marg) %>% dplyr::select(L1:L3, S2, cond)
      theta_L_star = sapply(hlist, function(h){
        hi = apply(cond.prob %>% dplyr::select(-S2, -cond),1,h)
        hi*cond.prob$cond
      }) %>% as.data.frame(.)
      cnt.pt %>% left_join(.,cbind(cond.prob, theta_L_star) %>% as.data.frame(.) %>% group_by(L1,L2) %>% summarise_at(vars(contains("h")), funs(sum(.))))
    }else if (pts == "S3"){
      cnt.pt = est.C %>% dplyr::select(L1,L3, pts) %>% group_by(L1,L3) %>% summarise_at(vars(contains("S")), funs(sum(.)))
      cond.prob = est.C%>% left_join(., est.C %>% group_by(L1,L3) %>% summarise_at(vars(contains("pijk_PS3")), funs(sum(.)))%>% setNames(c("L1","L3","marg"))) %>% 
        group_by(L1,L3) %>% mutate(cond = pijk_PS3/marg) %>% dplyr::select(L1:L3, S3, cond)
      theta_L_star = sapply(hlist, function(h){
        hi = apply(cond.prob %>% dplyr::select(-S3, -cond),1,h)
        hi*cond.prob$cond
      }) %>% as.data.frame(.)
      cnt.pt %>% left_join(.,cbind(cond.prob, theta_L_star) %>% as.data.frame(.) %>% group_by(L1,L3) %>% summarise_at(vars(contains("h")), funs(sum(.))))
      
    }else if (pts == "S4"){
      cnt.pt =  est.C %>% dplyr::select(L1, pts) %>% group_by(L1) %>% summarise_at(vars(contains("S")), funs(sum(.)))
      cond.prob = est.C%>% left_join(., est.C %>% group_by(L1) %>% summarise_at(vars(contains("pijk_PS4")), funs(sum(.)))%>% setNames(c("L1","marg"))) %>% 
        group_by(L1) %>% mutate(cond = pijk_PS4/marg) %>% dplyr::select(L1:L3, S4, cond)
      theta_L_star = sapply(hlist, function(h){
        hi = apply(cond.prob %>% dplyr::select(-S4, -cond),1,h)
        hi*cond.prob$cond
      }) %>% as.data.frame(.)
      cnt.pt %>% left_join(.,cbind(cond.prob, theta_L_star) %>% as.data.frame(.) %>% group_by(L1) %>% summarise_at(vars(contains("h")), funs(sum(.))))
    }
  })
  
  names(m_LR_R_FI) = paste0("S",1:4)
  
  RA_FI= lapply(1:length(m_LR_R_FI), function(s){
    theta.s = m_LR_R_FI[[s]] %>%ungroup(.) %>% dplyr::select(h1:h6)
    cnt.s = m_LR_R_FI[[s]] %>% dplyr::select(contains("S")) %>% unlist(.)
    apply(theta.s*cnt.s,2,sum)
  })%>% Reduce("rbind",.) %>% apply(.,2,sum)
  
  RA_FI.mean = RA_FI/n
  
  # Step 2 EIF:
  C_R.cnt.pt.EIF = C_R.cnt.pt %>% group_by(pt,L1,L2,L3) %>% summarise_at(vars(count), sum) %>% ungroup(.)
  # EIF_hs = sapply(paste0("h",1:6), function(h){EIF(C_R.cnt.pt.EIF, m_LR_R, h)})
  EIF_FI_hs = sapply(paste0("h",1:6), function(h_par){
    lapply(1:nrow(C_R.cnt.pt.EIF), function(i){
      row.i = C_R.cnt.pt.EIF[i,]
      L_i = row.i %>% dplyr::select(L1:L3) %>% select_if(~ !any(is.na(.)))
      R_i = row.i$pt
      
      mLS2 = tryCatch(expr = {m_LR_R_FI$S2 %>% filter((L1 == L_i$L1) & (L2 == L_i$L2))%>% ungroup(.)%>% dplyr::select(h_par)%>% unlist(.)},error = function(e){0}) 
      OLS2 = tryCatch(expr = {Orlr$uniq.pt.S2 %>% filter(L1 == L_i$L1, L2 == L_i$L2) %>% dplyr::select(Sel.Odds)%>% unlist(.)},error = function(e){0}) 
      
      mLS3 = tryCatch(expr = {m_LR_R_FI$S3 %>% filter(L1 == L_i$L1, L3 == L_i$L3) %>% ungroup(.)%>% dplyr::select(h_par)%>% unlist(.)},error = function(e){0}) 
      OLS3 = tryCatch(expr = {Orlr$uniq.pt.S3 %>% filter(L1 == L_i$L1 & L3 == L_i$L3) %>% dplyr::select(Sel.Odds)%>% unlist(.)},error = function(e){0}) 
      
      mLS4 = tryCatch(expr = {m_LR_R_FI$S4 %>% filter(L1 == L_i$L1) %>% ungroup(.)%>% dplyr::select(h_par)%>% unlist(.)},error = function(e){0}) 
      OLS4 = tryCatch(expr = {Orlr$uniq.pt.S4 %>% filter(L1 == L_i$L1 ) %>% dplyr::select(Sel.Odds)%>% unlist(.)},error = function(e){0}) 
      
      EIF_110_110 = mLS2*((R_i=="S2") - (OLS2)*(R_i=="S1")) 
      EIF_101_101 = mLS3*((R_i=="S3") - (OLS3)*(R_i=="S1")) 
      EIF_110_100 = mLS2*((R_i=="S2") - (OLS2)*(R_i=="S1"))*OLS4 
      EIF_101_100 = mLS3*((R_i=="S3") - (OLS3)*(R_i=="S1"))*OLS4 
      EIF_100_100 = mLS4*((R_i == "S4")-OLS4*(R_i == "S2"|R_i == "S3"))
      
      EIF_L_R = EIF_110_100+EIF_101_101+EIF_110_100 +EIF_101_100 +EIF_100_100
      
      return(unlist(EIF_L_R))
    }) %>% Reduce("c",.)
    
  })
  
  EIF_FI_df =  as.data.frame(cbind(C_R.cnt.pt.EIF, EIF_FI_hs)) %>% 
    left_join(., IPW.res %>% dplyr::select(L1,L2,L3,pi_L)) %>%
    mutate_at(vars(pi_L), ~replace_na(.,0))
  
  L_FI_SEMI = sapply(paste0("h",1:6), function(h){
    theta.1d = apply(EIF_FI_df%>% dplyr::select(L1:L3), 1, eval(parse(text = h)))
    L_semi = EIF_FI_df%>% dplyr::select(pt:count, h, pi_L) %>% rename(., "hi" = h) %>%
      mutate(theta.1d = theta.1d) %>% 
      mutate(theta.1d = ifelse(pt == "S1",theta.1d, NA)) %>%
      mutate(semi_i = ifelse(!is.na(theta.1d),(theta.1d/pi_L)+hi,hi), cnt_semi_i = count*semi_i)
    sum(L_semi$cnt_semi_i)
  })
  SEMI_FI.mean = L_FI_SEMI /n
  names(SEMI_FI.mean) = paste0("SEMI.FI_", c("Y2","Y3", "Y2Y3", "Y1Y2Y3", "Y10Y2Y3", "Y10Y20Y30"))
  names(RA_FI.mean)  = paste0("RA.FI_", c("Y2","Y3", "Y2Y3", "Y1Y2Y3", "Y10Y2Y3", "Y10Y20Y30"))
  return(list(RA_FI.mean = RA_FI.mean,
              SEMI_FI.mean = SEMI_FI.mean,
              C_R.cnt.pt.EIF = C_R.cnt.pt.EIF, EIF_FI_hs = EIF_FI_hs, EIF_FI_df = EIF_FI_df, L_FI_SEMI = L_FI_SEMI, m_LR_R_FI = m_LR_R_FI))
  
  
}


