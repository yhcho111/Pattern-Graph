h1 = function(L){L["L2"]}
h2 = function(L){L["L3"]}
h3 = function(L){L["L2"]*L["L3"]}
h4 = function(L){L["L1"]*L["L2"]*L["L3"]}
h5 = function(L){(1-L["L1"])*L["L2"]*L["L3"]}
h6 = function(L){(1-L["L1"])*(1-L["L2"])*(1-L["L3"])}
hlist = list(h1 = h1, h2=h2, h3=h3, h4=h4, h5 = h5, h6 = h6)


RA_SEMI.imp = function(C_R, short.pop, S1, S2, S3, S4, Orlr, n, IPW.res, C_R.cnt.pt, hlist){
  #Step1 : MI to calculate MI for m(lr,r)
  # Ste1-1: prepare the inputs for donor 
  Donor.for.S2 = C_R %>% dplyr::select(L1,L2,L3,S1) %>% group_by(L1,L2) %>% 
                 summarise_at(vars(contains('S')), sum) %>% dplyr::rename(marg.S1 = S1) %>%  #marg.S1 = f(L1,L2|R=111) count
                 left_join(., C_R %>% dplyr::select(L1:L3, S1)) %>% 
                 dplyr::mutate(est.pijk = S1/marg.S1) %>% dplyr::select(-marg.S1,-S1)
  Donor.for.S3 = C_R %>% dplyr::select(L1,L2,L3,S1) %>% group_by(L1,L3) %>% 
                 summarise_at(vars(contains('S')), sum) %>% dplyr::rename(marg.S1 = S1) %>%  #marg.S1 = f(L1,L3|R=111) count
                 left_join(., C_R %>% dplyr::select(L1:L3, S1)) %>% 
                 dplyr::mutate(est.pijk = S1/marg.S1) %>% dplyr::select(-marg.S1,-S1)
  
  #Donor probability
  
  core.set.S2 = C_R %>% dplyr::select(L1, L2, S2) %>% group_by(L1, L2) %>% summarise_at(vars(contains("s")), sum)  # real observed S2
  core.set.S3 = C_R %>% dplyr::select(L1, L3, S3) %>% group_by(L1, L3) %>% summarise_at(vars(contains("s")), sum)  # real observed S3 
  
  n_s2 = sum(unlist(core.set.S2$S2)); n_s3 =  sum(unlist(core.set.S3$S3));
  p_s2_s4 = C_R %>% dplyr::select(L1, S2) %>% group_by(L1) %>% summarise_at(vars(S2), sum) %>% dplyr::mutate(S2 = S2/sum(S2)) %>% dplyr::mutate(n_S2 = n_s2) #p(L1|R=110)
  p_s3_s4 = C_R %>% dplyr::select(L1, S3) %>% group_by(L1) %>% summarise_at(vars(S3), sum) %>% dplyr::mutate(S3 = S3/sum(S3)) %>% dplyr::mutate(n_S3 = n_s3) #p(L1|R=101)
  p_s23_s4 = p_s2_s4 %>% left_join(., p_s3_s4) %>% 
             dplyr::mutate(p_s2 = S2*n_S2, p_s3 = S3*n_S3, p_s2_s3 = p_s2+p_s3) %>% #p_S2 = p(L1|R=110)n_{110}
             dplyr::mutate( p_s2 = p_s2/p_s2_s3, p_s3 = p_s3/p_s2_s3) #p_s2 = P(S=110)  # path probability
  
  Donor.for.S4_S2 = core.set.S2 %>% group_by(L1) %>% summarise_at(vars(contains("S")), sum) %>% dplyr::rename(marg.S2 = S2) %>%  #marg.S2: marginal f(L1|R=110)
                    left_join(., core.set.S2) %>% dplyr::mutate(est.pijk = S2/ marg.S2) %>% dplyr::select(-marg.S2,-S2) # est.pijk = f(L2|L1, R=110)
  Donor.for.S4_S3 = core.set.S3 %>% group_by(L1) %>% summarise_at(vars(contains("S")), sum) %>% dplyr::rename(marg.S3 = S3) %>%  #marg.S3: marginal f(L1|R=101)
                    left_join(., core.set.S3)%>% dplyr::mutate(est.pijk = S3/ marg.S3) %>% dplyr::select(-marg.S3,-S3) # est.pijk = f(L3|L1, R=101)
  
  
  # Step1-2: do MI for each lr
  IMPUTED = lapply(1:nrow(short.pop), function(i){
    # print(i)
    obs.i = ungroup(short.pop[i,])
    pts = unlist(short.pop[i,"pt"])
    obs.val = obs.i[paste0("L",1:3)[eval(parse(text = pts))]] %>% ungroup(.)
    count.i = short.pop[i,"count"] %>% unlist(.)
    
    if (pts == "S2"|pts == "S3"){
      choose.donor = eval(parse(text = paste0("Donor.for.",pts))) %>% plyr::match_df(.,obs.val) %>% ungroup(.) %>%  # choose maching rows from donor
                     # dplyr::select(-paste0("L",1:4)[eval(parse(text = pts))]) %>% 
                     mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric) %>% 
                     dplyr::select(contains("L"), est.pijk) %>% #merge(., obs.val) %>% # merge obs_i with donor_i
                     dplyr::select(L1, L2, L3,est.pijk) 
      
      RA_i = sapply(hlist, function(h){apply(choose.donor%>% dplyr::select(L1:L3),1,h)}) %>% 
              as.data.frame(.) %>% 
              dplyr::mutate(est.pijk = choose.donor$est.pijk) %>% 
              dplyr::mutate(across(contains("h"), ~.*est.pijk)) %>% 
              dplyr::select(-est.pijk) %>% apply(.,2,sum) %>% t(.) %>%as.data.frame(.) %>% dplyr::mutate(cnt = count.i)  
      RA_i 
    } else if (pts == "S4"){
      # Donor Probability
      donor.prob = p_s23_s4 %>% match_df(., obs.val) %>% dplyr::select(L1, p_s2,p_s3) %>% gather(., donors, pijk, p_s2:p_s3) # Donors with probability S4 vs S5
      dn_res = c() # save h1:h7 calculated by each donor
      # 
      for(dn_i in 1:nrow(donor.prob)){
        dn = paste0("S",2:3)[dn_i]
        # Step1: imp S4 with donor S2 (S3): impute L2 (L3)
        imp.st1 = eval(parse(text = paste0("Donor.for.S4_",dn))) %>% match_df(obs.val) %>% ungroup(.) %>% 
          mutate_if(is.factor, as.character) %>%mutate_if(is.character, as.numeric) %>% 
          dplyr::select(contains("L")|contains("est")) %>% dplyr::rename(est.pijk1 = est.pijk)# first entry imputed & Keep probability
        
        #Step 2-0: This returns all possible imputed values with probabilities:
        imp.Final = lapply(1:nrow(imp.st1), function(k){
          obs.i.k = imp.st1[k,] %>% as.data.frame(.) 
          # Step2: imp S4 (S5) with S2 (S3): impute L2
          imp.st2= eval(parse(text = paste0("Donor.for.", dn))) %>% match_df(.,obs.i.k) %>% ungroup(.) %>%  # chose maching rows from donor
                  # dplyr::select(-paste0("L",1:4)[eval(parse(text = dn.st2))])%>% 
                  mutate_if(is.factor, as.character) %>%mutate_if(is.character, as.numeric) %>% 
                  dplyr::select(contains("L"), est.pijk)%>% dplyr::rename(est.pijk2 = est.pijk) %>% merge(., obs.i.k) %>% # merge obs_i wiht donor_i
                  dplyr::mutate(est.pijk = est.pijk1*est.pijk2)%>% dplyr::select(contains("L"), est.pijk)
          imp.st2}) %>% Reduce("rbind",.)  
        imp.Final
        
        RA_i = sapply(hlist, function(h){apply(imp.Final%>% dplyr::select(L1,L2,L3),1,h)}) %>% as.data.frame(.) %>% 
          dplyr::mutate(est.pijk = imp.Final$est.pijk) %>% 
          dplyr::mutate(across(contains("h"), ~.*est.pijk)) %>% 
          dplyr::select(-est.pijk) %>% apply(.,2,sum) %>% t(.) %>%as.data.frame(.) %>% dplyr::mutate(cnt = count.i)  
        res_i = cbind(donor.prob[dn_i,], as.vector(RA_i))
        dn_res = rbind(dn_res,res_i) # save h1:h7 imputed by each path
      }
      
      RA_i = dn_res %>% as.data.frame(.) %>% dplyr::mutate(across(contains("h"), ~.*pijk))  %>% 
        dplyr::select(contains("h")) %>% apply(.,2,sum) %>% t(.)%>% as.data.frame(.)%>% dplyr::mutate(cnt = count.i)
      RA_i 
      # # Donor Probability
      # donor.prob = p_s23_s4 %>% match_df(., obs.val) %>% dplyr::select(L1, p_s2,p_s3) %>% gather(., donors, pijk, p_s2:p_s3) # Donors with probability S4 vs S5
      # dn_res = c() # save h1:h7 calculated by each donor
      # # 
      # for(dn_i in 1:nrow(donor.prob)){
      #   dn = paste0("S",2:3)[dn_i]
      #   donor.prob.i = unlist(donor.prob[dn_i, "pijk"])
      #   # Step1: imp S4 with donor S2 (S3): impute L2 (L3)
      #   imp.st1 = eval(parse(text = paste0("Donor.for.S4_",dn))) %>% match_df(obs.val) %>% ungroup(.) %>% 
      #             mutate_if(is.factor, as.character) %>%mutate_if(is.character, as.numeric) %>% 
      #             dplyr::select(contains("L")|contains("est")) %>% dplyr::rename(est.pijk1 = est.pijk)# first entry imputed & Keep probability
      #         
      #   #Step 2-0: This returns all possible imputed values with probabilities:
      #   imp.Final = lapply(1:nrow(imp.st1), function(k){
      #               obs.i.k = imp.st1[k,] %>% as.data.frame(.) 
      #               # Step2: imp S4 (S5) with S2 (S3): impute L2
      #               imp.st2= eval(parse(text = paste0("Donor.for.", dn))) %>% match_df(.,obs.i.k) %>% ungroup(.) %>%  # chose maching rows from donor
      #                 # dplyr::select(-paste0("L",1:4)[eval(parse(text = dn.st2))])%>% 
      #                 mutate_if(is.factor, as.character) %>%mutate_if(is.character, as.numeric) %>% 
      #                 dplyr::select(contains("L"), est.pijk)%>% dplyr::rename(est.pijk2 = est.pijk) %>% merge(., obs.i.k) %>% # merge obs_i wiht donor_i
      #                 dplyr::mutate(est.pijk = est.pijk1*est.pijk2)%>% dplyr::select(L1,L2,L3, est.pijk)
      #               imp.st2}) %>% Reduce("rbind",.) %>% dplyr::mutate(est.pijk = est.pijk*donor.prob.i) 
      #             
      #   imp.Final
      #   dn_res = rbind(dn_res,as.matrix(imp.Final)) %>% as.data.frame(.) # save h1:h7 imputed by each path
      # }
      # RA_i = sapply(hlist, function(h){apply(dn_res%>% dplyr::select(L1,L2,L3),1,h)}) %>% as.data.frame(.) %>% 
      #   dplyr::mutate(est.pijk = dn_res$est.pijk) %>% 
      #   dplyr::mutate(across(contains("h"), ~.*est.pijk)) %>% 
      #   dplyr::select(-est.pijk) %>% apply(.,2,sum) %>% t(.) %>%as.data.frame(.) %>% dplyr::mutate(cnt = count.i)  
      # # res_i = cbind(donor.prob[dn_i,], as.vector(RA_i))
      # 
      # 
      # RA_i 
    }else if(pts == "S1"){
      RA_i = sapply(hlist, function(h){ apply(obs.i%>% dplyr::select(L1:L3),1,h)}) %>% t(.) %>% as.data.frame(.) %>% 
             setNames(paste0("h",1:6)) %>% dplyr::mutate(cnt = count.i)
      RA_i       
    }
    RA_i }) %>% Reduce("rbind",.)
  
  #   
  #   thetas = sapply(hlist, function(h){
  #     apply(IMP_oneset,1,h) 
  #   }) %>% apply(.,2,sum)
  #   thetas
  # }) %>% Reduce("rbind",.)
  MI_total = IMPUTED%>% as.data.frame(.) %>% dplyr::mutate(across(contains("h"), ~.*cnt)) %>% apply(.,2,sum) %>%t(.) %>% as.data.frame(.) %>% dplyr::select(-cnt)
  RA.mean = MI_total/n
  names(RA.mean) = paste0("RA_", c("Y2","Y3", "Y2Y3", "Y1Y2Y3", "Y10Y2Y3", "Y10Y20Y30"))
  
  
  ############## SEMI Parametric Estimator ##########################
  comp_set = cbind(short.pop,IMPUTED) %>% as.data.frame(.)
  
  m_LR_R = vector("list",4)
  names(m_LR_R) = paste0("S",1:4)
  
  for(s in names(m_LR_R)){
    m_LR_R[[s]] = comp_set%>% filter(pt == s) %>% dplyr::select(L1:L3,h1:h6)
  }
  
  SEMI.mean = sapply(paste0("h",1:6), function(h){
    comb.set =  short.pop %>%
                left_join(., Orlr$uniq.pt.S2 %>% dplyr::select(L1, L2,Sel.Odds) %>% dplyr::rename(.,O2 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
                left_join(., Orlr$uniq.pt.S3 %>% dplyr::select(L1, L3,Sel.Odds) %>% dplyr::rename(.,O3 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
                left_join(., Orlr$uniq.pt.S4 %>% dplyr::select(L1,Sel.Odds) %>% dplyr::rename(.,O4 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
                left_join(., m_LR_R$S1 %>% dplyr::select(L1,L2,L3,h) %>% dplyr::rename(., m1 = h)) %>%
                left_join(., m_LR_R$S2 %>% dplyr::select(L1,L2, h) %>% dplyr::rename(., m2 = h)) %>%
                left_join(., m_LR_R$S3 %>% dplyr::select(L1,L3, h) %>% dplyr::rename(., m3 = h)) %>%
                left_join(.,m_LR_R$S4 %>% dplyr::select(L1, h) %>% dplyr::rename(., m4 = h)) %>%
                left_join(., IPW.res%>% dplyr::select(L1:L3,pi_L)) %>%
                dplyr::mutate(across(O2:m4, replace_na, 0)) %>%
                dplyr::mutate(EIF_S2 = m2*((pt=="S2")-O2*(pt=="S1")), 
                              EIF_S3 = m3*((pt=="S3")-O3*(pt=="S1")), # for every observation, calculate EIF_sum
                              EIF_S4 = (m2*((pt=="S2")-O2*(pt=="S1"))*O4 + 
                                        m3*((pt=="S3")-O3*(pt=="S1"))*O4 +
                                        m4*((pt=="S4")-O4*(pt=="S2"|pt=="S3"))))%>% 
                dplyr::mutate(across(EIF_S2:EIF_S4, replace_na, 0))%>%
                dplyr::mutate(IPW.sum = ifelse(pt=="S1", m1/pi_L,0)) %>%
                dplyr::mutate(EIF = (IPW.sum+EIF_S2+EIF_S3+EIF_S4)) %>% dplyr::mutate(EIF.tot = EIF*count)
    
    sum(comb.set$EIF.tot)/n
  })
  
  
  names(SEMI.mean) = paste0("SEMI_", c("Y2","Y3", "Y2Y3", "Y1Y2Y3", "Y10Y2Y3", "Y10Y20Y30"))
  
  
  
 return(list(RA.mean = RA.mean,
             SEMI.mean = SEMI.mean,
             m_LR_R = m_LR_R ))
   
 
}


