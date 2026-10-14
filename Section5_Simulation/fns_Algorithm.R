## Step 0 : Model Parameter Estimation:
#EM algorithm:
Q_fn= function(ini.norm_pijk, C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI, n){
  
  core_samp = C_R.cnt.pt %>% filter(pt=="S1") %>%  mutate(t.pijk = ini.norm_pijk) # Cate,pt,count, L1:L3, pijk(t), current prob
  T.cate.pijk = core_samp %>% dplyr::select(cate, L1:L3, t.pijk) #
  oneiter_C_R = C_R # L1:L3,nc,S1:S4,cate unreal situation
  
  mis_samp = sapply(paste0("S",2:4), function(pts){ 
    T.IPW.res_FI = IPW.res_FI %>% dplyr::select(L1:L3,paste0("P_",pts)) %>% setNames(c(paste0("L",1:3), "P_pt")) # L1:L3, \hat{P}(R=pts|L)
    T.C_R = C_R %>% dplyr::select(cate, L1:L3, all_of(pts))  # Real information we only know marginal information
    T.mrg_cnt = T.C_R %>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% 
      group_by(across(contains("L"))) %>% summarise_at(vars(pts),sum) # Marginal count by observed variables (Real Situation)
    T.cnt = T.C_R %>% dplyr::select(-pts) %>% left_join(., T.mrg_cnt) %>% 
      left_join(., T.cate.pijk) %>% #f(L1,L2,L3|\hat{\theta}(t))
      left_join(., T.IPW.res_FI) %>% #P(R=pts|L1,L2, L3)
      mutate(pmis_obs = t.pijk*P_pt) # Calculate pmis_obs = P(L1=l1,L2=l2,L3=l3)*P(R=pts|L) numerator
    T.mrg.cnt = T.cnt %>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% # Remove unobserved variables
      group_by(across(contains("L"))) %>% summarise_at(vars(pmis_obs),sum) %>% # Summation wrt unobserved variables
      dplyr::rename(mrg = pmis_obs)  # calculate denominator: \sum_{unobs}P(L1=l1,L2=l2,L3=l3)*P(R=pts|L)
    T.Final.pt =  T.cnt %>% left_join(.,T.mrg.cnt) %>% dplyr::mutate(pmis_obs = pmis_obs/mrg, # p(l_{unobs}|l_{obs},R=r)
                                                                     new.cnt = eval(parse(text = paste0(pts,"*pmis_obs"))))
    T.Final.pt$new.cnt})
  oneiter_C_R[,paste0("S",2:4)] = mis_samp
  oneiter_C_R = oneiter_C_R %>% mutate(new.nc = S1+S2+S3+S4)
  
  
  #Maximize:
  Max_1iter_C_R = oneiter_C_R %>% dplyr::select(cate,L1:L3,new.nc) %>% dplyr::mutate(mle.nc = 1)
  # new.fit = (Iter_Prob_Fit(Max_1iter_C_R)) #sec)
  # new.fit/sum(new.fit)
  Max_1iter_C_R = Max_1iter_C_R %>% dplyr::mutate(across(c(L1,L2,L3), factor, levels = c(1,0)))
  fit = (loglm(new.nc ~ L1+L2+L3+L1:L2+L1:L3+L2:L3, data = Max_1iter_C_R, param = T, fitted =T))
  new.p.ijk = c(aperm(fit$fitted , aperm = c("L1","L2","L3")))/n
  
  new.Q.t = sum(oneiter_C_R$new.nc*log(new.p.ijk))
  prev.Q.t = sum(oneiter_C_R$new.nc*log(ini.norm_pijk))
  
  
  return(list( new.p.ijk = new.p.ijk,
               new.Q.t = new.Q.t,
               prev.Q.t = prev.Q.t))
}



EM_pijk = function(ini.norm_pijk ,C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI,n){
  
  old.pijk = ini.norm_pijk
  i = 0
  repeat{
    i = i+1
    # old.pijk = new.pijk
    Step1_EM = Q_fn(ini.norm_pijk = old.pijk,C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI,n)
    new.Q = Step1_EM$new.Q.t
    new.pijk = Step1_EM$new.p.ijk
    print(new.Q); print(Step1_EM$prev.Q.t); print(new.Q > Step1_EM$prev.Q.t);
    # print(new.pijk);sum(new.pijk)
    diff = sum((new.pijk-old.pijk)^2)
    if(diff> 1e-8){
      old.pijk = new.pijk
      old.Q = new.Q
    } else{break}
  }
  return(new.pijk)
}




## Selection Odds Orlr
#Inverse probability weighting:
OrLr= function(C_R){
  O_S2 = C_R %>% dplyr::select(L1,L2,S1,S2) #%>% mutate(S = S1+S2) 
  uniq.pt  = O_S2 %>% mutate(L1 = factor(L1), L2 = factor(L2)) %>% 
    group_by(L1,L2) %>% summarize_at(., vars(S2,S1),sum) %>% as.data.frame(.) #  L1 L2 S2 S1 
  logis.S2 = glm(cbind(S2,S1)~L1+L2, family="binomial", data=uniq.pt)
  coef.S2 = coef(logis.S2)
  Sel.Odds <- apply((uniq.pt %>% mutate(intercept = 1)  %>% 
                       mutate_if(is.factor, as.character) %>% 
                       mutate_if(is.character, as.numeric) %>% 
                       dplyr::select(intercept, L1,L2) %>% t(.))*coef.S2 ,2,sum) %>% exp(.)
  
  Sel.Odds.df <- uniq.pt %>% mutate(intercept = 1)  %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate_if(is.character, as.numeric) %>% 
    dplyr::select(intercept, L1,L2) 
  
  
  uniq.pt.S2 = uniq.pt %>% mutate(Sel.Odds = Sel.Odds)
  
  
  
  O_S3 = C_R %>% dplyr::select(L1,L3,S1,S3) #%>% mutate(S = S1+S2) 
  uniq.pt  = O_S3 %>% mutate(L1 = factor(L1), L3 = factor(L3)) %>% 
    group_by(L1,L3) %>% summarize_at(., vars(S3,S1),sum) %>% as.data.frame(.) #L1 L3 S3 S1
  logis.S3 = glm(cbind(S3,S1)~L1+L3, family="binomial", data=uniq.pt)
  coef.S3 = coef(logis.S3)
  Sel.Odds <- apply((uniq.pt %>% mutate(intercept = 1)  %>% 
                       mutate_if(is.factor, as.character) %>% 
                       mutate_if(is.character, as.numeric) %>% 
                       dplyr::select(intercept, L1,L3) %>% t(.))*coef.S3 ,2,sum) %>% exp(.)
  
  uniq.pt.S3 = uniq.pt %>% mutate(Sel.Odds = Sel.Odds)
  
  
  O_S4 = C_R %>%  mutate(S = S2+S3) %>% dplyr::select(-S1,-S2,-S3) # S: parenets of S4
  uniq.pt  = O_S4 %>% mutate(L1 = factor(L1)) %>% group_by(L1) %>% summarize_at(., vars(S4,S),sum) %>% as.data.frame(.)
  logis.S4 = glm(cbind(S4,S)~L1, family="binomial", data=uniq.pt)
  coef.S4 = coef(logis.S4)
  Sel.Odds <- apply((uniq.pt %>% mutate(intercept = 1)  %>%
                       mutate_if(is.factor, as.character) %>%
                       mutate_if(is.character, as.numeric) %>%
                       dplyr::select(intercept, L1) %>% t(.))*coef.S4 ,2,sum) %>% exp(.)
  uniq.pt.S4 = uniq.pt %>% mutate(Sel.Odds = Sel.Odds)
  # 
  
  return(list(uniq.pt.S2 = uniq.pt.S2
              ,uniq.pt.S3 = uniq.pt.S3
              , uniq.pt.S4 = uniq.pt.S4
              , coef.S2 = coef.S2
              , coef.S3 = coef.S3
              , coef.S4 = coef.S4)
  )
}

IPW = function(C_R, Orlr){
  core.samp = C_R %>% dplyr::select(L1,L2,L3,S1)
  AL.1 = lapply(1:nrow(core.samp), function(i){
    T.L = core.samp[i,]
    T.L <- T.L %>% dplyr::mutate(Q_S1.L = 1,
                                 Q_S2.L = Orlr$uniq.pt.S2 %>% filter(L1==T.L$L1 & L2==T.L$L2) %>% dplyr::select(Sel.Odds) %>% unlist(.), #O_{110} = Q_{110}
                                 Q_S3.L = Orlr$uniq.pt.S3 %>% filter(L1==T.L$L1 & L3==T.L$L3) %>% dplyr::select(Sel.Odds) %>% unlist(.), #O_{101} = Q_{101}
                                 O_S4.L = Orlr$uniq.pt.S4 %>% filter(L1==T.L$L1) %>% dplyr::select(Sel.Odds) %>% unlist(.), 
                                 Q_S4.L = O_S4.L*(Q_S2.L+Q_S3.L),
                                 pi_L = 1/sum(Q_S1.L+Q_S2.L+Q_S3.L+Q_S4.L))})%>% Reduce("rbind",.)
  return(AL.1)
}

## Fractional Imputation:
FI_imp.df = function(C_R, long.pop, S1, S2, S3,S4, M,IPW.res){
  
  core.set = C_R %>% dplyr::select(L1,L2,L3,S1)
  FI.wt.S4 = IPW.res %>% dplyr::select(L1:L3,Q_S4.L) 
  
  Donor.for.S2 = core.set %>% left_join(.,core.set %>% group_by(L1,L2) %>% summarize_at(vars(S1),sum) %>% dplyr::rename(mrg = S1)) %>% mutate(wt = S1/mrg) %>% dplyr::select(L1:L3,wt) # wt: f(L3|L1,L2,R=111)
  # Donor.for.S3 = core.set %>% left_join(.,core.set %>% group_by(L1,L3) %>% summarize_at(vars(S1),sum) %>%mutate(mrg = S1) %>% dplyr::select(-S1)) %>% mutate(wt = S1/mrg)%>% dplyr::select(L1:L3,wt)
  Donor.for.S4 = core.set %>% left_join(.,core.set %>% group_by(L1) %>% summarize_at(vars(S1),sum) %>% mutate(mrg = S1) %>% dplyr::select(-S1)) %>% mutate(mrg.prob = S1/mrg)%>% dplyr::select(L1:L3,mrg.prob) %>%
    left_join(., FI.wt.S4) %>% mutate(wt = Q_S4.L)
  
  long.pop$id = 1:nrow(long.pop)
  
  
  imp.S1 = long.pop %>% filter(pt == "S1") %>% mutate(wt = 1)
  
  imp.S2 = long.pop %>% filter(pt == "S2") %>% dplyr::select(-L3) %>% left_join(.,Donor.for.S2) %>% dplyr::select(pt,L1,L2,L3,id,wt)
  imp.S3 = long.pop %>% filter(pt == "S3") %>% dplyr::select(-L2) %>% left_join(.,Donor.for.S3)%>% dplyr::select(pt,L1,L2,L3,id,wt)
  
  mis.S4.df = long.pop %>% filter(pt == "S4") %>% dplyr::select(-L2,-L3)
  imp.S4 = lapply(1:nrow(mis.S4.df), function(i){
    mis.i = mis.S4.df[i,]
    sel.donor = Donor.for.S4 %>% filter(L1 == mis.i$L1)
    M.imp.i = sel.donor[rep(1:nrow(sel.donor),rmultinom(1,size = M, sel.donor$mrg.prob)),] %>% dplyr::select(L1:L3, wt)
    M.imp.i$wt = M.imp.i$wt/sum(M.imp.i$wt)
    mis.i %>% left_join(., M.imp.i)
  })%>% Reduce("rbind",.) %>% as.data.frame(.) %>% dplyr::select(pt,L1,L2,L3, id,wt)
  
  FI_dat = rbind(imp.S1, imp.S2, imp.S3, imp.S4)
  
  return( FI_dat)
}


## FI with EM algorithm:
FI_EM.df = function(est.pijk, long.pop, S1, S2, S3,S4, IPW.res_FI, n){
  new.IPW.res_FI  = IPW.res_FI; new.IPW.res_FI$est.pijk = est.pijk;
  long.pop$id = 1:nrow(long.pop)
  imp.S1 = long.pop %>% filter(pt == "S1") %>% dplyr::select(id,pt,L1:L3) %>%mutate(final_prob = 1)
  
  ref_prob_pts = lapply(paste0("S",2:4), function(pts){
    T.IPW.res_FI = new.IPW.res_FI %>% dplyr::select(L1:L3,paste0("P_",pts),est.pijk) %>% setNames(c(paste0("L",1:3), "P_pt","est.pijk")) %>%
      mutate(t.new.prob = P_pt*est.pijk)
    T.mrg.prob = T.IPW.res_FI%>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% 
      group_by(across(contains("L"))) %>% summarise_at(vars(t.new.prob),sum) %>% dplyr::rename(t.mrg.prob = t.new.prob) 
    T.pt.mrg.prob = T.IPW.res_FI %>% left_join(., T.mrg.prob)%>% mutate(final_prob = t.new.prob/t.mrg.prob) %>% mutate(pt = pts)
    T.pt.mrg.prob%>% dplyr::select(pt,L1:L3,final_prob)}) %>% Reduce("rbind",.)
  # names( ref_prob_pts ) = paste0("S",2:4)
  imp_df = lapply(paste0("S", 2:4), function(pts){
    # T.pt.df = long.pop %>% filter(pt == pts) %>% group_by_all() %>% summarise(COUNT = n())
    T.pt.df =  long.pop %>% filter(pt == pts) %>% 
      dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>%
      left_join(., ref_prob_pts %>% filter(pt == pts)) %>% 
      dplyr::select(id, pt, L1,L2,L3,final_prob)
    T.pt.df}) %>% Reduce("rbind",.) %>% as.data.frame(.)
  return(as.data.frame(rbind(imp.S1, imp_df)))
}

h1 = function(L){L["L3"]}
h2 = function(L){L["L1"]*L["L2"]*L["L3"]}
hlist = list(h1 = h1, h2=h2)

## Regression Adjustment
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
    }else if(pts == "S1"){
      RA_i = sapply(hlist, function(h){ apply(obs.i%>% dplyr::select(L1:L3),1,h)}) %>% t(.) %>% as.data.frame(.) %>% 
        setNames(paste0("h",1:6)) %>% dplyr::mutate(cnt = count.i)
      RA_i       
    }
    RA_i }) %>% Reduce("rbind",.)
  
  MI_total = IMPUTED%>% as.data.frame(.) %>% dplyr::mutate(across(contains("h"), ~.*cnt)) %>% apply(.,2,sum) %>%t(.) %>% as.data.frame(.) %>% dplyr::select(-cnt)
  RA.mean = MI_total/n
  names(RA.mean) = paste0("RA_", c("Y3", "Y1Y2Y3"))
  
  
  ############## SEMI Parametric Estimator ##########################
  comp_set = cbind(short.pop,IMPUTED) %>% as.data.frame(.)
  
  m_LR_R = vector("list",4)
  names(m_LR_R) = paste0("S",1:4)
  
  for(s in names(m_LR_R)){
    m_LR_R[[s]] = comp_set%>% filter(pt == s) %>% dplyr::select(L1:L3,h1:h6)
  }
  
  SEMI.mean = sapply(paste0("h",1:2), function(h){
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
  
  
  names(SEMI.mean) = paste0("SEMI_", c("Y3",  "Y1Y2Y3"))
  
  
  
  return(list(RA.mean = RA.mean,
              SEMI.mean = SEMI.mean,
              m_LR_R = m_LR_R ))
  
  
}

obs_lik = function(pars, C_R.cnt.pt){
  t.norm_pijk = c(pars, 1-sum(pars))
  pt.pijk = C_R.cnt.pt %>% filter(pt =="S1") %>% mutate(t.norm_pijk = t.norm_pijk)%>% dplyr::select(cate, t.norm_pijk)
  New.margin.pijk = C_R.cnt.pt %>% left_join(., pt.pijk)
  Marg.prob = lapply(paste0("S", 1:4), function(pts){
    # paste0("L",1:3)[eval(parse(text = pts))]
    marg.prob = New.margin.pijk %>% filter(pt == pts) %>% dplyr::select(c(paste0("L",1:3)[eval(parse(text = pts))],"t.norm_pijk")) %>%
      group_by(across(contains("L"))) %>% summarize_at(.,vars(t.norm_pijk),sum) %>%
      mutate(pt = pts) #%>% ungroup() %>% dplyr::select(pt, t.norm_pijk)
    New.margin.pijk %>% filter(pt == pts) %>% dplyr::select(-t.norm_pijk)%>% left_join(., marg.prob)
  }) %>% Reduce("rbind",.) %>% mutate(cnt.pijk = count*log(t.norm_pijk))
  # Marg.prob %>% filter(pt == "S1") %>% dplyr::select(t.norm_pijk) %>% sum(.)
  print(-sum(Marg.prob$cnt.pijk))
  return(-sum(Marg.prob$cnt.pijk))
}




SEMI_Kim.imp = function(C_R, short.pop, S1, S2, S4, Orlr, n, IPW.res, C_R.cnt.pt, hlist, est.pijk){
  #Step0 : p(R=r|L)
  est.mis.prob = IPW.res %>% mutate(P_S1 = pi_L, P_S2 = pi_L*Q_S2.L, P_S3 = pi_L*Q_S3.L, P_S4 = pi_L*Q_S4.L) %>% dplyr::select(L1:L3,P_S1:P_S4) # P(R=r|L)
  new.C_R = C_R %>% dplyr::mutate(est.pijk = est.pijk) %>% dplyr::select(L1:L3, est.pijk) #f(L|\hat{\theta})
  
  #Step1: f(y1,y2|x)
  p_y1y2_x = new.C_R %>% left_join(., new.C_R %>% group_by(L1) %>% summarize_at(vars(est.pijk), sum) %>% dplyr::rename(marg.est.pijk = est.pijk)) %>% #f(L1)
    dplyr::mutate(est.pijk = est.pijk/marg.est.pijk) #f(L)\f(L1) = f(L2,L3|L1)
  
  E_gL_x =  p_y1y2_x %>%
    cbind(., sapply(hlist, function(h){apply( p_y1y2_x  %>% dplyr::select(L1:L3),1,h)}) %>%  as.data.frame(.)) %>% #cbind(L, f(L2,L3|L1), g(L))
    dplyr::mutate(across(contains("h"), ~.*est.pijk)) %>% #g(L) *f(L2,L3|L1) 
    group_by(L1) %>% summarize_at(vars(contains("h")),sum)#E(g(L)|L1)
  
  #Step1: f(y2|x, y1)
  p_y2_xy1 = new.C_R %>% left_join(., new.C_R %>% group_by(L1, L2) %>% summarize_at(vars(est.pijk), sum) %>% dplyr::rename(marg.est.pijk = est.pijk)) %>% #f(L1,L2)
    dplyr::mutate(est.pijk = est.pijk/marg.est.pijk) #f(L)/f(L1,L2) = f(L3|L1,L2)
  
  E_gL_xy1 =  p_y2_xy1 %>%
    cbind(., sapply(hlist, function(h){apply( p_y2_xy1 %>% dplyr::select(L1:L3),1,h)}) %>%  as.data.frame(.)) %>% # cbind(L, f(L3|L1,L2), g(L))
    dplyr::mutate(across(contains("h"), ~.*est.pijk)) %>% # f(L3|L1,L2)*g(L)
    group_by(L1, L2) %>% summarize_at(vars(contains("h")),sum) #E(g(L)|L1,L2)
  
  #Step1: f(y1|x, y2)
  p_y1_xy2 = new.C_R %>% left_join(., new.C_R %>% group_by(L1, L3) %>% summarize_at(vars(est.pijk), sum) %>% dplyr::rename(marg.est.pijk = est.pijk)) %>% # f(L1, L3)
    dplyr::mutate(est.pijk = est.pijk/marg.est.pijk) #f(L2|L1,L3)
  
  E_gL_xy2 =  p_y1_xy2 %>%
    cbind(., sapply(hlist, function(h){apply( p_y1_xy2 %>% dplyr::select(L1:L3),1,h)}) %>%  as.data.frame(.)) %>%
    dplyr::mutate(across(contains("h"), ~.*est.pijk)) %>%
    group_by(L1, L3) %>% summarize_at(vars(contains("h")),sum) # E[g(L)|L1, L3]
  #Step E(g_i|x,y1,y2)=gi          
  E_gL_xy1y2= cbind(C_R%>% dplyr::select(L1:L3) , sapply(hlist, function(h){ apply(C_R%>% dplyr::select(L1:L3),1,h)})) %>% as.data.frame(.) %>% as_tibble(.) #
  ############## SEMI Parametric Estimator ##########################
  # comp_set = cbind(short.pop,IMPUTED) %>% as.data.frame(.)
  
  m_LR_R = vector("list",4)
  names(m_LR_R) = paste0("b", 1:4)
  m_LR_R[[1]] = E_gL_x
  m_LR_R[[2]] = E_gL_xy1
  m_LR_R[[3]] = E_gL_xy2
  m_LR_R[[4]] = E_gL_xy1y2
  
  ############  Calculate E(O_r|\bar{r})
  
  # E(P(R=111|x,Y2, y1)|X=x, Y1=y1)
  E_p11_xy1 = p_y2_xy1 %>% left_join(.,est.mis.prob %>% dplyr::select(L1:L3, P_S1) ) %>% 
    dplyr::mutate(e_p11_xy1 = est.pijk*P_S1) %>% #P(R=111|x,y1,y2)*f(y2|x,y1)
    group_by(L1,L2) %>% summarize_at(vars(e_p11_xy1), sum) #\sum_{y_{2}}P(R=111|x,y1,y2)*f(y2|x,y1)
  
  # p_y2_xy1 %>% left_join(.,est.mis.prob %>% dplyr::select(L1:L3, P_S1)) %>% 
  # left_join(., Orlr$uniq.pt.S2 %>% dplyr::select(L1, L2, Sel.Odds) %>% dplyr::rename(.,O2 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
  # left_join(., Orlr$uniq.pt.S3 %>% dplyr::select(L1, L3, Sel.Odds) %>% dplyr::rename(.,O3 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
  # left_join(., Orlr$uniq.pt.S4 %>% dplyr::select(L1,Sel.Odds) %>% dplyr::rename(.,O4 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>% 
  # mutate(p11_1 = (1/(1+O2*(1+O4)+O3*(1+O4)))*est.pijk) %>%
  # group_by(L1, L2) %>% summarize_at(vars(p11_1),sum)
  
  # E(P(R=111|x,Y1,y2)|X=x, Y2=y2)
  E_p11_xy2 = p_y1_xy2 %>%  left_join(.,est.mis.prob %>% dplyr::select(L1:L3, P_S1) ) %>%
    dplyr::mutate(e_p11_xy2 = est.pijk*P_S1) %>% # P(R=111|x,y1,y3)*f(y1|x,y2)
    group_by(L1,L3) %>% summarize_at(vars(e_p11_xy2), sum) #\sum_{y1}P(R=111|x,y1,y3)*f(y1|x,y2)
  
  
  
  SEMI.mean = sapply(paste0("h",1:6), function(h){
    comb.set = short.pop %>%
      left_join(., Orlr$uniq.pt.S2 %>% dplyr::select(L1, L2, Sel.Odds) %>% dplyr::rename(.,O2 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
      left_join(., Orlr$uniq.pt.S3 %>% dplyr::select(L1, L3, Sel.Odds) %>% dplyr::rename(.,O3 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
      left_join(., Orlr$uniq.pt.S4 %>% dplyr::select(L1,Sel.Odds) %>% dplyr::rename(.,O4 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
      left_join(., m_LR_R$b1 %>% dplyr::select(L1,h) %>% dplyr::rename(., b1 = h)) %>%
      left_join(., m_LR_R$b2 %>% dplyr::select(L1,L2,h) %>% dplyr::rename(., b2 = h)) %>%
      left_join(., m_LR_R$b3 %>% dplyr::select(L1,L3,h) %>% dplyr::rename(., b3 = h)) %>%
      left_join(., m_LR_R$b4 %>% dplyr::select(L1,L2,L3, h) %>% dplyr::rename(., b4 = h)) %>%
      left_join(., est.mis.prob %>% dplyr::select(L1:L3, P_S1)) %>%
      # dplyr::mutate(across(O2:b4, replace_na, 0)) %>%
      left_join(., E_p11_xy1) %>% left_join(., E_p11_xy2) %>%
      dplyr::mutate(EIF_S2 = ((pt=="S2")-O2*(pt=="S1"))*((b1*(1+O4))+((b2-b1)/(e_p11_xy1*(1+O2)))), # for every observation, calculate EIF_sum changed from e_p11_xy1 to e_p11_xy1*(1+O2)
                    EIF_S3 = ((pt=="S3")-O3*(pt=="S1"))*((b1*(1+O4))+((b3-b1)/(e_p11_xy2*(1+O3)))), #((I(R=101))-O3*(I(R=111)))*((m_x*(1+O4))+((m_xy2-m_x)/(pi_11*(1+O3))))
                    EIF_S4 = ((pt=="S4")-O4*(pt=="S2"|pt=="S3"))*b1)%>% 
      dplyr::mutate(across(EIF_S2:EIF_S4, replace_na, 0)) %>%
      dplyr::mutate(IPW.sum = ifelse(pt=="S1", b4/P_S1,0)) %>%
      dplyr::mutate(EIF = (IPW.sum+EIF_S2+EIF_S3+EIF_S4)) %>% mutate(EIF.tot = EIF*count)
    
    sum(comb.set$EIF.tot)/n
  })
  
  SEMI.mean = unlist(SEMI.mean)
  names(SEMI.mean) = paste0("SEMI2_", c("Y2","Y3", "Y2Y3", "Y1Y2Y3", "Y10Y2Y3", "Y10Y20Y30"))
  
  SEMI.mean
  
  
  ####### created 04182023
  # E(1/P(R=111|x,Y2, y1)|X=x, Y1=y1)
  E_p11_inv_xy1 = p_y2_xy1 %>% left_join(., est.mis.prob %>% dplyr::select(L1:L3, P_S1) ) %>%
    dplyr::mutate(e_p11_inv_xy1 = est.pijk/P_S1) %>% # f(y2|x,y1)/P(R=111|x, y1, y2)
    group_by(L1,L2) %>% summarize_at(vars(e_p11_inv_xy1), sum) #\sum_{y_{2}} f(y2|x,y1)/P(R=111|x, y1, y2) = E[1/P(R=111|x, y1, y2)|x,y1]
  
  # E(1/P(R=111|x,Y1,y2)|X=x, Y2=y2)
  E_p11_inv_xy2 = p_y1_xy2 %>%  left_join(.,est.mis.prob %>% dplyr::select(L1:L3, P_S1) ) %>%
    dplyr::mutate(e_p11_inv_xy2 = est.pijk/P_S1) %>% # f(y1|x,y2)/P(R=111|x, y1, y2)
    group_by(L1,L3) %>% summarize_at(vars(e_p11_inv_xy2), sum) #\sum_{y_{1}}f(y1|x,y2)/P(R=111|x, y1, y2) = E[1/P(R=111|x,y1,y2)|x,y2]
  
  SEMI.mean2 = sapply(paste0("h",1:6), function(h){
    comb.set2 = short.pop %>%
      left_join(., Orlr$uniq.pt.S2 %>% dplyr::select(L1, L2, Sel.Odds) %>% dplyr::rename(.,O2 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
      left_join(., Orlr$uniq.pt.S3 %>% dplyr::select(L1, L3, Sel.Odds) %>% dplyr::rename(.,O3 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
      left_join(., Orlr$uniq.pt.S4 %>% dplyr::select(L1,Sel.Odds) %>% dplyr::rename(.,O4 = Sel.Odds) %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric)) %>%
      left_join(., m_LR_R$b1 %>% dplyr::select(L1,h) %>% dplyr::rename(., b1 = h)) %>%
      left_join(., m_LR_R$b2 %>% dplyr::select(L1,L2,h) %>% dplyr::rename(., b2 = h)) %>%
      left_join(., m_LR_R$b3 %>% dplyr::select(L1,L3,h) %>% dplyr::rename(., b3 = h)) %>%
      left_join(., m_LR_R$b4 %>% dplyr::select(L1,L2,L3, h) %>% dplyr::rename(., b4 = h)) %>%
      left_join(., est.mis.prob %>% dplyr::select(L1:L3, P_S1)) %>%
      left_join(., E_p11_inv_xy1) %>% left_join(., E_p11_inv_xy2) %>%
      dplyr::mutate(EIF1 = b1, EIF2 =(b2-b1)*e_p11_inv_xy1/(1+O2), EIF3 = (b3-b1)*e_p11_inv_xy2/(1+O3), EIF4 = (b4-b3-b2+b1)/P_S1)%>%
      # left_join(., E_p11_xy1) %>% left_join(., E_p11_xy2) %>%
      # dplyr::mutate(EIF1 = b1, EIF2 =(b2-b1)/((1+O2)*e_p11_xy1), EIF3 = (b3-b1)*e_p11_xy2/((1+O3)*e_p11_xy2), EIF4 = (b4-b3-b2+b1)/P_S1)%>%
      dplyr::mutate(across(EIF1:EIF4, replace_na, 0))%>%
      # dplyr::mutate(IPW.sum = ifelse(pt=="S1", b4/P_S1,0)) %>%
      dplyr::mutate(EIF = (EIF1+EIF2+EIF3+EIF4)) %>% mutate(EIF.tot = EIF*count)
    
    sum(comb.set2$EIF.tot)/n
  })
  
  SEMI.mean2 = unlist(SEMI.mean2)
  names(SEMI.mean2) = paste0("SEMI3_", c("Y2","Y3", "Y2Y3", "Y1Y2Y3", "Y10Y2Y3", "Y10Y20Y30"))
  SEMI.mean
  SEMI.mean2
  
  
  
  return(list(SEMI.mean = SEMI.mean, SEMI.mean2 = SEMI.mean2))
  
  
}


############# For the saturated model:

Q_fn.ST= function(ini.norm_pijk, C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI, n){
  
  core_samp = C_R.cnt.pt %>% filter(pt=="S1") %>%  mutate(t.pijk = ini.norm_pijk) # Cate,pt,count, L1:L3, pijk(t), current prob
  T.cate.pijk = core_samp %>% dplyr::select(cate, L1:L3, t.pijk) #
  oneiter_C_R = C_R # L1:L3,nc,S1:S4,cate unreal situation
  
  mis_samp = sapply(paste0("S",2:4), function(pts){ 
    T.IPW.res_FI = IPW.res_FI %>% dplyr::select(L1:L3,paste0("P_",pts)) %>% setNames(c(paste0("L",1:3), "P_pt")) # L1:L3, \hat{P}(R=pts|L)
    T.C_R = C_R %>% dplyr::select(cate, L1:L3, all_of(pts))  # Real information we only know marginal information
    T.mrg_cnt = T.C_R %>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% 
      group_by(across(contains("L"))) %>% summarise_at(vars(pts),sum) # Marginal count by observed variables (Real Situation)
    T.cnt = T.C_R %>% dplyr::select(-pts) %>% left_join(., T.mrg_cnt) %>% 
      left_join(., T.cate.pijk) %>% #f(L1,L2,L3|\hat{\theta}(t))
      left_join(., T.IPW.res_FI) %>% #P(R=pts|L1,L2, L3)
      mutate(pmis_obs = t.pijk*P_pt) # Calculate pmis_obs = P(L1=l1,L2=l2,L3=l3)*P(R=pts|L) numerator
    T.mrg.cnt = T.cnt %>% dplyr::select(-paste0("L",1:3)[!eval(parse(text = pts))]) %>% # Remove unobserved variables
      group_by(across(contains("L"))) %>% summarise_at(vars(pmis_obs),sum) %>% # Summation wrt unobserved variables
      dplyr::rename(mrg = pmis_obs)  # calculate denominator: \sum_{unobs}P(L1=l1,L2=l2,L3=l3)*P(R=pts|L)
    T.Final.pt =  T.cnt %>% left_join(.,T.mrg.cnt) %>% dplyr::mutate(pmis_obs = pmis_obs/mrg, # p(l_{unobs}|l_{obs},R=r)
                                                                     new.cnt = eval(parse(text = paste0(pts,"*pmis_obs"))))
    T.Final.pt$new.cnt})
  oneiter_C_R[,paste0("S",2:4)] = mis_samp
  oneiter_C_R = oneiter_C_R %>% mutate(new.nc = S1+S2+S3+S4)
  new.p.ijk = oneiter_C_R$new.nc/n
  new.Q.t = sum(oneiter_C_R$new.nc*log(new.p.ijk))
  prev.Q.t = sum(oneiter_C_R$new.nc*log(ini.norm_pijk))
  return(list( new.p.ijk = new.p.ijk,
               new.Q.t = new.Q.t,
               prev.Q.t = prev.Q.t))
}



EM_pijk.ST = function(ini.norm_pijk ,C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI,n){
  
  old.pijk = ini.norm_pijk
  i = 0
  repeat{
    i = i+1
    # old.pijk = new.pijk
    Step1_EM = Q_fn.ST(ini.norm_pijk = old.pijk,C_R, C_R.cnt.pt, S1, S2, S3, S4, IPW.res_FI,n)
    new.Q = Step1_EM$new.Q.t
    new.pijk = Step1_EM$new.p.ijk
    print(new.Q); print(Step1_EM$prev.Q.t); print(new.Q > Step1_EM$prev.Q.t);
    # print(new.pijk);sum(new.pijk)
    diff = sum((new.pijk-old.pijk)^2)
    if(diff> 1e-8){
      # if(!(all((abs(new.pijk-old.pijk))<1e-8))){
      # if (i<200){
      # if((abs((new.Q-old.Q)/old.Q))>10e-8){
      # if(all(abs((new.pijk-old.pijk)/new.pijk)>1e-8)){
      old.pijk = new.pijk
      old.Q = new.Q
    } else{break}
  }
  return(new.pijk)
}


######## Variance Estimation
est.FI.EM.b = function(est.pijk, coef.S2, coef.S3, coef.S4, C_R, n, MODEL){
  S1 = c(T,T,T); S2 = c(T,T,F); S3 = c(T,F,T); S4 = c(T,F,F);
  
  
  # Data generate
  C_R.nc = rmultinom(1,size=n, est.pijk)
  C_R.b = C_R %>% dplyr::select(L1:L3, cate) %>% mutate(nc = C_R.nc)
  
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
      T.nc = .["nc"] %>% unlist(.)
      T.probc = .[c("P_S1","P_S2","P_S3","P_S4")] %>% unlist(.)
      rmultinom(1,size=T.nc, T.probc)}) %>% as.data.frame %>% t(.) %>% as.data.frame(.) %>%
    setNames(paste0("S",1:4))# %>% mutate(nc = S1+S2+S3+S4)
  nc_by_mis.b = cbind(pic_by_mis.b, nc_by_mis.b) %>% as.data.frame(.) %>% dplyr::select(L1:L3,nc, S1:S4)
  C_R.b  = nc_by_mis.b %>% left_join(., C_R.b%>%dplyr::select(L1:L3, cate))
  
  # Create Estimates:
  Orlr.b = OrLr(C_R.b)
  IPW.res.b = IPW(C_R.b, Orlr.b)
  
  IPW.res_FI.b = IPW.res.b %>% 
    mutate(P_S2 = Q_S2.L*pi_L, P_S3 = Q_S3.L*pi_L, P_S4 = Q_S4.L*pi_L, P_S1 = pi_L, lat = P_S1 + P_S2 + P_S3 + P_S4) %>%
    dplyr::select(L1:L3,P_S1, P_S2, P_S3, P_S4)
  C_R.pt.b = C_R.b %>% dplyr::select(cate, L1, L2, L3)
  C_R.cnt.b = C_R.b %>% dplyr::select(cate, S1:S4) %>% gather(., pt, count, S1:S4)
  C_R.cnt.pt.b = C_R.cnt.b %>% left_join(., C_R.pt.b) %>% 
    mutate(L3 = ifelse(pt == "S2"|pt == "S4",NA,L3),
           L2 = ifelse(pt == "S3"|pt == "S4",NA,L2))
  long.pop.b = C_R.cnt.pt.b[rep(1:nrow(C_R.cnt.pt.b), times = C_R.cnt.pt.b$count),] %>% dplyr::select(pt, L1:L3)
  
  
  ini.norm_pijk = rep(1/8,8)#C_R$S1/sum(C_R$S1) 
  
  if(MODEL!="ST"){
    est.pijk.b = EM_pijk(ini.norm_pijk ,C_R.b, C_R.cnt.pt.b, S1, S2, S3, S4, IPW.res_FI.b,n)
    
    FI.EM.b = FI_EM.df(est.pijk.b, long.pop.b, S1, S2, S3,S4, IPW.res_FI.b, n)
    FI.EM.sum.b = FI.EM.b %>% mutate(FI.EM.b_Y2 = (L2)*final_prob,
                                     FI.EM.b_Y3 =L3*final_prob,
                                     FI.EM.b_Y2Y3 =(L2*L3)*final_prob,
                                     FI.EM.b_Y1Y2Y3 =(L1*L2*L3)*final_prob,
                                     FI.EM.b_Y10Y2Y3= ((1-L1)*L2*L3)*final_prob,
                                     FI.EM.b_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*final_prob)%>%
      dplyr::select(contains("FI.EM.b")) %>% apply(.,2,sum)
    FI.EM.mean.b = FI.EM.sum.b/n
    
    FI.EM.means = FI.EM.mean.b
  }else{
    est.pijk.ST.b = EM_pijk.ST(ini.norm_pijk ,C_R.b, C_R.cnt.pt.b, S1, S2, S3, S4, IPW.res_FI.b,n)
    FI.EM.ST.b = FI_EM.df(est.pijk.ST.b, long.pop.b, S1, S2, S3,S4, IPW.res_FI.b, n)
    FI.EM.sum.ST.b = FI.EM.ST.b %>% mutate(FI.EM.ST.b_Y2 = (L2)*final_prob,
                                           FI.EM.ST.b_Y3 =L3*final_prob,
                                           FI.EM.ST.b_Y2Y3 =(L2*L3)*final_prob,
                                           FI.EM.ST.b_Y1Y2Y3 =(L1*L2*L3)*final_prob,
                                           FI.EM.ST.b_Y10Y2Y3= ((1-L1)*L2*L3)*final_prob,
                                           FI.EM.ST.b_Y10Y20Y30= ((1-L1)*(1-L2)*(1-L3))*final_prob)%>%
      dplyr::select(contains("FI.EM.ST.b")) %>% apply(.,2,sum)
    FI.EM.ST.mean.b = FI.EM.sum.ST.b/n
    FI.EM.means = FI.EM.ST.mean.b
  }
  
  return(FI.EM.means)
  
}

