## Step 0 : Model Parameter Estimation:
#EM algorithm:
EM_initials = function(C_R){
  
  core.S1 = C_R %>% filter(pts == "S1")
  ini.m1 = mean(core.S1$L1);
  ini.m2 = mean(core.S1$L2);
  ini.s11 = var(core.S1$L1);ini.s22 =var(core.S1$L2);
  ini.s12 = mean((core.S1$L1)*(core.S1$L2)) - mean(core.S1$L1)*mean(core.S1$L2)
  ini.rho <- ini.s12/sqrt(ini.s11*ini.s22)
  inis = paste0("ini.",c("m1","m2","s11","s22","rho")) %>% sapply(., function(.){eval(parse(text = .))})
  return(inis)
}


# Define the log-likelihood function
loglik <- function(par, imp.SS) {
  m1 <- par["ini.m1"]; m2 <- par["ini.m2"];
  sigma1 <- exp(par["ini.s11"]); sigma2 <- exp(par["ini.s22"]); rho <- tanh(par["ini.rho"]);
  
  log1 <- -0.5*log((1-rho^2)*sigma1^2*sigma2^2) 
  
  log2s <- imp.SS %>% dplyr::mutate(s1 = (L1-m1)^2/(sigma1^2), 
                                    s2 = (L2-m2)^2/(sigma2^2),
                                    s12 = -2*rho*(L1-m1)*(L2-m2)/(sigma1*sigma2),
                                    log2 = (log1 -(0.5/(1-rho^2))*(s1+s2+s12))*final_wt.t)
  ll <- -sum(log2s$log2)
  
  
  return(ll)
}

## Selection Odds Orlr
#Inverse probability weighting:
OrLr= function(C_R){
  O_S2 = C_R %>%dplyr::select(ID,L1,S1,S2) %>%filter(S1==1|S2==1) %>% mutate(S1 = factor(S1), S2 = factor(S2))
  logis.S2 = glm(S2~L1, family="binomial", data=O_S2)
  coef.S2 = coef(logis.S2)
  Sel.Odds <- apply((O_S2 %>% mutate(intercept = 1)  %>% 
                       dplyr::select(intercept, L1) %>% t(.))*coef.S2 ,2,sum) %>% exp(.)
  # plot(Sel.Odds/(1+Sel.Odds), fitted.values(logis.S2)); abline(0,1,col="red")
  uniq.pt.S2 = O_S2 %>% mutate(Sel.Odds = Sel.Odds)
  
  
  O_S3 = C_R %>%dplyr::select(ID,L1,L2,S1,S3) %>% filter(S1==1|S3==1) %>% mutate(S1 = factor(S1), S3 = factor(S3))
  logis.S3 = glm(S3~L2, family="binomial", data=O_S3)
  coef.S3 = coef(logis.S3)
  Sel.Odds <- apply((O_S3 %>% mutate(intercept = 1)  %>% 
                       mutate_if(is.factor, as.character) %>% 
                       mutate_if(is.character, as.numeric) %>% 
                       dplyr::select(intercept, L2) %>% t(.))*coef.S3 ,2,sum) %>% exp(.)
  
  uniq.pt.S3 = O_S3 %>% mutate(Sel.Odds = Sel.Odds)
  
  O_S4 = C_R %>% filter(S1!=1)%>% dplyr::select(ID,S2,S3,S4) %>% mutate(S4_PA = S3+S2) %>% 
    dplyr::select(-S3,-S2)  %>% mutate(S4 = factor(S4), S4_PA = factor(S4_PA))
  logis.S4 = glm(S4~1, family="binomial", data=O_S4)
  coef.S4 = coef(logis.S4)
  Sel.Odds <- apply((O_S4 %>% mutate(intercept = 1)  %>% 
                       mutate_if(is.factor, as.character) %>% 
                       mutate_if(is.character, as.numeric) %>% 
                       dplyr::select(intercept) %>% t(.))*coef.S4 ,2,sum) %>% exp(.)
  uniq.pt.S4 = O_S4 %>% mutate(Sel.Odds = Sel.Odds)
  
  return(list(uniq.pt.S2 = uniq.pt.S2
              , uniq.pt.S3 = uniq.pt.S3 
              , uniq.pt.S4 = uniq.pt.S4
              , coef.S2 = coef.S2
              , coef.S3 = coef.S3
              , coef.S4 = coef.S4))
}

#--------------------------------------------------------
# IPW
#--------------------------------------------------------
IPW = function(C_R, Orlr){
  core.samp = C_R %>% dplyr::select(L1,L2,S1,ID) %>% filter(S1==1)
  AL.1 = lapply(1:nrow(core.samp), function(i){
    T.L = core.samp[i,]
    T.L <- T.L %>% mutate(Q_S1.L = 1,
                          Q_S2.L = Orlr$uniq.pt.S2 %>% filter(L1==T.L$L1 ) %>% dplyr::select(Sel.Odds) %>% unlist(.),
                          Q_S3.L = Orlr$uniq.pt.S3 %>% filter(L2==T.L$L2 ) %>% dplyr::select(Sel.Odds) %>% unlist(.),
                          Q_S4.L = Orlr$uniq.pt.S4$Sel.Odds[1] %>% unlist(.), 
                          Q_S4.L = Q_S4.L *(Q_S2.L+Q_S3.L),
                          pi_L = 1/sum(Q_S1.L+Q_S2.L+Q_S3.L+Q_S4.L))})%>% Reduce("rbind",.)
  return(AL.1)
}


#--------------------------------------------------------
# FI
#--------------------------------------------------------
impY_hm_piL_kim = function(ini.norm_pijk, C_R, coef.S2, coef.S3, coef.S4, S1, S2, S3, S4, M){
  
  # Generate imputed values only one!!!!!!
  
  #Simulate theta for the proposal density, where f_{proposal} = K^{-1}\sum_{k=1}^{K}
  ini_trans = ini.norm_pijk
  ini_trans["ini.s11"] = log(sqrt(abs(ini_trans["ini.s11"])));
  ini_trans["ini.s22"] = log(sqrt(abs(ini_trans["ini.s22"])));
  ini_trans["ini.rho"] = atanh(ini_trans["ini.rho"])
  imp.SS = C_R %>% filter(S1==1) %>% mutate(final_wt.t = 1)
  news_optim <- optim(ini_trans, loglik, hessian = T, imp.SS = imp.SS, method = "BFGS")
  news <- news_optim$par
  news["ini.s11"] = exp(news["ini.s11"])^2
  news["ini.s22"] = exp(news["ini.s22"])^2
  news["ini.rho"] = tanh(news["ini.rho"])
  news.hess <- news_optim$hessian
  news.obsinfo <- solve(news.hess)
  
  theta.gen.M <- mvrnorm(M, mu = news, Sigma = news.obsinfo) %>% as.data.frame(.) %>% mutate(sel.row = row_number()) %>% 
    mutate(ini.rho = replace(ini.rho, ini.rho>1, max(ini.rho[ini.rho<1])),
           ini.s11 = replace(ini.s11, ini.s11<0, min(ini.s11[ini.s11>0])),
           ini.s22 = replace(ini.s22, ini.s22<0, min(ini.s22[ini.s22>0])))
  
  
  # Step 2: Simulate Each Variable:
  # complete S1: R = (T, T)
  imp.S1 = C_R %>% filter(pts == "S1") %>% dplyr::select(ID, L1:L2, pts)  %>% dplyr::select(ID, L1:L2, pts)
  
  # S2: R = (T, F), impute L2 from the proposal den:
  mis.S2 = C_R %>% filter(pts == "S2") 
  
  imp.S2 = sapply(1:M, function(m){
    sel.row.theta = sample(1:nrow(theta.gen.M), nrow(mis.S2), replace = T) # SAMPLE from the Proposal density:
    samp.iter = mis.S2 %>% dplyr::select(ID, L1) %>% dplyr::mutate(sel.row = sel.row.theta) %>% 
      left_join(., theta.gen.M) %>% 
      mutate(cond.m = ini.m2 + sqrt(ini.s22/ini.s11)*ini.rho*(L1-ini.m1), #E(L2|L1;theta(k)) 
             cond.v = (1-ini.rho^2)*ini.s22)                            # var(L2|L1;theta(k))
    rnorm(nrow(mis.S2), mean = samp.iter$cond.m, sd = sqrt(samp.iter$cond.v))}) %>% 
    as.data.frame(.) %>% setNames(1:M) %>% mutate(ID = mis.S2$ID) %>% gather(., key = "iter", value= "L2", '1':'100') %>% 
    left_join(., mis.S2 %>% dplyr::select(ID, L1)) %>% 
    mutate(Q_S2 = exp(coef.S2[1]+(coef.S2[2]*L1)), 
           Q_S3 = exp(coef.S3[1]+(coef.S3[2]*L2)),
           Q_S4 = exp(coef.S4[1])*(Q_S3+Q_S2), 
           pi_L = Q_S2 /(1+Q_S2+Q_S3+Q_S4)) #P(R=10|L1,L2*)
  
  imp.S2.hm = lapply(1:nrow(imp.S2), function(i){ # Calculate proposal density
    L1_i = unlist(imp.S2[i,"L1"]); L2_i = unlist(imp.S2[i,"L2"]);
    a = theta.gen.M %>% mutate(cond.m = ini.m2 + sqrt(ini.s22/ini.s11)*ini.rho*(L1_i-ini.m1), #E(L2|L1;theta.k)
                               cond.v = (1-ini.rho^2)*ini.s22, # Var(L2|L1;theta.k)
                               hm = dnorm(L2_i, mean= cond.m, sd = sqrt(cond.v))) %>%
      dplyr::select(hm) %>% unlist(.) %>% mean(.)
    a}) %>% Reduce("c",.)
  
  imp.S2$hm = imp.S2.hm
  
  # pattern S3, imputed values L2:
  mis.S3 = C_R %>% filter(pts == "S3") 
  
  imp.S3 = sapply(1:M, function(m){
    sel.row.theta = sample(1:nrow(theta.gen.M), nrow(mis.S3), replace = T) # SAMPLE from the proposal density
    samp.iter = mis.S3 %>% dplyr::select(ID, L2) %>% dplyr::mutate(sel.row = sel.row.theta) %>% 
      left_join(., theta.gen.M) %>% 
      mutate(cond.m =ini.m1 + sqrt(ini.s11/ini.s22)*ini.rho*(L2-ini.m2), #E[L1|L2; theta.k]
             cond.v = (1-ini.rho^2)*ini.s11) # var(L2|L1; theta.k)
    rnorm(nrow(mis.S3), mean = samp.iter$cond.m, sd = sqrt(samp.iter$cond.v))}) %>% 
    as.data.frame(.) %>% setNames(1:M) %>% mutate(ID = mis.S3$ID) %>% gather(., key = "iter", value= "L1", '1':'100') %>% 
    left_join(., mis.S3 %>% dplyr::select(ID, L2)) %>% 
    mutate(Q_S2 = exp(coef.S2[1]+(coef.S2[2]*L1)), 
           Q_S3 = exp(coef.S3[1]+(coef.S3[2]*L2)),
           Q_S4 = exp(coef.S4[1])*(Q_S3+Q_S2), 
           pi_L = Q_S3 /(1+Q_S2+Q_S3+Q_S4)) # P(R=01|L1*,L2)
  
  imp.S3.hm = lapply(1:nrow(imp.S3), function(i){
    L1_i = unlist(imp.S3[i,"L1"]); L2_i = unlist(imp.S3[i,"L2"]);
    theta.gen.M %>% mutate(cond.m =ini.m1 + sqrt(ini.s11/ini.s22)*ini.rho*(L2_i-ini.m2), #E[L1|L2; theta.k]
                           cond.v = (1-ini.rho^2)*ini.s11, #E[L1|L2;theta.k]
                           hm = dnorm(L1_i, mean= cond.m, sd = sqrt(cond.v))) %>% dplyr::select(hm) %>% unlist(.) %>% mean(.)}) %>%
    Reduce("c",.)
  imp.S3$hm = imp.S3.hm
  
  # pattern S4, imputed values:
  mis.S4 = C_R %>% filter(pts == "S4") 
  
  imp.S4 = lapply(1:M, function(m){
    sel.row.theta = sample(1:nrow(theta.gen.M), nrow(mis.S4), replace = T) #SAMPLE from the proposal mixture density
    samp.iter.temp = mis.S4 %>% dplyr::select(ID) %>% dplyr::mutate(sel.row = sel.row.theta) %>% left_join(., theta.gen.M) 
    samp.iter.temp$L1_Z = rnorm(nrow(mis.S4)); samp.iter.temp$L2_Z = rnorm(nrow(mis.S4)); 
    samp.iter = samp.iter.temp %>%
      mutate(L1 = ini.m1 +sqrt(ini.s11)*L1_Z, #simulate (L1*,L2*) from f(L1*,L2*;theta.k)
             L2 = ini.m2 + sqrt(ini.s22)*(ini.rho*L1_Z + sqrt(1 - ini.rho^2)*L2_Z),
             iter = m) %>% dplyr::select(ID,L1,L2, iter)
    samp.iter #u = m1 + s1 * X and v = m2 + s2 * (rho X + sqrt(1 - rho^2) Y), u,v ~iid N(0,1)
  }) %>% Reduce("rbind",.) %>% as.data.frame(.) %>%
    mutate(Q_S2 = exp(coef.S2[1]+(coef.S2[2]*L1)), 
           Q_S3 = exp(coef.S3[1]+(coef.S3[2]*L2)),
           Q_S4 = exp(coef.S4[1])*(Q_S3+Q_S2), 
           pi_L = Q_S4 /(1+Q_S2+Q_S3+Q_S4)) #P(R=00|L1*,L2*)
  
  imp.S4.hm = lapply(1:nrow(imp.S4), function(i){
    L1_i = unlist(imp.S4[i,"L1"]); L2_i = unlist(imp.S4[i,"L2"]);
    theta.gen.M %>% mutate(cond.m = ini.m1 + sqrt(ini.s11/ini.s22)*ini.rho*(L2_i-ini.m2), # E(L1|L2)
                           cond.v = (1-ini.rho^2)*ini.s11,#Var[L1|L2]
                           hm = dnorm(L1_i, mean = cond.m, sd = sqrt(cond.v))*dnorm(L2_i, mean = ini.m2, sd = sqrt(ini.s22))) %>% #f(L1*|L2*)*f(L2;theta.k)
      dplyr::select(hm) %>% unlist(.) %>% mean(.)}) %>% Reduce("c",.) #K^{-1}\sum_{k=1}^{K}f(L1*|L2*;\theta.k)*f(L2;theta.k)
  
  imp.S4$hm = imp.S4.hm
  
  return(list(imputed.S1 = imp.S1,
              imputed.S2  = imp.S2,
              imputed.S3 = imp.S3,
              imputed.S4  = imp.S4))
}

one.iter.EM.kim = function(ini.norm_pijk, imp.Y.vals, C_R, M){
  #current model parameters:
  for(i in 1:length(ini.norm_pijk)){assign(names(ini.norm_pijk)[i], ini.norm_pijk[i])}
  ini.sigma = matrix(c(ini.s11, ini.rho*sqrt(ini.s11*ini.s22), ini.rho*sqrt(ini.s11*ini.s22), ini.s22), ncol=2) #Cov(L1,L2)
  
  # E-step by each pattern, Calculate fractional weights:
  # S2 : L2 missing
  imp.S2 = imp.Y.vals$imputed.S2 %>%
    mutate(cond.m = ini.m2 + sqrt(ini.s22/ini.s11)*ini.rho*(L1-ini.m1),
           cond.v = (1-ini.rho^2)*ini.s22,
           ran.wt.t = dnorm(L2, mean = cond.m, sd = sqrt(cond.v)), # current f(L2*|L1;theta_t)
           new_wt.t = pi_L*ran.wt.t/hm) # Fractional weights proportionality {new_wt.t = p(R=10|L)*f(L2|L1;theta_t)/hm}
  
  imp.S2.SS = imp.S2 %>% left_join(., imp.S2 %>% group_by(ID)%>% summarise_at(vars(contains("_wt.t")), sum) %>% setNames(c("ID","sum_new.wt.t"))) %>% # normalized fractional weights
    mutate(final_wt.t = new_wt.t/sum_new.wt.t) %>% dplyr::select(ID, L1:L2, final_wt.t)
  
  # S3 :L1 missing
  imp.S3 = imp.Y.vals$imputed.S3 %>%
    mutate(cond.m =ini.m1 + sqrt(ini.s11/ini.s22)*ini.rho*(L2-ini.m2), # E[L1|L2; theta.t]
           cond.v = (1-ini.rho^2)*ini.s11, # Var[L1|L2; theta.t]
           ran.wt.t = dnorm(L1, mean = cond.m, sd = sqrt(cond.v)), # current f(L1|L2;theta_t)
           new_wt.t = pi_L*ran.wt.t/hm) # Fractional weights proportionality {new_wt.t = p(R=01|L)*f(L1|L2; theta_t)/hm}
  
  
  imp.S3.SS = imp.S3 %>% left_join(., imp.S3 %>% group_by(ID)%>% summarise_at(vars(contains("_wt.t")), sum) %>% setNames(c("ID","sum_new.wt.t")) ) %>%
    mutate(final_wt.t = new_wt.t/sum_new.wt.t) %>% dplyr::select(ID, L1:L2, final_wt.t)
  
  # S4: L2, L3 missing
  imp.S4 = imp.Y.vals$imputed.S4 %>% mutate(cond.m = ini.m1 + sqrt(ini.s11/ini.s22)*ini.rho*(L2-ini.m2), # E(L1|L2*;theta.t)
                                            cond.v = (1-ini.rho^2)*ini.s11, # Var(L1|L2*;theta.t)
                                            ran.wt.t = dnorm(L1, mean = cond.m, sd = sqrt(cond.v))*dnorm(L2, mean = ini.m2, sd = sqrt(ini.s22)), #f(L1*|L2*;theta.t)f(L2*;theta.t)
                                            new_wt.t = pi_L*ran.wt.t/hm)
  
  # 
  imp.S4.SS = imp.S4 %>% left_join(., imp.S4 %>% group_by(ID)%>% summarise_at(vars(contains("_wt.t")), sum) %>% setNames(c("ID","sum_new.wt.t")) ) %>%
    mutate(final_wt.t = new_wt.t/sum_new.wt.t) %>% dplyr::select(ID, L1:L2, final_wt.t)
  
  ## S1 : fully observed:
  imp.S1.SS = C_R %>% filter(pts == "S1") %>% dplyr::select(ID, L1:L2) %>% dplyr::mutate(final_wt.t=1)
  
  # lapply(1:4, function(i){eval(parse(text = paste0("head(imp.S",i,".SS)")))})
  # Calculate E-step: 
  imp.SS = rbind(imp.S1.SS, imp.S2.SS, imp.S3.SS, imp.S4.SS) %>% as.data.frame(.) %>% arrange(ID) #%>% 
  
  # 
  # # M-step:
  ini_trans = ini.norm_pijk
  ini_trans["ini.s11"] = log(sqrt(abs(ini_trans["ini.s11"])));
  ini_trans["ini.s22"] = log(sqrt(abs(ini_trans["ini.s22"])));
  ini_trans["ini.rho"] = atanh(ini_trans["ini.rho"])
  prev.ll = loglik(ini_trans, imp.SS)
  prev.grad = grad(loglik, ini_trans,imp.SS = imp.SS)
  # print(prev.grad)
  # 
  news_optim <- optim(ini_trans, loglik,imp.SS = imp.SS, method = "BFGS")
  news <- news_optim$par
  new.ll = loglik(news, imp.SS)
  new.grad = grad(loglik, news,imp.SS = imp.SS)
  # print(new.grad)
  news["ini.s11"] = exp(news["ini.s11"])^2
  news["ini.s22"] = exp(news["ini.s22"])^2
  news["ini.rho"] = tanh(news["ini.rho"])
  # new.ll = loglik(u_trans, imp.SS)
  
  # 
  news.names = paste0("ini.",c("m1","m2","s11","s22","rho")) 
  names(news) = news.names
  
  
  return(list(news = news,
              new.ll = new.ll,
              prev.ll = prev.ll,
              prev.grad = prev.grad,
              new.grad = new.grad))
  
}


EM_pars.kim = function(ini.norm_pijk, imp.Y.vals.yang, C_R, M){
  
  old.pijk = ini.norm_pijk
  i = 0
  repeat{
    i = i+1
    # old.pijk = new.pijk
    Step1_EM =one.iter.EM.kim (ini.norm_pijk = old.pijk, imp.Y.vals.yang, C_R, M) 
    # new.Q = Step1_EM$new.ll
    new.pijk = Step1_EM$news
    print(new.pijk); print(Step1_EM$new.ll );
    # print(new.Q); print(Step1_EM$prev.Q.t); 
    print(Step1_EM$new.ll < Step1_EM$prev.ll); # because want to find the minimizer 
    # print(new.pijk);sum(new.pijk)
    diff = sum(abs(new.pijk-old.pijk))
    if(diff>1e-8){
      old.pijk = new.pijk
    } else{break}
  }
  return(new.pijk)
}

FI_EM.df.kim = function(est.pijk.kim, imp.Y.vals.kim, coef.S2, coef.S3, coef.S4, S1, S2, S3, S4,M){
  
  for(i in 1:length(est.pijk.kim)){assign(names(est.pijk.kim)[i], est.pijk.kim[i])}
  ini.s12 = ini.rho*sqrt(ini.s11*ini.s22)
  ini.mus = c(ini.m1, ini.m2); ini.sigma = matrix(c(ini.s11, ini.s12, ini.s12, ini.s22), ncol=2)
  
  imp_df = lapply(paste0("S",1:4), function(pts){
    
    imp.Y.vals_pts = eval(parse(text = paste0("imp.Y.vals.kim$imputed.",pts)))
    
    if(pts == "S2"){ #L2 missing
      imp.Y.update = imp.Y.vals_pts %>%
        mutate(cond.m = ini.m2 + (ini.s12*(L1-ini.m1)/ini.s11), #E[L2|L1; theta.t]
               cond.v = ini.s22-((ini.s12^2)/ini.s11), #Var[E2|L1; theta.t]
               ran.wt.t = dnorm(L2, mean = cond.m, sd = sqrt(cond.v)), #  f(L2*|L1,theta_t)
               new_wt.t = pi_L*ran.wt.t/hm) %>% #P(R=10|L1,L2*) f(L2*|L1,theta_t)/hm
        left_join(., {.} %>% group_by(ID)%>% summarise_at(vars(contains("_wt.t")), sum) %>% setNames(c("ID","sum_new.wt.t"))) %>%
        mutate(final_wt.t = new_wt.t/sum_new.wt.t)%>% dplyr::select(ID,L1,L2,final_wt.t)
    }else if (pts == "S3"){ #L1 missing
      imp.Y.update = imp.Y.vals_pts %>%
        mutate(cond.m = ini.m1 + (ini.s12*(L2-ini.m2)/ini.s22), # E[L1|L2; theta.t]
               cond.v = ini.s11-((ini.s12^2)/ini.s22), # Var[L1|L2; theta.t]
               ran.wt.t = dnorm(L1, mean = cond.m, sd = sqrt(cond.v)), #f(L1*|L2;theta_t)
               new_wt.t = pi_L*ran.wt.t/hm) %>% # P(R=01|L1*,L2)*f(L1*|L2;theta.t)/hm
        left_join(., {.} %>% group_by(ID)%>% summarise_at(vars(contains("_wt.t")), sum) %>% setNames(c("ID","sum_new.wt.t")) ) %>%
        mutate(final_wt.t = new_wt.t/sum_new.wt.t)%>% dplyr::select(ID,L1,L2,final_wt.t)
    }else if (pts == "S4"){
      imp.Y.update = imp.Y.vals_pts %>% mutate(cond.m = ini.m1 + (ini.s12*(L2-ini.m2)/ini.s22), # E(L1|L2; theta.t)
                                               cond.v = ini.s11-((ini.s12^2)/ini.s22),# Var(L1|L2; theta.t)
                                               ran.wt.t = dnorm(L1, mean = cond.m, sd = sqrt(cond.v))*dnorm(L2, mean = ini.m2, sd = sqrt(ini.s22)), # f(L1*|L2*; theta.t)f(L2*;theta.t)
                                               new_wt.t = pi_L*ran.wt.t/hm) %>% # P(R=00|L1*,L2*)f(L1*,L2*;theta.t)/hm
        left_join(., {.} %>% group_by(ID)%>% summarise_at(vars(contains("_wt.t")), sum) %>% setNames(c("ID","sum_new.wt.t")) ) %>%
        mutate(final_wt.t = new_wt.t/sum_new.wt.t)%>% dplyr::select(ID,L1,L2,final_wt.t)
    }else {
      imp.Y.update = imp.Y.vals_pts %>% mutate(final_wt.t = 1) %>% dplyr::select(ID, L1,L2, final_wt.t)
    }
    
    return(imp.Y.update) }) %>% Reduce("rbind",.) %>% as.data.frame(.)
  
  
  return(imp_df )
}

#--------------------------------------------------------
# Chen RA
#--------------------------------------------------------
RA_imp.df = function(C_R, S1, S2, S3,S4, M){
  
  core.set = C_R %>% filter(pts == "S1")
  Donor.for.S2 = lm(L2 ~ L1, data = core.set); coef.DS2 = coef(Donor.for.S2); sig2.DS2 = sum(residuals(Donor.for.S2)^2)/nrow(core.set)  # E[L2|L1, R=11]
  Donor.for.S3 = lm(L1 ~  L2, data = core.set); coef.DS3 = coef(Donor.for.S3); sig2.DS3 = sum(residuals(Donor.for.S3)^2)/nrow(core.set) # E[L1|L2, R=11]
  
  core_S2 = C_R %>% filter(pts == "S2")
  Donor.for.S4_S2 = lm(L1 ~1,  data = core_S2); coef.DS4_S2 = coef(Donor.for.S4_S2); sig2.DS4_S2 = sum(residuals(Donor.for.S4_S2)^2)/nrow(core_S2) # E[L1|R = 10]
  
  core_S3 = C_R %>% filter(pts == "S3")
  Donor.for.S4_S3 = lm(L2 ~ 1, data = core_S3); coef.DS4_S3 = coef(Donor.for.S4_S3); sig2.DS4_S3 = sum(residuals(Donor.for.S4_S3)^2)/nrow(core_S3) # E[L2|R = 01]
  
  
  # 
  RA_dat = lapply(paste0("S",1:4), function(pt){
    if(pt == "S1"){
      imp.df =  C_R %>% filter(pts == pt) %>% dplyr::select(ID, L1, L2) #%>% mutate(wt = 1)
    }else if (pt == "S2"){
      imp.df = C_R %>% filter(pts == pt)%>% dplyr::select(ID, L1) %>%
        mutate(L2= coef.DS2[1]+coef.DS2[2]*L1) # L2* = E[L2|L1_i, R=11]
    }else if (pt == "S3"){
      imp.df = C_R %>% filter(pts == pt)%>% dplyr::select(ID, L2) %>%
        mutate(L1 = coef.DS3[1]+coef.DS3[2]*L2) # L1* = E[L1|L2_i, R=11]
    }else if (pt == "S4"){
      mis.S4 = C_R %>% filter(pts == pt) %>% dplyr::select(ID, L1, L2)
      
      n_S3 = nrow(core_S3); n_S2 = nrow(core_S2); p_S2 = n_S2/(n_S2+n_S3); p_S3 = 1-p_S2; #path probability: p_S2 =  P(00->10->00) p_S3 = P(00->01->11)
      imp.df = lapply(1:nrow(mis.S4), function(i){
        mis.i = mis.S4[i,]
        
        mis.i; mis.i.path1 = mis.i;mis.i.path2 = mis.i;
        mis.i.path1$L1 = coef.DS4_S2[1] #L1.1* = E[L1|R=10]
        mis.i.path1$L2 = coef.DS2[1]+(coef.DS2[2]*coef.DS4_S2[1]) # L2.1* = E[L2|L1*,R=11]
        
        mis.i.path2$L2 = coef.DS4_S3[1] # L2.2* = E[L2|R=01]
        mis.i.path2$L1 = coef.DS3[1]+(coef.DS3[2]* mis.i.path2$L2) # L1.2* = E[L1|L2*,R=11]
        mis.i$L1 = p_S2*mis.i.path1$L1+(1-p_S2)*mis.i.path2$L1 # L1* = p_s2*L1.1*+p_s3*L1.2*
        mis.i$L2 = p_S2*mis.i.path1$L2+(1-p_S2)*mis.i.path2$L2 # L2* = p_s2*L2.1*+p_s3*L2.2*
        mis.i}) %>% Reduce("rbind",.)%>% as.data.frame(.)
    }
    return(imp.df)}) %>% Reduce("rbind",.) %>% as.data.frame(.) %>% arrange(ID)
  
  RA_mean = RA_dat %>%
    mutate(RA_L1 = L1,
           RA_L2 = L2,
           RA_L1L2 = (L1+L2)) %>% group_by(ID) %>% summarise_at(vars(contains("RA_")), mean) %>% dplyr::select(-ID) %>%apply(.,2,mean)
  
  
  
  return( list(RA_mean = RA_mean,
               RA_dat = RA_dat) )
}

#--------------------------------------------------------
# Chen AIPW
#--------------------------------------------------------

Semi_imp= function(C_R, RA_dat,coef.S2, coef.S3, coef.S4, M, n, IPW.res, ini.norm_pijk){
  est.C = C_R
  # Step1: RA m_LR_R
  m_LR_R <- RA_dat %>% 
    mutate(h1 = L1 , h2 = L2, h3 = (L1+L2)) %>% 
    group_by(ID) %>% summarise_at(vars(contains("h")), mean) 
  
  m_Ll_Rr =  est.C %>% dplyr::select(ID,pts, L1, L2)%>% left_join(., m_LR_R) %>% 
    mutate(O_LS2 = exp(coef.S2[1]+coef.S2[2]*L1),
           O_LS3 = exp(coef.S3[1]+coef.S3[2]*L2),
           O_LS4 = exp(coef.S4[1]))
  # for(i in 1:length(ini.norm_pijk)){assign(names(ini.norm_pijk)[i],ini.norm_pijk[i])}
  # ini.mus = c(ini.m1, ini.m2); ini.sigma = matrix(c(ini.s11, ini.s12, ini.s12, ini.s22), ncol=2)
  # 
  core.set = C_R %>% filter(pts == "S1")
  Donor.for.S2 = lm(L2 ~ L1, data = core.set); coef.DS2 = coef(Donor.for.S2); sig2.DS2 = sum(residuals(Donor.for.S2)^2)/nrow(core.set)
  Donor.for.S3 = lm(L1 ~  L2, data = core.set); coef.DS3 = coef(Donor.for.S3); sig2.DS3 = sum(residuals(Donor.for.S3)^2)/nrow(core.set)
  
  core_S2 = C_R %>% filter(pts == "S2")
  Donor.for.S4_S2 = lm(L1 ~1,  data = core_S2); coef.DS4_S2 = coef(Donor.for.S4_S2); sig2.DS4_S2 = sum(residuals(Donor.for.S4_S2)^2)/nrow(core_S2)
  
  core_S3 = C_R %>% filter(pts == "S3")
  Donor.for.S4_S3 = lm(L2 ~ 1, data = core_S3); coef.DS4_S3 = coef(Donor.for.S4_S3); sig2.DS4_S3 = sum(residuals(Donor.for.S4_S3)^2)/nrow(core_S3)
  
  n_S3 = nrow(core_S3); n_S2 = nrow(core_S2); p_S2 = n_S2/(n_S2+n_S3); p_S3 = 1-p_S2;
  
  EIF_m_Ll_Rr = lapply(paste0("S",1:4), function(pt){
    
    if(pt == "S1"){
      
      fil_pts =  m_Ll_Rr %>% filter(pts == pt) #E[g|L,R=111]=g
      
      imp.df_S2 = fil_pts %>% mutate(L2 = NA) %>% dplyr::select(ID, L1,L2) %>%  #E[g|L1,R=10]
        dplyr::mutate(L2= coef.DS2[1]+coef.DS2[2]*L1) %>%
        mutate(mL_S2_h1 = L1,
               mL_S2_h2 = L2,
               mL_S2_h3 = (L1+L2)) %>% dplyr::select(ID, mL_S2_h1, mL_S2_h2, mL_S2_h3)
      
      imp.df_S3 = fil_pts %>% mutate(L1 = NA) %>% dplyr::select(ID, L1,L2) %>%  #e[g|L2,R=01]
        dplyr::mutate(L1 = coef.DS3[1]+coef.DS3[2]*L2) %>% 
        mutate(mL_S3_h1 = L1,
               mL_S3_h2 = L2,
               mL_S3_h3 = (L1+L2)) %>% dplyr::select(ID, mL_S3_h1, mL_S3_h2, mL_S3_h3)
      
      
      EIF_df_S1 = fil_pts %>% rename(mL_S1_h1 = h1, mL_S1_h2 = h2, mL_S1_h3 = h3) %>% 
        left_join(.,imp.df_S2) %>% left_join(.,imp.df_S3) %>% 
        # left_join(., imp.df_S4) %>%
        mutate(EIF_h1 = -mL_S2_h1*O_LS2*(1+O_LS4) -mL_S3_h1*O_LS3*(1+O_LS4),
               EIF_h2 = -mL_S2_h2*O_LS2*(1+O_LS4) -mL_S3_h2*O_LS3*(1+O_LS4),
               EIF_h3 = -mL_S2_h3*O_LS2*(1+O_LS4) -mL_S3_h3*O_LS3*(1+O_LS4)) %>%
        dplyr::select(ID, pts, L1,L2,EIF_h1, EIF_h2, EIF_h3)
      
      EIF_df = EIF_df_S1
      
    }else if(pt == "S2"){
      fil_S2 = m_Ll_Rr %>% filter(pts == pt) 
      
      mis_S2_S4 = fil_S2 %>% mutate(L1 = NA) %>% dplyr::select(ID, L1,L2)
      
      imp.df_S2_S4  = lapply(1:nrow(mis_S2_S4), function(i){ #E[g|R=00]
        mis.i = mis_S2_S4[i,]
        mis.i; mis.i.path1 = mis.i;mis.i.path2 = mis.i; #path1: 00->10->11
        mis.i.path1$L1 = coef.DS4_S2[1]
        mis.i.path1$L2 = coef.DS2[1]+(coef.DS2[2]*mis.i.path1$L1)
        
        mis.i.path2$L2 = coef.DS4_S3[1] # path2: 00->01->11
        mis.i.path2$L1 = coef.DS3[1]+(coef.DS3[2]*mis.i.path2$L2)
        mis.i$L1 = p_S2*mis.i.path1$L1+(1-p_S2)*mis.i.path2$L1
        mis.i$L2 = p_S2*mis.i.path1$L2+(1-p_S2)*mis.i.path2$L2
        mis.i}) %>% Reduce("rbind",.)%>% as.data.frame(.) %>%
        mutate(mL_S4_h1 = L1,
               mL_S4_h2 = L2,
               mL_S4_h3 = (L1+L2)) %>% dplyr::select(ID, mL_S4_h1, mL_S4_h2, mL_S4_h3)
      
      
      EIF_df_S2 = fil_S2 %>% mutate(mL_S2_h1 = h1, mL_S2_h2 = h2, mL_S2_h3 = h3) %>% 
        left_join(.,imp.df_S2_S4) %>%
        mutate(EIF_h1 = mL_S2_h1*(1+O_LS4) -mL_S4_h1*O_LS4,
               EIF_h2 = mL_S2_h2*(1+O_LS4) -mL_S4_h2*O_LS4,
               EIF_h3 = mL_S2_h3*(1+O_LS4) -mL_S4_h3*O_LS4) %>%
        dplyr::select(ID, pts, L1,L2, EIF_h1, EIF_h2, EIF_h3)
      
      
      EIF_df = EIF_df_S2
      
    }else if (pt == "S3"){
      fil_S3 = m_Ll_Rr %>% filter(pts == pt) 
      
      mis_S3_S4 = fil_S3 %>% mutate(L1 = NA) %>% dplyr::select(ID, L1, L2)
      
      imp.df_S3_S4 =lapply(1:nrow(mis_S3_S4), function(i){ #E[g|R=00]
        mis.i = mis_S3_S4[i,]
        mis.i; mis.i.path1 = mis.i;mis.i.path2 = mis.i;
        mis.i.path1$L1 = coef.DS4_S2[1]
        mis.i.path1$L2 = coef.DS2[1]+(coef.DS2[2]*mis.i.path1$L1)
        
        mis.i.path2$L2 = coef.DS4_S3[1]
        mis.i.path2$L1 = coef.DS3[1]+(coef.DS3[2]*mis.i.path2$L2)
        mis.i$L1 = p_S2*mis.i.path1$L1+(1-p_S2)*mis.i.path2$L1
        mis.i$L2 = p_S2*mis.i.path1$L2+(1-p_S2)*mis.i.path2$L2
        mis.i}) %>% Reduce("rbind",.)%>% as.data.frame(.) %>%
        mutate(mL_S4_h1 = L1,
               mL_S4_h2 = L2,
               mL_S4_h3 = (L1+L2)) %>% dplyr::select(ID, mL_S4_h1, mL_S4_h2, mL_S4_h3)
      
      
      EIF_df_S3 = fil_S3 %>% mutate(mL_S3_h1 = h1, mL_S3_h2 = h2, mL_S3_h3 = h3) %>% 
        left_join(.,imp.df_S3_S4) %>%
        mutate(EIF_h1 = mL_S3_h1*(1+O_LS4) -mL_S4_h1*O_LS4,
               EIF_h2 = mL_S3_h2*(1+O_LS4) -mL_S4_h2*O_LS4,
               EIF_h3 = mL_S3_h3*(1+O_LS4) -mL_S4_h3*O_LS4) %>%
        dplyr::select(ID, pts, L1,L2,EIF_h1, EIF_h2, EIF_h3)
      EIF_df = EIF_df_S3
    }else if(pt == "S4"){
      EIF_df_S4 =  m_Ll_Rr %>% filter(pts == pt) %>% 
        dplyr::select(ID, pts, L1,L2,h1, h2, h3) %>% 
        setNames(c("ID","pts", "L1","L2","EIF_h1", "EIF_h2", "EIF_h3"))
      
      EIF_df = EIF_df_S4
    }
    
    EIF_df
  }) %>% Reduce("rbind",.) %>% as.data.frame(.) %>%arrange(ID)
  
  SEMI_sum =  EIF_m_Ll_Rr %>% 
    left_join(., IPW.res %>% mutate(h1 = L1/pi_L, h2 = L2/pi_L, h3 = (L1+L2)/pi_L) %>% dplyr::select(ID, h1:h3)) %>%
    mutate(h1 = ifelse(is.na(h1)==T, 0,h1),
           h2 = ifelse(is.na(h2)==T, 0,h2),
           h3 = ifelse(is.na(h3)==T, 0,h3)) %>% 
    mutate(semi_sum_h1 = h1+EIF_h1,
           semi_sum_h2 = h2+EIF_h2,
           semi_sum_h3 = h3+EIF_h3) %>% 
    dplyr::select(semi_sum_h1, semi_sum_h2, semi_sum_h3) %>% apply(.,2,sum)
  SEMI.mean = SEMI_sum/n 
  names(SEMI.mean) = paste0("SEMI_",c("L1","L2","L1L2"))
  
  return(SEMI.mean= SEMI.mean)
}

#--------------------------------------------------------
# Bootstrp Variance Estimation
#--------------------------------------------------------
est.FI.EM.b = function(est.pijk, coef.S2, coef.S3, coef.S4, C_R, n, MODEL){
  #Create Sample n= 1000 ->nc:
  s12.b = est.pijk["ini.rho"]*sqrt(est.pijk["ini.s11"]*est.pijk["ini.s22"])
  T.mean.b = c(est.pijk["ini.m1"], est.pijk["ini.m2"])
  T.Sigma.b = matrix(c(est.pijk["ini.s11"], s12.b, s12.b, est.pijk["ini.s22"]), ncol = 2)
  samp.b = mvrnorm(n = n, mu = T.mean.b, Sigma = T.Sigma.b) %>% as.data.frame(.) %>% setNames(paste0("L",1:2)) %>% mutate(ID = row_number()) 
  
  # Create Missing:
  pic_by_mis.b = samp.b %>%
    mutate(O_S2 = exp(coef.S2[1] + coef.S2[2]*L1), 
           O_S3 = exp(coef.S3[1] + coef.S3[2]*L2), 
           O_S4 = exp(coef.S4[1]) ) %>% # Add intercept: more realistic model assumption
    mutate(Q_S1 = 1, Q_S2 = O_S2, Q_S3 = O_S3, Q_S4 = O_S4*(Q_S2+Q_S3)) %>%
    mutate(pi_L = 1/(Q_S1+Q_S2+Q_S3+Q_S4)) %>%
    mutate(P_S1 = pi_L, P_S2 = Q_S2*pi_L, P_S3 = Q_S3*pi_L, P_S4 = Q_S4*pi_L) ## create P(R=r|L) This is real probability for creating missingness
  
  nc_by_mis.b = pic_by_mis.b %>% 
    apply(.,1,function(.){
      T.probc = .[c("P_S1","P_S2","P_S3","P_S4")] %>% unlist(.)
      rmultinom(1,size=1, T.probc)}) %>% as.data.frame(.) %>% t(.) %>% as.data.frame(.) %>%
    setNames(paste0("S",1:4))# %>% mutate(nc = S1+S2+S3+S4)
  
  nc_by_mis.b = cbind(pic_by_mis.b, nc_by_mis.b) %>% as.data.frame(.) %>% dplyr::select(ID, L1:L2,S1:S4) 
  
  ## Dataframe with missingness:
  C_R.b = nc_by_mis.b%>%
    mutate(L2 = ifelse(S2==1 |S4 ==1, NA,L2),
           L1 = ifelse(S3==1 |S4 ==1, NA,L1)) %>% 
    mutate(pts = case_when((S1==1) ~ "S1",
                           (S2==1) ~ "S2",
                           (S3==1) ~ "S3",
                           (S4==1)  ~ "S4")) #%>% 
  #mutate(ID = row_number())
  
  S1 = c(T,T); S2 = c(T,F); S3 = c(F,T);S4 = c(F,F);
  
  
  #Inverse probability weighting:
  Orlr.b = OrLr(C_R.b)
  # for(i in 1:length(Orlr)){assign(names(Orlr)[[i]], Orlr[[i]])} #"coef.S2"    "coef.S3"    "coef.S4"
  for(i in 1:length(Orlr.b)){assign(names(Orlr.b)[[i]], Orlr.b[[i]])} #"coef.S2"    "coef.S3"    "coef.S4"
  
  IPW.res.b = IPW(C_R.b, Orlr.b)
  
  M = 100 # number of imputation
  
  #Fractional Imputation:
  # #EM algorithm:
  ini.norm_pijk = EM_initials(C_R.b) 
 
    imp.Y.vals.kim.b = impY_hm_piL_kim(ini.norm_pijk, C_R.b, coef.S2, coef.S3, coef.S4, S1, S2, S3, S4,M)
    est.pijk.kim.b = EM_pars.kim(ini.norm_pijk, imp.Y.vals.kim.b, C_R.b, M);
    FI.EM.kim.b = FI_EM.df.kim(est.pijk.kim.b, imp.Y.vals.kim.b, coef.S2, coef.S3, coef.S4, S1, S2, S3, S4,M);
    FI.EM.sum.kim.b = FI.EM.kim.b%>% mutate(FI.EM.K.b_L1 = (L1)*final_wt.t,
                                            FI.EM.K.b_L2 =L2*final_wt.t,
                                            FI.EM.K.b_L1L2 =(L1+L2)*final_wt.t)%>%
      dplyr::select(contains("FI.EM.K.b")) %>% apply(.,2,sum)
    FI.EM.mean.kim.b = FI.EM.sum.kim.b/n
    FI.EM.means.b = FI.EM.mean.kim.b
  
    return(FI.EM.means.b)
  
}
