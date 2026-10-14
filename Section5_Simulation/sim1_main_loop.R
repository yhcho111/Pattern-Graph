rm(list=ls(all=TRUE))
source("fns_Algorithm.R")

required_files[!file.exists(required_files)]

required_pkgs <- c("lme4","sampling","data.table","dplyr","doParallel",
                   "fastGHQuad","merDeriv","MASS","numDeriv","Matrix",
                   "tidyr","reldist","plyr")
all(sapply(required_pkgs, require, character.only=TRUE))

## read model parameters
lambda  <- Scenario[s,"lambda"];  a_i    <- Scenario[s,"a_i"]
b_j     <- Scenario[s,"b_j"];     c_k    <- Scenario[s,"c_k"]
ab_ij   <- Scenario[s,"ab_ij"];   ac_ik  <- Scenario[s,"ac_ik"]
bc_jk   <- Scenario[s,"bc_jk"];   abc_ijk <- Scenario[s,"abc_ijk"]  # 추가

o_s2_0  <- Scenario[s,"o_s2_0"]; o_s2_L1 <- Scenario[s,"o_s2_L1"]
o_s2_L2 <- Scenario[s,"o_s2_L2"]
o_s3_0  <- Scenario[s,"o_s3_0"]; o_s3_L1 <- Scenario[s,"o_s3_L1"]
o_s3_L3 <- Scenario[s,"o_s3_L3"]
o_s4_0  <- Scenario[s,"o_s4_0"]; o_s4_L1 <- Scenario[s,"o_s4_L1"]

N <- Scenario[s,"N"]; n <- Scenario[s,"n"]

L1 <- c(-1,1); L2 <- c(-1,1); L3 <- c(-1,1)

## DGP: saturated 
C <- expand.grid(L1, L2, L3) %>%
  setNames(paste0("L",1:3)) %>%
  mutate(
    lambda  = lambda,
    a_i     = L1 * a_i,
    b_j     = L2 * b_j,
    c_k     = L3 * c_k,
    ab_ij   = L1 * L2 * ab_ij,
    ac_ik   = L1 * L3 * ac_ik,
    bc_jk   = L2 * L3 * bc_jk,
    abc_ijk = L1 * L2 * L3 * abc_ijk,
    pijk    = exp(lambda + a_i + b_j + c_k +
                    ab_ij + ac_ik + bc_jk + abc_ijk)
  ) %>%
  mutate(pijk = pijk/sum(pijk)) %>%
  mutate(L1=ifelse(L1==-1,0,L1),
         L2=ifelse(L2==-1,0,L2),
         L3=ifelse(L3==-1,0,L3)) %>%
  arrange(desc(L1), desc(L2), desc(L3))

C$nc <- rmultinom(1, size=n, C$pijk)

#------------------------------------------------------
#
#------------------------------------------------------

t.mean <- C %>%
  mutate(t_Y3        = L3*pijk,
         t_Y1Y2Y3    = L1*L2*L3*pijk) %>%
  dplyr::select(contains("t")) %>% apply(.,2,sum)

## Missing mechanism
pic_by_mis <- C %>%
  mutate(O_S2 = exp(o_s2_0 + o_s2_L1*L1 + o_s2_L2*L2),
         O_S3 = exp(o_s3_0 + o_s3_L1*L1 + o_s3_L3*L3),
         O_S4 = exp(o_s4_0 + o_s4_L1*L1)) %>%
  mutate(Q_S1=1, Q_S2=O_S2, Q_S3=O_S3, Q_S4=O_S4*(Q_S2+Q_S3)) %>%
  mutate(pi_L=1/(Q_S1+Q_S2+Q_S3+Q_S4)) %>%
  mutate(P_S1=pi_L, P_S2=Q_S2*pi_L, P_S3=Q_S3*pi_L, P_S4=Q_S4*pi_L)

nc_by_mis <- pic_by_mis %>%
  apply(.,1,function(.){
    T.nc    <- .["nc"] %>% unlist()
    T.probc <- .[c("P_S1","P_S2","P_S3","P_S4")] %>% unlist()
    rmultinom(1, size=T.nc, T.probc)
  }) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
  setNames(paste0("S",1:4))

nc_by_mis <- cbind(pic_by_mis, nc_by_mis) %>%
  as.data.frame() %>% dplyr::select(L1:L3, nc, S1:S4)

C_R <- C %>% left_join(.,nc_by_mis) %>%
  mutate(cate=paste0("c",1:nrow(.))) %>%
  dplyr::select(L1:L3, nc, S1:S4, cate)

S1<-c(T,T,T); S2<-c(T,T,F); S3<-c(T,F,T); S4<-c(T,F,F)

#------------------------------------------------------
# complete-case
#------------------------------------------------------

samp.m.sum <- C_R %>%
  mutate(samp_Y3        = L3*nc,
         samp_Y1Y2Y3    = L1*L2*L3*nc) %>%
  dplyr::select(contains("samp_")) %>% apply(.,2,sum)
samp.mean <- samp.m.sum/n

C_R.pt     <- C_R %>% dplyr::select(cate, L1, L2, L3)
C_R.cnt    <- C_R %>% dplyr::select(cate, S1:S4) %>% gather(.,pt,count,S1:S4)
C_R.cnt.pt <- C_R.cnt %>% left_join(.,C_R.pt) %>%
  mutate(L3=ifelse(pt=="S2"|pt=="S4",NA,L3),
         L2=ifelse(pt=="S3"|pt=="S4",NA,L2))

#------------------------------------------------------
# IPW
#------------------------------------------------------

## Estiamte phi
Orlr    <- OrLr(C_R)
for(i in 1:length(Orlr)) assign(names(Orlr)[[i]], Orlr[[i]])
IPW.res <- IPW(C_R, Orlr)

IPW.sum <- IPW.res %>%
  mutate(IPW_Y3        = L3*S1/pi_L,
         IPW_Y1Y2Y3    = L1*L2*L3*S1/pi_L) %>%
  dplyr::select(contains("IPW_")) %>% apply(.,2,sum)
IPW.denom <- IPW.res %>% mutate(ipw.denom=S1/pi_L) %>%
  dplyr::select(ipw.denom) %>% sum() %>% unlist()
IPW.mean  <- IPW.sum/IPW.denom

IPW.res_FI <- IPW.res %>%
  mutate(P_S2=Q_S2.L*pi_L, P_S3=Q_S3.L*pi_L,
         P_S4=Q_S4.L*pi_L, P_S1=pi_L,
         lat=P_S1+P_S2+P_S3+P_S4) %>%
  dplyr::select(L1:L3, P_S1, P_S2, P_S3, P_S4)

long.pop  <- C_R.cnt.pt[rep(1:nrow(C_R.cnt.pt),times=C_R.cnt.pt$count),] %>%
  dplyr::select(pt, L1:L3)
short.pop <- long.pop %>% group_by(pt,L1,L2,L3) %>% dplyr::summarize(count=n())

#------------------------------------------------------
# RA & AIPW of Chen
#------------------------------------------------------
SEMI_RA  <- RA_SEMI.imp(C_R, short.pop, S1, S2, S3, S4,
                        Orlr, n, IPW.res, C_R.cnt.pt, hlist)
SEMI.mean <- SEMI_RA$SEMI.mean
RA.mean   <- SEMI_RA$RA.mean %>% unlist()
m_LR_R    <- SEMI_RA$m_LR_R

#------------------------------------------------------
# FI _ two-way
#------------------------------------------------------

ini.norm_pijk <- rep(1/8,8)
est.pijk      <- EM_pijk(ini.norm_pijk, C_R, C_R.cnt.pt,
                         S1, S2, S3, S4, IPW.res_FI, n)
FI.EM         <- FI_EM.df(est.pijk, long.pop, S1, S2, S3, S4, IPW.res_FI, n)
FI.EM.sum     <- FI.EM %>%
  mutate(FI.EM_Y3        = L3*final_prob,
         FI.EM_Y1Y2Y3    = L1*L2*L3*final_prob) %>%
  dplyr::select(contains("FI.EM")) %>% apply(.,2,sum)
FI.EM.mean <- FI.EM.sum/n

#------------------------------------------------------
# FI.ST
#------------------------------------------------------
est.pijk.ST   <- EM_pijk.ST(ini.norm_pijk, C_R, C_R.cnt.pt,
                            S1, S2, S3, S4, IPW.res_FI, n)
FI.EM.ST      <- FI_EM.df(est.pijk.ST, long.pop, S1, S2, S3, S4, IPW.res_FI, n)
FI.EM.sum.ST  <- FI.EM.ST %>%
  mutate(FI.EM.ST_Y3        = L3*final_prob,
         FI.EM.ST_Y1Y2Y3    = L1*L2*L3*final_prob) %>%
  dplyr::select(contains("FI.EM.ST")) %>% apply(.,2,sum)
FI.EM.ST.mean <- FI.EM.sum.ST/n

#------------------------------------------------------
# Bootstrap variance (Section 3.2): B replicates
#------------------------------------------------------
B <- 100
FI.EM.b  <- t(sapply(1:B, function(b)
  est.FI.EM.b(est.pijk,    coef.S2, coef.S3, coef.S4, C_R, n, MODEL = "TW")))
FI.EM.ST.b <- t(sapply(1:B, function(b)
  est.FI.EM.b(est.pijk.ST, coef.S2, coef.S3, coef.S4, C_R, n, MODEL = "ST")))

FI.EM.sd    <- apply(FI.EM.b[,    c("FI.EM.b_Y3",    "FI.EM.b_Y1Y2Y3")],    2, sd)
FI.EM.ST.sd <- apply(FI.EM.ST.b[, c("FI.EM.ST.b_Y3", "FI.EM.ST.b_Y1Y2Y3")], 2, sd)
