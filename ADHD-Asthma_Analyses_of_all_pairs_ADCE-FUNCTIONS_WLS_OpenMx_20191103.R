##########################################################################
###
### Functions to be sourced
###
###
###
##########################################################################

vlen <- 4
dlen <- 4
tlen <- vlen+dlen

##########################################################################
### Correlations
corrFun <- function( dat=dat ){
### Model
baseMod <- mxModel( 'Base' ,
# Manifest endogenous means
  mxMatrix('Full', vlen , 1 , free=F , values=0,
    dimnames=list(varNames,NA) , name='eMean' ),

# Manifest endogenous thresholds
  mxMatrix('Full', 1 , 4 ,free=T,values=thrstart,labels=paste(thrLab,'_F',sep=''),
    dimnames=list(NA,varNames) , name='ThrF' ),
  mxMatrix('Full', 1 , 4 ,free=T,values=thrstart,labels=paste(thrLab,'_H',sep=''),
    dimnames=list(NA,varNames) , name='ThrH' ),

# Loading from latent endogenous to manifest
  mxMatrix('Full',4,8,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(c(rep(T,2),rep(F,2)),2) , rep(c(rep(F,2),rep(T,2)),2) ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(c(rep(.001,2),rep(0,2)),2) , rep(c(rep(0,2),rep(.001,2)),2) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','na','na',
              'betaByr_asthma','betaByr_adhd','na','na',
              'na','na','betaSex_asthma','betaSex_adhd',
              'na','na','betaByr_asthma','betaByr_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r') ) ,name='LambdaY' ) ,

# Loading from latent endogenous to manifest MZ
  mxMatrix('Full',4,6,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(T,8)  ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(.001,8) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','betaSex_asthma','betaSex_adhd',
              'betaByr_asthma','betaByr_adhd','betaByr_asthma','betaByr_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb') ) ,name='LambdaYMZ' ) ,

# Loading from latent endogenous to manifest DZ
  mxMatrix('Full',4,7,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(T,2),rep(F,2),rep(T,4), rep(F,2),rep(T,2) ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(.001,2),rep(0,2),rep(.001,4) ,rep(0,2),rep(.001,2) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','na','na',
              'betaByr_asthma','betaByr_adhd','betaByr_asthma','betaByr_adhd',
              'na','na','betaSex_asthma','betaSex_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r') ) ,name='LambdaYDZ' ) ,


# Covariance of latent endogenous
  mxMatrix('Symm', 8 , 8 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r')) , name= 'Psi' ) ,
  mxMatrix('Symm', 6 , 6 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb')) , name= 'PsiMZ' ) ,
  mxMatrix('Symm', 7 , 7 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r')) , name= 'PsiDZ' ) ,

# Covariance of latent exogenuous
  mxMatrix('Stand', 2 , 2 , byrow=F , ubound=.99 , lbound=-.99 ,
    free=  c(T),labels='rPH' ,values=.3, name='PWTN' ),
# MZ
  mxMatrix('Symm', 2 , 2 , byrow=F , ubound=.99 , lbound=-.99 ,
    free=T,labels=c('ICCmzAs','CTCTmz','ICCmzAd') ,values=.3, name='PBTWmz' ),
# DZ
  mxMatrix('Symm', 2 , 2 , byrow=F , ubound=.99 , lbound=-.99 ,
    free=  c(T),labels=c('ICCdzAs','CTCTDz','ICCdzAd') ,values=.3, name='PBTWdz' ),
# FU
  mxMatrix('Symm', 2 , 2 , byrow=F , ubound=.99 , lbound=-.99 ,
    free=  c(T),labels=c('ICCfuAs','CTCTfu','ICCfuAd') ,values=.3, name='PBTWfu' ),
# MH
  mxMatrix('Symm', 2 , 2 , byrow=F , ubound=.99 , lbound=-.99 ,
    free=  c(T),labels=c('ICCmhAs','CTCTmh','ICCmhAd') ,values=.3, name='PBTWmh' ),
# PH
  mxMatrix('Symm', 2 , 2 , byrow=F , ubound=.99 , lbound=-.99 ,
    free=  c(T),labels=c('ICCphAs','CTCTph','ICCphAd') ,values=.3, name='PBTWph' ),

# Put together
  mxAlgebra( rbind(cbind(PWTN,PBTWmz),cbind(PBTWmz,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMZ' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWdz),cbind(PBTWdz,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiDZ' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWfu),cbind(PBTWfu,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiFU' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWmh),cbind(PBTWmh,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMH' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWph),cbind(PBTWph,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiPH' ),

# Loadings from latent exogenous to endogenous
  mxMatrix('Diag',4,4,free=F,values=1,name='GaDiag'),
  mxMatrix('Zero',4,4,name='Ga0'),
  mxMatrix('Zero',2,4,name='Ga0MZ'),
  mxMatrix('Zero',3,4,name='Ga0DZ'),
  mxAlgebra(rbind( GaDiag , Ga0 ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='Ga' ),
  mxAlgebra(rbind( GaDiag , Ga0MZ ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='GaMZ' ),
  mxAlgebra(rbind( GaDiag , Ga0DZ ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='GaDZ' ),

# Latent variable equation
  mxMatrix('Diag',8,8,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r')) , name='Be' ),
  mxMatrix('Diag',6,6,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb')) , name='BeMZ' ),
  mxMatrix('Diag',7,7,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r')) , name='BeDZ' ),

# Residual Variance of manifest variable
  mxMatrix('Zero',4,4,free=F,
    dimnames=list(varNames,varNames) , name='ThetaE' ),

# Matrices for non-loading onto latent Exogenous
  mxMatrix('Full',nrow=4,ncol=1,free=F,values=0,dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),NA) , name='Ka' ),
  mxMatrix('Full',nrow=0,ncol=4,free=F,dimnames=list(NULL,c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) ,name='Lx' ),
  mxMatrix('Full',nrow=0,ncol=0,free=F,dimnames=list(NULL,NULL) , name='Td'),
  mxMatrix('Full',nrow=0,ncol=1,free=F,dimnames=list(NULL,NULL) , name='Tx' ),
  mxMatrix('Full',nrow=0,ncol=4,free=F,dimnames=list(NULL,varNames) , name='Th')
)

# MZ
mzMod <- mxModel( 'MZ' ,
# Data
  mxData( dat[dat$fam==1,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',6,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaYMZ', PS='Base.PsiMZ',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.BeMZ', 
    LX='Base.Lx',GA='Base.GaMZ',  PH='Base.PhiMZ',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# DZ
dzMod <- mxModel( 'DZ' ,
# Data
  mxData(  dat[dat$fam==2,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',7,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaYDZ', PS='Base.PsiDZ',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.BeDZ', 
    LX='Base.Lx',GA='Base.GaDZ',  PH='Base.PhiDZ',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# Full
fuMod <- mxModel( 'FU' ,
# Data
  mxData(  dat[dat$fam==3,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiFU',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# MH
mhMod <- mxModel( 'MH' ,
# Data
  mxData(  dat[dat$fam==4,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiMH',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrH',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# PH
phMod <- mxModel( 'PH' ,
# Data
  mxData(  dat[dat$fam==5,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiPH',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrH',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

#
modCorr <- mxModel( 'Correlations' ,
  baseMod , mzMod , dzMod , fuMod , mhMod , phMod ,
  mxFitFunctionMultigroup( c( 'MZ','DZ','FU','MH','PH' ) )
)
# Fit model
modCorrFit <- mxRun( modCorr , silent=T )

### Output
# Correlations
out <-c(
mxEval( Base.PWTN[1,2] , modCorrFit) ,
mxEval( Base.PBTWmz[1,1] , modCorrFit) ,
mxEval( Base.PBTWdz[1,1] , modCorrFit) ,
mxEval( Base.PBTWfu[1,1] , modCorrFit) ,
mxEval( Base.PBTWmh[1,1] , modCorrFit) ,
mxEval( Base.PBTWph[1,1] , modCorrFit) ,

mxEval( Base.PBTWmz[2,2] , modCorrFit) ,
mxEval( Base.PBTWdz[2,2] , modCorrFit) ,
mxEval( Base.PBTWfu[2,2] , modCorrFit) ,
mxEval( Base.PBTWmh[2,2] , modCorrFit) ,
mxEval( Base.PBTWph[2,2] , modCorrFit) ,

mxEval( Base.PBTWmz[1,2] , modCorrFit) ,
mxEval( Base.PBTWdz[1,2] , modCorrFit) ,
mxEval( Base.PBTWfu[1,2] , modCorrFit) ,
mxEval( Base.PBTWmh[1,2] , modCorrFit), 
mxEval( Base.PBTWph[1,2] , modCorrFit) 
)
out
}
##########################################################################


##########################################################################
### ADCE-model
adceFun <- function( dat=dat ){
### Model
baseMod <- mxModel( 'Base' ,
# Manifest endogenous means
  mxMatrix('Full', vlen , 1 , free=F , values=0,
    dimnames=list(varNames,NA) , name='eMean' ),

# Manifest endogenous thresholds
  mxMatrix('Full', 1 , 4 ,free=T,values=thrstart,labels=paste(thrLab,'_F',sep=''),
    dimnames=list(NA,varNames) , name='ThrF' ),
  mxMatrix('Full', 1 , 4 ,free=T,values=thrstart,labels=paste(thrLab,'_H',sep=''),
    dimnames=list(NA,varNames) , name='ThrH' ),

# Loading from latent endogenous to manifest
  mxMatrix('Full',4,8,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(c(rep(T,2),rep(F,2)),2) , rep(c(rep(F,2),rep(T,2)),2) ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(c(rep(.001,2),rep(0,2)),2) , rep(c(rep(0,2),rep(.001,2)),2) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','na','na',
              'betaByr_asthma','betaByr_adhd','na','na',
              'na','na','betaSex_asthma','betaSex_adhd',
              'na','na','betaByr_asthma','betaByr_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r') ) ,name='LambdaY' ) ,

# Loading from latent endogenous to manifest MZ
  mxMatrix('Full',4,6,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(T,8)  ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(.001,8) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','betaSex_asthma','betaSex_adhd',
              'betaByr_asthma','betaByr_adhd','betaByr_asthma','betaByr_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb') ) ,name='LambdaYMZ' ) ,

# Loading from latent endogenous to manifest DZ
  mxMatrix('Full',4,7,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(T,2),rep(F,2),rep(T,4), rep(F,2),rep(T,2) ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(.001,2),rep(0,2),rep(.001,4) ,rep(0,2),rep(.001,2) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','na','na',
              'betaByr_asthma','betaByr_adhd','betaByr_asthma','betaByr_adhd',
              'na','na','betaSex_asthma','betaSex_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r') ) ,name='LambdaYDZ' ) ,


# Covariance of latent endogenous
  mxMatrix('Symm', 8 , 8 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r')) , name= 'Psi' ) ,
  mxMatrix('Symm', 6 , 6 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb')) , name= 'PsiMZ' ) ,
  mxMatrix('Symm', 7 , 7 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r')) , name= 'PsiDZ' ) ,

# Covariance of latent exogenuous
# Put together
  mxAlgebra( rbind(cbind(PWTN,PBTWmz),cbind(PBTWmz,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMZ' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWdz),cbind(PBTWdz,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiDZ' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWfu),cbind(PBTWfu,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiFU' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWmh),cbind(PBTWmh,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMH' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWph),cbind(PBTWph,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiPH' ),


# Paths for A, C and E
  mxMatrix('Stand', 2 , 2  ,free=T,values=0,lbound=-.99,ubound=.99,
    labels=c('a12') , name='Ac' ),
  mxMatrix('Stand', 2 , 2  ,free=T,values=0,lbound=-.99,ubound=.99,
    labels=c('d12') , name='Dc' ),
  mxMatrix('Stand', 2 , 2  ,free=T,values=0,lbound=-.99,ubound=.99,
    labels=c('c12') , name='Cc' ),
  mxMatrix('Stand', 2 , 2  ,free=T,values=0,lbound=-.99,ubound=.99,
    labels=c('e12') , name='Ec' ),
  mxMatrix('Diag', 2 , 2  ,free=T,values=.6,lbound=0,ubound=.99,
    labels=c('a11','a22'),name='Av'),
  mxMatrix('Diag', 2 , 2  ,free=T,values=.5,lbound=0,ubound=.99,
    labels=c('d11','d22'),name='Dv'),
  mxMatrix('Diag', 2 , 2  ,free=T,values=.4,lbound=0,ubound=.99,
    labels=c('c11','c22'),name='Cv'),
  mxAlgebra( Av%*%Ac%*%Av , name='A2' ),
  mxAlgebra( Dv%*%Dc%*%Dv , name='D2' ),
  mxAlgebra( Cv%*%Cc%*%Cv , name='C2' ),
# Algebra to make E
  mxMatrix('Full' , 2 , 1  ,free=F,values=1,name='Ones21'),
  mxAlgebra( vec2diag( sqrt( Ones21-diag2vec(A2+D2+C2) ) ) , name='Ev' ),
  mxAlgebra( Ev%*%Ec%*%Ev , name='E2' ),
# Expected covariance matrices
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , A2+D2+C2 ) , cbind( A2+D2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMZ'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .5%x%A2+.25%x%D2+C2 ) , cbind( .5%x%A2+.25%x%D2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiDZ'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .5%x%A2+.25%x%D2+C2 ) , cbind( .5%x%A2+.25%x%D2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiFU'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .25%x%A2+C2 ) , cbind( .25%x%A2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMH'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .25%x%A2 ) , cbind( .25%x%A2 , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiPH'  ),
  
# Loadings from latent exogenous to endogenous
  mxMatrix('Diag',4,4,free=F,values=1,name='GaDiag'),
  mxMatrix('Zero',4,4,name='Ga0'),
  mxMatrix('Zero',2,4,name='Ga0MZ'),
  mxMatrix('Zero',3,4,name='Ga0DZ'),
  mxAlgebra(rbind( GaDiag , Ga0 ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='Ga' ),
  mxAlgebra(rbind( GaDiag , Ga0MZ ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='GaMZ' ),
  mxAlgebra(rbind( GaDiag , Ga0DZ ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='GaDZ' ),

# Latent variable equation
  mxMatrix('Diag',8,8,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r')) , name='Be' ),
  mxMatrix('Diag',6,6,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb')) , name='BeMZ' ),
  mxMatrix('Diag',7,7,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r')) , name='BeDZ' ),

# Residual Variance of manifest variable
  mxMatrix('Zero',4,4,free=F,
    dimnames=list(varNames,varNames) , name='ThetaE' ),

# Matrices for non-loading onto latent Exogenous
  mxMatrix('Full',nrow=4,ncol=1,free=F,values=0,dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),NA) , name='Ka' ),
  mxMatrix('Full',nrow=0,ncol=4,free=F,dimnames=list(NULL,c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) ,name='Lx' ),
  mxMatrix('Full',nrow=0,ncol=0,free=F,dimnames=list(NULL,NULL) , name='Td'),
  mxMatrix('Full',nrow=0,ncol=1,free=F,dimnames=list(NULL,NULL) , name='Tx' ),
  mxMatrix('Full',nrow=0,ncol=4,free=F,dimnames=list(NULL,varNames) , name='Th')
)

# MZ
mzMod <- mxModel( 'MZ' ,
# Data
  mxData( dat[dat$fam==1,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',6,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaYMZ', PS='Base.PsiMZ',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.BeMZ', 
    LX='Base.Lx',GA='Base.GaMZ',  PH='Base.PhiMZ',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# DZ
dzMod <- mxModel( 'DZ' ,
# Data
  mxData( dat[dat$fam==2,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',7,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaYDZ', PS='Base.PsiDZ',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.BeDZ', 
    LX='Base.Lx',GA='Base.GaDZ',  PH='Base.PhiDZ',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# Full
fuMod <- mxModel( 'FU' ,
# Data
  mxData( dat[dat$fam==3,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiFU',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# MH
mhMod <- mxModel( 'MH' ,
# Data
  mxData( dat[dat$fam==4,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiMH',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrH',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# PH
phMod <- mxModel( 'PH' ,
# Data
  mxData( dat[dat$fam==5,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiPH',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrH',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

#
modADCE <- mxModel( 'ADCE' ,
  baseMod , mzMod , dzMod , fuMod , mhMod , phMod ,
  mxFitFunctionMultigroup( c( 'MZ','DZ','FU','MH','PH' ) )
)
# Fit model
modADCEFit <- mxRun( modADCE , silent=T )


### Output
out <- c(
# Phenotypic correlation
mxEval( (Base.A2+Base.D2+Base.C2+Base.E2)[1,2] , modADCEFit ),

# Univariate heritabilities
mxEval( Base.A2[1,1] , modADCEFit ),
mxEval( Base.A2[2,2] , modADCEFit ),
mxEval( Base.D2[1,1] , modADCEFit ),
mxEval( Base.D2[2,2] , modADCEFit ),
mxEval( Base.A2[1,1]+Base.D2[1,1] , modADCEFit ),
mxEval( Base.A2[2,2]+Base.D2[2,2] , modADCEFit ),

# Univariate C
mxEval( Base.C2[1,1] , modADCEFit ),
mxEval( Base.C2[2,2] , modADCEFit ),

# Univariate E
mxEval( Base.E2[1,1] , modADCEFit ),
mxEval( Base.E2[2,2] , modADCEFit ),

# ADCE-correlations
mxEval( Base.Ac[1,2] , modADCEFit ),
mxEval( Base.Dc[1,2] , modADCEFit ),
mxEval( cov2cor(Base.A2+Base.D2)[1,2] , modADCEFit ),
mxEval( Base.Cc[1,2] , modADCEFit ),
mxEval( Base.Ec[1,2] , modADCEFit ),

# ADCE-explained covariance
mxEval( Base.A2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modADCEFit ),
mxEval( Base.D2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modADCEFit ),
mxEval( (Base.A2+Base.D2)[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modADCEFit ),
mxEval( Base.C2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modADCEFit ),
mxEval( Base.E2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modADCEFit )
)
return( out )
break()
}

##########################################################################




##########################################################################
### ACE-model


aceFun <- function( dat=dat ){
### Model
baseMod <- mxModel( 'Base' ,
# Manifest endogenous means
  mxMatrix('Full', vlen , 1 , free=F , values=0,
    dimnames=list(varNames,NA) , name='eMean' ),

# Manifest endogenous thresholds
  mxMatrix('Full', 1 , 4 ,free=T,values=thrstart,labels=paste(thrLab,'_F',sep=''),
    dimnames=list(NA,varNames) , name='ThrF' ),
  mxMatrix('Full', 1 , 4 ,free=T,values=thrstart,labels=paste(thrLab,'_H',sep=''),
    dimnames=list(NA,varNames) , name='ThrH' ),

# Loading from latent endogenous to manifest
  mxMatrix('Full',4,8,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(c(rep(T,2),rep(F,2)),2) , rep(c(rep(F,2),rep(T,2)),2) ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(c(rep(.001,2),rep(0,2)),2) , rep(c(rep(0,2),rep(.001,2)),2) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','na','na',
              'betaByr_asthma','betaByr_adhd','na','na',
              'na','na','betaSex_asthma','betaSex_adhd',
              'na','na','betaByr_asthma','betaByr_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r') ) ,name='LambdaY' ) ,

# Loading from latent endogenous to manifest MZ
  mxMatrix('Full',4,6,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(T,8)  ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(.001,8) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','betaSex_asthma','betaSex_adhd',
              'betaByr_asthma','betaByr_adhd','betaByr_asthma','betaByr_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb') ) ,name='LambdaYMZ' ) ,

# Loading from latent endogenous to manifest DZ
  mxMatrix('Full',4,7,
    free=  c( rep(F,0)      , rep(F,1) , rep(F,3) ,
              rep(F,1)      , rep(F,1) , rep(F,2) ,
              rep(F,2)      , rep(F,1) , rep(F,1) ,
              rep(F,3)      , rep(F,1) , rep(F,0) ,
              rep(T,2),rep(F,2),rep(T,4), rep(F,2),rep(T,2) ) , 
    values=c( rep(0,0)      , rep(1,1) , rep(0,3) ,
              rep(0,1)      , rep(1,1) , rep(0,2) ,
              rep(0,2)      , rep(1,1) , rep(0,1) ,
              rep(0,3)      , rep(1,1) , rep(0,0) ,
              rep(.001,2),rep(0,2),rep(.001,4) ,rep(0,2),rep(.001,2) ) , 
    labels=c( rep('Z',0)      , rep('O',1) , rep('Z',3) ,
              rep('Z',1)      , rep('O',1) , rep('Z',2) ,
              rep('Z',2)      , rep('O',1) , rep('Z',1) ,
              rep('Z',3)      , rep('O',1) , rep('Z',0) ,
              'betaSex_asthma','betaSex_adhd','na','na',
              'betaByr_asthma','betaByr_adhd','betaByr_asthma','betaByr_adhd',
              'na','na','betaSex_asthma','betaSex_adhd'),
    dimnames=list(varNames,c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r') ) ,name='LambdaYDZ' ) ,


# Covariance of latent endogenous
  mxMatrix('Symm', 8 , 8 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r')) , name= 'Psi' ) ,
  mxMatrix('Symm', 6 , 6 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb')) , name= 'PsiMZ' ) ,
  mxMatrix('Symm', 7 , 7 , free=F ,values=0,
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r')) , name= 'PsiDZ' ) ,

# Covariance of latent exogenuous
# Put together
  mxAlgebra( rbind(cbind(PWTN,PBTWmz),cbind(PBTWmz,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMZ' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWdz),cbind(PBTWdz,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiDZ' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWfu),cbind(PBTWfu,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiFU' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWmh),cbind(PBTWmh,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMH' ),
  mxAlgebra( rbind(cbind(PWTN,PBTWph),cbind(PBTWph,PWTN)) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiPH' ),


# Paths for A, C and E
  mxMatrix('Stand', 2 , 2  ,free=T,values=0,lbound=-.99,ubound=.99,
    labels=c('a12') , name='Ac' ),
  mxMatrix('Stand', 2 , 2  ,free=F,values=0,lbound=-.99,ubound=.99,
    labels=c('d12') , name='Dc' ),
  mxMatrix('Stand', 2 , 2  ,free=T,values=0,lbound=-.99,ubound=.99,
    labels=c('c12') , name='Cc' ),
  mxMatrix('Stand', 2 , 2  ,free=T,values=0,lbound=-.99,ubound=.99,
    labels=c('e12') , name='Ec' ),
  mxMatrix('Diag', 2 , 2  ,free=T,values=.6,lbound=0,ubound=.99,
    labels=c('a11','a22'),name='Av'),
  mxMatrix('Diag', 2 , 2  ,free=F,values=0,lbound=0,ubound=.99,
    labels=c('d11','d22'),name='Dv'),
  mxMatrix('Diag', 2 , 2  ,free=T,values=.4,lbound=0,ubound=.99,
    labels=c('c11','c22'),name='Cv'),
  mxAlgebra( Av%*%Ac%*%Av , name='A2' ),
  mxAlgebra( Dv%*%Dc%*%Dv , name='D2' ),
  mxAlgebra( Cv%*%Cc%*%Cv , name='C2' ),
# Algebra to make E
  mxMatrix('Full' , 2 , 1  ,free=F,values=1,name='Ones21'),
  mxAlgebra( vec2diag( sqrt( Ones21-diag2vec(A2+D2+C2) ) ) , name='Ev' ),
  mxAlgebra( Ev%*%Ec%*%Ev , name='E2' ),
# Expected covariance matrices
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , A2+D2+C2 ) , cbind( A2+D2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMZ'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .5%x%A2+.25%x%D2+C2 ) , cbind( .5%x%A2+.25%x%D2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiDZ'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .5%x%A2+.25%x%D2+C2 ) , cbind( .5%x%A2+.25%x%D2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiFU'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .25%x%A2+C2 ) , cbind( .25%x%A2+C2  , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiMH'  ),
  mxAlgebra( rbind(cbind( A2+D2+C2+E2 , .25%x%A2 ) , cbind( .25%x%A2 , A2+D2+C2+E2 ) ) ,
    dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='PhiPH'  ),
  
# Loadings from latent exogenous to endogenous
  mxMatrix('Diag',4,4,free=F,values=1,name='GaDiag'),
  mxMatrix('Zero',4,4,name='Ga0'),
  mxMatrix('Zero',2,4,name='Ga0MZ'),
  mxMatrix('Zero',3,4,name='Ga0DZ'),
  mxAlgebra(rbind( GaDiag , Ga0 ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='Ga' ),
  mxAlgebra(rbind( GaDiag , Ga0MZ ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='GaMZ' ),
  mxAlgebra(rbind( GaDiag , Ga0DZ ) ,
    dimnames=list( c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r') , c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) , name='GaDZ' ),

# Latent variable equation
  mxMatrix('Diag',8,8,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r')) , name='Be' ),
  mxMatrix('Diag',6,6,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb')) , name='BeMZ' ),
  mxMatrix('Diag',7,7,free=F,values=0, 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'),c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r')) , name='BeDZ' ),

# Residual Variance of manifest variable
  mxMatrix('Zero',4,4,free=F,
    dimnames=list(varNames,varNames) , name='ThetaE' ),

# Matrices for non-loading onto latent Exogenous
  mxMatrix('Full',nrow=4,ncol=1,free=F,values=0,dimnames=list(c(paste0('Phi1_',1:2),paste0('Phi2_',1:2)),NA) , name='Ka' ),
  mxMatrix('Full',nrow=0,ncol=4,free=F,dimnames=list(NULL,c(paste0('Phi1_',1:2),paste0('Phi2_',1:2))) ,name='Lx' ),
  mxMatrix('Full',nrow=0,ncol=0,free=F,dimnames=list(NULL,NULL) , name='Td'),
  mxMatrix('Full',nrow=0,ncol=1,free=F,dimnames=list(NULL,NULL) , name='Tx' ),
  mxMatrix('Full',nrow=0,ncol=4,free=F,dimnames=list(NULL,varNames) , name='Th')
)

# MZ
mzMod <- mxModel( 'MZ' ,
# Data
  mxData( dat[dat$fam==1,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',6,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaYMZ', PS='Base.PsiMZ',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.BeMZ', 
    LX='Base.Lx',GA='Base.GaMZ',  PH='Base.PhiMZ',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# DZ
dzMod <- mxModel( 'DZ' ,
# Data
  mxData( dat[dat$fam==2,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',7,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaYDZ', PS='Base.PsiDZ',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.BeDZ', 
    LX='Base.Lx',GA='Base.GaDZ',  PH='Base.PhiDZ',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# Full
fuMod <- mxModel( 'FU' ,
# Data
  mxData( dat[dat$fam==3,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiFU',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrF',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# MH
mhMod <- mxModel( 'MH' ,
# Data
  mxData( dat[dat$fam==4,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiMH',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrH',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

# PH
phMod <- mxModel( 'PH' ,
# Data
  mxData( dat[dat$fam==5,] , type='raw' ),

# Endogenous latent variable means
  mxMatrix('Full',8,1,free=F,values=0,
    labels=c(rep(paste0('am',1:2),2),'data.sex','data.yob','data.sex_rel','data.yob_rel') , 
    dimnames=list(c(paste0('Psi1_',1:2),paste0('Psi2_',1:2),'Sb','Bb','Sb_r','Bb_r'), NA) , name='Alpha' ),

# Expectation in LISREL setup
  mxExpectationLISREL(LY='Base.LambdaY', PS='Base.Psi',TE='Base.ThetaE',TY='Base.eMean', AL = 'Alpha',  BE='Base.Be', 
    LX='Base.Lx',GA='Base.Ga',  PH='Base.PhiPH',  TD='Base.Td', TH='Base.Th',TX = 'Base.Tx', KA = 'Base.Ka' ,
    thresholds='Base.ThrH',dimnames=varNames , threshnames=varNames ),
  
  mxFitFunctionWLS( type='WLS' ) 
)

#
modACE <- mxModel( 'ACE' ,
  baseMod , mzMod , dzMod , fuMod , mhMod , phMod ,
  mxFitFunctionMultigroup( c( 'MZ','DZ','FU','MH','PH' ) )
)
# Fit model
modACEFit <- mxRun( modACE , silent=T )

### Output
out <- c(
# Phenotypic correlation
mxEval( (Base.A2+Base.D2+Base.C2+Base.E2)[1,2] , modACEFit ),

# Univariate heritabilities
mxEval( Base.A2[1,1] , modACEFit ),
mxEval( Base.A2[2,2] , modACEFit ),
mxEval( Base.D2[1,1] , modACEFit ),
mxEval( Base.D2[2,2] , modACEFit ),
mxEval( Base.A2[1,1]+Base.D2[1,1] , modACEFit ),
mxEval( Base.A2[2,2]+Base.D2[2,2] , modACEFit ),

# Univariate C
mxEval( Base.C2[1,1] , modACEFit ),
mxEval( Base.C2[2,2] , modACEFit ),

# Univariate E
mxEval( Base.E2[1,1] , modACEFit ),
mxEval( Base.E2[2,2] , modACEFit ),

# ACE-correlations
mxEval( Base.Ac[1,2] , modACEFit ),
mxEval( Base.Dc[1,2] , modACEFit ),
mxEval( cov2cor(Base.A2+Base.D2)[1,2] , modACEFit ),
mxEval( Base.Cc[1,2] , modACEFit ),
mxEval( Base.Ec[1,2] , modACEFit ),

# ACE-explained covariance
mxEval( Base.A2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modACEFit ),
mxEval( Base.D2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modACEFit ),
mxEval( (Base.A2+Base.D2)[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modACEFit ),
mxEval( Base.C2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modACEFit ),
mxEval( Base.E2[1,2]/((Base.A2+Base.D2+Base.C2+Base.E2)[1,2]) , modACEFit )
)
out
}



##########################################################################
##########################################################################
########################## END OF FILE ###################################
##########################################################################
##########################################################################