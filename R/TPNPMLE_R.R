#' @title Penalized Non-Parametric Maximum-Likelihood Estimation (PNPMLEs) for Cohort Samplings with Time Matching under Cox's Regression Model
#' @description The function utilizes a self-consistency iterative algorithm to calculate PNPMLEs by adding penalty function for cohort samplings with time matching under Cox's regression model. In addition to compute PNPMLEs, it can also estimate asymptotic varance, as described in Wang et al. (2019+). The Cox's regression model is \deqn{\lambda(t|z)=\lambda_{0}(t)\exp(z^T\beta).}
#' @param data The description is the same as the statement of \code{TNPMLE} function.
#' @param iteration1 The number of iteration for computing (P)NPMLEs.
#' @param iteration2 The number of iteration for computing profile likelihoods which are used to estimate asymptotic variance.
#' @param converge The description is the same as the statement of \code{TNPMLE} function.
#' @param penalty The choice of penalty, it can be SCAD, HARD or LASSO.
#' @param penaltytuning The tuning parameter for penalty function, it is a sequence of numeric vector.
#' @param fold The \code{fold} information for cross validation. Without loss of generality, we note that \code{fold} value have to be bigger than one (>1) and cohort size is divisible by \code{fold} value. However, if cohort size is not able to be divided, we are going to partition off cohort into several suitable parts according to \code{fold} value automaticly for cross-validation. 
#' @param cut The cut point. When \eqn{\hat{\beta}_j} is smaller than the cut point, we set \eqn{\hat{\beta}_j} be zero, i.e. remove the corresponding covariate from our model to do variable selection.
#' @param seed The seed of the random number generator to obtain reproducible results.
#' @return Returns a list with components
#' @return \item{num}{The numbers of case and observed subjects.}
#' @return \item{iloop}{The final number of iteration for computing PNPMLEs.}
#' @return \item{diff}{The sup-norm distance between the last two iterations of the estimates of the relative risk coefficients.}
#' @return \item{cvl}{The cross-validated profile log-likelihood.}
#' @return \item{tuning}{The suitable tuning parameter, such that the maximum of cross-validated profile log-likelihood is attained.}
#' @return \item{likelihood}{The log likelihood value of PNPMLEs.}
#' @return \item{pnpmle}{The estimated regression coefficients with their corresponding estimated standard errors and p-values.}
#' @return \item{Lpnpmle}{The estimated cumulative baseline hazards function.}
#' @return \item{Ppnpmle}{The empirical distribution of covariates which are missing for unobserved subjects.}
#' @return \item{elements}{The description is the same as the statement of \code{TNPMLE} function.}
#' @return \item{Adata}{The description is the same as the statement of \code{TNPMLE} function.}
#' @references Wang JH, Pan CH, Chang IS*, and Hsiung CA (2019) Penalized full likelihood approach to variable selection for Cox's regression model under nested case-control sampling. published in Lifetime Data Analysis <doi:10.1007/s10985-019-09475-z>.
#' @note The missing value (NA) in the DATA is not allowed in this version.
#' @seealso See \code{\link{TNPMLE}}.
#' @import MASS splines survival stats
#' @export
#' @examples
#' set.seed(100)
#' library(splines)
#' library(survival)
#' library(MASS)
#' beta=c(1,0)
#' lambda=0.3
#' cohort=100
#' covariate=2+length(beta)
#' z=matrix(rnorm(cohort*length(beta)),nrow=cohort)
#' rate=1/(runif(cohort,1,3)*exp(z%*%beta))
#' c=rexp(cohort,rate)
#' u=-log(runif(cohort,0,1))/(lambda*exp(z%*%beta))
#' time=apply(cbind(u,c),1,min)
#' status=(u<=c)+0
#' casenum=sum(status)
#' odata=cbind(time,status,z)
#' odata=data.frame(odata)
#' a=order(status)
#' data=matrix(0,cohort,covariate)
#' data=data.frame(data)
#' for (i in 1:cohort){
#' data[i,]=odata[a[cohort-i+1],]
#' }
#' ncc=matrix(0,cohort,covariate)
#' ncc=data.frame(data)
#' aa=order(data[1:casenum,1])
#' for (i in 1:casenum){
#' ncc[i,]=data[aa[i],]
#' }
#' control=1
#' q=matrix(0,casenum,control)
#' for (i in 1:casenum){
#' k=c(1:cohort)
#' k=k[-(1:i)]
#' sumsc=sum(ncc[i,1]<ncc[,1][(i+1):cohort])
#' if (sumsc==0) {
#' 			q[i,]=c(1)
#' } else {
#' 			q[i,]=sample(k[ncc[i,1]<ncc[,1][(i+1):cohort]],control)
#' }
#' }
#' cacon=c(q,1:casenum)
#' k=c(1:cohort)
#' owf=k[-cacon]
#' wt=k[-owf]
#' owt=k[-wt]
#' ncct=matrix(0,cohort,covariate)
#' ncct=data.frame(ncct)
#' for (i in 1:length(wt)){
#' ncct[i,]=ncc[wt[i],]
#' }
#' for (i in 1:length(owt)){
#' ncct[length(wt)+i,]=ncc[owt[i],]
#' }
#' d=length(wt)+1
#' ncct[d:cohort,3:covariate]=-9
#' TPNPMLEtest=TPNPMLE(ncct,100,30,0,"SCAD",seq(0.10,0.13,0.005),2,1e-05,1)

###################################################################
#              CREATE R PACKAGE for NCCPNPMLE with PENALTY       #
###################################################################


TPNPMLE=function(data,iteration1,iteration2,converge,penalty,penaltytuning,fold,cut,seed){

ARRANGEDATA=function(data){
cohort=dim(data)[1]
covariate=dim(data)[2]
casenum=sum(data[,2])
nonobsnum=sum(data[,covariate]==-9)
nonobsindex=which(data[,covariate]==-9)
obsnum=sum(data[,covariate]!=-9)
obsindex=which(data[,covariate]!=-9)
case1=matrix(0,cohort,covariate)
case1=data.frame(case1)
for (i in 1:obsnum){
case1[i,]=data[obsindex[i],]
}
for (i in (obsnum+1):cohort){
case1[i,]=data[nonobsindex[i-obsnum],]
}
case2=matrix(0,cohort,covariate)
case2=data.frame(case2)
case2=case1
a=order(case1[1:obsnum,2])
for (i in 1:obsnum){
case2[i,]=case1[a[obsnum-i+1],]
}
case3=matrix(0,cohort,covariate)
case3=data.frame(case3)
case3=case2
aa=order(case2[,1][1:casenum])
for (i in 1:casenum){
case3[i,]=case2[aa[i],]
}
ARRANGEDATAlist=list(arrangedata=case3)
return(ARRANGEDATAlist)
}

#################################################################
#################################################################

ncct=data.frame(ARRANGEDATA(data)$arrangedata)
expand=function(yy) c(yy%*%t(yy))
pi=1e-08
casenum=sum(ncct[,2])
cohort=dim(ncct)[1]
covariate=dim(ncct)[2]
obsnum=sum(ncct[,covariate]!=-9)
nonobsnum=cohort-obsnum
d=obsnum+1
numbeta=covariate-2

w=matrix(0,obsnum,numbeta)
wextend=matrix(0,obsnum,(numbeta*numbeta))
for (j in 1:obsnum){
w[j,]=as.numeric(ncct[j,3:covariate])
wextend[j,]=matrix(w[j,] %*% t(w[j,]),1,(numbeta*numbeta))
}
windex=(!duplicated(w))
numw=sum(windex)
w=as.matrix(w[windex,])
wextend=as.matrix(wextend[windex,])

zf=matrix(0,cohort,numbeta)
zfextend=matrix(0,cohort,numbeta*numbeta)
for (j in 1:obsnum){
zf[j,]=as.numeric(ncct[j,3:covariate])
zfextend[j,]=matrix(zf[j,] %*% t(zf[j,]),1,(numbeta*numbeta))
}

casetime=ncct[,1][1:casenum]
                       
                       ######################################################################
                       #                          Cross validation                          #
                       ######################################################################
set.seed(seed)
acv=sample(sample(sample(cohort)))
if (cohort/fold==round(cohort/fold)) division=1 else division=0
if (division==1) {
basesize=cohort/fold
subcohort=basesize*(fold-1)
vindex=matrix(0,fold,subcohort)
vindex[1,]=c(1:subcohort)
vindex[2,]=c((basesize+1):cohort)
if (fold>2) {
for (i in 2:(fold-1)){
vindex[(i+1),]=c((basesize*i+1):cohort,1:(basesize*(i-1)))
}
} 
cvchhat=matrix(0,subcohort,fold)
} else {
basesize=ceiling(cohort/fold)
subcohort=basesize*(fold-1)
vindex=matrix(0,fold,subcohort)
extra=basesize-(cohort-subcohort)
vindex[1,]=c(1:subcohort)
vindex[2,]=c((basesize+1):cohort,rep(0,extra))
if (fold>2){
for (i in 2:(fold-1)){
vindex[(i+1),]=c((basesize*i+1):cohort,1:(basesize*(i-1)),rep(0,extra))
}
}
cvchhat=matrix(0,subcohort,fold)
}
cvbhat=matrix(0,numbeta,fold)

for (v in 1:fold){ 
subcohort=sum(vindex[v,]!=0)
subncct=ncct[sort(acv[vindex[v,]]),]
subwt=sum((subncct[,covariate]!=-9))
subowt=subcohort-subwt
subcasenum=sum(subncct[,2])
subd=subwt+1

subw=matrix(0,subwt,numbeta)
subwextend=matrix(0,subwt,(numbeta*numbeta))
for (j in 1:subwt){
subw[j,]=as.numeric(subncct[j,3:covariate])
subwextend[j,]=matrix(subw[j,] %*% t(subw[j,]),1,(numbeta*numbeta))
}
subwindex=(!duplicated(subw))
numsubw=sum(subwindex)
subw=as.matrix(subw[subwindex,])
subwextend=as.matrix(subwextend[subwindex,])

subzf=matrix(0,subcohort,numbeta)
subzfextend=matrix(0,subcohort,numbeta*numbeta)
for (j in 1:subwt){
subzf[j,]=as.numeric(subncct[j,3:covariate])
subzfextend[j,]=matrix(subzf[j,] %*% t(subzf[j,]),1,(numbeta*numbeta))
}

subcasetime=subncct[,1][1:subcasenum]

############################################################################
############################################################################ initial value

subncct=as.matrix(subncct)
covariatedata=matrix(0,subwt,numbeta)
covariatedata=subncct[1:subwt,3:covariate]
subncct=data.frame(subncct)
us=coxph(Surv(subncct[1:subwt,1],subncct[1:subwt,2])~covariatedata)

inibeta=matrix(c(us$coef),numbeta,1)
inip=c(rep(1/numsubw,numsubw)) 

subcf=c(rep(0,subcasenum))
new=matrix(exp(subzf[1:subwt,]%*%inibeta),subwt,1)
subcf=1/(((t(matrix(rep(subncct[1:subwt,1],subcasenum),subwt,subcasenum))>=subncct[1:subcasenum,1])+0)%*%new)

inich=c(rep(0,subcohort))
inich[1:subcasenum]=cumsum(subcf)
subhh=apply(t(matrix(rep(subncct[1:subcasenum,1],(subcohort-subcasenum)),subcasenum,(subcohort-subcasenum)))<=subncct[(subcasenum+1):subcohort,1],1,sum)+1
inich[(subcasenum+1):subcohort]=c(0,inich[1:subcasenum])[subhh]

##############################################################
subr3=matrix(1,1,subwt)
subr4=matrix(1,1,subowt)

alpha=matrix(0,subcohort,numsubw)
if (numsubw==subwt) {
 diag(alpha[1:numsubw,1:numsubw])=1
} else {
for (i in 1:subwt){
for (j in 1:numsubw){
if (sum(subzf[i,]==subw[j,])==numbeta) alpha[i,j]=1 else alpha[i,j]=0
}
}
}
subtmC=(t(matrix(rep(subncct[,1],subcasenum),subcohort,subcasenum))>=subcasetime)+0

for (kkk in 1:iteration1) {

wb=matrix(exp(exp(subw%*%inibeta)%o%(-inich[subd:subcohort])),numsubw,subowt)
alpha[subd:subcohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

########################################################

bhat=matrix(0,numbeta,1)
order1a=t(subr3%*%(subzf*c(subncct[,2]-inich*exp(subzf%*%inibeta)))[1:subwt,])
order1b=t(subr4%*%((alpha[subd:subcohort,]*-inich[subd:subcohort])%*%(subw*c(exp(subw%*%inibeta)))))
order1=order1a+order1b
order2a=matrix(subr3%*%((inich*exp(subzf%*%inibeta))[1:subwt,1]*subzfextend[1:subwt,]),numbeta,numbeta)
order2b=matrix(subr4%*%((alpha[subd:subcohort,]*inich[subd:subcohort])%*%(subwextend*c(exp(subw%*%inibeta)))),numbeta,numbeta)
o2c11=t(wb)%*%(subwextend*c(exp(subw%*%inibeta)^2*phat))
o2c12=t(wb)%*%phat
o2c21=matrix(t(wb)%*%(subw*c(exp(subw%*%inibeta)*phat)),subowt,numbeta)
o2c1=o2c11*c(o2c12)*-inich[subd:subcohort]

if (numbeta==1){
o2c2=c((o2c21^2)*inich[subd:subcohort])
} else {
o2c2=t(apply(o2c21,1,expand))*inich[subd:subcohort]
}

order2c=matrix(subr4%*%(inich[subd:subcohort]*((o2c1+o2c2)/c(o2c12^2))),numbeta,numbeta)
order2=-(order2a+order2b+order2c)
bhat=inibeta-ginv(order2)%*%order1

#################################################################
diff=max(abs(bhat-inibeta))
if (diff<converge) break
#################################################################

sigma21=matrix(0,subcohort,1)
sigma21[1:subwt,1]=(exp(subzf%*%bhat))[1:subwt]
sigma21[subd:subcohort,]=alpha[subd:subcohort,]%*%exp(subw%*%bhat)
f=1/(subtmC%*%sigma21)
chhat=c(rep(0,subcohort))
chhat[1:subcasenum]=cumsum(f)
chhat[(subcasenum+1):subcohort]=c(0,chhat[1:subcasenum])[subhh]

inibeta=bhat
inip=phat
inich=chhat
}
subbhatnpmle=bhat
subchhatnpmle=chhat
subphatnpmle=phat

cvbhat[,v]=subbhatnpmle
extra=dim(cvchhat)[1]-subcohort
cvchhat[,v]=c(subchhatnpmle,rep(0,extra))
}

###################################################################################
################################################################################### CV

nccnpmle=TNPMLE(data,iteration1=300,iteration2=30,converge=0)
bhatnpmle=matrix(c(nccnpmle$npmle),numbeta,1)
phatnpmle=nccnpmle$Pnpmle
chhatnpmle=nccnpmle$Lambdanpmle

tuning=c(rep(0,numbeta))
cv=c(rep(0,length(penaltytuning)))

for (m in 1:length(penaltytuning)){
tuning=penaltytuning[m]*c(rep(1,numbeta))
subbhatpnpmle=matrix(0,numbeta,fold)

for (v in 1:fold){ 
subcohort=sum(vindex[v,]!=0)
subncct=ncct[sort(acv[vindex[v,]]),]
subwt=sum((subncct[,covariate]!=-9))
subowt=subcohort-subwt
subcasenum=sum(subncct[,2])
subd=subwt+1

subw=matrix(0,subwt,numbeta)
subwextend=matrix(0,subwt,(numbeta*numbeta))
for (j in 1:subwt){
subw[j,]=as.numeric(subncct[j,3:covariate])
subwextend[j,]=matrix(subw[j,] %*% t(subw[j,]),1,(numbeta*numbeta))
}
subwindex=(!duplicated(subw))
numsubw=sum(subwindex)
subw=as.matrix(subw[subwindex,])
subwextend=as.matrix(subwextend[subwindex,])

subzf=matrix(0,subcohort,numbeta)
subzfextend=matrix(0,subcohort,numbeta*numbeta)
for (j in 1:subwt){
subzf[j,]=as.numeric(subncct[j,3:covariate])
subzfextend[j,]=matrix(subzf[j,] %*% t(subzf[j,]),1,(numbeta*numbeta))
}

subcasetime=subncct[,1][1:subcasenum]

###########################################################################################

inibeta=matrix(c(cvbhat[,v]),numbeta,1) 
inip=c(rep(1/numsubw,numsubw)) 
inich=c(cvchhat[,v][vindex[v,]!=0])   

if (penalty=="HARD"){
error=(pi/(2*(subcohort)*max(2*tuning)))*min(abs(inibeta))
} else {
error=(pi/(2*(subcohort)*max(tuning)))*min(abs(inibeta))
}

subr3=matrix(1,1,subwt)
subr4=matrix(1,1,subowt)

alpha=matrix(0,subcohort,numsubw)
if (numsubw==subwt) {
 diag(alpha[1:numsubw,1:numsubw])=1
} else {
for (i in 1:subwt){
for (j in 1:numsubw){
if (sum(subzf[i,]==subw[j,])==numbeta) alpha[i,j]=1 else alpha[i,j]=0
}
}
}

subhh=apply(t(matrix(rep(subncct[1:subcasenum,1],(subcohort-subcasenum)),subcasenum,(subcohort-subcasenum)))<=subncct[(subcasenum+1):subcohort,1],1,sum)+1
subtmC=(t(matrix(rep(subncct[,1],subcasenum),subcohort,subcasenum))>=subcasetime)+0

for (kkk in 1:iteration1) {
wb=matrix(exp(exp(subw%*%inibeta)%o%(-inich[subd:subcohort])),numsubw,subowt)
alpha[subd:subcohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

#########################################################

bhat=matrix(0,numbeta,1)
order1a=t(subr3%*%(subzf*c(subncct[,2]-inich*exp(subzf%*%inibeta)))[1:subwt,])
order1b=t(subr4%*%((alpha[subd:subcohort,]*-inich[subd:subcohort])%*%(subw*c(exp(subw%*%inibeta)))))
order1=order1a+order1b
order2a=matrix(subr3%*%((inich*exp(subzf%*%inibeta))[1:subwt,1]*subzfextend[1:subwt,]),numbeta,numbeta)
order2b=matrix(subr4%*%((alpha[subd:subcohort,]*inich[subd:subcohort])%*%(subwextend*c(exp(subw%*%inibeta)))),numbeta,numbeta)
o2c11=t(wb)%*%(subwextend*c(exp(subw%*%inibeta)^2*phat))
o2c12=t(wb)%*%phat
o2c21=matrix(t(wb)%*%(subw*c(exp(subw%*%inibeta)*phat)),subowt,numbeta)
o2c1=o2c11*c(o2c12)*-inich[subd:subcohort]

if (numbeta==1){
o2c2=c((o2c21^2)*inich[subd:subcohort])
} else {
o2c2=t(apply(o2c21,1,expand))*inich[subd:subcohort]
}

order2c=matrix(subr4%*%(inich[subd:subcohort]*((o2c1+o2c2)/c(o2c12^2))),numbeta,numbeta)
order2=-(order2a+order2b+order2c)

if (penalty=="SCAD"){
index=matrix(0,numbeta,1)
for (i in 1:numbeta){
value=3.7*tuning[i]-abs(inibeta[i,])
if (value>0) index[i,]=value else index[i,]=0
}
deri=tuning*(ifelse(abs(inibeta)<=tuning,1,0))+(ifelse(abs(inibeta)>tuning,1,0))*(index/(3.7-1))
} else if (penalty=="HARD") {
deri=(-2*abs(inibeta)+2*tuning)*(ifelse(abs(inibeta)<=tuning,1,0))
} else {
deri=tuning
}

abe=c(deri/(abs(inibeta)+error))
abeta=diag(abe,numbeta)

bhat=inibeta-ginv(order2-subcohort*abeta)%*%(order1-subcohort*abeta%*%inibeta)

#################################################################
diff=max(abs(bhat-inibeta))
if (diff<converge) break
#################################################################

sigma21=matrix(0,subcohort,1)
sigma21[1:subwt,1]=(exp(subzf%*%bhat))[1:subwt]
sigma21[subd:subcohort,]=alpha[subd:subcohort,]%*%exp(subw%*%bhat)
f=1/(subtmC%*%sigma21)
chhat=c(rep(0,subcohort))
chhat[1:subcasenum]=cumsum(f)
chhat[(subcasenum+1):subcohort]=c(0,chhat[1:subcasenum])[subhh]

inibeta=bhat
inip=phat
inich=chhat
}

###########################################################
for (i in 1:numbeta){
if (abs(bhat[i,])<=cut) bhat[i,]=0 else bhat[i,]=bhat[i,]
}
subbhatpnpmle[,v]=bhat
}

cvbhatpnpmle=matrix(0,numbeta,cohort)
for (v in 1:fold){
subcohort=sum(vindex[v,]!=0)
dex=sort(acv[-vindex[v,]])
for (i in 1:(cohort-subcohort)){
cvbhatpnpmle[,dex[i]]=subbhatpnpmle[,v]
}
}

###########################################################profile likelihood for CV

r3=matrix(1,1,obsnum)
r4=matrix(1,1,nonobsnum)
wr3=matrix(1,1,numw)

alpha=matrix(0,cohort,numw)
if (numw==obsnum) {
 diag(alpha[1:numw,1:numw])=1
} else {
for (i in 1:obsnum){
for (j in 1:numw){
if (sum(zf[i,]==w[j,])==numbeta) alpha[i,j]=1 else alpha[i,j]=0
}
}
}

inip=phatnpmle 
inich=chhatnpmle
hh=apply(t(matrix(rep(ncct[1:casenum,1],(cohort-casenum)),casenum,(cohort-casenum)))<=ncct[(casenum+1):cohort,1],1,sum)+1
tmC=(t(matrix(rep(ncct[,1],casenum),cohort,casenum))>=casetime)+0

for (kkk in 1:iteration2){

for (i in d:cohort){
alpha2=wr3%*%(exp(-inich[i]*exp(w %*% cvbhatpnpmle[,i]))*inip)
alpha1=exp(-inich[i]*exp(w %*% cvbhatpnpmle[,i]))*inip
alpha1=c(alpha1)
alpha[i,]=alpha1/c(alpha2)
}
phat=apply(alpha,2,mean)

sigma21=matrix(0,cohort,1)
for (i in 1:obsnum){
sigma21[i,1]=exp(zf[i,]%*%cvbhatpnpmle[,i])
}
for (i in d:cohort){
ee=wr3%*%(exp(w %*% cvbhatpnpmle[,i])*alpha[i,])
sigma21[i,]=c(ee)
}
f=c(1/(tmC%*%sigma21))
chhat=c(rep(0,cohort))
chhat[1:casenum]=cumsum(f)
chhat[(casenum+1):cohort]=c(0,chhat[1:casenum])[hh]

inip=phat
inich=chhat
}

##########################################################

baseline=diff(c(0,chhat[1:casenum]))
la=0
for (i in 1:obsnum){
la=la+log(as.numeric(exp(-chhat[i]*exp(zf[i,] %*% cvbhatpnpmle[,i]))))
}
lb=0
for (i in 1:casenum){
lb=lb+log(as.numeric(baseline[i]*exp(zf[i,] %*% cvbhatpnpmle[,i])))
}
lc=sum(log(alpha[1:obsnum,]%*%(phat)))
ld=0
for (i in d:cohort){
ee=as.numeric(wr3%*%(exp(-chhat[i]*exp(w %*% cvbhatpnpmle[,i]))*phat))
ld=ld+log(ee)
}
cv[m]=la+lb+lc+ld
}
suit=order(-cv)[1]
tuningover=penaltytuning[suit]*c(rep(1,numbeta))

###############################################################compute full data based on SCAD given optimal tuning parameter

inibeta=bhatnpmle
inip=phatnpmle 
inich=chhatnpmle

alpha=matrix(0,cohort,numw)
if (numw==obsnum) {
 diag(alpha[1:numw,1:numw])=1
} else {
for (i in 1:obsnum){
for (j in 1:numw){
if (sum(zf[i,]==w[j,])==numbeta) alpha[i,j]=1 else alpha[i,j]=0
}
}
}

if (penalty=="HARD"){
error=(pi/(2*(cohort)*max(2*tuningover)))*min(abs(inibeta))
} else {
error=(pi/(2*(cohort)*max(tuningover)))*min(abs(inibeta))
}

for (iloop in 1:iteration1) {

wb=matrix(exp(exp(w%*%inibeta)%o%(-inich[d:cohort])),numw,nonobsnum)
alpha[d:cohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

########################################################

bhat=matrix(0,numbeta,1)
order1a=t(r3%*%(zf*c(ncct[,2]-inich*exp(zf%*%inibeta)))[1:obsnum,])
order1b=t(r4%*%((alpha[d:cohort,]*-inich[d:cohort])%*%(w*c(exp(w%*%inibeta)))))
order1=order1a+order1b
order2a=matrix(r3%*%((inich*exp(zf%*%inibeta))[1:obsnum,1]*zfextend[1:obsnum,]),numbeta,numbeta)
order2b=matrix(r4%*%((alpha[d:cohort,]*inich[d:cohort])%*%(wextend*c(exp(w%*%inibeta)))),numbeta,numbeta)
o2c11=t(wb)%*%(wextend*c(exp(w%*%inibeta)^2*phat))
o2c12=t(wb)%*%phat
o2c21=matrix(t(wb)%*%(w*c(exp(w%*%inibeta)*phat)),nonobsnum,numbeta)
o2c1=o2c11*c(o2c12)*-inich[d:cohort]

if (numbeta==1){
o2c2=c((o2c21^2)*inich[d:cohort])
} else {
o2c2=t(apply(o2c21,1,expand))*inich[d:cohort]
}

order2c=matrix(r4%*%(inich[d:cohort]*((o2c1+o2c2)/c(o2c12^2))),numbeta,numbeta)
order2=-(order2a+order2b+order2c)

if (penalty=="SCAD"){
index=matrix(0,numbeta,1)
for (i in 1:numbeta){
value=3.7*tuningover[i]-abs(inibeta[i,])
if (value>0) index[i,]=value else index[i,]=0
}
deri=tuningover*(ifelse(abs(inibeta)<=tuningover,1,0))+(ifelse(abs(inibeta)>tuningover,1,0))*(index/(3.7-1))
} else if (penalty=="HARD"){
deri=(-2*abs(inibeta)+2*tuningover)*(ifelse(abs(inibeta)<=tuningover,1,0))
} else {
deri=tuningover
}

abe=c(deri/(abs(inibeta)+error))
abeta=diag(abe,numbeta)
bhat=inibeta-ginv(order2-cohort*abeta)%*%(order1-cohort*abeta%*%inibeta)

#################################################################
diff=max(abs(bhat-inibeta))
if (diff<converge) break
##########################################################################

sigma21=matrix(0,cohort,1)
sigma21[1:obsnum,1]=(exp(zf%*%bhat))[1:obsnum]
sigma21[d:cohort,]=alpha[d:cohort,]%*%exp(w%*%bhat)
f=c(1/(tmC%*%sigma21))
chhat=c(rep(0,cohort))
chhat[1:casenum]=cumsum(f)
chhat[(casenum+1):cohort]=c(0,chhat[1:casenum])[hh]

inibeta=bhat
inip=phat
inich=chhat
}

##########################################################

for (i in 1:numbeta){
if (abs(bhat[i,])<=cut) bhat[i,]=0 else bhat[i,]=bhat[i,]
}

###########################################################################

wb=matrix(exp(exp(w%*%bhat)%o%(-inich[d:cohort])),numw,nonobsnum)
alpha[d:cohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

########################################################################

sigma21=matrix(0,cohort,1)
sigma21[1:obsnum,1]=(exp(zf%*%bhat))[1:obsnum]
sigma21[d:cohort,]=alpha[d:cohort,]%*%exp(w%*%bhat)
f=c(1/(tmC%*%sigma21))
chhat=c(rep(0,cohort))
chhat[1:casenum]=cumsum(f)
chhat[(casenum+1):cohort]=c(0,chhat[1:casenum])[hh]

chhatpnpmle=chhat
phatpnpmle=phat
bhatpnpmle=bhat

############################################## compute asymptotic variance;likelihood4
##############################################

baseline=diff(c(0,chhatpnpmle[1:casenum]))
la=sum((-chhatpnpmle*exp(zf%*%bhatpnpmle))[1:obsnum,])
lb=sum(log(baseline*exp(zf[ncct[,2]==1,]%*%bhatpnpmle)))
lc=sum(log(alpha[1:obsnum,]%*%(phatpnpmle)))
ld=sum(log(t(matrix(exp(exp(w%*%bhatpnpmle)%o%(-chhatpnpmle[d:cohort])),numw,nonobsnum))%*%phatpnpmle))

likelihood4=la+lb+lc+ld
aic=-2*likelihood4+2*sum(bhatpnpmle[,1]!=0)
bic=-2*likelihood4+log(cohort)*sum(bhatpnpmle[,1]!=0)

###################################################################

nons=sum(bhatpnpmle[,1]!=0)

if (nons!=0) {

ws=as.matrix(w[,bhatpnpmle[,1]!=0])
zfs=as.matrix(zf[,bhatpnpmle[,1]!=0])
rebhatpnpmle=matrix(c(bhatpnpmle[bhatpnpmle[,1]!=0]),nons,1)

likelihood1=matrix(0,nons,nons)
ga=1/sqrt(cohort)

for (row in 1:nons){
for (col in 1:nons){
if (row>col) next

########################################################################likelihood1

inip=phatpnpmle
inich=chhatpnpmle

erow=c(rep(0,nons))
ecol=c(rep(0,nons))
erow[row]=1
ecol[col]=1

rebhatpnpmle1=rebhatpnpmle+ga*(erow+ecol)

for (kkk in 1:iteration2) {

wb=matrix(exp(exp(ws%*%rebhatpnpmle1)%o%(-inich[d:cohort])),numw,nonobsnum)
alpha[d:cohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

sigma21=matrix(0,cohort,1)
sigma21[1:obsnum,1]=(exp(zfs%*%rebhatpnpmle1))[1:obsnum]
sigma21[d:cohort,]=alpha[d:cohort,]%*%exp(ws%*%rebhatpnpmle1)
f=c(1/(tmC%*%sigma21))
chhat=c(rep(0,casenum))
chhat[1:casenum]=cumsum(f)
chhat[(casenum+1):cohort]=c(0,chhat[1:casenum])[hh]

inip=phat
inich=chhat
}

##############################################likelihood1

baseline=diff(c(0,chhat[1:casenum]))
la=sum((-chhat*exp(zfs%*%rebhatpnpmle1))[1:obsnum,])
lb=sum(log(baseline*exp(zfs[ncct[,2]==1,]%*%rebhatpnpmle1)))
lc=sum(log(alpha[1:obsnum,]%*%(phat)))
ld=sum(log(t(matrix(exp(exp(ws%*%rebhatpnpmle1)%o%(-chhat[d:cohort])),numw,nonobsnum))%*%phat))
likelihood1[row,col]=la+lb+lc+ld
}
}

for (row in 1:nons){
for (col in 1:nons){
if (row<=col) next
likelihood1[row,col]=likelihood1[col,row]
}
}

######################################################likelihood2
likelihood2=c(rep(0,nons))
for (row in 1:nons){

inip=phatpnpmle
inich=chhatpnpmle
erow=c(rep(0,nons))
erow[row]=1
rebhatpnpmle2=rebhatpnpmle+ga*erow

for (kkk in 1:iteration2){
wb=matrix(exp(exp(ws%*%rebhatpnpmle2)%o%(-inich[d:cohort])),numw,nonobsnum)
alpha[d:cohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

sigma21=matrix(0,cohort,1)
sigma21[1:obsnum,1]=(exp(zfs%*%rebhatpnpmle2))[1:obsnum]
sigma21[d:cohort,]=alpha[d:cohort,]%*%exp(ws%*%rebhatpnpmle2)
f=c(1/(tmC%*%sigma21))
chhat=c(rep(0,casenum))
chhat[1:casenum]=cumsum(f)
chhat[(casenum+1):cohort]=c(0,chhat[1:casenum])[hh]
inip=phat
inich=chhat
}

##############################################likelihood2

baseline=diff(c(0,chhat[1:casenum]))
la=sum((-chhat*exp(zfs%*%rebhatpnpmle2))[1:obsnum,])
lb=sum(log(baseline*exp(zfs[ncct[,2]==1,]%*%rebhatpnpmle2)))
lc=sum(log(alpha[1:obsnum,]%*%(phat)))
ld=sum(log(t(matrix(exp(exp(ws%*%rebhatpnpmle2)%o%(-chhat[d:cohort])),numw,nonobsnum))%*%phat))
likelihood2[row]=la+lb+lc+ld
}

sigma=matrix(0,nons,nons)
sigma=(-likelihood1+matrix(likelihood2,nons,nons)+t(matrix(likelihood2,nons,nons))-likelihood4)/(cohort*ga^2)

asysepnpmle=c(rep(0,numbeta))
asysepnpmle[bhatpnpmle[,1]!=0]=sqrt(diag(solve(sigma)/cohort))
wald=(bhatpnpmle[bhatpnpmle[,1]!=0]/asysepnpmle[bhatpnpmle[,1]!=0])^2
pvalue=c(rep(0,numbeta))
pvalue[bhatpnpmle[,1]!=0]=1-pchisq(wald,1)

pnpmletable=cbind(bhatpnpmle,asysepnpmle,pvalue)
colnames(pnpmletable)=c("PNPMLE","ASYSE","PVALUE")
numtable=rbind(casenum,obsnum)
plotelements=cbind(ncct[1:casenum,1],chhatpnpmle[1:casenum],exp(-chhatpnpmle[1:casenum]))

} else {

asysepnpmle=c(rep(0,numbeta))
pvalue=c(rep(0,numbeta))
pnpmletable=cbind(bhatpnpmle,asysepnpmle,pvalue)
colnames(pnpmletable)=c("PNPMLE","ASYSE","PVALUE")
numtable=rbind(casenum,obsnum)
likelihood=cbind(likelihood4,aic,bic)
colnames(likelihood)=c("log.likelihood","AIC","BIC")
plotelements=cbind(ncct[1:casenum,1],chhatpnpmle[1:casenum],exp(-chhatpnpmle[1:casenum]))

}

TPNPMLElist=list(num=numtable,iloop=iloop,diff=diff,cvl=cv,tuning=penaltytuning[suit],likelihood=likelihood4,pnpmle=pnpmletable,Lambdapnpmle=chhatpnpmle,Ppnpmle=phatpnpmle,elements=plotelements,adata=ncct)
return(TPNPMLElist)
}



