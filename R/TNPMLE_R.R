#' @title Non-Parametric Maximum-Likelihood Estimation for Cohort Samplings with Time Matching under Cox's Regression Model
#' @description The function utilizes a self-consistency iterative algorithm to calculate NPMLEs for cohort samplings with time matching under Cox's regression model. In addition to compute NPMLEs, it can also estimate asymptotic varance, as described in Wang et al. (2019). The Cox's regression model is \deqn{\lambda(t|z)=\lambda_{0}(t)\exp(z^T\beta).}
#' @param data The \eqn{N \times P} matrix of data. There are \eqn{N} individuals in matrix, with one individual in each row. The \eqn{P} columns orderly included the observable times which are time-to-event or censoring times and without ties at the event times, the status is a binary variable with 1 indicating the event has occured and 0 indicating (right) censoring, and the (\eqn{P-2}) covariates which only observed for some individuals. Note that the covariate of those unobserved individuals are denoted by \eqn{-9}, not missing value (NA) and the observed covariates values are not the same as \eqn{-9}.
#' @param iteration1 The number of iteration for computing NPMLEs.
#' @param iteration2 The number of iteration for computing profile likelihoods which are used to estimate asymptotic variance.
#' @param converge The parameter influence the convergence of the algorithm, if the sup-norm of \eqn{\left(\hat{\beta}_{(k)}-\hat{\beta}_{(k-1)}\right)} is smaller than the thresholding value, we then declare the estimates are converge, stop computing estimates, otherwise the number of iteration for computing estimates is the \code{iteration1} term.
#' @return Returns a list with components
#' @return \item{num}{The numbers of case and observed subjects.}
#' @return \item{iloop}{The final number of iteration for computing NPMLEs.}
#' @return \item{diff}{The sup-norm distance between the last two iterations of the estimates of the relative risk coefficients.}
#' @return \item{likelihood}{The log likelihood value of NPMLEs.}
#' @return \item{npmle}{The estimated regression coefficients with their corresponding estimated standard errors and p-values.}
#' @return \item{Lnpmle}{The estimated cumulative baseline hazards function.}
#' @return \item{Pnpmle}{The empirical distribution of covariates which are missing for unobserved subjects.}
#' @return \item{elements}{A list which is used to plot cumulative baseline hazards function and baseline survival function. The \eqn{n \times 3} matrix of data, \eqn{n} is the total number of case and the \eqn{3} columns orderly included the order observed time of case, the estimated cumulative baseline hazards function and estimated baseline survival function.}
#' @return \item{Adata}{Arranging original data to let our analysis performed conveniently. There are three steps for this arrangement, the 1st step divides original data into observed and unobserved groups, then put them on top and bottom, respectively; the 2nd step divides the observed data of 1st step into case and control groups; the final step order the case data of 2nd step by observed time from low to high.}
#' @references Wang JH, Pan CH, Chang IS*, and Hsiung CA (2019) Penalized full likelihood approach to variable selection for Cox's regression model under nested case-control sampling. published in Lifetime Data Analysis <doi:10.1007/s10985-019-09475-z>.
#' @note The missing value (NA) in the DATA is not allowed in this version.
#' @seealso See \code{\link{TPNPMLE}}.
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
#' TNPMLEtest=TNPMLE(data=ncct,iteration1=100,iteration2=30,converge=0)

                ##################################################################
                #              CREATE NCCNPMLE FUNCTION for R PACKAGE            #
                ##################################################################


TNPMLE=function(data,iteration1,iteration2,converge){

#####################################################################
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

ncct=data.frame(ARRANGEDATA(data)$arrangedata)
casenum=sum(ncct[,2])
cohort=dim(ncct)[1]
covariate=dim(ncct)[2]
obsnum=sum(ncct[,covariate]!=-9)
nonobsnum=cohort-obsnum
d=obsnum+1
numbeta=covariate-2

############################################ 
#produce W1,...,WJ                         #  
############################################

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

                                ################################################################ 
                                #                            NPMLE                             #
                                ################################################################ 
################################
#initial value                 #
################################

ncct=as.matrix(ncct)
covariatedata=matrix(0,obsnum,numbeta)
covariatedata=ncct[1:obsnum,3:covariate]
ncct=data.frame(ncct)
inibeta=matrix(0,numbeta,1)
inip=c(rep(1/numw,numw)) 
inich=c(rep(0,cohort))

hh=apply(t(matrix(rep(ncct[1:casenum,1],(cohort-casenum)),casenum,(cohort-casenum)))<=ncct[(casenum+1):cohort,1],1,sum)+1
r3=matrix(1,1,obsnum)
r4=matrix(1,1,nonobsnum)
expand=function(yy) c(yy%*%t(yy))

alpha=matrix(0,cohort,numw)
if (numw==obsnum) {
 diag(alpha[1:numw,1:numw])=1
} else {
for (i in 1:obsnum){
for (j in 1:numw){
if (sum(abs(zf[i,]-w[j,]))<1e-10) alpha[i,j]=1 else alpha[i,j]=0
}
}
}

tmC=(t(matrix(rep(ncct[,1],casenum),cohort,casenum))>=casetime)+0

for (iloop in 1:iteration1){
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
bhat=inibeta-ginv(order2)%*%order1

#################################################################
diff=max(abs(bhat-inibeta))
if (diff<converge) break
#################################################################

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

bhatnpmle=bhat
chhatnpmle=chhat
phatnpmle=phat

########################################################
#compute asymptotic variance;likelihood4               #
########################################################

baseline=diff(c(0,chhatnpmle[1:casenum]))
la=sum((-chhatnpmle*exp(zf%*%bhatnpmle))[1:obsnum,])
lb=sum(log(baseline*exp(zf[ncct[,2]==1,]%*%bhatnpmle)))
lc=sum(log(alpha[1:obsnum,]%*%(phatnpmle)))
ld=sum(log(t(matrix(exp(exp(w%*%bhatnpmle)%o%(-chhatnpmle[d:cohort])),numw,nonobsnum))%*%phatnpmle))
likelihood4=la+lb+lc+ld
#aic=-2*likelihood4+2*numbeta
#bic=-2*likelihood4+log(cohort)*numbeta

###################################################################

ga=1/sqrt(cohort)
likelihood1=matrix(0,numbeta,numbeta)
for (row in 1:numbeta){
for (col in 1:numbeta){
if (row>col) next

####################################################################likelihood1
inip=phatnpmle
inich=chhatnpmle

erow=c(rep(0,numbeta))
ecol=c(rep(0,numbeta))
erow[row]=1
ecol[col]=1

bhatnpmle1=bhatnpmle+ga*(erow+ecol)

for (kkk in 1:iteration2) {
wb=matrix(exp(exp(w%*%bhatnpmle1)%o%(-inich[d:cohort])),numw,nonobsnum)
alpha[d:cohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

sigma21=matrix(0,cohort,1)
sigma21[1:obsnum,1]=(exp(zf%*%bhatnpmle1))[1:obsnum]
sigma21[d:cohort,]=alpha[d:cohort,]%*%exp(w%*%bhatnpmle1)
f=c(1/(tmC%*%sigma21))
chhat=c(rep(0,casenum))
chhat[1:casenum]=cumsum(f)
chhat[(casenum+1):cohort]=c(0,chhat[1:casenum])[hh]

inip=phat
inich=chhat
}

##############################################likelihood1

baseline=diff(c(0,chhat[1:casenum]))
la=sum((-chhat*exp(zf%*%bhatnpmle1))[1:obsnum,])
lb=sum(log(baseline*exp(zf[ncct[,2]==1,]%*%bhatnpmle1)))
lc=sum(log(alpha[1:obsnum,]%*%phat))
ld=sum(log(t(matrix(exp(exp(w%*%bhatnpmle1)%o%(-chhat[d:cohort])),numw,nonobsnum))%*%phat))
likelihood1[row,col]=la+lb+lc+ld
}
}

for (row in 1:numbeta){
for (col in 1:numbeta){
if (row<=col) next
likelihood1[row,col]=likelihood1[col,row]
}
}

######################################################likelihood2

likelihood2=c(rep(0,numbeta))
for (row in 1:numbeta){

inip=phatnpmle
inich=chhatnpmle 
erow=c(rep(0,numbeta))
erow[row]=1
bhatnpmle2=bhatnpmle+ga*erow

for (kkk in 1:iteration2){
wb=matrix(exp(exp(w%*%bhatnpmle2)%o%(-inich[d:cohort])),numw,nonobsnum)
alpha[d:cohort,]=t(wb*inip)*c(1/(t(wb)%*%inip))
phat=apply(alpha,2,mean)

sigma21=matrix(0,cohort,1)
sigma21[1:obsnum,1]=(exp(zf%*%bhatnpmle2))[1:obsnum]
sigma21[d:cohort,]=alpha[d:cohort,]%*%exp(w%*%bhatnpmle2)
f=c(1/(tmC%*%sigma21))
chhat=c(rep(0,casenum))
chhat[1:casenum]=cumsum(f)
chhat[(casenum+1):cohort]=c(0,chhat[1:casenum])[hh]
inip=phat
inich=chhat
}

##############################################likelihood2

baseline=diff(c(0,chhat[1:casenum]))
la=sum((-chhat*exp(zf%*%bhatnpmle2))[1:obsnum,])
lb=sum(log(baseline*exp(zf[ncct[,2]==1,]%*%bhatnpmle2)))
lc=sum(log(alpha[1:obsnum,]%*%(phat)))
ld=sum(log(t(matrix(exp(exp(w%*%bhatnpmle2)%o%(-chhat[d:cohort])),numw,nonobsnum))%*%phat))
likelihood2[row]=la+lb+lc+ld
}

sigma=matrix(0,numbeta,numbeta)
sigma=(-likelihood1+matrix(likelihood2,numbeta,numbeta)+t(matrix(likelihood2,numbeta,numbeta))-likelihood4)/(cohort*ga^2)

asyse=sqrt(diag(solve(sigma)/cohort))
wald=(bhatnpmle/asyse)^2
pvalue=1-pchisq(wald,1)
npmletable=cbind(bhatnpmle,asyse,pvalue)
colnames(npmletable)=c("NPMLE","ASYSE","PVALUE")
numtable=rbind(casenum,obsnum)
plotelements=cbind(ncct[1:casenum,1],chhatnpmle[1:casenum],exp(-chhatnpmle[1:casenum]))
TNPMLElist=list(num=numtable,iloop=iloop,diff=diff,likelihood=likelihood4,npmle=npmletable,Lambdanpmle=chhatnpmle,Pnpmle=phatnpmle,elements=plotelements,adata=ncct)
return(TNPMLElist)
}

