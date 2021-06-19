1:50000
1:100

proc.time()
stat_all<-rep(0,10000)
for(i in 1:10000)
{
	stat_all[i]<-sum(sample(1:50000,2000)<=100)
}
proc.time()

hist(stat_all,xlim=c(0,42))
abline(v=40)
abline(v=6)

p_estimated<-sum(stat_all>=8)/length(stat_all)
p_estimated
p_estimated<0.05



p value 

The probability that the null is corrected (conditional to the observed data)
Monte Carlo Approach 




#####
p<-0.1
aa<-c()
for(i in 1:1000)
{
	aa<-c(aa,sum(runif(1000,0,1)<p))
}

bb<-rbinom(1000,1000,p)

par(mfcol=c(1,2))
hist(aa)
hist(bb)

stat_all<-rep(0,50000)
for(i in 1:50000)
{
	stat_all[i]<-sum(sample(1:50000,2000)<=100)
}


par(mfcol=c(1,2))
hist(stat_all)
hist(rhyper(50000, m=2000, n=50000, k=100))


plot(dnorm((-40:40)/10,0,01))
plot(pnorm((-40:40)/10,0,01)) p(X<input|X~the distribution) 
pnorm(1,0,1)

pnorm
qnorm
rnorm


unif

1-phyper(10, m=2000, n=50000, k=100)


