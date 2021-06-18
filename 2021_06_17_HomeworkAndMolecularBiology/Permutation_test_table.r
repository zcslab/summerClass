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