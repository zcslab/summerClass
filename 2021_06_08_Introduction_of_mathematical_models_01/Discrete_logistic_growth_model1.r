Discrete-time Logistic Population Growth Simulations

#Set number of generations to plot
PlotGen <- 200

###Set parameters for simulation 1###

# set initial population size
N0 <- 10

#Set R
R <- 2
#Set K
K <- 1000


###Run simulation###

# initialize vector to hold results 
PopSize <- N0 

# create variable to hold the current population size
PopNow <- N0 

# calculate population sizes and append to popsize


PopSize <- N0 
Pop_T <- N0 
print(PopSize)
for(i in 1:10) { 
        Pop_T<-PopSize[i]
	  Pop_Tplus1<- Pop_T+ Pop_T*R*(1-Pop_T/K)    #discrete logistic
        if (Pop_Tplus1< 0) {Pop_Tplus1<- 0}
        PopSize <- c(PopSize,Pop_Tplus1)                #add result to vector	
	  print(PopSize)
}


###Plot results###

tvals <- 1:length(PopSize) #create a vector of time values equal to length of results vector

plot(tvals[1:PlotGen], PopSize[1:PlotGen],type="o",col="red", 
xlab="Generation",ylab="Population size",pch=16,cex=.75)
