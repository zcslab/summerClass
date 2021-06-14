#Multispecies Logistic Population Growth Simulations

#Set number of generations to plot
PlotGen <-4000

###Set parameters for Competitors 1 and 2 and Predator###

# set initial population sizes
CompA0 <- 10
CompB0 <- 10
Pred0 <- 2

#Set r
rA <- .01
rB <- .01

#Set K
KA <- 500
KB <- 500

#Set competition coefficients
alpha <- 0.9   #effect of B on A
beta <- 0.9    #effect of A on B

#Set predator-prey parameters
PredEff <- 0.001   #parameter "a" in the Lotka-Volterra predator-prey model
PredConv <- 0.5  #parameter "f" in the Lotka-Volterra predator-prey model
PredStarve <- 0.01  #parameter "q" in the Lotka-Volterra predator-prey model
PredPref <- 0.5   #predator's preference: from 0 (eats only A) to 0.5 (eats A, B equally) to 1 (eats only B)

###Run simulation###

# initialize vector to hold results 
CompAPopSize <- CompA0 
CompBPopSize <- CompB0
PredPopSize <- Pred0

# create variable to hold the current population size
CompAPopNow <- CompA0 
CompBPopNow <- CompB0
PredPopNow <- Pred0

# calculate a pair of predation preference adjustments based on PredPref
PredPrefA <- 1 - PredPref
PredPrefB <- PredPref

# calculate population sizes and append to popsize
for(i in 1:PlotGen) { 
	  # implement competition for Competitor A
        CompAPopNow <- CompAPopNow + CompAPopNow*rA*(1-CompAPopNow/KA - alpha*CompBPopNow/KA)    #discrete logistic
        if (CompAPopNow < 0) {CompAPopNow <- 0}

	  # implement predation losses for Competitor A
  	  PredEatsA <- PredPrefA*PredEff*CompAPopNow*PredPopNow
	  CompAPopNow <- CompAPopNow - PredEatsA
        if (CompAPopNow < 0) {CompAPopNow <- 0}

	  # implement competition for Competitor B
        CompBPopNow <- CompBPopNow + CompBPopNow*rB*(1-CompBPopNow/KB - beta*CompAPopNow/KB)    #discrete logistic
        if (CompBPopNow < 0) {CompBPopNow <- 0}

	  # implement predation losses for Competitor B
  	  PredEatsB <- PredPrefB*PredEff*CompBPopNow*PredPopNow
	  CompBPopNow <- CompBPopNow - PredEatsB
        if (CompBPopNow < 0) {CompBPopNow <- 0}

	  # implement predator population growth
	  PredPopNow <-PredPopNow + PredConv*(PredEatsA + PredEatsB) - PredStarve*PredPopNow
	  if (PredPopNow < 0) {PredPopNow <- 0}

        CompAPopSize <- c(CompAPopSize,CompAPopNow)                #add results to vectors
        CompBPopSize <- c(CompBPopSize,CompBPopNow)                
	  PredPopSize <- c(PredPopSize,PredPopNow)
}

###Plot results###
#find maximum value for y axis
maxes <- c(max(CompAPopSize), max(CompBPopSize), max(PredPopSize))
ymax <-max(maxes)
YAxisLength <- max(pretty(ymax*1.05))

tvals <- 1:PlotGen #create a vector of time values equal to length of results vector

#par(yaxp=c(0, YAxisLength, 5))
plot(tvals[1:PlotGen], CompAPopSize[1:PlotGen],type="o",col="red", 
xlab="Generation",ylab="Population size",pch=16,cex=.75, ylim=c(0, YAxisLength))
points(tvals[1:PlotGen], CompBPopSize[1:PlotGen], col="blue")
points(tvals[1:PlotGen], PredPopSize[1:PlotGen], col= "black")


