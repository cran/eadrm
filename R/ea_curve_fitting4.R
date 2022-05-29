#' Fits a dose-response curve using an evolutionary algorithm
#'
#' Uses an evolutionary algorithm to fit a set of four possible
#' dose-response models. The evolutionary parameters and stopping
#' rules can be customized by the user.
#'
#' @param obs A vector of response values (y-values).
#' @param xvals A vector of doses (x-values).
#' @param model Type of dose-response model to fit. Possible values
#' include "h3", "h4", and "h5" (corresponding to 3-parameter,
#' 4-parameter, and 5-parameter log-logistic models, respectively),
#' "e" (corresponding to an exponential model) and "all" (which allows
#' the procedure to evaluate all four types of models). Defaults to "h4".
#' @param pop.size The number of initial potential solutions.
#' Defaults to 1000.
#' @param stable.pop.size This quantity is divided by the number of
#' tournaments to calculate the number of children in each generation.
#' Defaults to 200.
#' @param num.tournaments The number of tournaments in each generation.
#' Defaults to 20.
#' @param tournament.size The number of players (i.e., models to
#' consider) in each tournament. Defaults to 10.
#' @param max.generations The maximum number of generations. If this
#' number is reached, the algorithm immediately terminates. Defaults to
#' 500.
#' @param stop.generations The algorithm will also terminate if there is
#' no improvement in fitness in stop.generation generations. Defaults to
#' 50.
#' @param epsilon If three successive new models produce an improvement
#' of less than epsilon in fitness, the procedure will terminate.
#' Defaults to 0.1.
#' @section Details:
#' The procedure will initially generate pop.size possible solutions.
#' The fitness for each solution will be calculated. In each generation,
#' a series of tournaments are performed. The children of the surviving
#' models from the previous generation are mutated. The model with the
#' best fitness "wins" each tournament and survives to the next
#' generation. The procedure continues until the maximum number of
#' generations is reached or the fitness fails to improve substantially
#' over a sufficient number of generations.
#' @return An object of class eadrm, which is a list containing the
#' following elements:
#' \describe{
#' \item{Model:}{Specifies the type of model (i.e., "h3", "h4", "h5", or
#' "e")}
#' \item{R2:}{Fitness for the final model}
#' \item{params:}{A vector of coefficients for the final model}
#' \item{xvals:}{The original x values (concentrations) for the model}
#' \item{yvals:}{The original y values (responses) for the model}
#' }
#' @references
#' Ma, J., Bair, E., Motsinger-Reif, A. "Nonlinear Dose-response Modeling
#' of High-Throughput Screening Data Using an Evolutionary Algorithm",
#' Dose Response 18(2):1559325820926734 (2020).
#' @importFrom stats runif
#' @export
#' @examples
#' ea.fit <- eadrm(CarboA$y, CarboA$x)
eadrm <- function(obs, xvals, model='h4', pop.size=1000,
                  stable.pop.size=200, num.tournaments=20,
                  tournament.size=10, max.generations=500,
                  stop.generations=50, epsilon=0.1)
{
  ## Initial Hillslope Estimates ##
  Top <- min(obs)
  Bottom <- max(obs)
  num.conc <- length(xvals)
  EC50 <- EC50.start(obs,xvals)
  W <- 2.5
  f<-0.5
  SST <- sum((obs - mean(obs))^2)
  num.children <- (stable.pop.size / num.tournaments) - 1     ### default 200/20 - 1 = 9
  hs.params <- c(Top,Bottom,EC50,W)  ## the third hs.params is EC50!!!


  Top2<-max(obs)
  Bottom2<-max(obs)
  hs3.params <- c(Top2,EC50,W)
  hs5.params <- c(Top2,Bottom2,EC50,W,f)
  #print('Starting values')
  #print('Top, Bottom, EC50, W')
  #print(hs.params)


  ## Initital Exponential Estimates ##
  B = 1
  k = 2.5
  exp.params <- c(B,k)

  ## Initial Population ##
  population <- initial.population(pop.size,obs,xvals,hs.params,hs3.params, hs5.params,exp.params,model)     ### default pop size = 1000

  ## Variable used for convergence measuring ##
  fitness.history <- NULL
    ###  Being Generational Simulation ###
    allTimeBestFitness <- c(-400,-300,-200,-100)           #####*******************#####
    allTimeBestFit <- NULL               #####*******************#####
    generation <- 0
    last.step <- 0

    while (max(allTimeBestFitness[2:4]-allTimeBestFitness[1:3])>epsilon &
    	   generation<max.generations & last.step<stop.generations)               ### iterate until no improvement or max.generations is reached
    {
      generation <- generation + 1
      last.step <- last.step + 1
      ###  Tournament Selection ###
      max.int  <- length(population)
      next.generation <- NULL
      tournament.winners <- list()
      numWinners <- 0
      for(tourn in 1:num.tournaments)                   ### default num.tournaments = 20
      {
        tournament <- list()
        max.fitness <- -Inf
        winner <- -Inf
        ### selects 'player' with highest r^2 from 10 random individuals from population
        for(player in 1:tournament.size)                ### default tournament.size = 10
        {
          player.index <- as.integer(runif(1)*max.int+1)
          current.player <- population[[player.index]]
          tournament[[player]] <- current.player
          current.fitness <- as.numeric(current.player[2])
          if(current.fitness > max.fitness)
          {
            max.fitness <- current.fitness
            winner <- player.index
          }
        }
        tournament.winners[[numWinners <- numWinners+1]] <- population[[winner]]  ### list of 20 tournament winners
      }

      next.gen <- list()
      gen.count <- 0
      c <- 0.2*((max.generations - generation)/max.generations)  ## Set mutation decary rate
      fitnesses <- NULL
      n <- length(tournament.winners)
      for(i in 1:n)
      {
      	fitnesses <- c(fitnesses, tournament.winners[[i]][2])
      }
      ranks <- n - rank(fitnesses) + 1
      for(i in 1:n)
      {
        current <- tournament.winners[[i]]
        ##Allow tournament winner to have children##
        child.count <- 1
        #while(child.count <= num.children)
        totalChildren <- 2*num.children*(ranks[i]/num.tournaments)
        while(child.count <= totalChildren)
        {
           current.params <- as.numeric(current[3:length(current)])
           child <- new.individual(current.params, c)
           pred.vals <- NULL
           if(current[1] == 'h4') pred.vals <- predict.Hillslope.values(current.params,xvals)
           if(current[1] == 'h3') pred.vals <- predict.Hillslope3.values(current.params,xvals)
           if(current[1] == 'h5') pred.vals <- predict.Hillslope5.values(current.params,xvals)
           if(current[1] == 'e') pred.vals <- predict.exponential.values(current.params,xvals)
           child.fitness <- fitness(current[1],obs,pred.vals)
           next.gen[[gen.count <- gen.count+1]] <- c(current[1],child.fitness,child)
           child.count <- child.count + 1
         }
         ### Add parent back in###
         next.gen[[gen.count <- gen.count+1]] <- current
      }

      a <- 0                                       ### 'a' is the sum of r^2 among tournament winners
      for(i in 1:length(tournament.winners))
      {
        a <- a + as.numeric(tournament.winners[[i]][2])
      }
      mean.fitness <-  a / length(tournament.winners)
      fitness.history <- rbind(fitness.history,mean.fitness)
      population <- next.gen



      winner.max <- as.numeric(population[[1]][2])
      winner.row <- 1
      for(i in 2:length(population)){
      	thisFitness <- as.numeric(population[[i]][2])
        if(thisFitness > winner.max){
        	winner.max <- thisFitness
        	winner.row <- i
        }
      }


      winner <- population[[winner.row]]
      if(winner.max > allTimeBestFitness[4]){
         allTimeBestFitness <- c(allTimeBestFitness[2:4], winner.max)
         allTimeBestFit <- winner
	 last.step <- 0
      }

    }  # end generation loop
    #print(population[[1]][2])
#    print(allTimeBestFitness)
#    print(generation)

    params <- NULL
    if(winner[1] == "h4")
    {
      params <- rbind(c(winner[3:6]))
      colnames(params) <- c("EMAX","EMIN","EC50","W")
    } else if(winner[1] == "h3")
      {
        params <- rbind(c(winner[3:5]))
        colnames(params) <- c("EMAX","EC50","W")
    }else if(winner[1] == "h5")
    {
      params <- rbind(c(winner[3:7]))
      colnames(params) <- c("EMAX","EMIN","EC50","W","f")
    }else if(winner[1] == "e")
    {
      params <- rbind(c(winner[3:4]))
      colnames(params) <- c("B","K")
    }
  cur.out <- list(Model=winner[1],R2=winner[2],params=params,xvals=xvals,yvals=as.vector(obs))
  class(cur.out) <- "eadrm"
	return(cur.out)

     ##### Down to here #########

}

## Creates an initial population of specified size seeded by initial paramter estimates ##
initial.population <- function(size,obs,xvals,hs.params,hs3.params, hs5.params,exp.params,model)
{
  population <- list()
  SST <- sum((obs - mean(obs))^2)
  count <- 0
  if(model == 'h4')
  {
    ##Generate initial population##
    temp.pop <- NULL
    for(i in 1:size)
    {
      individual <- new.initial.individual(hs.params)               ### an 'individual' is just a set of parameters
      predicted.vals <- predict.Hillslope.values(individual,xvals)
      R.squared <- fitness('h4',obs,predicted.vals)
      new.member <- c('h4',R.squared,individual)                     ### each member of population is a model, an r^2, and an 'individual'
      population[[count <- count+1]] <- new.member
    }
  } else if(model == 'h3')
    {
      ##Generate initial population##
      for(i in 1:size)
      {
        individual <- new.initial.individual(hs3.params)
        predicted.vals <- predict.Hillslope3.values(individual,xvals)
        R.squared <- fitness('h3',obs,predicted.vals)
        new.member <- c('h3',R.squared,individual)
        population[[count <- count+1]] <- new.member
      }
    } else if(model == 'h5')
    {
      ##Generate initial population##
      for(i in 1:size)
      {
        individual <- new.initial.individual(hs5.params)
        predicted.vals <- predict.Hillslope5.values(individual,xvals)
        R.squared <- fitness('h5',obs,predicted.vals)
        new.member <- c('h5',R.squared,individual)
        population[[count <- count+1]] <- new.member
      }
    } else if(model == 'e')
    {
      ##Generate initial population##
      for(i in 1:size)
      {
        individual <- new.initial.individual(exp.params)
        predicted.vals <- predict.exponential.values(individual,xvals)
        R.squared <- fitness('e',obs,predicted.vals)
        new.member <- c('e',R.squared,individual)
        population[[count <- count+1]] <- new.member
      }
    } else if(model == 'all')
      {
        ##Generate initial population##
        for(i in 1:(size/4))
        {
          individual <- new.initial.individual(hs.params)
          predicted.vals <- predict.Hillslope.values(individual,xvals)
          R.squared <- fitness('h4',obs,predicted.vals)
          new.member <- c('h4',R.squared,individual)
          population[[count <- count+1]] <- new.member
        }
      for(i in 1:(size/4))
      {
        individual <- new.initial.individual(hs3.params)
        predicted.vals <- predict.Hillslope3.values(individual,xvals)
        R.squared <- fitness('h3',obs,predicted.vals)
        new.member <- c('h3',R.squared,individual)
        population[[count <- count+1]] <- new.member
      }
      for(i in 1:(size/4))
      {
        individual <- new.initial.individual(hs5.params)
        predicted.vals <- predict.Hillslope5.values(individual,xvals)
        R.squared <- fitness('h5',obs,predicted.vals)
        new.member <- c('h5',R.squared,individual)
        population[[count <- count+1]] <- new.member
      }
        for(i in 1:(size/4))
        {
          individual <- new.initial.individual(exp.params)
          predicted.vals <- predict.exponential.values(individual,xvals)
          R.squared <- fitness('e',obs,predicted.vals)
          new.member <- c('e',R.squared,individual)
          population[[count <- count+1]] <- new.member
        }
      }
  return(population)
}

## Individual in initial population. Allowed to Mutate up to 100% from parental seed ##
new.initial.individual <- function(params)
{
    newParams <- NULL
    for(i in 1:length(params))
    {
      ran <- as.integer(runif(1)*100+1)
      if(ran%%2 == 1) ran <- -ran
      mutated <- params[i] + (params[i] * (ran/100))
      newParams <- c(newParams, mutated)
    }

    ## make sign of exponential parameter negative to fit inhibitory responses ##
    ran <- as.integer(runif(1)*100+1)
    if(ran%%2 == 1) newParams[length(newParams)] <- newParams[length(newParams)] * -1
    return(newParams)
}

#### Begin utility funtions ####

## Predict values based on HillSlope model and xvals ##
predict.Hillslope.values <- function(params,xvals)
{
    y.vals <- params[1] -
        (params[1]-params[2])/(1+(xvals/params[3])^(params[4]))
    return(y.vals)
}

##Predict values based on 3-p logistic model
predict.Hillslope3.values <- function(params,xvals)
{

    y.vals <- (params[1])/(1+exp(params[3]*(log(xvals)-log(params[2]))))
    return(y.vals)
}

##Predict values based on 5-p logistic model

predict.Hillslope5.values <- function(params,xvals)
{
    y.vals <- params[2]+(params[1]-params[2])/
        ((1+exp(params[4]*(log(xvals)-log(params[3]))))^params[5])
    return(y.vals)
}


## Predict values based on exponential model and xvals ##
predict.exponential.values <- function(params,xvals)
{
    y.vals <- params[1]*exp(params[2]*xvals)
    return(y.vals)
}


new.individual <- function(params, c)
{
	newParams <- params*(1 + runif(length(params), -c, c))
    return(newParams)
}


## Calculates fitness of indivdual ##
fitness <- function(model,obs,predicted.vals)
{
  n1<-length(obs)
  if(model == 'h4') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+4*log(n1))
  if(model == 'h3') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+3*log(n1))
  if(model == 'h5') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+5*log(n1))
  if(model == 'e') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+2*log(n1))
  return(aiccminus)
}


## Get start value for EC50, which is the xval in the "middle" ##
EC50.start <- function(obs,xvals)
{
  obs <- obs[order(xvals)]
  xvals <- sort(xvals)
  if(length(xvals) %% 2 == 0)
  {
      first <- length(xvals) / 2
      last  <- first + 1
      EC50  <- 10^( (log10(xvals[first]) + log10(xvals[last]))/2 )
  } else {
      EC50 <- xvals[round(length(xvals) / 2)]
  }

  return(EC50)
}
