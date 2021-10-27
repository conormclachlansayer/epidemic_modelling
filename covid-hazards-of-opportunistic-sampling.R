#### README
#### Section 1 contains a covid simulation definition, a util for running the
####           simulation k times (whilst also computing some interesting
####           statistics), and also some plotting utils
#### Section 2 performs the simulations (also where simulation parameters are defined), produces visualizations,
####           and comments on how we can compare our results to inferences made from the ZOE dataset

##### Section 1

covid_model <- function(params) {
  ### FUNCTION: Simulates an epidemic for num_days for a population of size n
  ### Input: params list containing n: population, lambda: overall viral infectivity parameter,
  ### num_days: number of days to run the model, ni: initially infected
  ### Output: No of daily new infections for pop n, for 0.1% random sample 
  ### and 'cautious' lowest 10% of betas (as number and proportion of overall population)
  
  n <- params$n ; lambda <- params$lambda ; num_days <- params$num_days ; ni <- params$ni
  
  beta <- rlnorm(n, 0, 0.5) # Initialise individual betas
  beta <- beta/mean(beta) # standardising betas to have mean 1 over population
  
  # Notation: (S)usceptible = 0, (E)xposed = 1, (I)nfected = 2, and (R)ecovered/transition to serious disease = 3
  pop <- rep(0,n) # Initialise susceptible state for entire population
  
  pop[sample(1:n, ni)] <- 1 # randomly expose ni people in population
  
  lowest_beta <- beta <= quantile(beta, probs = 0.1) # find cautious 10%
  
  rand_0.1p <- sample(n, n*0.001) # picking random 0.1% sample of population
  
  e_to_i <- 1/3 # probability person moves from state exposed to infected
  i_to_r <- 1/5 # probability person moves from state infected to recovered
  
  n_I <- n_beta <- n_rand <- rep(0,num_days)
  
  for (i in 2:num_days){ # simulation of model over days 2 to t
    
    #### Start of day - initialize probabilities
    
    u <- runif(n) # generate chance of moving to next state in day i
    s_to_e <- lambda * beta * sum(beta[pop==2]) # generate individual prob of moving from S to E (from day i to j)
    
    #### Moving between states during day
    
    # newly infecteds are just the people who will move from state E to I
    
    newly_infected <- pop==1&u<e_to_i
    pop[pop==2&u<i_to_r] <- 3 # move people from I to R
    pop[newly_infected] <- 2  # move people from E to I
    pop[pop==0&u<s_to_e] <- 1 # move people from S to E
    
    #### End of day counts
    
    n_I[i] <- sum(newly_infected)
    n_beta[i] <- sum(newly_infected[lowest_beta])
    n_rand[i] <- sum(newly_infected[rand_0.1p])
    
  }
  
  list(new_inf = n_I,
       new_inf_betas = n_beta,
       new_inf_sample = n_rand,
       new_inf_std = n_I/n, # new infection rate of entire population
       new_inf_betas_std = n_beta/sum(beta <= quantile(beta, probs = 0.1)), # new infection rate of bottom 10% beta
       new_inf_sample_std = n_rand/(0.001*n)  # new infection rate of 0.1% random sample
  )
}

simulate_epidemic_k_times <- function(k, params){
  ### FUNCTION Simulates epidemic model k times
  ### Input: k; number of times to simulate model, params; same as covid_model params
  ### Output: No of daily new infections for pop n, for 0.1% random sample 
  ### and 'cautious' lowest 10% of betas (as number and proportion of overall population)
  ### A Vector of peak days across population, cautious 10% and random sample
  ### A Vector of number of days epidemic affected >1% of population, cautious 10% and random sample
  
  epi_list<- list() # to record k simulations
  peak_all <- peak_beta <- peak_rand <- rep(0,k) # to store times of epidemic peak over the k simulations
  num_bad_days_all <- num_bad_days_beta <- num_bad_days_rand <- rep(0, k) # number of days where the epidemic affected 1% of sub population
  
  for (i in 1:k){
    model <- covid_model(params)
    epi_list[[i]] <- model
    
    peak_all[i] <- which(model$new_inf == max(model$new_inf))[1] # finding peak day for all pop, beta pop and random pop
    peak_beta[i] <- which(model$new_inf_betas == max(model$new_inf_betas))[1]
    peak_rand[i] <- which(model$new_inf_sample == max(model$new_inf_sample))[1]
    
    num_bad_days_all[i] <- sum(model$new_inf_std > 0.01) # calculating length of epidemic (as defined above)
    num_bad_days_beta[i] <- sum(model$new_inf_betas_std > 0.01)
    num_bad_days_rand[i] <- sum(model$new_inf_sample_std > 0.01)
  }
  epi_list <- c(epi_list, list(peak_all, peak_beta, peak_rand, num_bad_days_all, num_bad_days_beta, num_bad_days_rand))
  names(epi_list) <- c(1:k, "peak_all", "peak_beta", "peak_rand", "num_bad_days_all", "num_bad_days_beta", "num_bad_days_rand")
  
  epi_list
}

plot_covid_model_graphs <- function(epi, labels = TRUE) {
  ### FUNCTION plot_covid_model_graphs creates plots showing daily infections over t days across the whole population, the cautious 10% and a 0.1% random sample
  ### Input: epi: covid_model output
  ### Output: A standardized plot comparing daily infections across the whole population, the cautious 10% and the 0.1% random sample
  
  # converting to %
  pop <- epi$new_inf_std*100
  beta <- epi$new_inf_betas_std*100
  rand <- epi$new_inf_sample_std*100
  
  # get peak day of epidemic amongst sub populations
  # if there are multiple peaks, just take the first
  peak_pop <- which(pop == max(pop))[1]
  peak_beta <- which(beta == max(beta))[1]
  peak_rand <- which(rand == max(rand))[1]
  
  plot(pop, # plotting points for daily population infections (in %)
       ylim=c(0,max(pop, beta, rand)*1.1), # limit of y axis: largest of max infection plus room for visual purposes
       xlab=if(labels){"Day of Epidemic"}else{""},ylab=if(labels){"Daily Incidence (% of sub population)"}else{""}) # setting axes title's
  
  points(beta, col="chocolate3"); points(rand, col= "blue") # plotting daily infections for cautious 10%
  lines(pop); lines(beta, col="chocolate3"); lines(rand, col= "blue") # and random 0.1% sample (in %)
  
  if(labels){ # option to turn off title and text on plots, default is on
    title(main = "Daily new infections among whole Population, a 0.1% random sample, and the cautious 10%")
    
    text(y = max(pop)*1.03, x = peak_pop, paste("Day ", peak_pop), cex = 1.2) # highlight day of peak infections
    text(y = max(beta)*1.03, x = peak_beta, paste("Day ", peak_beta), cex = 1.2, col="chocolate3")
    text(y = max(rand)*1.05, x = peak_rand, paste("Day ", peak_rand), cex = 1.2, col= "blue")
    
    text(y = max(pop)*1.03, x = peak_pop*.5, paste("Population incidence peak: ", round(max(pop),2),"%"), cex = 1.2) # displaying peak daily infections
    text(y = max(beta)*1.03, x = peak_beta*.5, paste("Cautious 10% incidence peak: ", round(max(beta),2),"%"), cex = 1.2, col="chocolate3")
    text(y = max(rand)*1.05, x = peak_rand*.5, paste("Random 0.1% incidence peak: ", round(max(rand),2),"%"), cex = 1.2, col= "blue")
  }
}

plot_boxplots <- function(entire, beta, rand, title) {
  ### FUNCTION create box plots on a single page for our 3 sub populations of interest
  ### Input: entire: vector of data for entire population, beta:   vector of data for cautious 10%
  # rand:   vector of data for random 0.1% sample
  
  par(mfcol=c(1,3), oma=c(2, 0, 2, 0))
  ymin <- min(entire, beta, rand)
  ymax <- max(entire, beta, rand)
  ylim <- c(ymin, ymax)
  boxplot(entire, xlab="Entire Population",  ylim=ylim)
  boxplot(beta,   xlab="Cautious 10%",       ylim=ylim, col = "chocolate3")
  boxplot(rand,   xlab="0.1% Random Sample", ylim=ylim, col = "blue")
  title(main=title, outer=TRUE)
}

##### Section 2: Analysis of Epidemic model and Results
n <- 5.5e+6
params <- list(n = n, lambda = 0.4/n, num_days = 150, ni = 10)

# Run the model once, plotting graph which shows how daily infection 
# trajectories differ between our three sub populations of interest

epi <- covid_model(params)
par(mfcol=c(1,1),mar=c(4,4,1,1))
plot_covid_model_graphs(epi)

# Run the model another 10 times in order to visualise the variability between simulations.

epi_list <- simulate_epidemic_k_times(k = 10, params)

# here we produce 10 plots similar to the single run,
par(mfcol=c(5,2),mar=c(4,4,1,1))
for (i in 1:10){
  plot_covid_model_graphs(epi_list[[i]], labels = FALSE)
}
title(main="Epidemic over 10 simulations: Daily Incidence (% of sub population) vs Day of Epidemic", outer=TRUE)

plot_boxplots(epi_list[["peak_all"]], epi_list[["peak_beta"]], epi_list[["peak_rand"]], #show how the peak days were distributed,
              title="Variation of epidemic peak day among each sub population")
# boxplots show that epidemic peak occurred later for cautious 10% compared to rest (i.e. a lag in timing)

plot_boxplots(epi_list[["num_bad_days_all"]], epi_list[["num_bad_days_beta"]], #show how the length of the epidemics were distributed.
              epi_list[["num_bad_days_rand"]], # defined length of epidemic as number of days more than 1 percent of pop were affected
              title="Variation in the number of days where the epidemic affected more than 1% of each sub population")
# boxplots show much shorter epidemics for cautious 10% compared to rest

