### Modelling pSIA + oSIA strategy ###
rm(list=ls())	
# Load required libraries
require(SimInf)
library(tidyverse)
library(dplyr)
library(data.table)

setwd("/Users/idcvmauz/Dropbox/PolioPreventingOutbreaks/16_EconRisk/polio_siminf/ES")

# Define parameters
beta <- 0.375
gamma <- 0.125
r <- beta/gamma
caseinfection <- 1/200
I0 <- 0
reps <- 10 # number of stochastic simulations
pop2 <- 500000
pop1 <- 1000000-pop2
ri_coverage <- 0.5
prop_Rv <- ri_coverage # RI coverage
prop_S0 <- (1-prop_Rv)-0.003 #subtract just so S1, S2 and S3 compartments aren't empty at start of model
prop_S1 <- 0.001 # can't have these compartments empty, so assign a v nominal number of individuals
prop_S2 <- 0.001
prop_S3 <- 0.001
sia_coverage <- 0.5
ipv_coverage <- ri_coverage
ve <- 0.5 # vaccine efficacy for bOPV
doses <- 3
seroconversion <- (1 -(1-ve)^doses)
birth_rate1 <- 5*10^-4 * pop1 
death_rate1 <- birth_rate1
birth_rate2 <- 5*10^-4 * pop2 
death_rate2 <- birth_rate2
births_S <- round(birth_rate1 * (1-ri_coverage), digits=0)
births_S3 <- round(birth_rate1 * ri_coverage, digits=0)
births2_S <- round(birth_rate2 * (1-ri_coverage), digits=0)
births2_S3 <- round(birth_rate2 * ri_coverage, digits=0)

# Generate times for random importations of infection
set.seed(123) # Set seed for reproducibility
generate_two_numbers <- function(start, end) {
  sample(start:end, 2)
}

random_importations <- c(
  generate_two_numbers(2, 366),
  generate_two_numbers(367, 731),
  generate_two_numbers(732, 1096),
  generate_two_numbers(1097, 1461),
  generate_two_numbers(1462, 1825)
)

# Timeliness assumptions
current_time <- 0
time_horizon <- 5 * 365 # 5 years in days
interval <- 35

# Mobility between populations A and B
mobility1 <- 0.0001 # from A to B
mobility2 <- 0.001 # from B to A

# Define surveillance sensitivities in both nodes
afp_sens1 <- 0.05
es_sens1 <- 0
afp_sens2 <- afp_sens1
es_sens2 <- afp_sens1

# Define the 'select' matrix
Et <- matrix(c(1, rep(0,17), # births R0
               rep(0,8),1,rep(0,9), # births S3
               rep(1,11),0,1,1,0,0,1,0, # deaths
               1, rep(0,17), # dynamic importation
               1,0,1,0,1,0,1,0,1,rep(0, 9), # vaccination
               0,1,0,1,0,1,0,1,rep(0,10), # NO seroconversion
               0,1,0,rep(0,15), # seroconversionS1
               0,0,0,1,0,rep(0,13), # seroconversionS2
               0,0,0,0,0,1,0,rep(0,11), # seroconversionS3
               1, rep(0,17), # infection SS
               rep(0,8),1,0,0,0,0,0,0,0,0,0,# seroconversionS3r 
               rep(0,7),1,rep(0,10), # seroconversionS4
               rep(0,9),1,rep(0,8), # seroMAX
               rep(0,8),1,rep(0,9)), # ipv
             nrow = 18, ncol = 14, 
             dimnames = list(c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v","S3r", "S4", "Rv", 
                               "Ii", "I", "C", "ES", "Icum", "Rn", "Ccum"), 
                             c("birthsS0", "birthsS3", "deaths", "importations", "vaccination",
                               "no_seroconversion", "seroconversionS1", "seroconversionS2", "seroconversionS3", "infection",
                               "seroconversionS3r", "seroconversionS4", "seroMAX", "ipv")))

# Define the 'shift' matrix
Nt <- matrix(c(1,0,1,0,1,0,1,0,1,rep(0, 9), # vaccination
               0,1,0,1,0,1,0,2,rep(0,10), # NO seroconversion
               0,9,0,7,0,5,0,2,0,1,rep(0,8), # seroconversion
               12, rep(0,17), # importation
               rep(0,8),2,rep(0,9), # ipv vaccination
               rep(0,18)), # migration between 
             nrow = 18, ncol = 6,
             dimnames = list(c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v","S3r","S4", "Rv", 
                               "Ii", "I", "C", "ES", "Icum", "Rn", "Ccum"), 
                             c("vaccination", "no_seroconversion", "seroconversion", "importation", "ipv", "migration")))


create.StocSIRobr0 <- function(I0, beta, gamma, caseinfection, f_time, initial_state, events) {
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", 
                    "Rv", "Ii", "I", "C", "ES", "Icum", "Rn", "Ccum")
  
  transitions <- c(
    "S0 -> (1-caseinfection)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I + Icum",
    "S0 -> caseinfection*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
    "S1 -> (1-caseinfection)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I + Icum",
    "S1 -> caseinfection*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
    "S2 -> (1-caseinfection)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I + Icum",
    "S2 -> caseinfection*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
    "S3 -> (1-caseinfection)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I + Icum",
    "S3 -> caseinfection*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
    "S3r -> (1-caseinfection)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I + Icum",
    "S3r -> caseinfection*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
    "I -> gamma*I -> Rn",
    "C -> gamma*C -> Rn"
  )
  
  tspan <- seq(from = 1, to = f_time, by = 1) 
  
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = events, 
                E = Et, N = Nt,
                gdata = c(beta = beta, gamma = gamma, caseinfection = caseinfection), 
                u0 = initial_state, 
                tspan = tspan))
}

prepare_events <- function(current_time, sims_with_cases_node1, sims_with_cases_node2, importations_df, interval) {
  days <- (current_time + 1):(current_time + interval)  # Ensure days are within the current interval
  
  # Define events for node 1
  birthsS_events <- data.frame(event = "enter", time = rep(days, each = 1), node = 1, dest = 0, n = births_S, proportion = 0, select = 1, shift = 6)
  birthsS3_events <- data.frame(event = "enter", time = rep(days, each = 1), node = 1, dest = 0, n = births_S3, proportion = 0, select = 2, shift = 6)
  deaths_events <- data.frame(event = "exit", time = rep(days, each = 1), node = 1, dest = 0, n = death_rate1, proportion = 0, select = 3, shift = 6)
  
  # Define events for node 2
  births2S_events <- data.frame(event = "enter", time = rep(days, each = 1), node = 2, dest = 0, n = births2_S, proportion = 0, select = 1, shift = 6)
  births2S3_events <- data.frame(event = "enter", time = rep(days, each = 1), node = 2, dest = 0, n = births2_S3, proportion = 0, select = 2, shift = 6)
  deaths2_events <- data.frame(event = "exit", time = rep(days, each = 1), node = 2, dest = 0, n = death_rate2, proportion = 0, select = 3, shift = 6)
  
  # Movement events between node 1 and node 2
  move_1_to_2 <- data.frame(event = "intTrans", time = rep(days, each = 1), node = 1, dest = 2, n = 0, proportion = mobility1, select = 1, shift = 6)
  move_2_to_1 <- data.frame(event = "intTrans", time = rep(days, each = 1), node = 2, dest = 1, n = 0, proportion = mobility2, select = 1, shift = 6)
  
  # Define events for both nodes 1 and 2
  ipv_events <- data.frame(event = "intTrans", time = rep(days, each = 2), node = rep(c(1, 2), times = length(days)), dest = 0, n = 0, proportion = ipv_coverage, select = 14, shift = 5)
  seroconversionS3r_events <- data.frame(event = "intTrans", time = rep(days, each = 2), node = rep(c(1, 2), times = length(days)), dest = 0, n = 0, proportion = (1 -(1-ve)^3), select = 11, shift = 5)
  importation_events <- data.frame(event = "intTrans", time = rep(random_importations, each = 2), node = rep(c(1, 2), times = length(random_importations)), dest = 0, n = 1, proportion = 0, select = 4, shift = 4)
  
  # # Preventative campaigns can occur in both nodes
  # psia_annual_events <- data.frame(event = "intTrans", time = rep(c(1, 366, 731, 1096, 1461) + current_time, each = 2), node = rep(c(1, 2), times = 5), dest = 0, n = 0, proportion = sia_coverage, select = 5, shift = 1)
  # psia_no_seroconversion_events <- data.frame(event = "intTrans", time = rep(c(3, 368, 733, 1098, 1463) + current_time, each = 2), node = rep(c(1, 2), times = 5), dest = 0, n = 0, proportion = 1, select = 6, shift = 2)
  # psia_seroconversionS1_events <- data.frame(event = "intTrans", time = rep(c(2, 367, 732, 1097, 1462) + current_time, each = 2), node = rep(c(1, 2), times = 5), dest = 0, n = 0, proportion = (1 -(1-ve)^1), select = 7, shift = 3)
  # psia_seroconversionS2_events <- data.frame(event = "intTrans", time = rep(c(2, 367, 732, 1097, 1462) + current_time, each = 2), node = rep(c(1, 2), times = 5), dest = 0, n = 0, proportion = (1 -(1-ve)^2), select = 8, shift = 3)
  # psia_seroconversionS3_events <- data.frame(event = "intTrans", time = rep(c(2, 367, 732, 1097, 1462) + current_time, each = 2), node = rep(c(1, 2), times = 5), dest = 0, n = 0, proportion = (1 -(1-ve)^3), select = 9, shift = 3)
  # psia_seroconversionS4_events <- data.frame(event = "intTrans", time = rep(c(2, 367, 732, 1097, 1462) + current_time, each = 2), node = rep(c(1, 2), times = 5), dest = 0, n = 0, proportion = 1, select = 12, shift = 3)
  # psia_seroMAX_events <- data.frame(event = "intTrans", time = rep(c(3, 368, 733, 1098, 1463) + current_time, each = 2), node = rep(c(1, 2), times = 5), dest = 0, n = 0, proportion = 1, select = 13, shift = 3)
  # 
  all_events <- rbind(birthsS_events, birthsS3_events, deaths_events,
                      births2S_events, births2S3_events, deaths2_events,
                      move_1_to_2, move_2_to_1,
                      ipv_events, seroconversionS3r_events, importation_events)
                      # psia_annual_events, psia_no_seroconversion_events, psia_seroconversionS1_events, psia_seroconversionS2_events,
                      # psia_seroconversionS3_events, psia_seroconversionS4_events, psia_seroMAX_events)
  
  # Outbreak response events based on cases or infections in either node
  if (!is.null(c(sims_with_cases_node1, sims_with_cases_node2)) && length(c(sims_with_cases_node1, sims_with_cases_node2)) > 0) {
    obr_events_node1 <- data.frame(event = "intTrans", time = rep(current_time + interval, each = 1), node = 1, dest = 0, n = 0, proportion = sia_coverage, select = 5, shift = 1)
    no_seroconversion_events_node1 <- data.frame(event = "intTrans", time = rep(current_time + interval + 2, each = 1), node = 1, dest = 0, n = 0, proportion = 1, select = 6, shift = 2)
    seroconversionS1_events_node1 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 1, dest = 0, n = 0, proportion = (1 -(1-ve)^1), select = 7, shift = 3)
    seroconversionS2_events_node1 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 1, dest = 0, n = 0, proportion = (1 -(1-ve)^2), select = 8, shift = 3)
    seroconversionS3_events_node1 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 1, dest = 0, n = 0, proportion = (1 -(1-ve)^3), select = 9, shift = 3)
    seroconversionS4_events_node1 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 1, dest = 0, n = 0, proportion = 1, select = 12, shift = 3)
    seroMAX_events_node1 <- data.frame(event = "intTrans", time = rep(current_time + interval + 2, each = 1), node = 1, dest = 0, n = 0, proportion = 1, select = 13, shift = 3)
    
    obr_events_node2 <- data.frame(event = "intTrans", time = rep(current_time + interval, each = 1), node = 2, dest = 0, n = 0, proportion = sia_coverage, select = 5, shift = 1)
    no_seroconversion_events_node2 <- data.frame(event = "intTrans", time = rep(current_time + interval + 2, each = 1), node = 2, dest = 0, n = 0, proportion = 1, select = 6, shift = 2)
    seroconversionS1_events_node2 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 2, dest = 0, n = 0, proportion = (1 -(1-ve)^1), select = 7, shift = 3)
    seroconversionS2_events_node2 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 2, dest = 0, n = 0, proportion = (1 -(1-ve)^2), select = 8, shift = 3)
    seroconversionS3_events_node2 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 2, dest = 0, n = 0, proportion = (1 -(1-ve)^3), select = 9, shift = 3)
    seroconversionS4_events_node2 <- data.frame(event = "intTrans", time = rep(current_time + interval + 1, each = 1), node = 2, dest = 0, n = 0, proportion = 1, select = 12, shift = 3)
    seroMAX_events_node2 <- data.frame(event = "intTrans", time = rep(current_time + interval + 2, each = 1), node = 2, dest = 0, n = 0, proportion = 1, select = 13, shift = 3)
    
    all_events <- rbind(all_events, 
                        obr_events_node1, no_seroconversion_events_node1, seroconversionS1_events_node1, 
                        seroconversionS2_events_node1, seroconversionS3_events_node1, seroconversionS4_events_node1, seroMAX_events_node1,
                        obr_events_node2, no_seroconversion_events_node2, seroconversionS1_events_node2, 
                        seroconversionS2_events_node2, seroconversionS3_events_node2, seroconversionS4_events_node2, seroMAX_events_node2)
  }
  
  return(all_events)
}

# Initialize lists to store results for each node
model_results_node1 <- list()
model_results_node2 <- list()
# Initialize lists to track oSIA events for each node
all_osias_node1 <- list()
all_osias_node2 <- list()
# Initialize empty vectors at the beginning of your main loop
selected_sims_node1 <- integer(0)
selected_sims_node2 <- integer(0)

# Set up the main simulation loop for the specified number of repetitions
for (rep in 1:reps) {
  cat("Running repetition:", rep, "\n")

  initial_state <- data.frame(
    node = c(1, 2),  # Nodes 1 and 2 represent populations A and B
    S0 = c(((prop_S0 * pop1)-1), ((prop_S0 * pop2)-1)), #subtract 1 for initial infection
    S0v = c(0, 0),
    S1 = c(prop_S1 * pop1, prop_S1 * pop2),
    S1v = c(0, 0),
    S2 = c(prop_S2 * pop1, prop_S2 * pop2),
    S2v = c(0, 0),
    S3 = c(prop_S3 * pop1, prop_S3 * pop2),
    S3v = c(0, 0),
    S3r = c(0, 0),
    S4 = c(0, 0),
    Rv = c(prop_Rv * pop1, prop_Rv * pop2),
    Ii = c(0, 0),
    I = c(1, 1),
    C = c(0, 0),
    ES = c(0, 0),
    Icum = c(0, 0),
    Rn = c(0, 0),
    Ccum = c(0, 0)
  )
  
  # Initialize the model
  events <- prepare_events(current_time = current_time, sims_with_cases_node1 = selected_sims_node1, sims_with_cases_node2 = selected_sims_node2, importations_df = importations_df, interval = interval)
  model <- create.StocSIRobr0(I0 = I0, beta = beta, gamma = gamma, caseinfection = caseinfection, f_time = interval, initial_state = initial_state, events = events)
  
  current_time <- 0
  model_result_node1 <- list()
  model_result_node2 <- list()
  
  while (current_time < time_horizon) {
    next_time <- min(current_time + interval, time_horizon)
    f_time <- next_time - current_time
    
    # Ensure f_time is positive and does not exceed the time horizon
    if (f_time <= 0) break
    
    # Adjust the interval if it goes beyond the time_horizon
    if (next_time > time_horizon) {
      f_time <- time_horizon - current_time
    }
    
    results <- run(model, threads = 1)
    
    # Access the state of compartments for node 1 and node 2 separately
    state_node1 <- trajectory(results)[trajectory(results)$node == 1, ]
    state_node2 <- trajectory(results)[trajectory(results)$node == 2, ]
    
    # Adjust time in state to be continuous
    state_node1$time <- state_node1$time + current_time
    state_node2$time <- state_node2$time + current_time
    
    # Save the results
    model_result_node1[[length(model_result_node1) + 1]] <- state_node1
    model_result_node2[[length(model_result_node2) + 1]] <- state_node2
    
    # Collect simulations requiring oSIA based on 'C' and 'I' compartments
    # Node 1
    sims_with_cases_AFP_node1 <- which(state_node1[nrow(state_node1), "C"] > 0)
    if (length(sims_with_cases_AFP_node1) > 0) {
      num_selected_sims_AFP_node1 <- ceiling(length(sims_with_cases_AFP_node1) * afp_sens1)
      selected_sims_AFP_node1 <- sample(sims_with_cases_AFP_node1, num_selected_sims_AFP_node1)
    } else {
      selected_sims_AFP_node1 <- integer(0)
    }
    
    sims_with_infections_ES_node1 <- which(state_node1[nrow(state_node1), "I"] >= 2)
    if (length(sims_with_infections_ES_node1) > 0) {
      valid_sims_ES_node1 <- sapply(sims_with_infections_ES_node1, function(sim) {
        infection_times_node1 <- which(state_node1[, "I"] == sim)
        if (length(infection_times_node1) < 2) {
          return(FALSE)
        }
        any(diff(infection_times_node1) >= 31)
      })
      sims_with_detected_infections_node1 <- sims_with_infections_ES_node1[valid_sims_ES_node1]
      num_selected_sims_ES_node1 <- ceiling(length(sims_with_detected_infections_node1) * es_sens1)
      selected_sims_ES_node1 <- sample(sims_with_detected_infections_node1, num_selected_sims_ES_node1)
    } else {
      selected_sims_ES_node1 <- integer(0)
    }
    
    selected_sims_node1 <- unique(c(selected_sims_AFP_node1, selected_sims_ES_node1))
    
    if (length(selected_sims_node1) > 0) {
      all_osias_node1[[length(all_osias_node1) + 1]] <- data.frame(simulation = rep, time = current_time)
    }
    
    # Node 2
    sims_with_cases_AFP_node2 <- which(state_node2[nrow(state_node2), "C"] > 0)
    if (length(sims_with_cases_AFP_node2) > 0) {
      num_selected_sims_AFP_node2 <- ceiling(length(sims_with_cases_AFP_node2) * afp_sens2)
      selected_sims_AFP_node2 <- sample(sims_with_cases_AFP_node2, num_selected_sims_AFP_node2)
    } else {
      selected_sims_AFP_node2 <- integer(0)
    }
    
    sims_with_infections_ES_node2 <- which(state_node2[nrow(state_node2), "I"] >= 2)
    if (length(sims_with_infections_ES_node2) > 0) {
      valid_sims_ES_node2 <- sapply(sims_with_infections_ES_node2, function(sim) {
        infection_times_node2 <- which(state_node2[, "I"] == sim)
        if (length(infection_times_node2) < 2) {
          return(FALSE)
        }
        any(diff(infection_times_node2) >= 31)
      })
      sims_with_detected_infections_node2 <- sims_with_infections_ES_node2[valid_sims_ES_node2]
      num_selected_sims_ES_node2 <- ceiling(length(sims_with_detected_infections_node2) * es_sens2)
      selected_sims_ES_node2 <- sample(sims_with_detected_infections_node2, num_selected_sims_ES_node2)
    } else {
      selected_sims_ES_node2 <- integer(0)
    }
    
    selected_sims_node2 <- unique(c(selected_sims_AFP_node2, selected_sims_ES_node2))
    
    if (length(selected_sims_node2) > 0) {
      all_osias_node2[[length(all_osias_node2) + 1]] <- data.frame(simulation = rep, time = current_time)
    }
    
    # Update the current time
    current_time <- next_time
    
    # Prepare the events for the next interval
    events <- prepare_events(current_time = current_time, sims_with_cases_node1 = selected_sims_node1, sims_with_cases_node2 = selected_sims_node2, importations_df = importations_df, interval = interval)
    initial_state_node1 <- state_node1[nrow(state_node1), -1] # Exclude the time column for initial_state
    initial_state_node2 <- state_node2[nrow(state_node2), -1] # Exclude the time column for initial_state
    
    model_node1 <- create.StocSIRobr0(I0 = I0, beta = beta, gamma = gamma, caseinfection = caseinfection, f_time = f_time, initial_state = initial_state_node1, events = events)
    model_node2 <- create.StocSIRobr0(I0 = I0, beta = beta, gamma = gamma, caseinfection = caseinfection, f_time = f_time, initial_state = initial_state_node2, events = events)
  }
  
  # Combine the results for the current repetition
  model_results_node1[[rep]] <- do.call(rbind, model_result_node1)
  model_results_node2[[rep]] <- do.call(rbind, model_result_node2)
}

# Combine all repetitions into a single dataframe for each node
final_results_node1 <- do.call(rbind, model_results_node1)
final_results_node2 <- do.call(rbind, model_results_node2)

# Subset the data to keep only the first 5 years
final_results_node1 <- final_results_node1[final_results_node1$time <= 1825, ]
final_results_node2 <- final_results_node2[final_results_node2$time <= 1825, ]
final_results_node1$simulation <- rep(1:reps, each=1825)
final_results_node1$node <- 1
final_results_node2$simulation <- rep(1:reps, each=1825)
final_results_node2$node <- 2

setwd("/Users/idcvmauz/Dropbox/PolioPreventingOutbreaks/16_EconRisk/polio_siminf/ES/Variable_catchment/Catchment_500000")

# Save the final results to CSV files
output_filename_node1 <- paste0("popA_outputs_sens_", round(afp_sens1 * 100), ".csv")
output_filename_node2 <- paste0("popB_outputs_sens_", round(es_sens2 * 100), ".csv")
write.csv(final_results_node1, output_filename_node1, row.names = FALSE)
write.csv(final_results_node2, output_filename_node2, row.names = FALSE)

# Save the oSIA occurrences to separate CSV files for each node
if (length(all_osias_node1) > 0) {
  oSIA_df_node1 <- do.call(rbind, all_osias_node1)
  output_filename_osia_node1 <- paste0("popA_simulations_requiring_oSIA_", round(afp_sens1 * 100), ".csv")
  write.csv(oSIA_df_node1, output_filename_osia_node1, row.names = FALSE)
}

if (length(all_osias_node2) > 0) {
  oSIA_df_node2 <- do.call(rbind, all_osias_node2)
  output_filename_osia_node2 <- paste0("popB_simulations_requiring_oSIA_", round(es_sens2 * 100), ".csv")
  write.csv(oSIA_df_node2, output_filename_osia_node2, row.names = FALSE)
}

    
    
### SUMMARY STATS ###

# Summarize the total infections and cases over the entire time horizon 

# Calculate total Ccum and Icum for each simulation
summary_stats_A <- final_results_node1 %>%
  group_by(simulation) %>%
  summarize(
    total_c = max(Ccum, na.rm = TRUE),
    total_i = max(Icum, na.rm = TRUE)
  )

# Calculate mean and CIs for Ccum and Icum across all simulations
final_summary_A <- summary_stats_A %>%
  summarize(
    mean_cases = round(mean(total_c, na.rm = TRUE),0),
    cases_lwr = round(quantile(total_c, 0.05, na.rm = TRUE), 0),
    cases_upr = round(quantile(total_c, 0.95, na.rm = TRUE), 0),
    
    mean_infections = round(mean(total_i, na.rm = TRUE),0),
    infections_lwr = round(quantile(total_i, 0.05, na.rm = TRUE), 0),
    infections_upr = round(quantile(total_i, 0.95, na.rm = TRUE), 0)
  )

# Calculate total Ccum and Icum for each simulation
summary_stats_B <- final_results_node2 %>%
  group_by(simulation) %>%
  summarize(
    total_c = max(Ccum, na.rm = TRUE),
    total_i = max(Icum, na.rm = TRUE)
  )

# Calculate mean and CIs for Ccum and Icum across all simulations
final_summary_B <- summary_stats_B %>%
  summarize(
    mean_cases = round(mean(total_c, na.rm = TRUE),0),
    cases_lwr = round(quantile(total_c, 0.05, na.rm = TRUE), 0),
    cases_upr = round(quantile(total_c, 0.95, na.rm = TRUE), 0),
    
    mean_infections = round(mean(total_i, na.rm = TRUE),0),
    infections_lwr = round(quantile(total_i, 0.05, na.rm = TRUE), 0),
    infections_upr = round(quantile(total_i, 0.95, na.rm = TRUE), 0)
  )


outbreak_counts_A <- oSIA_df_node1 %>%
  group_by(simulation) %>%
  summarise(num_outbreaks = n())
round(mean(outbreak_counts_A$num_outbreaks),0)

outbreak_counts_B <- oSIA_df_node2 %>%
  group_by(simulation) %>%
  summarise(num_outbreaks = n())
round(mean(outbreak_counts_B$num_outbreaks),0)

