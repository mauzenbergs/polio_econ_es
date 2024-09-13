This repo includes code for running the previously published polio SimInf model with parameters for surveillance 

Important changes to this model structure include: 
- Parameters for the proportion of infections that are detected via AFP and environmental surveillance (ES)
- For infections detected via ES, in order for an outbreak to be declared (as per standard operating procedures),
  the isolates must be detected more than one month apart and not be genetically linked, the model accounts for this
- The model allows for migration between two separate nodes (one node is a population within an ES catchment area
  the other population is outside the ES catchment area)
- The model also allows for different timeliness assumptions (time from sample collection through to outbreak response)
- The model takes inventory of the simulations that result in an outbreak (at least one paralytic case or two separate ES detections meeting criteria for ES),
  then depending on the proportion of infections that are actually detected (referred to here as surveillance sensitivity) the model triggers an oSIA in certain sims

How to run LSHTM’s polio SimInf model: A tutorial

This document will go line by line through the R script: “ES_oSIA.R” which is the code for running the transmission model that accounts for non-polio AFP and environmental surveillance assumptions, alongside outbreak response SIAs. 

Line(s) of code	and description

1-9	Load in required libraries and working directory

10-38	Define parameters – the parameters as they are set up assume there are two different populations (each with different population sizes and birth and death rates). 

40-52	Generate random times (days) that an importation of infection occurs. We want this to be a random event, so generating a random series of days (at a rate of 2 importations per year) is needed each time you run the model.

54-57	Timeliness assumptions – how many days between the detection of an infection or a case will an oSIA occur. Vary this parameter to model imperfect timeliness (i.e. interval = 74 days). We assume a 5 year time horizon to align with GPEI strategic plans, but this can be varied for other analyses.
59-61	Define the rate of movement between populations in the model

63-67	Define what proportion of infections and cases will actually be detected (and trigger an oSIA). In the analysis, we assume that AFP = ES sensitivity, but this can be varied here.

69-102	To run the SimInf model, model compartments and events must be listed within a ‘select’ and a ‘shift’ matrix. The select matrix selects which compartments are involved in each of the events and the shift matrix dictates how individuals shift between compartments (i.e. if shift = 3 the individual will move three compartments in the matrix).

104-132	This is the function for the underlying SimInf model. It specifies the transitions of the compartments and assumptions for initial conditions, transmission parameters, and scheduled events.

134-198	Prepare events function: this function includes all of the events that will occur in the model. Because we have two different populations in this model example, we define events separately for the two populations (nodes). We assume that every day of the model births and deaths occur as well as movement between nodes. If pSIAs were desired, lines 156-163 include coding for annual pSIAs. Lines 172-188 include code for oSIAs only in the simulations that triggered an outbreak response – it is necessary to have an explicit event for everything that happens within the vaccination process – i.e. getting vaccinated, not seroconverting, or potential for seroconversion based on previous doses received (S1, S2, S3, S4). These oSIA events are specified for each of the populations. Note that two ‘events’ cannot occur on the same day in the SimInf structure, this is why the vaccination event happens on day t and seroconversion happens on day t+1.

200-208	Initialise lists to store the results for each simulation

210-234	Set up the main simulation loop based on the number of repetitions desired (i.e. 10,000) and define the initial state assumptions (here we assume that based on RI coverage a certain proportion of the population is already immune). We introduce 1 infection at the start of the model.

Main simulation loop	It is important to understand the limitations of the SimInf package. For example, the model framework can only run scheduled in pre-specified simulations or nodes. This is challenging because for oSIAs, we do not know ahead of time which simulations will result in an outbrea. Therefore, we need to stop and start the model continuously (based on the defined time interval, i.e. 35 days). We then take inventory of which simulations had outbreak and we tell the model to conduct oSIA events only in these simulations. For this reason, the main simulation loop stores model outputs for each compartment during each time interval for which the model is stopped. It also records the simulations that required an oSIA to keep track of the number of oSIAs required.

244-269	Checks to make sure the model runs for only the specified time horizon, the model runs continuously and the model results are stored each time the model is stopped and started

270-331	For each of the nodes (corresponding to populations A and B in this example), we must specify the conditions for an outbreak. In these lines you will see that sims_with_infections_ES must have at least 2 infections that are at least 31 days apart. This detection also depends on ‘es_sens.’ For the sims_with_cases_AFP, the simulations must have at least 1 case, and the proportion of the cases that are detected depends on ‘afp_sens.’ If simulations meet the criteria for an outbreak, then an oSIA is triggered only in these simulations. We store this data through the entire model to keep track of which simulations required an oSIA and when. We then loop through this entire process for the number of specified ‘reps.’

350-381	We then combine all data stored across all simulations into one output file for saving. We store compartmental data separately to oSIA occurrences.

385-end	Examples of summary statistics that can be run once the model finishes.
