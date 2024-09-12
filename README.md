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
