
source("summer_model.R")

# static functions for use later

add_vtp_latency_parameters = function(.parameters, change_time_unit=365.25) {
  # function to add the latency parameters estimated by Ragonnet et al from our paper in Epidemics to the existing parameter dictionary
  vtp_latency_parameters = c(early_progression = 1.1e-3, stabilisation = 1.0e-2, late_progression = 5.5e-6)
  for (parameter in names(vtp_latency_parameters)) {
    vtp_latency_parameters[[parameter]] <- vtp_latency_parameters[[parameter]] * change_time_unit
  }
  .parameters <- c(.parameters, vtp_latency_parameters)
}

get_age_specific_latency_parameters = function(.parameters, change_time_unit=365.25) {
  # get the age-specific latency parameters estimated by Ragonnet et al
  age_stratified_parameters = list(early_progression=list("0W"=6.6e-3, "5W"=2.7e-3, "15W"=2.7e-4),
                                   stabilisation=list("0W"=1.2e-2, "5W"=1.2e-2, "15W"=5.4e-3),
                                   late_progression=list("0W"=1.9e-11, "5W"=6.4e-6, "15W"=3.3e-6))
  return(age_stratified_parameters)
}

get_all_age_specific_latency_parameters = function(change_time_unit=365.25) {
  # collate all the latency parameters together from the previous function
  age_specific_parameters <- get_age_specific_latency_parameters()
  revised_parameters <- list()
  for (parameter in names(age_specific_parameters)) {
    for (stratum in names(age_specific_parameters[[parameter]])) {
      if (stratum != "overwrite") {
        revised_parameters[[parameter]][[stratum]] <- age_specific_parameters[[parameter]][[stratum]] * change_time_unit
      }
    }
  }
  return(revised_parameters)
}


add_standard_latency_flows = function(.flows) {
  # adds our standard latency flows to the list of flows to be implemented in the model
  .flows <- c(.flows, list(c("standard_flows", "early_progression", "early_latent", "infectious"),
                           c("standard_flows", "stabilisation", "early_latent", "late_latent"),
                           c("standard_flows", "late_progression", "late_latent", "infectious")))
}

# main section of code to run
case_fatality_rate <- 0.4
untreated_disease_duration <- 3
parameters <- list(beta = 10.0,
                   recovery = case_fatality_rate / untreated_disease_duration,
                   infect_death = (1.0 - case_fatality_rate) / untreated_disease_duration,
                   universal_death_rate = 1.0 / 50.0)
parameters <- add_vtp_latency_parameters(parameters)
times <- seq(0, 200)
flows <- list(c("infection_frequency", "beta", "susceptible", "early_latent"),
              c("infection_frequency", "beta", "recovered", "early_latent"),
              c("standard_flows", "recovery", "infectious", "recovered"),
              c("compartment_death", "infect_death", "infectious"))

flows <- add_standard_latency_flows(flows)

tb_model <- StratifiedModel$new(
  times, c("susceptible", "early_latent", "late_latent", "infectious", "recovered"), c(infectious=1e-3),
  parameters, flows, birth_approach="replace_deaths")


tb_model$stratify("age", c(5, 15), c(), get_age_specific_latency_parameters(), verbose=FALSE)

tb_model$run_model()

interpreter <- ModelInterpreter$new(tb_model)
interpreter$plot_compartment("infectious")

# print(tb_model$transition_flows)
# print(tb_model$parameters)



