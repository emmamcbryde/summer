
source("summer_model.R")


create_arbitrary_time_variant_function = function(time) {
  365 / 13 * exp(-time)
}


# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

sir_model <- StratifiedModel$new(seq(from=0, to=60/365, by=1/365),
                          c("susceptible", "infectious", "recovered"),
                          c("infectious"=0.001),
                          list(beta=400, recovery=365/13, infect_death=1),
                          list(c("standard_flows", "recovery", "infectious", "recovered"),
                               c("infection_density", "beta", "susceptible", "infectious"),
                               c("compartment_death", "infect_death", "infectious")),
                          verbose=FALSE)
sir_model$stratify("hiv", c("negative", "positive"), c(),
                    list(recovery=list("negative"=0.7, "positive"=0.5),
                         infect_death=list("negative"=0.5)),
                  list("negative"=0.6, "positive"=0.4), verbose = FALSE)
# 
# sir_model$stratify("age", c(3, 2, 1), c(),
#                     list(recovery=list("2"=1, "3"=2),
#                          recoveryXhiv_positive=list("1"=2, "2"=365/13*.5, "3"=1, overwrite=c("2")),
#                          universal_death_rate=list("1"=1, "2"=2, "3"=3)), verbose=TRUE)

# sir_model$stratify("age", c(3, 2, 1), c(), verbose = TRUE)
# sir_model$add_time_variant("recovery", create_arbitrary_time_variant_function)

sir_model$run_model()

interpreter <- ModelInterpreter$new(sir_model)
interpreter$plot_compartment("infectious")
# interpreter$create_flowchart()

