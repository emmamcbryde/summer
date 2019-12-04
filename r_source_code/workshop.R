source("summer_model.R")

my_model <- StratifiedModel$new(seq(from=0, to=60/365, by=1/365),
                                 c("susceptible", "infectious", "recovered"),
                                 c("infectious"=0.001),
                                 list(beta=400, recovery=365/13, infect_death=1),
                                 list(c("standard_flows", "recovery", "infectious", "recovered"),
                                      c("infection_density", "beta", "susceptible", "infectious"),
                                      c("compartment_death", "infect_death", "infectious")),
                                 report_progress=FALSE)

#my_model$stratify("hiv", c("negative", "positive"), c(),
#                   list(recovery=list("negative"=0.7, "positive"=0.5),
#                        infect_death=list("negative"=0.5)),
#                   list("negative"=0.6, "positive"=0.4), report = FALSE)

my_model$run_model()

interpreter <- ModelInterpreter$new(my_model)
interpreter$plot_compartment("infectious")
# interpreter$create_flowchart()

