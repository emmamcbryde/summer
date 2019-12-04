
# SUMMER
# Scalable
# Universal
# Mathematical
# Model
# for Epidemics
# in R

library(deSolve)
library(R6)
# library(tidyverse)
# library(DiagrammeR)
# library(DiagrammeRsvg)
library(rsvg)
library(stringr)



# library(RSQLite)
# this file contains the main model builder function, that is intended to be agnostic
# to the type and complexity of model that the user wants to build
# all instructions for what the characteristics of the model are are separated out to a file that calls/sources this one

# outstanding tasks
# faster calculation of parameters to avoid repeatedly multiplying the non-time-variant quantities together at each time step
# extend functionality - automatic age stratification, heterogeneous mixing, different infectivity for different compartments, stochastic implementation

# global parameter for checking recursion depth
depth=0


# ------------------------------
# string manipulation functions
# -------------------------------

create_stratum_name = function(stratification_name, stratum_name, joining_string="X") {
  # generate a name string to represent a particular stratum within a requested stratification
  
  # :param stratification_name: str
  # the "stratification" or rationale for implementing the current stratification process
  # :param stratum_name: str
  # name of the stratum within the stratification
  # :param joining_string: logical
  # the character to add to the front to indicate that this string is the extension of the existing one
  # in SUMMER, capitals are reserved for non-user-requested strings, in this case "X" is used as the default
  # :return: str
  # the composite string for the stratification
  
  startum_name = paste(stratification_name, "_", stratum_name, sep="")
  ifelse(!is.na(joining_string), paste(joining_string, stratification_name, "_", stratum_name, sep = ""), startum_name)
}


create_stratified_name = function(stem, stratification_name, stratum_name) {
  # generate a standardised stratified compartment name
  
  # :param stem: str
  # the previous stem to the compartment or parameter name that needs to be extended
  # :param stratification_name: str
  # the "stratification" or rationale for implementing the current stratification process
  # :param stratum_name: str
  # name of the stratum within the stratification
  # :return: str
  # the composite name with the standardised stratification name added on to the old stem
  
  paste(stem, create_stratum_name(stratification_name, stratum_name), sep = "")
}



extract_x_positions = function(parameter, joining_string="X") {
  # find the positions within a string which are X and return as list, including length of list
  
  # :param parameter: str
  # the string for interrogation
  # :param joining_string: str
  # the string of interest whose character positions need to be found
  # :return: list
  # list of all the indices for where the X character occurs in the string, along with the total length of the list
  
  result = c()
  param = strsplit(parameter, "")[[1]]
  for (loc in seq(nchar(parameter))){
    if (as.character(param[loc]) == as.character(joining_string))
      result = c(result, loc)}
  result = c(result, length(param))
  result
}


extract_reversed_x_positions = function(parameter) {
  # find the positions within a string which are X and return as list reversed, including length of list
  
  # :params and return: see extract_x_positions above
  rev(extract_x_positions(parameter))
}



find_stem = function(stratified_string, joining_string="X") {
  # find the stem of the compartment name as the text leading up to the first occurrence of the joining string
  # (usually "X")
  # should run slightly faster than using find_name_components
  # 
  # :param stratified_string: str
  # the stratified string for the compartment or parameter name
  # :param joining_string: str
  # the string of interest whose character position will be used to truncate the stratified string
  # :return: int
  # the point at which the first occurrence of the joining string occurs
  find_name_components(stratified_string)[1]
}

find_stratum_index_from_string = function(compartment, stratification, remove_stratification_name=TRUE){
  # finds the stratum which the compartment (or parameter) name falls in when provided with the compartment name and the
  # name of the stratification of interest
  # for example, if the compartment name was infectiousXhiv_positiveXdiabetes_none and the stratification of interest
  # provided through the stratification argument was hiv, then
  # :param compartment: str
  # name of the compartment or parameter to be interrogated
  # :param stratification: str
  # the stratification of interest
  # :param remove_stratification_name: bool
  # whether to remove the stratification name and its trailing _ from the string to return
  # :return: str
  # the name of the stratum within which the compartment falls
  stratum_name <- c()
  for (name in unlist(find_name_components(compartment))){
    if (grepl(stratification, name)){
      stratum_name <- c(stratum_name, name)
    }
  }
  
  stratum <- ifelse(remove_stratification_name, unlist(strsplit(as.character(stratum_name[1]),"_"))[2], stratum_name[1]) 
  stratum
}

add_w_to_param_names = function(parameter_dict){
  # add a "W" string to the end of the parameter name to indicate that we should over-write up the chain 
  # :param parameter_dict: dict
  # the dictionary before the adjustments
  # :return: dict
  # same dictionary but with the "W" string added to each of the keys
  
}

find_name_components = function(compartment_or_parameter, joining_string="X"){
  # extract all the components of a stratified compartment or parameter name
  # :param compartment: str
  # name of the compartment or parameter to be interrogated
  # :return:
  #  list of the extracted compartment components
  
  component_list = unlist(strsplit(compartment_or_parameter, joining_string))
  component_list 
}

# ----------------------------------
# basic data manipulation functions
# ----------------------------------

increment_list_by_index = function(list_to_increment, index_to_increment, increment_value){
  # general method to increment the odes by a value specified as an argument
  
  # :param list_to_increment: list
  # the list to be incremented
  # :param index_to_increment: int
  # the index of the list that needs to be incremented
  # :param increment_value:
  #  the value to increment the list by
  # :return: list
  # the list after it has been incremented
  
  #if(index_to_increment != NA){
  #print(increment_value)
  if (!(identical(increment_value, numeric(0)) | identical(increment_value, integer(0))))
    list_to_increment[index_to_increment] = list_to_increment[index_to_increment] + increment_value
  
  
  
  list_to_increment
}


normalise_list = function(value_list) {
  # simple function to normalise the values from a list
  
  # :param value_dict: dict
  # dictionary whose values will be adjusted
  # :return: dict
  # same dictionary after values have been normalised to the total of the original values
  
  total <- sum(as.numeric(value_list))
  if (total != 1) {
    value_list <- lapply(value_list, function(value) value / total)
  }
  value_list
}

change_parameter_unit = function(parameter_dict, multiplier){
  # used to adapt the parameters according their unit - for example, could be used for models that are running in time
  # steps that are different from the time step assumed by the input parameter
  # 
  # :param parameter_dict: dict
  # dictionary whose values need to be adjusted
  # :param multiplier: float
  # multiplier
  # :return: dict
  # dictionary with values multiplied by the multiplier argument
  
}


element_list_multiplication = function(list_1, list_2){
  # multiply elements of two lists to return another list with the same dimensions
  
  # :param list_1: list
  # first list of numeric values to be multiplied
  # :param list_2: list
  # first list of numeric values to be multiplied
  # :return: list
  # list populated with multiplied values
  list = list_1 * list_2
  list
}

element_list_division = function(list_1, list_2){
  # divide elements of two lists to return another list with the same dimensions
  
  # :param list_1: list
  # first list of numeric values for numerators of division
  # :param list_2: list
  # first list of numeric values for denominators of division
  # :return: list
  # list populated with divided values
  list = list_1/list_2
  list
}

# ----------------------------------------
# functions for use as inputs to the model

# considering moving to different file
# ----------------------------------------

create_step_function_from_dict = function(input_dict){
  # create a step function out of dictionary with numeric keys and values, where the keys determine the values of the
  # independent variable at which the steps between the output values occur
  
  # :param input_dict: dict
  # dictionary in standard format as described above
  # :return: function
  # the function constructed from input data
  
}


# find the trailing text for the stratum of the compartment
find_stratum = function(stratified_string) {
  if (grepl("X", stratified_string)) {
    stratum <- substr(stratified_string, gregexpr(pattern="X", stratified_string)[[1]][1], 100)
  }
  else {
    ""
  }
}

# function to chop off the last piece of the character vector after the last occurrence of X (i.e. the last stratification)
remove_last_stratification = function(stratified_string) {
  x_positions <- gregexpr(pattern = "X", stratified_string)[[1]]
  substr(stratified_string, 1, x_positions[length(x_positions)] - 1)
}





# string is cleaned up for presentation
capitalise_compartment_name = function(compartment) {
  compartment_capitalised <- compartment %>% str_replace_all('X', ' ') %>%
    str_replace_all('_', ' ') %>% str_to_title()
}




# general method to increment the odes by a value specified as an argument
increment_compartment = function(ode_equations, compartment_number, increment) {
  ode_equations[compartment_number] <- ode_equations[compartment_number] + increment
  ode_equations
}



create_multiplicative_function = function(multiplier) {
  
  multiplicative_function = function(value, time) {
    # print(paste("multiplier = ", multiplier, " value= ", value))
    multiplier * unlist(unname(value))
  }
  multiplicative_function
  
}



create_function_of_function = function(outer_function, inner_function, max_depth=4) {
  function_to_return = function(depth) {
    if (depth > max_depth){
      return(1)
    }
    
    depth = depth + 1
    
    return(outer_function(inner_function()))
  }
  return(function_to_return(depth+1))
}



EpiModel <- R6Class(
  "EpiModel",
  
  #   general epidemiological model for constructing compartment-based models, typically of infectious disease
  #   transmission. See README.md for full description of purpose and approach of this model.
  #
  #   :attribute times: list
  #       time steps at which outputs are to be evaluated
  #   :attribute compartment_types: list
  #       strings representing the compartments of the model
  #   :attribute initial_conditions: dict
  #       keys are compartment types, values are starting population values for each compartment
  #       note that not all compartment_types must be included as keys here
  #   :attribute parameters: dict
  #       constant parameter values
  #   :attribute requested_flows: list of dicts in standard format
  #       list with each element being a model flow, with fixed key names according to the type of flow implemented
  #   :attribute initial_conditions_to_total: bool
  #       whether to add the initial conditions up to a certain total if this value hasn't yet been reached through
  #       the initial_conditions argument
  #   :attribute infectious_compartment: str
  #       name of the infectious compartment for calculation of intercompartmental infection flows
  #   :attribute birth_approach: str
  #       approach to allowing entry flows into the model, must be add_crude_birth_rate, replace_deaths or no_births
  #   :attribute verbose: bool
  #       whether to output progress in model construction as this process proceeds
  #   :attribute reporting_sigfigs: int
  #       number of significant figures to output to when reporting progress
  #   :attribute entry_compartment: str
  #       name of the compartment that births come in to
  #   :attribute starting_population: numeric
  #       value for the total starting population to be supplemented to if initial_conditions_to_total requested
  #   :attribute starting_compartment: str
  #       optional name of the compartment to add population recruitment to
  #   :attribute equilibrium_stopping_tolerance: float
  #       value at which relative changes in compartment size trigger stopping when equilibrium reached
  #   :attribute integration_type: str
  #       integration approach for numeric solution to odes, must be odeint or solveivp currently
  
  # attributes that shouldn't be changed by the user
  public = list(
    
    # attributes that are fed in as inputs (so defaults can be set as arguments to the initialisation method)
    times = NULL,
    compartment_types = NULL,
    compartment_names = NULL,
    initial_conditions = NULL,
    parameters = NULL,
    initial_conditions_to_total = NULL,
    infectious_compartment = NULL,
    birth_approach = NULL,
    verbose = NULL,
    reporting_sigfigs = NULL,
    entry_compartment = NULL,
    starting_population = NULL,
    starting_compartment = NULL,
    equilibrium_stopping_tolerance = NULL,
    output_connections = NULL,
    tracked_quantities = NULL,
    
    # other required attributes
    compartment_values = c(),
    unstratified_flows = data.frame(),
    transition_flows = data.frame(),
    death_flows = data.frame(),
    time_variants = list(),
    outputs = NULL,
    mixing_categories = list(),
    infectiousness_adjustments = c(),
    derived_outputs = list(),
    flows_to_track = list(),
    infectious_indices = list(),
    infectiousness_multipliers = c(),
    infectious_compartments = c(),
    
    # __________
    # most general methods
    
    find_parameter_value = function(parameter_name, time) {
      #   find the value of a parameter with time-variant values trumping constant ones
      #
      #   :param parameter_name: string for the name of the parameter of interest
      #   :param time: model integration time (if needed)
      #   :return: parameter value, whether constant or time variant
      if (parameter_name %in% names(self$time_variants)) {
        return(self$time_variants[[parameter_name]](time))
      }
      else {
        return(self$parameters[[parameter_name]])
      }
    },
    
    output_to_user = function(comment) {
      #   short function to save the if statement in every call to output some information, may be adapted later and was
      #   more important to the R version of the repository
      #
      #   :param: comment: string for the comment to be displayed to the user
      if (self$verbose & is.character(comment)) {
        writeLines(comment)
      }
      else if (self$verbose) {
        print(comment)
      }
    },
    
    # __________
    # model construction methods
    
    # initialise basic model characteristics from inputs and check appropriately requested
    initialize = function(times, compartment_types, initial_conditions, parameters, requested_flows,
                          initial_conditions_to_total=TRUE, infectious_compartment="infectious", birth_approach="no_births",
                          verbose=FALSE, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                          starting_compartment="", equilibrium_stopping_tolerance=1e6, output_connections=list(), tracked_quantities=list()) {
      #   construction method to create a basic (and at this stage unstratified) compartmental model, including checking
      #   that the arguments have been provided correctly (in a separate method called here)
      #
      #   :params: all arguments essentially become object attributes and are described in the first main docstring to
      #       this object class
      
      # ensure requests are fed in correctly
      self$check_and_report_attributes(
        times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
        infectious_compartment, birth_approach, verbose, reporting_sigfigs, entry_compartment, starting_population,
        starting_compartment, equilibrium_stopping_tolerance, output_connections, tracked_quantities)
      
      # convert input arguments to model attributes
      for (attribute_to_assign in c(
        "times", "compartment_types", "initial_conditions", "parameters", "infectious_compartment",
        "birth_approach", "verbose", "reporting_sigfigs", "entry_compartment", "starting_population",
        "starting_compartment", "infectious_compartment", "equilibrium_stopping_tolerance", "output_connections",
        "tracked_quantities")) {
        self[[attribute_to_assign]] <- get(attribute_to_assign)
      }
      
      # set initial conditions and implement flows
      self$set_initial_conditions(initial_conditions_to_total)
      
      # implement unstratified flows
      self$implement_flows(requested_flows)
      
      # add any missing quantities that will be needed
      self$add_default_quantities()
      
      # find the compartments that are infectious
      #self$infectious_indices <- as.logical(grepl(self$infectious_compartment , self$compartment_names))
      #self$infectious_indices_int <- which(self$infectious_indices)
    },
    
    
    
    
    check_and_report_attributes = function(
      .times, .compartment_types, .initial_conditions, .parameters, .requested_flows, .initial_conditions_to_total,
      .infectious_compartment, .birth_approach, .verbose, .reporting_sigfigs, .entry_compartment, .starting_population,
      .starting_compartment, .equilibrium_stopping_tolerance, .output_connections, .tracked_quantities) {
      #   check all input data have been requested correctly
      #
      #   :parameters: all parameters have come directly from the construction (initialize) method unchanged and have been
      #       renamed with a preceding . character
      
      # check that variables are of the expected type
      for (expected_numeric_variable in c(".times", ".reporting_sigfigs", ".starting_population", ".equilibrium_stopping_tolerance")) {
        if (!is.numeric(get(expected_numeric_variable))) {
          stop(c("expected numeric input for ", expected_numeric_variable))
        }
      }
      for (expected_string_variable in c(".compartment_types", ".infectious_compartment", ".birth_approach", ".entry_compartment",
                                         ".starting_compartment")) {
        if (!is.character(get(expected_string_variable))) {
          stop(c("expected character input for ", expected_string_variable))
        }
      }
      for (expected_boolean_variable in c(".initial_conditions_to_total", ".verbose")) {
        if (!is.logical(get(expected_boolean_variable))) {
          stop(c("expected boolean for ", expected_boolean_variable))
        }
      }
      for (expected_list_variable in c(".requested_flows", ".output_connections", ".tracked_quantities")) {
        if (!is.list(get(expected_list_variable))) {
          stop(c("expected list for ", expected_list_variable))
        }
      }
      
      # check some specific requirements
      if (!.infectious_compartment %in% .compartment_types) {
        stop("infectious compartment name is not one of the listed compartment types")
      }
      if (!.birth_approach %in% c("add_crude_birth_rate", "replace_deaths", "no_births")) {
        stop("requested birth approach unavailable")
      }
      if (is.unsorted(self$times)) {
        self$output_to_user("requested integration times are not sorted, now sorting")
        self$times <- sort(self$times)
      }
      
      # report on characteristics of inputs
      if (.verbose) {
        writeLines(paste("\nintegrating times are from ", round(.times[1], .reporting_sigfigs),
                         " to ", round(tail(.times, 1), .reporting_sigfigs), " time units are always arbitrary", sep=""))
        writeLines("\nunstratified requested initial conditions are:")
        for (compartment in names(.initial_conditions)) {
          writeLines(paste(compartment, ": ",
                           as.character(round(as.numeric(.initial_conditions[compartment]), .reporting_sigfigs)), sep=""))
        }
        writeLines("\nunstratified parameter values are:")
        for (parameter in names(.parameters)) {
          writeLines(paste(parameter, ": ", as.character(round(as.numeric(.parameters[parameter]), .reporting_sigfigs)), sep=""))
        }
        writeLines(paste("\ninfectious compartment is called:", .infectious_compartment))
        writeLines(paste("\nbirth approach is:", .birth_approach))
      }
    },
    
    set_initial_conditions = function(.initial_conditions_to_total) {
      #   set starting compartment values
      #
      #   :param .initial_conditions_to_total: logical
      #       unchanged from argument to initialize
      
      # keep copy of the compartment types for when the compartment names are stratified later
      self$compartment_names <- self$compartment_types
      
      # start from making sure all compartments are set to zero values
      self$compartment_values <- rep(0, length(self$compartment_names))
      names(self$compartment_values) <- self$compartment_names
      
      # set starting values of unstratified compartments to requested value
      for (compartment in names(self$initial_conditions)) {
        if (compartment %in% self$compartment_types) {
          self$compartment_values[compartment] <- self$initial_conditions[compartment]
        }
        else {
          stop(paste("compartment", compartment, "requested in initial conditions not found in model compartment types"))
        }
      }
      
      # sum to a total value if requested
      if (.initial_conditions_to_total) {
        self$sum_initial_compartments_to_total()
      }
    },
    
    sum_initial_compartments_to_total = function() {
      # make initial conditions sum to a certain value
      
      compartment <- self$find_remainder_compartment()
      remaining_population <- self$starting_population - Reduce("+", self$compartment_values)
      if (remaining_population < 0) {
        stop("total of requested compartment values is greater than the requested starting population")
      }
      self$output_to_user(paste("requested that total population sum to", self$starting_population))
      self$output_to_user(paste("remaining population of ", as.character(round(remaining_population, self$reporting_sigfigs)),
                                " allocated to ", compartment, " compartment", sep=""))
      self$compartment_values[compartment] <- remaining_population
    },
    
    find_remainder_compartment = function() {
      #   find the compartment to put the remaining population that hasn't been assigned yet when summing to total
      #
      #   :return: char
      #       name of the compartment to assign the remaining population size to
      
      if (nchar(self$starting_compartment) > 0 & !self$starting_compartment %in% self$compartment_types) {
        stop("starting compartment to populate with initial values not found in available compartments")
      }
      else if (nchar(self$starting_compartment) > 0) {
        return(self$starting_compartment)
      }
      else {
        self$output_to_user(paste(
          "no default starting compartment requested for unallocated population,",
          "so will be allocated to entry compartment", self$entry_compartment))
        return(self$entry_compartment)
      }
    },
    
    implement_flows = function(.requested_flows) {
      #   add all flows to create data frames from input lists
      #
      #   :param _requested_flows: list
      #     unchanged from argument to initialize
      
      for (flow in .requested_flows) {
        
        # check flow requested correctly
        if(!flow[2] %in% names(self$parameters)) {
          stop("flow parameter not found in parameter list")
        }
        if(!flow[3] %in% self$compartment_types) {
          stop("from compartment name not found in compartment types")
        }
        if(length(flow) > 3 & !flow[4] %in% self$compartment_types) {
          stop("to compartment name not found in compartment types")
        }
        
        # add flow to appropriate data frame
        if (flow[1] == "compartment_death") {
          self$add_death_flow(flow)
        }
        else {
          self$add_transition_flow(flow)
        }
        # add any tracked quantities that will be needed for calculating flow rates during integration
        if (grepl("infection", flow[1])) {
          self$tracked_quantities$infectious_population <- 0
        }
        if (flow[1] == "infection_frequency") {
          self$tracked_quantities$total_population <- 0
        }
        
      }
    },
    
    add_default_quantities = function() {
      #   add parameters and tracked quantities that weren't requested but will be needed
      
      # universal death rate
      if (!"universal_death_rate" %in% names(self$parameters)) {
        self$parameters$universal_death_rate <- 0
      }
      
      # birth approach-specific parameters
      if (self$birth_approach == "add_crude_birth_rate" & !"crude_birth_rate" %in% names(self$parameters)) {
        self$parameters <- c(self$parameters, c(crude_birth_rate = 0))
      }
      else if (self$birth_approach == "replace_deaths") {
        self$tracked_quantities$total_deaths <- 0
      }
      
      # for each derived output to be recorded, initialise a tracked quantities key to zero
      for (output in names(self$output_connections)) {
        self$tracked_quantities[[output]] <- 0
      }
      
      # parameters essential for latter stratification if called
      self$parameters$entry_fractions <- 1
    },
    
    add_transition_flow = function(.flow) {
      
      #   add a flow (row) to the data frame storing the flows
      if ('force_index' %in% colnames(self$transition_flows) | 'strain' %in% colnames(self$transition_flows))
        self$transition_flows <-
          rbind(self$transition_flows, data.frame(type=.flow[1], parameter=as.character(.flow[2]), from=.flow[3], to=.flow[4],
                                                  implement=length(self$stratifications), strain=NA, force_index=1, stringsAsFactors=FALSE))
      else
        self$transition_flows <-
          rbind(self$transition_flows, data.frame(type=.flow[1], parameter=as.character(.flow[2]), from=.flow[3], to=.flow[4],
                                                  implement=length(self$stratifications), stringsAsFactors=FALSE))
    },
    
    add_death_flow = function(.flow) {
      #   similarly for compartment-specific death flows
      self$death_flows <- rbind(self$death_flows, data.frame(type=.flow[1], parameter=as.character(.flow[2]), from=.flow[3],
                                                             implement=0, stringsAsFactors=FALSE))
    },
    
    find_all_infectious_indices = function(){
      # find all the compartment names that begin with one of the requested infectious compartments
      
      # :return: list
      # booleans for whether each compartment is infectious or not
      infectious_compartment_index = c()
      i_comp = 1
      for (comp in self$compartment_names)
        i_comp = i_comp + 1
      if (find_stem(comp) %in% self$infectious_compartment)
        infectious_compartment_index = c(infectious_compartment_index, i_comp)
      return(infectious_compartment_index)
    },
    
    prepare_to_run = function(){
      # primarily for use in the stratified version when over-written
      # here just find all of the compartments that are infectious and prepare some list indices to speed integration 
      self$infectious_indices = self$find_all_infectious_indices()
      self$transition_indices_to_implement = self$find_transition_indices_to_implement()
      self$death_indices_to_implement = self$find_death_indices_to_implement()
    },
    
    # __________
    # methods for model running
    
    
    
    run_model = function () {
      #   main function to integrate model odes, called externally in the master running script
      self$output_to_user("now integrating")
      self$prepare_to_run()
      self$prepare_stratified_parameter_calculations()
      
      #   unlike python version, only one integration solver included so far
      self$outputs <- as.data.frame(lsodar(self$compartment_values, self$times, self$make_model_function(),
                                           rootfunc = self$set_stopping_conditions()))
      self$output_to_user("integration complete")
    },
    
    prepare_stratified_parameter_calculations = function() {
      #   for use in the stratified version only
    },
    
    make_model_function = function() {
      #   create derivative function, different approach to python because only one solver implemented here
      epi_model_function <- function(time, compartment_values, parameters) {
        
        # update all the emergent model quantities needed for integration
        self$update_tracked_quantities(compartment_values)
        
        # apply flows
        self$apply_all_flow_types_to_odes(rep(0, length(compartment_values)), compartment_values, time)
      }
    },
    
    set_stopping_conditions = function() {
      #   if requested to stop when equilibrium is reached
      
      # never stop because impossible to find root
      if (is.null(self$equilibrium_stopping_tolerance)) {
        epi_model_function <- function(time, compartment_values, parameters) {1}
      }
      
      # stop once largest net compartment change falls below specified threshold
      else {
        epi_model_function <- function(time, compartment_values, parameters) {
          self$update_tracked_quantities(compartment_values)
          net_flows <- self$apply_all_flow_types_to_odes(rep(0, length(compartment_values)), compartment_values, time)
          max(abs(unlist(net_flows))) - self$equilibrium_stopping_tolerance
        }
      }
    },
    
    apply_all_flow_types_to_odes = function(.ode_equations, .compartment_values, .time) {
      #   apply all flow types to a vector of zeros (note deaths must come before births in case births replace deaths)
      #
      #   :param .ode_equations: numeric vector
      #       comes in as a list of zeros with length equal to that of the compartments
      #   :param .compartment_values: numeric vector
      #       working values of the compartment sizes
      #   :param .time: numeric scalar value
      #       current integration time
      #   :return: ode equations as list
      #       updated ode equations in same format but with all flows implemented
      .ode_equations <- self$apply_transition_flows(.ode_equations, .compartment_values, .time)
      if (nrow(self$death_flows) > 0) {
        .ode_equations <- self$apply_compartment_death_flows(.ode_equations, .compartment_values, .time)
      }
      .ode_equations <- self$apply_universal_death_flow(.ode_equations, .compartment_values, .time)
      .ode_equations <- self$apply_birth_rate(.ode_equations, .compartment_values)
      list(.ode_equations)
    },
    
    apply_transition_flows = function(.ode_equations, .compartment_values, .time) {
      #   apply fixed or infection-related intercompartmental transition flows to odes
      #
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      for (n_flow in which(self$transition_flows$implement == length(self$strata))) {
        #for (n_flow in  self$find_transition_indices_to_implement()) {
        
        
        # find adjusted parameter value
        parameter_value <- self$get_parameter_value(self$transition_flows$parameter[n_flow], .time)
        # print(self$transition_flows$parameter[n_flow])
        #print('parameter_value')
        #print(parameter_value)
        
        # find from compartment and "infectious population" (which is 1 for standard flows)
        infectious_population <- self$find_infectious_multiplier(n_flow)
        # print(paste('infectious_population for flow ', n_flow) )
        # print(infectious_population)
        
        # calculate the flow and apply to the odes
        from_compartment <- match(self$transition_flows$from[n_flow], names(self$compartment_values))
        net_flow <- unlist(unname(parameter_value)) * .compartment_values[[from_compartment]] * infectious_population
        # print(paste('net_flow',net_flow) )
        
        .ode_equations <- increment_list_by_index(.ode_equations, from_compartment, -net_flow)
        .ode_equations <- increment_list_by_index(.ode_equations, match(self$transition_flows$to[n_flow], names(self$compartment_values)), net_flow)
        
        # track any quantities dependent on flow rates
        self$track_derived_outputs(flow, net_flow)
      }
      
      # add another element to the derived outputs vector
      self$extend_derived_outputs(.time)
      
      # return flow rates
      .ode_equations
    },
    
    find_transition_indices_to_implement = function(){
      #for over-writing in stratified version, here just returns the indices of all the transition flows, as they all
      #need to be implemented
      
      #:return: list
      #integers for all the rows of the transition matrix
      c(which(self$transition_flows$implement == length(self$strata)))
    },
    
    find_death_indices_to_implement = function(){
      #for over-writing in stratified version, here just returns the indices of all the transition flows, as they all
      #need to be implemented
      
      #:return: list
      #integers for all the rows of the transition matrix
      c(seq(nrow(self$death_flows)))
    },
    
    find_net_transition_flow = function(n_flow, .time, .compartment_values){
      # common code to finding transition flows during and after integration packaged into single function
      
      # :param n_flow: int
      # row of interest in transition flow dataframe
      # :param _time: float
      # time step, which may be time during integration or post-integration time point of interest
      # :param _compartment_values: list
      # list of current compartment sizes
      # :return: float
      # net transition between the two compartments being considered
      
      # find adjusted parameter value
      parameter_value = self$get_parameter_value(self$transition_flows$parameter[n_flow], .time)
      
      # find from compartment and the "infectious population" (which equals one for non-infection-related flows)
      infectious_population = self$find_infectious_multiplier(n_flow)
      
      # find the index of the origin or from compartment
      from_compartment = self$compartment_names$index(self$transition_flows$origin[n_flow])
      
      if(self$transition_flows$type[n_flow] == "customised_flows"){
        return(parameter_value * self$customised_flow_functions[n_flow](n_flow))
      } else {
        parameter_value * .compartment_values[from_compartment] * infectious_population
      }
    },
    
    track_derived_outputs = function(.flow, .net_flow) {
      #   calculate derived quantities to be tracked, which are stored as the self.derived_outputs dictionary for the
      #   current working time step
      #
      #   :param .flow: row of dataframe
      #       row of the dataframe representing the flow being considered in the preceding method
      #   :param .net_flow: numeric scalar value
      #       previously calculated magnitude of the transition flow
      for (output_type in names(self$output_connections)) {
        if (grepl(self$output_connections[[output_type]]["from"], .flow$from) &
            grepl(self$output_connections[[output_type]]["to"], .flow$to)){
          self$tracked_quantities[[output_type]] <- self$tracked_quantities[[output_type]] + as.numeric(.net_flow)
        }
      }
    },
    
    extend_derived_outputs = function(.time) {
      #   add the derived quantities being tracked to the end of the tracking vector, taking the self.derived_outputs
      #   dictionary for a single time point and updating the derived outputs dictionary of lists for all time points
      #
      #   :param .time: numeric scalar value
      #       current time in integration process
      self$derived_outputs$times <- c(self$derived_outputs$times, .time)
      for (output_type in names(self$output_connections)) {
        self$derived_outputs[[output_type]] <- c(self$derived_outputs[[output_type]], self$tracked_quantities[[output_type]])
      }
    },
    
    apply_compartment_death_flows = function(.ode_equations, .compartment_values, .time) {
      #   equivalent method to for transition flows above, but for deaths
      #
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      for (f in which(self$death_flows$implement == length(self$strata))) {
        flow <- self$death_flows[f,]
        parameter_value <- self$get_parameter_value(flow$parameter, .time)
        from_compartment <- match(flow$from, names(self$compartment_values))
        net_flow <- unlist(unname(parameter_value)) * .compartment_values[[from_compartment]]
        .ode_equations <- increment_list_by_index(.ode_equations, from_compartment, -net_flow)
        if ("total_deaths" %in% names(self$tracked_quantities)) {
          self$tracked_quantities$total_deaths <- self$tracked_quantities$total_deaths + net_flow
        }
      }
      .ode_equations
    },
    
    apply_universal_death_flow = function(.ode_equations, .compartment_values, .time) {
      #   apply the population-wide death rate to all compartments
      #
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      n_comp = 1
      for (compartment in names(self$compartment_values)) {
        net_flow <- unlist(unname(self$get_compartment_death_rate(compartment, .time))) * .compartment_values[[n_comp]]
        .ode_equations <- increment_list_by_index(.ode_equations, n_comp, -net_flow)
        
        # track deaths in case births are meant to replace deaths
        if ("total_deaths" %in% names(self$tracked_quantities)) {
          self$tracked_quantities$total_deaths <- self$tracked_quantities$total_deaths + net_flow
        }
        n_comp = n_comp + 1
      }
      .ode_equations
    },
    
    get_compartment_death_rate = function(.compartment, .time){
      # find the universal or population-wide death rate for a particular compartment
      
      # :param _compartment: str
      # name of the compartment
      # :param _time: float
      # current integration time
      # :return: float
      # death rate
      self$get_parameter_value(paste("universal_death_rateX", .compartment, sep='' ), .time)
    },
    
    apply_birth_rate = function(.ode_equations, .compartment_values) {
      #   apply a birth rate to the entry compartments
      #
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      .ode_equations <- increment_list_by_index(.ode_equations, match(self$entry_compartment, names(self$compartment_values)),
                                                self$find_total_births(.compartment_values))
    },
    
    find_total_births = function (.compartment_values) {
      #   work out the total births to apply dependent on the approach requested
      #
      #   :param _compartment_values:
      #       as for preceding methods
      #   :return: numeric scalar value
      #       total rate of births to be implemented in the model
      if (self$birth_approach == "add_crude_birth_rate") {
        return(self$parameters$crude_birth_rate * sum(.compartment_values))
      }
      else if (self$birth_approach == "replace_deaths") {
        return(self$tracked_quantities$total_deaths)
      }
      else {
        return(0)
      }
    },
    
    find_infectious_multiplier = function(n_flow) {
      #   find the multiplier to account for the infectious population in dynamic flows
      #
      #   :param flow_type: char
      #       type of flow, as per the standard naming approach to flow types for the dataframes flow attribute
      #   :return:
      #       the total infectious quantity, whether that be the number or proportion of infectious persons
      #       needs to return as one for flows that are not transmission dynamic infectiousness flows
      
      
      if (self$transition_flows$type[n_flow] == "infection_density") {
        self$infectious_population <- self$infectious_population
      }
      else if (self$transition_flows$type[n_flow] == "infection_frequency") {
        self$infectious_population <- self$infectious_population / self$infectious_denominators
      }
      else {
        return(1) 
      }
      
    },
    
    update_tracked_quantities = function(.compartment_values) {
      #   update quantities that emerge during model running (not pre-defined functions of time)
      #
      #   :param _compartment_values:
      #       as for preceding methods
      for (quantity in names(self$tracked_quantities)) {
        self$tracked_quantities[[quantity]] <- 0
        if (quantity == "infectious_population") {
          self$find_infectious_population(.compartment_values)
        }
        else if (quantity == "total_population") {
          self$tracked_quantities$total_population <- sum(.compartment_values)
        }
      }
    },
    
    find_infectious_population = function(compartment_values) {
      #   calculations to find the effective infectious population
      #
      #   :param _compartment_values:
      #       as for preceding methods
      for (compartment in which(names(self$compartment_values) == self$infectious_compartment)) {
        self$infectious_population <-
          self$infectious_population + compartment_values[compartment]
      }
      self$infectious_denominators <- sum(compartment_values)
    },
    
    
    get_single_parameter_component = function(.parameter, .time){
      # find the value of a parameter with time-variant values trumping constant ones
      
      # param _parameter: str
      # string for the name of the parameter of interest
      # :param time: float
      # model integration time (if needed)
      # :return: float
      # parameter value, whether constant or time variant
      if(.parameter %in% self$time_variants)
        return(self$time_variants$parameter[.time])
      else
        return(self$parameters[.parameter])
      
    },
    
    calculate_outputs_post_integration = function(){
      # find outputs for each requested time point, rather than at the time points that the model integration steps
      # occurred at, which are arbitrary and determined by the integration routine used
      
      self$derived_outputs_shadow["times"] = self$times
      for (output in self$output_connections){
        self$derived_outputs_shadow[output] = rep(0., length(self$times))
        transition_indices = self$find_output_transition_indices(output)
        n_time = 1
        for (time in self$times){
          self$restore_past_state(time)
          for (n_flow in transition_indices){
            net_flow = self$find_net_transition_flow(n_flow, time, self$compartment_values)
            self.derived_outputs_shadow[[output]][n_time] = self.derived_outputs_shadow[[output]][n_time] + net_flow
          }
          n_time = n_time + 1 
        }
        
      }
      
    },
    
    restore_past_state = function(time){
      # return compartment values and tracked quantities to the values current at a particular time during model
      # integration from the returned outputs structure
      # the model need not have been evaluated at this particular time point
      
      # :param time: float
      # time point to go back to
      
      self$compartment_values = self$outputs[self.times.index(time)]
      self$update_tracked_quantities(self$compartment_values)
    },
    
    find_output_transition_indices = function(output){
      # find the transition indices that are relevant to a particular output evaluation request from the
      # output_connections dictionary created from the user's request
      
      # :param output: str
      # name of the output of interest
      # :return: list
      # integers referencing the transition flows relevant to this output connection
      row_list = c()
      for (row in seq(length(self$transition_flows))){
        if ((self$transition_flows$implement[row] == length(self$all_stratifications)) &
            (find_stem(self$transition_flows$origin[row]) == self$output_connections[[output]]["origin"])&
            (find_stem(self$transition_flows$to[row]) == self$output_connections[[output]]["to"]) &
            (self$output_connections[[output]]["origin_condition"] %in% self$transition_flows$origin[row]) &
            (self$output_connections[[output]]["to_condition"] %in% self$transition_flows$to[row]))
          row_list = c(row_list,row)
      }
      row_list
    },
    
    get_parameter_value = function(.parameter, .time) {
      #   very simple, essentially place-holding, but need to split this out as a function in order to
      #   stratification later
      #
      #   :param .parameter: char
      #       parameter name
      #   :param .time: numeric scalar value
      #       current integration time
      #   :return: numeric scalar value
      #       parameter value
      self$get_single_parameter_component(.parameter, .time)
    }
  )
)


StratifiedModel <- R6Class(
  #   stratified version of the epidemiological model, inherits from EpiModel which is a concrete class and can run models
  #   independently (and could even include stratifications by using loops in a more traditional way to coding these
  #   models)
  #
  #   :attribute all_stratifications: character vector
  #       all the stratification names implemented so far
  #   :attribute removed_compartments: characeter vector
  #       all unstratified compartments that have been removed through the stratification process
  #   :attribute overwrite_parameters: character vector
  #       any parameters that are intended as absolute values to be applied to that stratum and not multipliers for the
  #       unstratified parameter further up the tree
  #   :attribute compartment_types_to_stratify
  #       see check_compartment_request
  #   :attribute heterogeneous_infectiousness
  #   :attribute infectiousness_adjustments
  #   :attribute parameter_components
  
  inherit = EpiModel,
  public = list(
    strata = c(),
    removed_compartments = c(),
    overwrite_parameters = c(),
    heterogeneous_infectiousness = FALSE,
    compartment_types_to_stratify = c(),
    infectious_populations = list(),
    infectious_denominators = NULL,
    strains = NULL,
    mixing_categories = c(),
    all_stratifications = NULL,
    infectiousness_adjustments = NULL,
    final_parameter_functions = list(),
    parameter_components = list(),
    mortality_components = list(),
    parameter_functions = list(),
    adaptation_functions = list(),
    mapped_adaptation_functions = list(),
    mixing_numerator_indices = c(),
    mixing_denominator_indices = c(),
    infectiousness_levels = NULL,
    infectious_indices = list(),
    infectious_compartments = NULL,
    infectiousness_multipliers = NULL,
    heterogeneous_mixing = FALSE,
    mixing_matrix = NULL,
    available_death_rates=c(""),
    strain_mixing_elements = data.frame(),
    strain_mixing_multipliers = data.frame(),
    overwrite_character = NULL,
    overwrite_key = c("W", "overwrite"),
    
    
    # __________
    # general methods
    
    add_compartment = function(new_compartment_name, new_compartment_value) {
      #   add a compartment by specifying its name and value to take
      #
      #   :param new_compartment_name: char
      #       name of the new compartment to be created
      #   :param new_compartment_value: numeric scalar value
      #       initial value to be assigned to the new compartment before integration
      self$compartment_names <- c(self$compartment_names, new_compartment_name)
      self$compartment_values[new_compartment_name] <- new_compartment_value
      self$output_to_user(paste("adding compartment:", new_compartment_name))
    },
    
    remove_compartment = function(compartment) {
      #   remove a compartment by taking the element out of the compartment_names and compartment_values attributes
      #   store name of removed compartment in removed_compartments attribute
      #
      #   :param compartment_name: char
      #       name of compartment to be removed
      self$removed_compartments <- c(self$removed_compartments, compartment)
      self$compartment_values <- self$compartment_values[names(self$compartment_values) != compartment]
      self$compartment_names <- self$compartment_names[self$compartment_names != compartment]
      self$output_to_user(paste("removing compartment:", compartment))
    },
    
    # __________
    # main master method for model stratification
    
    stratify = function(stratification_name, strata_request, compartment_types_to_stratify, requested_proportions=list(),
                        entry_proportions=list(), adjustment_requests=list(),  infectiousness_adjustments=c(), mixing_matrix=NULL, verbose=TRUE) {
      #   calls to initial preparation, checks and methods that stratify the various aspects of the model
      #
      #   :param stratification_name:
      #       see prepare_and_check_stratification
      #   :param strata_request:
      #       see find_strata_names_from_input
      #   :param compartment_types_to_stratify:
      #       see check_compartment_request
      #   :param adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      #   :param requested_proportions:
      #       see prepare_starting_proportions
      #   :param infectiousness_adjustments:
      #
      #   :param verbose: logical
      #       whether to report on progress, note that this can be changed at this stage from what was requested at
      #       the original unstratified model construction
      
      # check inputs correctly specified
      check_results <- self$prepare_and_check_stratification(
        stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, verbose)
      strata_names <- check_results[[1]]
      adjustment_requests <- check_results[[2]]
      
      # work out ageing flows (comes first so that the compartment names are still in the unstratified form)
      if (stratification_name == "age") {
        self$set_ageing_rates(strata_names)
      }
      
      # stratify the compartments
      requested_proportions <- self$prepare_starting_proportions(strata_names, requested_proportions)
      self$stratify_compartments(stratification_name, strata_names, requested_proportions)
      
      # stratify the flows
      self$stratify_transition_flows(stratification_name, strata_names, adjustment_requests)
      self$stratify_entry_flows(stratification_name, strata_names, requested_proportions)
      if (nrow(self$death_flows) > 0) {
        self$stratify_death_flows(stratification_name, strata_names, adjustment_requests)
      }
      self$stratify_universal_death_rate(stratification_name, strata_names, adjustment_requests)
      
      # check submitted mixing matrix and combine with existing matrix, if any
      self$prepare_mixing_matrix(mixing_matrix, stratification_name, strata_names)
      
      # implement heterogeneous mixing across multiple population groups
      self$prepare_implement_mixing()
      
      # if a multi-strain model
      if (stratification_name == "strain"){
        self$strains = strata_names
      }
      
      # heterogeneous infectiousness adjustments
      #self$apply_heterogeneous_infectiousness(stratification_name, strata_request, infectiousness_adjustments)
      self$prepare_infectiousness_levels(stratification_name, strata_names, infectiousness_adjustments)
    },
    
    # __________
    # pre-integration methods
    
    prepare_and_check_stratification = function(.stratification_name, .strata_request, .compartment_types_to_stratify,
                                                .adjustment_requests, .verbose) {
      #   initial preparation and checks
      #
      #   :param .stratification_name: char
      #       the name of the stratification - i.e. the reason for implementing this type of stratification
      #   :param .strata_request:
      #       see find_strata_names_from_input
      #   :param .compartment_types_to_stratify:
      #       see check_compartment_request
      #   :param .adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      #   :param .verbose:
      #       see stratify
      #   :return:
      #       .strata_names: list
      #           revised version of _strata_request after adaptation to class requirements
      #       .adjustment_requests:
      #           revised version of _adjustment_requests after adaptation to class requirements
      self$verbose <- .verbose
      
      # checks and reporting for age stratification and general starting message otherwise
      if (.stratification_name == "age") {
        .strata_request <- self$check_age_stratification(.strata_request, .compartment_types_to_stratify)
      }
      else {
        self$output_to_user(paste("\nimplementing stratification for:", .stratification_name))
      }
      
      # make sure all stratification names are characters
      if (!is.character(.stratification_name)) {
        .stratification_name <- as.character(.stratification_name)
        self$output_to_user(paste("converting stratification name", .stratification_name, "to character"))
      }
      
      # record stratification as attribute to model, find the names to apply strata and check compartment and parameter requests
      self$strata <- c(self$strata, .stratification_name)
      .strata_names <- self$find_strata_names_from_input(.strata_request)
      .adjustment_requests <- self$incorporate_alternative_overwrite_approach(.adjustment_requests)
      self$check_compartment_request(.compartment_types_to_stratify)
      .adjustment_requests <- self$check_parameter_adjustment_requests(.adjustment_requests, .strata_names)
      list(.strata_names, .adjustment_requests)
    },
    
    check_age_stratification = function(.strata_request, .compartment_types_to_stratify) {
      #   check that request meets the requirements for stratification by age
      #
      #   :parameters: all parameters have come directly from the stratification (stratify) method unchanged and have been
      #       renamed with a preceding . character
      #   :return: .strata_request: character vector
      #       revised names of the strata tiers to be implemented
      self$output_to_user(paste("implementing age-specific stratification with specific behaviour"))
      if (length(.compartment_types_to_stratify) > 0) {
        stop("requested age stratification, but compartment request should be passed as empty vector in order to apply to all compartments")
      }
      else if (!is.numeric(.strata_request)) {
        stop("inputs for age strata breakpoints are not numeric")
      }
      else if ("age" %in% self$strata) {
        stop("requested stratification by age, but this has specific behaviour and can only be applied once")
      }
      if (is.unsorted(.strata_request)) {
        .strata_request <- sort(.strata_request)
        self$output_to_user(paste("requested strata for age stratification not ordered, so have been sorted to:",
                                  paste(rep(", ", length(.strata_request)), .strata_request, collapse=""), sep=""))
      }
      if (!0 %in% .strata_request) {
        self$output_to_user(paste("adding age stratum called '0' as not requested, to represent those aged less than", .strata_request[1]))
        .strata_request <- c(0, .strata_request)
      }
      .strata_request
    },
    
    find_strata_names_from_input = function(.strata_request) {
      #   find the names of the strata to be implemented from a particular user request
      #
      #   :.strata_request: list or alternative format to be adapted
      #       strata requested in the format provided by the user (except for age, which is dealth with in the preceding
      #       method)
      #   :return: strata_names: character vector
      #       modified list of strata to be implemented in model
      if (length(.strata_request) == 1 & is.numeric(.strata_request)) {
        if (.strata_request %% 1 == 0 & .strata_request > 1) {
          strata_names <- seq(.strata_request)
          self$output_to_user(paste("integer provided strata labels for stratification, hence strata implemented are integers from 1 to", .strata_request))
        }
        else {
          stop("number passed as request for strata labels, but not an integer greater than one, unclear what to do, stratification failed")
        }
      }
      else if (length(.strata_request) > 0) {
        strata_names <- .strata_request
      }
      else {
        stop("requested to stratify, but no strata provided")
      }
      for (name in strata_names) {
        self$output_to_user(paste("adding stratum:", name))
      }
      as.character(strata_names)
    },
    
    check_compartment_request = function(.compartment_types_to_stratify) {
      #   check the requested compartments to be stratified has been requested correctly
      #
      #   :param .compartment_types_to_stratify: character vector
      #       the names of the compartment types that the requested stratification is intended to apply to
      
      # if vector of length zero passed, stratify all the compartment types in the model
      if (length(.compartment_types_to_stratify) == 0) {
        self$output_to_user("no compartment names specified for this stratification, so stratification applied to all model compartments")
        self$compartment_types_to_stratify <- self$compartment_types
      }
      
      # otherwise check all the requested compartments are available and implement the user request
      else if (length(setdiff(.compartment_types_to_stratify, self$compartment_types)) != 0) {
        stop("requested compartment or compartments to be stratified are not available in this model")
      }
      else {
        self$compartment_types_to_stratify <- .compartment_types_to_stratify
      }
    },
    
    incorporate_alternative_overwrite_approach = function(.adjustment_requests) {
      #   alternative approach to working out which parameters to overwrite
      #   can now put a capital W at the string's end to indicate that it is an overwrite parameter, as an alternative to
      #   submitting a separate dictionary key to represent the strata which need to be overwritten
      #
      #   :param .adjustment_requests: list
      #       user-submitted version of adjustment requests
      #   :return: revised_adjustments: list
      #       modified version of _adjustment_requests after working out whether any parameters began with W
      
      revised_adjustments <- list()
      for (parameter in names(.adjustment_requests)) {
        
        # accept the key representing the overwrite parameters
        revised_adjustments[[parameter]] <- list()
        if ("overwrite" %in% .adjustment_requests[[parameter]]) {
          revised_adjustments[[parameter]]$overwrite <- .adjustment_requests[[parameter]]$overwrite
        }
        else {
          revised_adjustments[[parameter]]$overwrite <- c()
        }
        
        # then loop through all the other keys of the user request
        for (stratum in names(.adjustment_requests[[parameter]])) {
          if (stratum == "overwrite") {
            next
          }
          
          # if the parameter ends in W, it is interpreted as an overwrite parameter and added to this key
          else if (substr(stratum, nchar(stratum), nchar(stratum)) == "W") {
            revised_adjustments[[parameter]][substr(stratum, 1, nchar(stratum) - 1)] <- .adjustment_requests[[parameter]][[stratum]]
            revised_adjustments[[parameter]]$overwrite <- c(revised_adjustments[[parameter]]$overwrite, substr(stratum, 1, nchar(stratum) - 1))
          }
          
          # otherwise just accept the parameter in its submitted form
          else {
            revised_adjustments[[parameter]][[stratum]] <- .adjustment_requests[[parameter]][[stratum]]
          }
        }
        .adjustment_requests[[parameter]] <- revised_adjustments[[parameter]]
      }
      revised_adjustments
    },
    
    check_parameter_adjustment_requests = function(.adjustment_requests, .strata_names) {
      #   check parameter adjustments have been requested appropriately and add parameter for any strata not referred to
      #
      #   :param .adjustment_requests: list
      #       version of the submitted adjustment_requests modified by incorporate_alternative_overwrite_approach
      #   :param .strata_names:
      #       see find_strata_names_from_input
      #   :return: .adjustment_requests
      #       modified version of _adjustment_requests after checking
      for (parameter in names(.adjustment_requests)) {
        
        # check all the requested strata for parameter adjustments were strata that were requested
        for (requested_stratum in names(.adjustment_requests[[parameter]])) {
          if (!requested_stratum %in% as.character(.strata_names) & requested_stratum != "overwrite") {
            stop(paste("stratum", requested_stratum, "requested in adjustments but unavailable"))
          }
        }
        
        # if any strata were not requested, assume a value of one
        for (stratum in as.character(.strata_names)) {
          if (!stratum %in% names(.adjustment_requests[[parameter]])) {
            .adjustment_requests[[parameter]][[stratum]] <- 1
            self$output_to_user(paste("no request made for adjustment to", parameter, "within stratum", stratum, "so using parent value by default"))
          }
        }
      }
      .adjustment_requests
    },
    
    prepare_starting_proportions = function(.strata_names, .requested_proportions) {
      #   prepare user inputs for starting proportions as needed
      #   must be specified with names that are strata being implemented during this stratification process
      #   note this applies to initial conditions and to entry flows
      #
      #   :param .strata_names:
      #       see find_strata_names_from_input
      #   :param .requested_proportions: list
      #       list with names for the stratum to assign starting population to and values the proportions to assign
      #   :return: list
      #       revised list of starting proportions after cleaning
      
      # assume an equal proportion of the total for the compartment if not otherwise specified
      for (stratum in .strata_names) {
        if (!stratum %in% names(.requested_proportions)) {
          starting_proportion <- 1 / length(.strata_names)
          .requested_proportions[as.character(stratum)] <- starting_proportion
          self$output_to_user(paste("no starting proportion requested for stratum", stratum,
                                    "so allocated", round(starting_proportion, self$reporting_sigfigs), "of total"))
        }
      }
      
      # normalise the dictionary before return, in case adding the missing groups as equal proportions exceeds one
      normalise_list(.requested_proportions)
    },
    
    stratify_compartments = function(.stratification_name, .strata_names, .requested_proportions) {
      #   stratify the model compartments, which affects the compartment_names and the compartment_values attributes
      #
      #   :param .stratification_name:
      #       see prepare_and_check_stratification
      #   :param .strata_names:
      #       see find_strata_names_from_input
      #   :param .requested_proportions:
      #       see prepare_starting_proportions
      
      # find the existing compartments that need stratification
      for (compartment in names(self$compartment_values)) {
        if (find_stem(compartment) %in% self$compartment_types_to_stratify) {
          
          # add and remove compartments
          for (stratum in .strata_names) {
            self$add_compartment(create_stratified_name(compartment, .stratification_name, stratum),
                                 self$compartment_values[[compartment]] * as.numeric(.requested_proportions[as.character(stratum)]))
          }
          self$remove_compartment(compartment)
        }
      }
    },
    
    stratify_transition_flows = function(.stratification_name, .strata_names, .adjustment_requests) {
      #   stratify flows depending on whether inflow, outflow or both need replication, using call to add_stratified_flows
      #   method below
      #
      #   :param .stratification_name:
      #       see prepare_and_check_stratification
      #   :param .strata_names:
      #       see find_strata_names_from_input
      #   :param .adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      for (n_flow in which((self$transition_flows$implement == length(self$strata) - 1))) {
        self$add_stratified_flows(n_flow, .stratification_name, .strata_names,
                                  find_stem(self$transition_flows$from[n_flow]) %in% self$compartment_types_to_stratify,
                                  find_stem(self$transition_flows$to[n_flow]) %in% self$compartment_types_to_stratify,
                                  .adjustment_requests)
      }
      self$output_to_user("stratified transition flows matrix:")
      self$output_to_user(self$transition_flows)
    },
    
    stratify_entry_flows = function(.stratification_name, .strata_names, .requested_proportions) {
      #   stratify entry/recruitment/birth flows according to requested entry proportion adjustments
      #   note this applies to initial conditions and to entry flows
      #
      #   :param .stratification_name:
      #       see prepare_and_check_stratification
      #   :param .strata_names:
      #       see find_strata_names_from_input
      #   :param .requested_proportions:
      #       see prepare_starting_proportions
      #   :return:
      #       normalised list of the compartments that the new entry flows should come in to
      entry_fractions <- list()
      if (self$entry_compartment %in% self$compartment_types_to_stratify) {
        for (stratum in .strata_names) {
          entry_fraction_name <- create_stratified_name("entry_fraction", .stratification_name, stratum)
          
          # specific behaviour for age stratification
          if (.stratification_name == "age" & as.character(stratum) == "0") {
            entry_fractions[entry_fraction_name] <- 1
            next
          }
          else if (.stratification_name == "age") {
            entry_fractions[entry_fraction_name] <- 0
            next
          }
          
          # where a request has been submitted
          else if (stratum %in% names(.requested_proportions)) {
            entry_fractions[entry_fraction_name] <- .requested_proportions[[stratum]]
            self$output_to_user(paste("assigning specified proportion of starting population to", stratum))
          }
          
          # otherwise if no request made
          else {
            entry_fractions[entry_fraction_name] <- 1 / length(.strata_names)
            self$output_to_user(paste("assuming", as.character(entry_fractions[entry_fraction_name]), "of starting population to be assigned to", stratum, "stratum by default"))
          }
        }
        self$parameters <- c(normalise_list(entry_fractions), self$parameters)
      }
    },
    
    stratify_death_flows = function(.stratification_name, .strata_names, .adjustment_requests) {
      #   add compartment-specific death flows to death_flows data frame attribute
      #
      #   :param .stratification_name:
      #       see prepare_and_check_stratification
      #   :param .strata_names:
      #        see find_strata_names_from_input
      #   :param .adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      for (n_flow in which(self$death_flows$implement == length(self$strata) - 1)) {
        if (find_stem(self$death_flows$from[n_flow]) %in% self$compartment_types_to_stratify) {
          for (stratum in .strata_names) {
            
            # get stratified parameter name if requested to stratify, otherwise use the unstratified one
            parameter_name <- self$add_adjusted_parameter(self$death_flows$parameter[n_flow], .stratification_name, stratum, .adjustment_requests)
            if (is.null(parameter_name)) {
              parameter_name <- self$death_flows$parameter[n_flow]
            }
            
            if('force_index' %in% colnames(self$transition_flows)){
              print('check')
              
            }
            
            
            # add the stratified flow
            self$death_flows <- rbind(self$death_flows,
                                      data.frame(type=self$death_flows$type[n_flow],
                                                 parameter=parameter_name,
                                                 from=create_stratified_name(self$death_flows$from[n_flow], .stratification_name, stratum),
                                                 implement=length(self$strata), stringsAsFactors=FALSE))
          }
        }
      }
    },
    
    stratify_universal_death_rate = function(.stratification_name, .strata_names, .adjustment_requests) {
      #   stratify the approach to universal, population-wide deaths (which can be made to vary by stratum)
      #   adjust each parameter that refers to the universal death rate according to user request
      #
      #   :param .stratification_name:
      #       see prepare_and_check_stratification
      #   :param .strata_names:
      #       see find_strata_names_from_input
      #   :param .adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      
      if(!(.stratification_name %in% self$full_stratification_list) & ("universal_death_rate" %in% .adjustment_requests) )
        stop("universal death rate can only be stratified when applied to all compartment types")
      else if (!(.stratification_name %in% self$full_stratification_list))
        print("universal death rate not adjusted as stratification not applied to all compartments")
      else
        print('check')
      
      if ("universal_death_rate" %in% self$time_variants)
        self$adaptation_functions["universal_death_rateX"] = self$time_variants["universal_death_rate"]
      else
        self$adaptation_functions[["universal_death_rateX"]] = (function(time){self$parameters["universal_death_rate"]})
      
      for (parameter in names(self$parameters)) {
        if (find_stem(parameter) == "universal_death_rate") {
          for (stratum in .strata_names) {
            self$add_adjusted_parameter(parameter, .stratification_name, stratum, .adjustment_requests)
          }
        }
      }
    },
    
    add_adjusted_parameter = function(.unadjusted_parameter, .stratification_name, .stratum, .adjustment_requests) {
      #   find the adjustment request that is relevant to a particular unadjusted parameter and stratum
      #   otherwise allow return of None
      #
      #   :param .unadjusted_parameter: char
      #       name of the unadjusted parameter value
      #   :param .stratification_name:
      #       see prepare_and_check_stratification
      #   :param .stratum:
      #       stratum being considered by the method calling this method
      #   :param .adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      #   :return: parameter_adjustment_name: str or None
      #       if returned as None, assumption will be that the original, unstratified parameter should be used
      #       otherwise create a new parameter name and value and store away in the appropriate model structure
      self$output_to_user(paste("\tmodifying", .unadjusted_parameter, "for", .stratum, "stratum of", .stratification_name))
      
      # find the adjustment requests that are extensions of the base parameter type being considered
      for (parameter_request in names(.adjustment_requests)) {
        if (startsWith(.unadjusted_parameter, parameter_request)) {
          parameter_adjustment_name <- create_stratified_name(.unadjusted_parameter, .stratification_name, .stratum)
          self$output_to_user(paste("modifying", .unadjusted_parameter, "for", .stratum, "of", .stratification_name,
                                    "with new parameter called", parameter_adjustment_name))
          
          # implement user request if requested (note that otherwise parameter will now be left out and assumed to be 1 during integration)
          if (.stratum %in% names(.adjustment_requests[[parameter_request]])) {
            self$parameters[parameter_adjustment_name] <- .adjustment_requests[[parameter_request]][[as.character(.stratum)]]
          }
          
          # overwrite parameters higher up the tree by listing which ones to be overwritten
          if (.stratum %in% .adjustment_requests[[parameter_request]]) {
            self$overwrite_parameters <- c(self$overwrite_parameters, parameter_adjustment_name)
          }
          return(parameter_adjustment_name)
        }
      }
    },
    
    prepare_mixing_matrix = function(.mixing_matrix, .stratification_name, .strata_names){
      
      #    check that the mixing matrix has been correctly specified and call the other relevant functions
      #  :param _mixing_matrix: ndarray
      #  array, which must be square representing the mixing of the strata within this stratification
      #  :param _stratification_name: str
      #  the name of the stratification - i.e. the reason for implementing this type of stratification
      #  :param _strata_names: list
      #  see find_strata_names_from_input
      if (is.null(.mixing_matrix))
        return
      else if (class(.mixing_matrix) != 'matrix')
        stop("submitted mixing matrix is wrong data type")
      else if (length(dim(.mixing_matrix)) != 2)  
        stop("submitted mixing matrix is not two-dimensional")
      else if (dim(.mixing_matrix)[1] != dim(.mixing_matrix)[2] )
        stop("submitted mixing is not square")
      else if (dim(.mixing_matrix)[1] != length(.strata_names)) 
        stop("mixing matrix does not sized to number of strata being implemented")
      self$combine_new_mixing_matrix_with_existing(.mixing_matrix, .stratification_name, .strata_names)
      
    },
    
    combine_new_mixing_matrix_with_existing = function(.mixing_matrix, .stratification_name, .strata_names){
      # master mixing matrix function to take in a new mixing matrix and combine with the existing ones
      # :param _mixing_matrix: ndarray
      # array, which must be square representing the mixing of the strata within this stratification
      # :param _stratification_name: str
      # the name of the stratification - i.e. the reason for implementing this type of stratification
      # :param _strata_names: list
      # see find_strata_names_from_input
      
      # if no mixing matrix yet, just convert the existing one to a dataframe
      if(is.null(self$mixing_matrix)){
        self$mixing_categories <- paste(.stratification_name, '_', .strata_names, sep='' )
        self$mixing_matrix <- .mixing_matrix
        
      }
      # otherwise take the kronecker product to get the new mixing matrix
      else {
        new_mixing_cols <- NULL
        for (each_col in self$mixing_categories){
          new_col_names <- paste(each_col,'X', .stratification_name, '_', .strata_names, sep="") 
          new_mixing_cols <- c(new_mixing_cols, new_col_names )
        }
        if(!is.null(.mixing_matrix)){
          self$mixing_categories <- new_mixing_cols
          self$mixing_matrix <- kronecker(as.matrix(self$mixing_matrix), .mixing_matrix)}
      }
    },
    
    find_transition_indices_to_implement = function(back_one=0){
      # find all the indices of the transition flows that need to be stratified
      # separated out as very short method in order that it can over-ride the version in the unstratified EpiModel
      
      # :param back_one: int
      # number to subtract from self.all_stratification, which will be one if this method is being called after the
      # stratification has been added
      # :return: list
      # list of indices of the flows that need to be stratified
      
      return(c(which(self$transition_flows$implement == length(self$all_stratifications)-back_one)))
    },
    
    find_death_indices_to_implement = function(back_one=0){
      # find all the indices of the death flows that need to be stratified
      # separated out as very short method in order that it can over-ride the version in the unstratified EpiModel
      
      # :param back_one: int
      # number to subtract from self.all_stratification, which will be one if this method is being called after the
      # stratification has been added
      # :return: list
      # list of indices of the flows that need to be stratified
      
      return(c(which(self$death_flows$implement == length(self$all_stratifications)-back_one)))
    },
    
    prepare_to_run = function(){
      # methods that can be run prior to integration to save various function calls being made at every time step
      self$prepare_stratified_parameter_calculations()
      self$prepare_infectiousness_calculations()
      #self$transition_indices_to_implement = self$find_transition_indices_to_implement1()
      #self$death_indices_to_implement = self$find_death_indices_to_implement1()
    },
    
    prepare_implement_mixing = function(){
      # methods to be run if there is a mixing matrix being applied at all, regardless of whether one is being
      # introduced during this stratification process
      if (!is.null(self$mixing_matrix)){
        self$find_mixing_indices()
        self$add_force_indices_to_transitions()
      }
    },
    
    find_mixing_indices = function(){
      # find the indices for the infectious compartments relevant to the column of the mixing matrix
      
      self$mixing_numerator_indices = list()
      self$mixing_denominator_indices = list()
      for (from_stratum in self$mixing_categories){
        #self$mixing_numerator_indices[from_stratum] 
        # create a list of list
        n_comp = 0
        for (compartment in self$compartment_names){
          if(identical(tail(unlist(find_name_components(compartment)),-1),unlist(find_name_components(from_stratum)))){
            self$mixing_denominator_indices[[from_stratum]] <- c(self$mixing_denominator_indices[[from_stratum]], n_comp)
            if(grepl(self$infectious_compartment,compartment)){
              self$mixing_numerator_indices[[from_stratum]] <- c(self$mixing_numerator_indices[[from_stratum]], n_comp)
            }
          }
          n_comp = n_comp + 1  
        }
        
      }
      
    },
    
    find_mixing_denominators = function(){
      # for each mixing category, create a list of the compartment numbers that are relevant
      
      # :return mixing_indices: list
      # indices of the compartments that are applicable to a particular mixing category
      
      if (is.null(self$mixing_matrix))
        return
      else{ 
        mixing_indices = list()
        for (category in self$mixing_categories){
          i_comp = 1
          for (compartment in self$compartment_names){
            if (find_name_components(category) %in% find_name_components(compartment))
              mixing_indices[[category]] <- c(mixing_indices[[category]], i_comp)
            i_comp = i_comp + 1}
        } 
        return(mixing_indices)
      }
    },
    
    find_strain_mixing_multipliers = function(mixing_indices){
      # find the relevant indices to be used to calculate the force of infection contribution to each strain from each
      # mixing category as a list of indices - and separately find multipliers as a list of the same length for
      # their relative infectiousness extracted from self.infectiousness_multipliers
      for(strain in c(self$strains, "all_strains")){
        mixing_categories_list = c()
        if (is.null(self$mixing_matrix))
          mixing_categories_list = c("all_population")
        else
          mixing_categories_list = self$mixing_categories
        
        for (category in mixing_categories_list){
          for (idx in mixing_indices[[category]]){
            if (self$infectious_indices[[strain]][idx]){
              self$strain_mixing_elements = rbind(self$strain_mixing_elements, data.frame(strain=strain, category=category, idx=idx))}
          }
          i_comp_list = self$strain_mixing_elements[ (self$strain_mixing_elements$strain==strain & self$strain_mixing_elements$category==category), "idx"]
          for(i_comp in i_comp_list ){
            self$strain_mixing_multipliers = rbind(self$strain_mixing_multipliers, data.frame( strain=strain, category=category, idx=self$infectiousness_multipliers[i_comp]))}
          
        }
      }
      
    },
    
    
    
    add_force_indices_to_transitions = function(){
      # find the indices from the force of infection vector to be applied for each infection flow and populate to the
      # force_index column of the flows frame 
      infection_flow_indices = c()
      # identify the indices of the infection-related flows to be implemented
      
      for (row in 1:nrow(self$transition_flows)){
        if(grepl('infection',self$transition_flows[row, 'type']) & (self$transition_flows[row,'implement'] != 0))
          infection_flow_indices <- c(infection_flow_indices, row)
      }
      
      self$transition_flows$force_index <- 1
      # loop through them and find the indices of the mixing matrix that will apply to that flow
      for (n_flow in infection_flow_indices){
        found = FALSE
        for (force_group in self$mixing_categories){
          #if(identical(tail(unlist(find_name_components(self$transition_flows$from[n_flow])), -1),
          #   find_name_components(force_group)))
          if (find_name_components(force_group) %in% find_name_components(self$transition_flows$from[n_flow]))
            self$transition_flows[n_flow, ]$force_index <- which(self$mixing_categories==force_group) + 1 #similar to ngroup
          if (found){
            print('mixing group found twice for transition flow number')
            
          }
          found = TRUE
          next
        }
      }
      
    },
    
    prepare_all_infectiousness_multipliers = function(){
      # find the infectiousness multipliers for each compartment being implemented in the model
      # start from assumption that each compartment is fully and equally infectious
      self$infectiousness_multipliers = c(rep(1., length(self$compartment_names)))
      
      # if infectiousness modification requested for the compartment type, multiply through by the current value
      n_comp = 1
      for (compartment in self$compartment_names){
        for (modifier in names(self$infectiousness_levels)){
          if (modifier %in% find_name_components(compartment))
            self$infectiousness_multipliers[n_comp] = self$infectiousness_multipliers[n_comp] * self$infectiousness_levels[modifier]
        }
        n_comp = n_comp + 1
      }
      
    },
    
    find_infectious_indices = function(){
      self$infectious_indices['all_strains'] <- list(as.logical(grepl(self$infectious_compartment , self$compartment_names)))
      if (!is.null(self$strains)){
        for (strain in self$strains){
          self$infectious_indices[strain]  <- list(as.logical(grepl(paste(self$infectious_compartment, '.*', strain,  sep="") , self$compartment_names)))
        }
      }
    },
    
    prepare_infectiousness_levels = function(.stratification_name, .strata_names, .infectiousness_adjustments){
      # store infectiousness adjustments as dictionary attribute to the model object, with first tier of keys the
      # stratification and second tier the strata to be modified
      
      # :param _stratification_name:
      #  see prepare_and_check_stratification
      # :param _strata_names:
      #  see find_strata_names_from_input
      # :param _infectiousness_adjustments: dict
      # requested adjustments to infectiousness for this stratification
      
      if (length(.infectiousness_adjustments) == 0) {
        self$output_to_user("heterogeneous infectiousness not requested for this stratification")
      }
      else if (!self$infectious_compartment %in% self$compartment_types_to_stratify) {
        stop("request for infectiousness adjustments passed, but stratification doesn't apply to the infectious compartment")
      }
      else {
        for (stratum in names(.infectiousness_adjustments)) {
          self$infectiousness_levels[[create_stratum_name( .stratification_name, stratum, joining_string="")]] <-
            .infectiousness_adjustments[[stratum]]
        }
      }
      self$find_infectious_indices()
      for(strain in names(self$infectious_indices)){
        self$infectious_compartments[[strain]] <- self$compartment_names[self$infectious_indices[[strain]]]
        self$infectiousness_multipliers[[strain]] <- rep(1,length(self$infectious_compartments[[strain]]))
        for(compartment_idx in seq_along(self$infectious_compartments[[strain]])){
          for (infectiousness_modifier in names(self$infectiousness_levels)){
            if (infectiousness_modifier %in% find_name_components(self$infectious_compartments[[strain]][compartment_idx])){
              self$infectiousness_multipliers[[strain]][compartment_idx] <- as.numeric(self$infectiousness_multipliers[[strain]][compartment_idx]) * as.numeric(self$infectiousness_levels[[infectiousness_modifier]])
            }
            
          }
        }
      }
      
    },
    
    prepare_infectiousness_calculations = function(){
      # aster method to run all the code concerned with preparation for force of infection calculations
      # infectiousness preparations
      self$prepare_all_infectiousness_multipliers()
      self$find_infectious_indices()
      
      # mixing preparations
      if (!is.null(self$mixing_matrix))
        self$add_force_indices_to_transitions()
      mixing_indices = self$find_mixing_denominators()
      
      # reconciling the strains and the mixing attributes together into one structure
      self$find_strain_mixing_multipliers(mixing_indices)
      
    },
    
    
    apply_heterogeneous_infectiousness = function(stratification_name, strata_request, infectiousness_adjustments) {
      #   work out infectiousness adjustments and set as model attributes
      #   this has not been fully documented, as we are intending to revise this to permit any approach to heterogeneous
      #   infectiousness or mixing assumptions
      if (length(infectiousness_adjustments) == 0) {
        self$output_to_user("heterogeneous infectiousness not requested for this stratification")
      }
      else if (!self$infectious_compartment %in% self$compartment_types_to_stratify) {
        stop("request for infectiousness adjustments passed, but stratification doesn't apply to the infectious compartment")
      }
      else {
        self$heterogeneous_infectiousness <- TRUE
        for (stratum in names(infectiousness_adjustments)) {
          if (!stratum %in% strata_request) {
            stop("stratum to have infectiousness modified not found within requested strata")
          }
          self$infectiousness_adjustments[[create_stratified_name("", stratification_name, stratum)]] <-
            infectiousness_adjustments[[stratum]]
        }
      }
    },
    
    set_ageing_rates = function(.strata_names) {
      #   set intercompartmental flows for ageing from one stratum to the next as the reciprocal of the width of the age
      #   bracket
      #
      #   :param .strata_names:
      #       see find_strata_names_from_input
      for (stratum_number in seq(length(.strata_names) - 1)) {
        start_age <- as.numeric(.strata_names[stratum_number])
        end_age <- as.numeric(.strata_names[stratum_number + 1])
        ageing_parameter_name <- paste("ageing", as.character(start_age), "to", as.character(end_age), sep="")
        ageing_rate <- 1 / (end_age - start_age)
        self$output_to_user(paste("ageing rate from age group", start_age, "to", end_age, "is", round(ageing_rate, self$reporting_sigfigs)))
        self$parameters[ageing_parameter_name] <- ageing_rate
        for (compartment in names(self$compartment_values)) {
          if('force_index' %in% colnames(self$transition_flows) | 'strain' %in% colnames(self$transition_flows))
            self$transition_flows <-
              rbind(self$transition_flows,
                    data.frame(type="standard_flows", parameter=ageing_parameter_name,
                               from=create_stratified_name(compartment, "age", start_age),
                               to=create_stratified_name(compartment, "age", end_age),
                               implement=length(self$strata), strain=NA, force_index=NA, stringsAsFactors=FALSE))
          else
            self$transition_flows <-
              rbind(self$transition_flows,
                    data.frame(type="standard_flows", parameter=ageing_parameter_name,
                               from=create_stratified_name(compartment, "age", start_age),
                               to=create_stratified_name(compartment, "age", end_age),
                               implement=length(self$strata), stringsAsFactors=FALSE))
        }
      }
    },
    
    add_stratified_flows = function(.n_flow, stratification_name, strata_names, stratify_from, stratify_to, adjustment_requests) {
      #   add additional stratified flow to the transition flow data frame attribute of the class
      #
      #   :param .n_flow: integer
      #       location of the unstratified flow within the transition flow attribute
      #   :param stratification_name:
      #       see prepare_and_check_stratification
      #   :param strata_names: character vector
      #       see find_strata_names_from_input
      #   :param stratify_from: logical
      #       whether to stratify the from/origin compartment
      #   :param stratify_to: logical
      #       whether to stratify the to/destination compartment
      #   :param _adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      parameter_name <- NULL
      if (stratify_from | stratify_to) {
        self$output_to_user(paste("for flow from", self$transition_flows$from[.n_flow], "to", self$transition_flows$to[.n_flow], "in stratification", stratification_name))
        
        # loop over each stratum in the requested stratification structure
        for (stratum in strata_names) {
          
          # find parameter name
          parameter_name <- self$add_adjusted_parameter(
            self$transition_flows$parameter[.n_flow], stratification_name, stratum, adjustment_requests)
          if (is.null(parameter_name)) {
            parameter_name <- self$sort_absent_parameter_request(stratification_name, strata_names, stratum, stratify_from, stratify_to, self$transition_flows$parameter[.n_flow])
          }
          self$output_to_user("parameter_name")
          self$output_to_user(parameter_name)
          # determine whether to and/or from compartments are stratified
          if (stratify_from) {
            from_compartment <- create_stratified_name(self$transition_flows$from[.n_flow], stratification_name, stratum)
          }
          else {
            from_compartment <- self$transition_flows$from[.n_flow]
          }
          if (stratify_to) {
            to_compartment <- create_stratified_name(self$transition_flows$to[.n_flow], stratification_name, stratum)
          }
          else {
            to_compartment <- self$transition_flows$to[.n_flow]
          }
          
          strain = 'None'
          if(stratification_name == "strain")
            strain = stratum
          else
            strain = self$transition_flows$strain[.n_flow]
          
          
          
          # add the new flow
          if('force_index' %in% colnames(self$transition_flows)){
            #self$transition_flows$strain = strain
            if(!('strain' %in% colnames(self$transition_flows)))
              self$transition_flows$strain = 'None'
            self$transition_flows <- rbind(self$transition_flows,data.frame(
              parameter=parameter_name, from=from_compartment, to=to_compartment, implement=length(self$strata), strain=strain, force_index=NA, type=self$transition_flows$type[.n_flow]))
          }
          else
            self$transition_flows <- rbind(self$transition_flows,data.frame(
              parameter=parameter_name, from=from_compartment, to=to_compartment, implement=length(self$strata), type=self$transition_flows$type[.n_flow]))  
          
        }
        
        # remove old flow
        self$transition_flows$implement[.n_flow] <- length(self$strata) - 1
      }
      else {
        retained_flow <- self$transition_flows[.n_flow,]
        retained_flow$implement <- retained_flow$implement + 1
        self$transition_flows <- rbind(self$transition_flows, retained_flow, stringsAsFactors=FALSE)
      }
    },
    
    sort_absent_parameter_request = function (.stratification_name, .strata_names, .stratum, .stratify_from, .stratify_to, unstratified_name) {
      #   work out what to do if a specific parameter adjustment has not been requested
      #
      #   :param .stratification_name:
      #       see prepare_and_check_stratification
      #   :param .strata_names:
      #       see find_strata_names_from_input
      #   :param .stratum:
      #   :param .stratify_from:
      #       seeadd_stratified_flows
      #   :param .stratify_to:
      #       see add_stratified_flows
      #   :param .n_flow: int
      #       index of the flow being dealt with
      #   :return: str
      #       parameter name for revised parameter than wasn't provided
      
      # default behaviour for parameters not requested is to split the parameter into equal parts from compartment not split but to compartment is
      if (!(.stratify_from) & .stratify_to) {
        self$output_to_user(paste("\tsplitting existing parameter value", unstratified_name, "into", length(.strata_names), "equal parts"))
        parameter_name <- create_stratified_name(unstratified_name, .stratification_name, .stratum)
        self$parameters[parameter_name] <- 1 / length(.strata_names)
        self$adaptation_functions[[parameter_name]] <- create_multiplicative_function(1.0/length(.strata_names))
        return(parameter_name)
      }
      
      # otherwise if no request, retain the existing parameter
      else {
        self$output_to_user(paste("\tretaining existing parameter value", unstratified_name))
        return(unstratified_name)
      }
      
    },
    
    prepare_stratified_parameter_calculations = function() {
      #   prior to integration commencing, work out what the components are 
      #   of each parameter being implemented
      
      # create list of all the parameters that we need to find the list of adjustments for
      parameters_to_adjust <- c()
      
      for (n_flow in seq(nrow(self$transition_flows))){
        if(self$transition_flows$implement[n_flow] == length(self$strata) &
           !(self$transition_flows$parameter[n_flow] %in% parameters_to_adjust))
          parameters_to_adjust <- c(parameters_to_adjust, self$transition_flows$parameter[n_flow])
      }
      for (n_flow in seq(nrow(self$death_flows))){
        if(self$death_flows$implement[n_flow] == length(self$strata) &
           !(self$death_flows$parameter[n_flow] %in% parameters_to_adjust))
          parameters_to_adjust <- c(parameters_to_adjust, self$death_flows$parameter[n_flow])
      }    
      
      
      # and adjust
      for (parameter in parameters_to_adjust){
        self$parameter_components[[parameter]] = self$find_transition_components(parameter)
        self$create_transition_functions(parameter, self$parameter_components[parameter])
      }
      
      # similarly for all model compartments
      for (compartment in self$compartment_names){
        self$mortality_components[[compartment]] = self$find_mortality_components(compartment)
        self$create_mortality_functions(compartment, self$mortality_components[compartment])
      }
    },
    
    find_mortality_components = function(.compartment){
      # find the sub-parameters for population-wide natural mortality that are relevant to a particular compartment
      # used in prepare_stratified_parameter_calculations for creating functions to find the mortality rate for each
      # compartment
      # similar to find_transition_components, except being applied by compartment rather than parameter
      
      # :param _compartment: str
      # name of the compartment of interest
      # :return: all_sub_parameters: list
      # list of all the mortality-related sub-parameters for the compartment of interest
      all_sub_parameters = c()
      compartments_strata = find_name_components(.compartment)
      compartments_strata = rev(compartments_strata)
      compartments_strata = c(compartments_strata, "")
      
      # loop through each stratification and adapt if the parametr is available
      for (stratum in compartments_strata){
        if (stratum %in% self$available_death_rates)
          all_sub_parameters = c(all_sub_parameters, paste('universal_death_rateX', stratum, sep=''))
        if (paste('universal_death_rateX', stratum, sep='') %in% self$overwrite_parameters)
          break
      }
      all_sub_parameters = rev(all_sub_parameters)
      all_sub_parameters
    },
    
    create_mortality_functions = function(.compartment, .sub_parameters) {
      # loop through all the components to the population-wide mortality and
      # create the recursive functions
      
      # :param _compartment: str
      # name of the compartment of interest
      # :param _sub_parameters: list
      # the names of the functions that need to update the upstream parameters
      # :return:
      
      self$final_parameter_functions[paste('universal_death_rateX', .compartment, sep='')] =
        self$adaptation_functions[.sub_parameters[[1]]]
      for (component in tail(unlist(unname(.sub_parameters)), -1)){
        # get the new function to act on the less stratified function (closer to the "tree-trunk")
        if (!(component %in% self$parameters))
          stop(paste("parameter component", component, " not found in parameters attribute"))
        else if (typeof(self$parameters[component]) == 'double')
          update_function = create_multiplicative_function(self$parameters[[component]])
        else if (typeof(self$parameters[component]) == 'character')
          update_function = create_time_variant_multiplicative_function(self$adaption_functions[component])
        else
          stop(paste("parameter component", component, "not in appropriate format"))
        
        
        
        # create the composite function
        self$final_parameter_functions[[paste('universal_death_rateX', .compartment, sep='')]] = create_function_of_function(upadate_function, 
                                                                                                                             self$final_parameter_function[[paste('universal_death_rateX', .compartment, sep='')]] )
      }
    },
    
    find_transition_components = function(.parameter){
      # finds each of the strings for the functions acting on the next function in the sequence
      
      # :param _parameter: str
      # full name of the parameter of interest
      
      sub_parameters = c()
      
      # work backwards to allow stopping for overwriting requests, then 
      # reverse in preparation for function creation 
      for(x_instance in extract_reversed_x_positions(.parameter)){
        component = substr(.parameter, 1, x_instance)
        if(substr(component, nchar(component), nchar(component) ) == 'X')
          component =  substr(component, 1, nchar(component)-1)
        sub_parameters = c(sub_parameters, component )
        if (component %in% self$overwrite_parameters)
          break
      }
      sub_parameters = rev(sub_parameters)
      sub_parameters
    },
    
    create_transition_functions = function(.parameter, .sub_parameters){
      # builds up each parameter to be implemented as a function, recursively creating an outer function that calls the
      # inner function
      
      # :param _parameter: str
      # full name of the parameter of interest
      # :param _sub_parameters: list
      # list of the strings representing the sub-parameters, including the base parameter as the stem and with all
      # of the relevant strata in the stratification sequence following
      update_function <- NA
      # print('---subparam---')
      # print(unlist(unname(.sub_parameters))[1])
      # print('---subparam---')
      
      
      self$final_parameter_functions[[.parameter]] = (function(time){self$parameters[unlist(unname(.sub_parameters))[1]]})
      
      for (component in tail(unlist(unname(.sub_parameters)), -1)){
        if (!(component %in% names(self$parameters))){
          stop(paste("parameter component", component, " not found in parameters attribute"))
        }
        else if ((typeof(self$parameters[[component]]) == 'double') | (typeof(self$parameters[[component]]) == 'integer')){
          update_function = create_multiplicative_function(unname(self$parameters[[component]]))}
        else if(typeof(self$parameters[[component]]) == 'character'){
          update_function = create_time_variant_multiplicative_function(self$adaptation_functions[component]) 
        }
        else {
          stop(paste("parameter component", component, "not in appropriate format"))
        }
        
        
        self$final_parameter_functions[[.parameter]] = create_function_of_function(update_function, self$final_parameter_functions[[.parameter]],   max_depth = length(self$strata) + 8 )
        
      }    
    },
    
    
    find_parameter_components = function(.parameter) {
      #   extract the components of the stratified parameter into a dictionary structure with values being a list of
      #   time-variant parameters, a list of constant parameters and the product of all the constant values applied
      #
      #   :param .parameter: str
      #       name of the parameter that we are tracking down the components of
      
      # work backwards through sub-strings of the parameter names from the full name to the name through to each X
      
      # start from base value as a function, time variable currently a place-holder
      return_starting_parameter_value = function() {
        self$parameters[[find_stem(.parameter)]]
      }
      self$parameter_functions[[.parameter]] = return_starting_parameter_value
      
      
      # find the parameter component to cycle through for recursive function calls, excluding the first one
      parameters_of_interest = c()
      # extract x positions except the last one
      x_positions <- head(extract_reversed_x_positions(.parameter), -1)
      
      for (x_instance in x_positions[2: length(x_positions)]) {
        component <- substr(.parameter, 1, x_instance - 1)
        parameters_of_interest <- c(parameters_of_interest, component)
        if (component %in% names(self$overwrite_parameters)) {
          break
        }
      }
      
      # cycle through applicable components and extend function recursively
      for (component in parameters_of_interest ) {
        if (component %in% names(self$parameters)) {
          self$parameter_functions[[.parameter]] <- create_function_of_function(
            self$find_parameter_adaptation(component), self$parameter_functions[[.parameter]] , max_depth = length(self$strata) + 4
          )
          if (component %in% names(self$overwrite_parameters)) {
            break
          }
          
        }
      }
      
      
      self$parameter_components[[.parameter]] <- list(time_variants = c(), constants = c(), constant_value = 1)
      for (x_instance in extract_reversed_x_positions(.parameter)) {
        component <- substr(.parameter, 1, x_instance - 1)
        is_time_variant <- component %in% self$time_variants
        if (component %in% self$overwrite_parameters & is_time_variant) {
          self$parameter_components[[.parameter]] <- list(time_variants = c(component), constants = c(), constant_value = 1)
          break
        }
        else if (component %in% self$overwrite_parameters & !is_time_variant) {
          self$parameter_components[[.parameter]] <- list(time_variants = c(), constants = c(component), constant_value = 1)
          break
        }
        else if (is_time_variant) {
          self$parameter_components[[.parameter]]$time_variants <- c(self$parameter_components[[.parameter]]$time_variants, component)
        }
        else if (component %in% names(self$parameters)) {
          self$parameter_components[[.parameter]]$constants <- c(component,
                                                                 self$parameter_components[[.parameter]]$constants)
        }
      }
      
      # pre-calculate the constant component by multiplying through all the constant values
      for (constant_parameter in self$parameter_components[[.parameter]]$constants) {
        self$parameter_components[[.parameter]]$constant_value <-
          self$parameter_components[[.parameter]]$constant_value * 
          self$parameters[[constant_parameter]]
      }
    },
    
    find_parameter_adaptation = function(.component) {
      #   finds the parameter function using its string in the case that it
      #   is a function or creates a multiplicative
      #
      #   function as the default if a numeric value is provided as the 
      #   default behaviour for this user request format
      #   :param .component: str
      #       name of the parameter component
      #   :return: function
      #       a function that can act within the recursive approach to 
      #       function creation described in the previous method
      parameter_value <- self$parameters[[.component]]
      if (is.character(parameter_value)) {
        return(self$mapped_adaptation_functions[[.component]])
      }
      if (is.numeric(parameter_value)) {
        return(create_multiplicative_function(parameter_value))
      }
    },
    
    # __________
    # methods to be called during the process of model running
    
    get_parameter_value = function(.parameter, .time) {
      #   using the approach specified in find_parameter_components calculate
      #   adjusted parameter value from pre-calculated
      #   product of constant components and time variants
      #
      #   :param .parameter: char
      #       name of the parameter whose value is needed
      #   :param .time: numeric scalar value
      #       time in model integration
      
      if (is.numeric(self$final_parameter_functions[[.parameter]]))
        return(self$final_parameter_functions[[.parameter]])
      param_func = self$final_parameter_functions[[.parameter]]
      return(param_func(time))
    },
    
    find_infectious_population = function(.compartment_values) {
      # find vectors for the total infectious populations and the total population that is needed in the case of
      # frequency-dependent transmission
      
      # :param _compartment_values: numpy array
      #     current values for the compartment sizes
      
      
      mixing_categories = c()
      if(is.na(self$mixing_matrix))
        mixing_categories = c('all_population')
      else
        mixing_categories = self$mixing_categories
      strain_list = c()
      if(is.na(self$strains))
        strain_list = c('all_strains')
      else
        strain_list = self$strains
      for (strain in strain_list){
        
        new_list = NA
        for (category in mixing_categories){
          elem_list = self$strain_mixing_elements[ (self$strain_mixing_elements$strain==strain & 
                                                      self$strain_mixing_elements$category==category), "idx"]
          mult_list = self$strain_mixing_multipliers[ (self$strain_mixing_multipliers$strain==strain & 
                                                         self$strain_mixing_multipliers$category==category), "idx"]
          sum_list  <- sum(element_list_multiplication(.compartment_values[elem_list], mult_list ))
          new_list = c(sum_list, tail(self$infectious_populations[[strain]], n=1))
          write(new_list, file="new_list.txt", append = TRUE)
        }
      }
      self$infectious_populations[[strain]] = new_list
      self$infectious_denominators = sum(.compartment_values)
    },
    
    find_infectious_multiplier = function(n_flow){
      # find the multiplier to account for the infectious population in dynamic flows
      
      # :param n_flow: int
      # index for the row of the transition_flows data frame
      # :return:
      #  the total infectious quantity, whether that is the number or proportion of infectious persons
      # needs to return as one for flows that are not transmission dynamic infectiousness flows
      
      if (grepl("infection", self$transition_flows$type[n_flow])){
        strain = "all_strains"
        if(is.na(self$strains)){
          strain = "all_strains"
        }else{
          strain = self$transition_flows$strain[n_flow]
        }
      } else
        return(1)
      mixing_elements = c()
      if(is.na(self$mixing_matrix))
        mixing_elements = c(1.0)
      else
        mixing_elements = c(self$mixing_matrix[self$transition_flows$force_index[n_flow]-1,]) #force indexsubtract by 2 to make
      
      if (grepl('_density', self$transition_flows$type[n_flow]))
        denominator = 1.0
      else{
        
        denominator = self$infectious_denominators
        
      }
      
      
      
      sum(element_list_multiplication(unname(unlist(self$infectious_populations[strain])),
                                      mixing_elements))/denominator
      
      
    },
    
    apply_birth_rate = function(.ode_equations, .compartment_values) {
      #   apply a population-wide death rate to all compartments
      #
      #   :parameters: all parameters have come directly from the apply_all_flow_types_to_odes method unchanged
      total_births = self$find_total_births(.compartment_values)
      
      # split the total births across entry compartments
      for (compartment in names(.compartment_values)) {
        if (find_stem(compartment) == self$entry_compartment) {
          
          # calculate adjustment to original stem entry rate
          entry_fraction <- 1
          x_positions <- extract_x_positions(compartment)
          if (!x_positions[1] == -1) {
            for (x_instance in seq(length(x_positions) - 1)) {
              if (identical(x_positions[x_instance], numeric(0)))
                break
              if (identical(x_positions[x_instance], integer(0)))
                break
              entry_fraction_param = paste("entry_fractionX",
                                           substr(compartment, x_positions[x_instance] + 1, 
                                                  x_positions[x_instance + 1]), sep="")
              if (entry_fraction_param %in% names(self$parameters))
                entry_fraction <- entry_fraction * self$parameters[[entry_fraction_param]]
            }
          }
          compartment_births <- entry_fraction * total_births
          .ode_equations <- increment_list_by_index(.ode_equations, match(compartment, 
                                                                          names(self$compartment_values)), compartment_births)
        }
      }
      .ode_equations
    }
  )
)


ModelInterpreter <- R6Class(
  "Interpreter",
  public = list(
    model = NULL,
    times = c(),
    outputs = NULL,
    compartment_types = list(),
    compartment_totals = data.frame(),
    compartment_capitalised = '',
    table_final_outputs = data.frame(),
    initialize = function(model) {
      self$model <- model
      self$outputs <- self$model$outputs
      self$times <- self$outputs$time
      self$compartment_types <- self$model$compartment_types
      self$find_compartment_totals()
      self$table_final_outputs <- cbind(self$outputs, self$compartment_totals)
      self$table_final_outputs <- self$table_final_outputs[!duplicated(as.list(self$table_final_outputs))]
      for (time_value in 1:length(self$table_final_outputs$time)) {
        self$table_final_outputs$time[[time_value]] <- self$table_final_outputs$time[[time_value]] * 365
      }
    },
    
    # sum all the compartment values of one type
    find_compartment_totals = function() {
      self$compartment_totals <-
        data.frame(matrix(NA, nrow=length(self$times), ncol=0))
      for (compartment_type in self$compartment_types) {
        self$compartment_totals[[compartment_type]] <- 0
        
        for (compartment in names(self$outputs)) {
          if (find_stem(compartment) == compartment_type) {
            self$compartment_totals[[compartment_type]] <-
              self$compartment_totals[[compartment_type]] + self$outputs[[compartment]]
          }
        }
      }
    },
    
    # allows the user to plot compartment/multiple compartments against time
    plot_function = function(compartment = c("susceptible", "infectious"), points = FALSE) {
      
      # final_data is populated with key/value columns that is put into ggplot
      final_data <- data.frame(interpreter$table_final_outputs %>%
                                 gather(Compartments, value, compartment))
      
      # ggplot is initiated
      plot <- ggplot(final_data, aes(time, value, color = Compartments)) + geom_line()
      
      # capitalise compartment name
      compartment_capitalised <- capitalise_compartment_name(compartment)
      
      # if points set as true, points signifying each day is placed
      if (points) {
        plot <- plot + geom_point(size = 0.4, color = 'red')
      }
      
      # labels are cleaned - depending on how many variables are used
      if (length(compartment) == 1) {
        plot <- plot + labs(title = paste(compartment_capitalised[[1]],
                                          "over time"),
                            x = "Time", y = compartment_capitalised[[1]])
      }
      else {
        plot <- plot + labs(title = paste("Compartment sizes over time"),
                            x = "Time", y = "Proportion of people")
      }
      
      # the legend is also cleaned up with necessary outputs
      plot <- plot + scale_color_discrete(breaks = compartment,
                                          labels = compartment_capitalised) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      
      # save function
      ggsave(paste('graph~',
                   paste(compartment, collapse = '~'),
                   ".pdf",
                   sep = ''))
    },
    
    # simple method to plot the size of a compartment
    plot_compartment = function(compartment) {
      plot(self$times, self$compartment_totals[[compartment]])
    },
    
    # creates a flowchart - stratified flowchart unless otherwise specified
    # parameters can also be presented unless otherwise specified
    create_flowchart = function(type = 'stratified', strata = NULL, parameters = TRUE) {
      
      
      #Pick type of input into the function, depending on whether the type of flowchart is
      if (type == 'stratified') {
        if (is.null(strata)){
          input_nodes <- names(self$model$compartment_values)
          new_df <- self$model$transition_flows
          type_of_flow <- new_df[which(new_df$implement == length(self$model$strata)),]
        }
        else{
          new_df <- self$model$transition_flows
          type_of_flow <- new_df[which(new_df$implement == strata),]
          input_nodes <- unique(type_of_flow$from,type_of_flow$to)
        }
        
      }
      else if (type == 'unstratified') {
        input_nodes <- self$model$compartment_types
        type_of_flow <- self$model$unstratified_flows
        
      }
      else {
        stop("Type needs to be either stratified or unstratified.")
      }
      
      #Inputs are sectioned according to the
      #stem value so colours can be added to each type.
      #broken_down_nodes list created
      
      broken_down_nodes = list()
      
      #broken_down_nodes is populated with different list for each stem value
      
      for (stem_value in 1:length(self$model$compartment_types)) {
        x_vector <- c()
        for (stem_type in 1:length(input_nodes)) {
          if (self$model$compartment_types[[stem_value]] == find_stem(input_nodes[[stem_type]])) {
            x_vector <- c(x_vector, input_nodes[[stem_type]])
          }
        }
        broken_down_nodes[[stem_value]] <- x_vector
      }
      
      #The colours of each stem value (compartment type) is created
      #The settings string is set
      
      settings <- ''
      
      #Settings is populated with the string necessary for grViz function
      
      for (list_different_nodes in 1:length(broken_down_nodes)) {
        nodes <- c()
        nodes <- paste(broken_down_nodes[[list_different_nodes]], collapse = ' ')
        settings <- paste(settings, 'node [shape = box,
                          fontname = Helvetica, style = filled, color =',
                          c('BlanchedAlmond', 'Grey',
                            'RosyBrown', 'LavenderBlush',
                            'Salmon', 'LightPink',
                            'PaleGreen', 'Thistle',
                            'Beige', 'PeachPuff',
                            'MintCream', 'AquaMarine',
                            'MistyRose', 'Tomato',
                            'Honeydew', 'LightCyan')[[sample(1:16,
                                                             1,
                                                             replace = FALSE,
                                                             prob = NULL)]],
                          ']', nodes)
      }
      
      #The pathways between nodes are set as empty
      
      connection_between_nodes <- ""
      connections <- ''
      
      #The pathway between nodes is populated from type_of_flow, as well as the parameters
      
      
      for (row in seq(nrow(type_of_flow))) {
        if (type_of_flow$implement[[row]]) {
          
          #Parameters are added or not added in to the flowchart depending on the setting
          
          if (parameters) {
            connections <- paste(' edge [label =',
                                 type_of_flow$parameter[[row]],
                                 ']')
          }
          else if (!parameters) {
            connections <- ''
          }
          connection_between_nodes <- paste(connection_between_nodes,
                                            connections,
                                            type_of_flow$from[[row]],
                                            "->", type_of_flow$to[[row]],
                                            ' ', sep = '')
        }}
      
      #The final string necessary for grViz is created here
      input_for_grViz <- ''
      input_for_grViz <- paste("digraph dot {
                               graph [layout = dot,
                               rankdir = LR]",
                               settings,
                               connection_between_nodes, '}')
      
      # 'X' are substituted for '_' and input into the function
      
      final_flowchart <- grViz(str_replace_all(input_for_grViz, "X", "_"))
      svg <- export_svg(final_flowchart)
      svg <- charToRaw(svg)
      rsvg_pdf(svg, paste('flowchart~', type , '.pdf', sep = ''))
    }
  )
  )

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
# hiv_mixing <- NULL
hiv_mixing <- matrix(1, ncol = 2, nrow = 2)
sir_model$stratify("hiv", c("negative", "positive"), c(), list("negative"=0.6),
                   list(recovery=list("negative"=0.7, "positive"=0.5),
                        infect_death=list("negative"=0.5),
                        entry_fraction=list("negative"=0.6, "positive"=0.4)),
                   #adjustment_requests = list(recovery = list("negative"=0.7)),
                   infectiousness_adjustments = c("positive"=0.5), mixing_matrix=hiv_mixing, verbose = FALSE)

sir_model$stratify("strain", c("sensitive", "resistant"),  c("infectious"), requested_proportions=c(), verbose = FALSE)
age_mixing <- NULL
# #age_mixing <- diag(4)
sir_model$stratify("age", c(1, 10, 3), c(), c(),
                   list(recovery=list("1"=0.5, "10"=0.8)),
                   #recoveryXhiv_positive=list("1"=2, "3"=365/13*.5, "10"=1, overwrite=c("2")),
                   #universal_death_rate=list("1"=1, "2"=2, "3"=3)),
                   infectiousness_adjustments=c("1"=0.8),
                   mixing_matrix = age_mixing, verbose=FALSE)

#sir_model$stratify("age", c(3, 2, 1), c(), verbose = TRUE)
# sir_model$add_time_variant("recovery", create_arbitrary_time_variant_function)

sir_model$run_model()

interpreter <- ModelInterpreter$new(sir_model)
# interpreter$create_flowchart()
interpreter$plot_compartment("infectious")



