
from summer_model import *
from db import InputDB
import matplotlib.pyplot
import os
import numpy
import scipy.integrate
import copy
from curve import scale_up_function


def provide_aggregated_latency_parameters():
    """
    function to add the latency parameters estimated by Ragonnet et al from our paper in Epidemics to the existing
    parameter dictionary
    """
    return {"early_progression": 1.1e-3, "stabilisation": 1.0e-2, "late_progression": 5.5e-6}


def provide_age_specific_latency_parameters():
    """
    simply list out all the latency progression parameters from Ragonnet et al as a dictionary
    """
    return {"early_progression": {0: 6.6e-3, 5: 2.7e-3, 15: 2.7e-4},
            "stabilisation": {0: 1.2e-2, 5: 1.2e-2, 15: 5.4e-3},
            "late_progression": {0: 1.9e-11, 5: 6.4e-6, 15: 3.3e-6}}


def get_adapted_age_parameters(
        age_breakpoints_, parameter_names=("early_progression", "stabilisation", "late_progression")):
    """
    get age-specific parameters adapted to any specification of age breakpoints
    """
    adapted_parameter_dict = {}
    for parameter in parameter_names:
        adapted_parameter_dict[parameter] = \
            add_w_to_param_names(
                change_parameter_unit(
                    get_parameter_dict_from_function(
                        create_step_function_from_dict(
                            provide_age_specific_latency_parameters()[parameter]), age_breakpoints_)))
    return adapted_parameter_dict


def add_standard_latency_flows(flows_):
    """
    adds our standard latency flows to the list of flows to be implemented in the model
    """
    flows_ += [
        {"type": "standard_flows", "parameter": "early_progression", "origin": "early_latent", "to": "infectious"},
        {"type": "standard_flows", "parameter": "stabilisation", "origin": "early_latent", "to": "late_latent"},
        {"type": "standard_flows", "parameter": "late_progression", "origin": "late_latent", "to": "infectious"}]
    return flows_


def convert_competing_proportion_to_rate(competing_flows):
    """
    convert a proportion to a rate dependent on the other flows coming out of a compartment
    """
    return lambda proportion: proportion * competing_flows / (1.0 - proportion)


def return_function_of_function(inner_function, outer_function):
    """
    general method to return a chained function from two functions
    """
    return lambda value: outer_function(inner_function(value))


if __name__ == "__main__":

    # set basic parameters, flows and times, except for latency flows and parameters, then functionally add latency
    case_fatality_rate = 0.4
    untreated_disease_duration = 3.0
    parameters = \
        {"beta": 10.0,
         "recovery": case_fatality_rate / untreated_disease_duration,
         "infect_death": (1.0 - case_fatality_rate) / untreated_disease_duration,
         "universal_death_rate": 1.0 / 50.0,
         "case_detection": 0.0}
    parameters.update(change_parameter_unit(provide_aggregated_latency_parameters()))

    times = numpy.linspace(1800., 2020.0, 201).tolist()
    flows = [{"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "early_latent"},
             {"type": "infection_frequency", "parameter": "beta", "origin": "recovered", "to": "early_latent"},
             {"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
             {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}]
    flows = add_standard_latency_flows(flows)

    tb_model = StratifiedModel(
        times, ["susceptible", "early_latent", "late_latent", "infectious", "recovered"], {"infectious": 1e-3},
        parameters, flows, birth_approach="replace_deaths")

    tb_model.add_transition_flow(
        {"type": "standard_flows", "parameter": "case_detection", "origin": "infectious", "to": "recovered"})

    # loading time-variant case detection rate
    input_database = InputDB(report=True)
    res = input_database.db_query("gtb_2015", column="c_cdr", is_filter="country", value="Mongolia")
    cdr_mongolia = res["c_cdr"].values / 1e2
    cdr_mongolia = numpy.concatenate(([0.0], cdr_mongolia))
    res = input_database.db_query("gtb_2015", column="year", is_filter="country", value="Mongolia")
    cdr_mongolia_year = res["year"].values
    cdr_mongolia_year = numpy.concatenate(([1950.], cdr_mongolia_year))
    cdr_scaleup = scale_up_function(cdr_mongolia_year, cdr_mongolia, smoothness=0.2, method=5)

    prop_to_rate = convert_competing_proportion_to_rate(1.0 / untreated_disease_duration)
    detect_rate = return_function_of_function(cdr_scaleup, prop_to_rate)

    tb_model.time_variants["case_detection"] = detect_rate

    age_breakpoints = [0, 5, 15]
    age_infectiousness = get_parameter_dict_from_function(logistic_scaling_function(15.0), age_breakpoints)
    #
    # x_values = numpy.linspace(0.0, 40.0, 1e3)
    # y_values = [logistic_scaling_function(15.0)(x) for x in x_values]
    # matplotlib.pyplot.plot(x_values, y_values)
    # matplotlib.pyplot.savefig("age_infectiousness_scaling")

    # cdr_values = [cdr_scaleup(x) for x in times]
    # matplotlib.pyplot.plot(times, cdr_values)
    # matplotlib.pyplot.savefig("mongolia_cdr_scaling")

    tb_model.stratify("age", age_breakpoints, [],
                      adjustment_requests=get_adapted_age_parameters(age_breakpoints),
                      infectiousness_adjustments=age_infectiousness,
                      report=False)
    tb_model.run_model()

    # get outputs - from here on the code is essentially rubbish - will be using PowerBI to view outputs
    infectious_population = tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_0")] + \
                            tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_5")] + \
                            tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_15")]

    matplotlib.pyplot.plot(times, infectious_population * 1e5)
    matplotlib.pyplot.xlim((1980., 2010.))
    matplotlib.pyplot.ylim((0.0, 500.0))
    matplotlib.pyplot.savefig("mongolia_cdr_output")


