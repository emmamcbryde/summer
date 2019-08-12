
from python_source_code.summer_model import *
from python_source_code.db import InputDB
import matplotlib.pyplot
import os
import numpy
import scipy.integrate
import copy
from python_source_code.curve import scale_up_function
import pandas as pd


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
        age_breakpoints, parameter_names=("early_progression", "stabilisation", "late_progression")):
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
                            provide_age_specific_latency_parameters()[parameter]), age_breakpoints), 365.251))
    return adapted_parameter_dict


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


def unpivot_outputs(model_object):
    """
    take outputs in the form they come out of the model object and convert them into a "long", "melted" or "unpiovted"
    format in order to more easily plug to PowerBI
    """
    output_dataframe = pd.DataFrame(model_object.outputs, columns=model_object.compartment_names)
    output_dataframe["times"] = model_object.times
    output_dataframe = output_dataframe.melt("times")
    for n_stratification in range(len(model_object.strata) + 1):
        column_name = "compartment" if n_stratification == 0 else model_object.strata[n_stratification - 1]
        output_dataframe[column_name] = \
            output_dataframe.apply(lambda row: row.variable.split("X")[n_stratification], axis=1)
        if n_stratification > 0:
            output_dataframe[column_name] = \
                output_dataframe.apply(lambda row: row[column_name].split("_")[1], axis=1)
    return output_dataframe.drop(columns="variable")


"""
standardised flow functions
"""


def add_standard_latency_flows(list_of_flows):
    """
    adds our standard latency flows to the list of flows to be implemented in the model

    :param list_of_flows: list
        existing flows for implementation in the model
    :return: list_of_flows: list
        list of flows updated to include the standard latency flows
    """
    list_of_flows += [
        {"type": "standard_flows", "parameter": "early_progression", "origin": "early_latent", "to": "infectious"},
        {"type": "standard_flows", "parameter": "stabilisation", "origin": "early_latent", "to": "late_latent"},
        {"type": "standard_flows", "parameter": "late_progression", "origin": "late_latent", "to": "infectious"}]
    return list_of_flows


def add_standard_natural_history_flows(list_of_flows):
    """
    adds our standard natural history to the list of flows to be implemented in the model

    :param list_of_flows: list
        existing flows for implementation in the model
    :return: list_of_flows: list
        list of flows updated to include the standard latency flows
    """
    list_of_flows += [
        {"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
        {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}]
    return list_of_flows


def add_standard_infection_flows(list_of_flows):
    """
    adds our standard infection processes to the list of flows to be implemented in the model

    :param list_of_flows: list
        existing flows for implementation in the model
    :return: list_of_flows: list
        list of flows updated to include the standard infection processes
    """
    list_of_flows += [
        {"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "early_latent"},
        {"type": "infection_frequency", "parameter": "beta", "origin": "recovered", "to": "early_latent"}]
    return list_of_flows


"""
main model construction functions
"""


def build_working_tb_model(tb_n_contact, cdr_adjustment=0.6, start_time=1800.):
    """
    current working tb model with some characteristics of mongolia applied at present
    """
    integration_times = numpy.linspace(start_time, 2020.0, 201).tolist()

    # set basic parameters, flows and times, then functionally add latency
    case_fatality_rate = 0.4
    untreated_disease_duration = 3.0
    parameters = \
        {"beta": tb_n_contact,
         "recovery": case_fatality_rate / untreated_disease_duration,
         "infect_death": (1.0 - case_fatality_rate) / untreated_disease_duration,
         "universal_death_rate": 1.0 / 50.0,
         "case_detection": 0.0}
    parameters.update(change_parameter_unit(provide_aggregated_latency_parameters(), 365.251))

    # sequentially add groups of flows
    flows = add_standard_infection_flows([])
    flows = add_standard_latency_flows(flows)
    flows = add_standard_natural_history_flows(flows)

    # compartments
    compartments = ["susceptible", "early_latent", "late_latent", "infectious", "recovered"]

    # define model
    _tb_model = \
        StratifiedModel(integration_times, compartments, {"infectious": 1e-3}, parameters, flows, birth_approach="replace_deaths")

    # add case detection process to basic model
    _tb_model.add_transition_flow(
        {"type": "standard_flows", "parameter": "case_detection", "origin": "infectious", "to": "recovered"})

    # create_flowchart(_tb_model, name="unstratified")

    # loading time-variant case detection rate
    input_database = InputDB()
    res = input_database.db_query("gtb_2015", column="c_cdr", is_filter="country", value="Mongolia")

    # add scaling case detection rate
    cdr_adjustment_factor = cdr_adjustment
    cdr_mongolia = res["c_cdr"].values / 1e2 * cdr_adjustment_factor
    cdr_mongolia = numpy.concatenate(([0.0], cdr_mongolia))
    res = input_database.db_query("gtb_2015", column="year", is_filter="country", value="Mongolia")
    cdr_mongolia_year = res["year"].values
    cdr_mongolia_year = numpy.concatenate(([1950.], cdr_mongolia_year))
    cdr_scaleup = scale_up_function(cdr_mongolia_year, cdr_mongolia, smoothness=0.2, method=5)
    prop_to_rate = convert_competing_proportion_to_rate(1.0 / untreated_disease_duration)
    detect_rate = return_function_of_function(cdr_scaleup, prop_to_rate)
    _tb_model.time_variants["case_detection"] = detect_rate

    # store scaling functions in database if required
    # function_dataframe = pd.DataFrame(times)
    # function_dataframe["cdr_values"] = [cdr_scaleup(t) for t in times]
    # store_database(function_dataframe, table_name="functions")

    # test strain stratification
    strain_only_model = copy.deepcopy(_tb_model)
    strain_only_model.stratify("strain", ["ds", "mdr"], ["early_latent", "late_latent", "infectious"], {},
                               verbose=False)
    create_flowchart(strain_only_model, name="stratified_by_strain")

    # test age stratification
    # age_only_model = copy.deepcopy(_tb_model)
    age_breakpoints = [0, 6, 13, 15]
    age_infectiousness = get_parameter_dict_from_function(logistic_scaling_function(15.0), age_breakpoints)
    age_params = get_adapted_age_parameters(age_breakpoints)
    # age_only_model.stratify("age", copy.deepcopy(age_breakpoints), [], {},
    #                         adjustment_requests=age_params,
    #                         infectiousness_adjustments=age_infectiousness,
    #                         verbose=False)
    # create_flowchart(age_only_model, name="stratified_by_age")

    # test organ stratification
    # organ_only_model = copy.deepcopy(_tb_model)
    # organ_only_model.stratify("smear",
    #                           ["smearpos", "smearneg", "extrapul"],
    #                           ["infectious"], adjustment_requests={}, verbose=False, requested_proportions={})
    # create_flowchart(organ_only_model, name="stratified_by_organ")

    # test risk stratification
    # risk_only_model = copy.deepcopy(_tb_model)
    # risk_only_model.stratify("risk",
    #                          ["urban", "urbanpoor", "ruralpoor"], [], requested_proportions={},
    #                          adjustment_requests={}, verbose=False)
    # create_flowchart(risk_only_model, name="stratified_by_risk")

    # _tb_model.stratify("strain", ["ds", "mdr"], ["early_latent", "late_latent", "infectious"], {},
    #                    verbose=False)
    # _tb_model.stratify("age", age_breakpoints, [], {},
    #                    adjustment_requests=age_params,
    #                    infectiousness_adjustments=age_infectiousness,
    #                    verbose=False)
    _tb_model.stratify("smear",
                       ["smearpos", "smearneg", "extrapul"],
                       ["infectious"],
                       infectiousness_adjustments={"smearneg": 0.24, "extrapul": 0.0},
                       verbose=False, requested_proportions={})
    # _tb_model.stratify("risk",
    #                    ["urban", "urbanpoor", "ruralpoor"], [], requested_proportions={},
    #                    adjustment_requests=[], verbose=False)
    return _tb_model


if __name__ == "__main__":

    tb_model = build_working_tb_model(40.0)

    # create_flowchart(tb_model, name="mongolia_flowchart")
    # tb_model.transition_flows.to_csv("transitions.csv")

    tb_model.run_model()

    # get outputs
    # infectious_population = tb_model.get_total_compartment_size(["infectious"])

    # print statements to enable crude manual calibration
    # time_2016 = [i for i in range(len(tb_model.times)) if tb_model.times[i] > 2016.][0]
    # print(time_2016)
    # print(infectious_population[time_2016] * 1e5)
    # print(cdr_mongolia)

    # output the results into a format that will be easily loadable into PowerBI
    # pbi_outputs = unpivot_outputs(tb_model)
    # store_database(pbi_outputs, table_name="pbi_outputs")

    # easy enough to output a flow diagram if needed:
    create_flowchart(tb_model)

    # output some basic quantities if not bothered with the PowerBI bells and whistles
    # tb_model.plot_compartment_size(["early_latent", "late_latent"])
    tb_model.plot_compartment_size(["infectious"], 1e5)

    # store outputs into database
    # tb_model.store_database()
    #
    # matplotlib.pyplot.plot(numpy.linspace(1800., 2020.0, 201).tolist(), infectious_population * 1e5)
    # matplotlib.pyplot.xlim((1980., 2020.))
    # matplotlib.pyplot.ylim((0.0, 1e3))
    matplotlib.pyplot.show()
    # matplotlib.pyplot.savefig("mongolia_cdr_output")



