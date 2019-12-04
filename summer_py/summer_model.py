import numpy
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot
import copy
import pandas as pd
from graphviz import Digraph
import os
import itertools
from sqlalchemy import create_engine


# use conda to install graphviz with "conda install python-graphviz" instead of manually setting path

"""
string manipulation functions
"""


def create_stratum_name(stratification_name, stratum_name, joining_string="X"):
    """
    generate a name string to represent a particular stratum within a requested stratification

    :param stratification_name: str
        the "stratification" or rationale for implementing the current stratification process
    :param stratum_name: str
        name of the stratum within the stratification
    :param joining_string: str
        the character to add to the front to indicate that this string is the extension of the existing one
        in SUMMER, capitals are reserved for non-user-requested strings, in this case "X" is used as the default
    :return: str
        the composite string for the stratification
    """
    return joining_string + "%s_%s" % (stratification_name, str(stratum_name))


def create_stratified_name(stem, stratification_name, stratum_name):
    """
    generate a standardised stratified compartment name

    :param stem: str
        the previous stem to the compartment or parameter name that needs to be extended
    :param stratification_name: str
        the "stratification" or rationale for implementing the current stratification process
    :param stratum_name: str
        name of the stratum within the stratification
    :return: str
        the composite name with the standardised stratification name added on to the old stem
    """
    return stem + create_stratum_name(stratification_name, stratum_name)


def extract_x_positions(parameter, joining_string="X"):
    """
    find the positions within a string which are X and return as list, including length of list

    :param parameter: str
        the string for interrogation
    :param joining_string: str
        the string of interest whose character positions need to be found
    :return: list
        list of all the indices for where the X character occurs in the string, along with the total length of the list
    """
    return [loc for loc in range(len(parameter)) if parameter[loc] == joining_string] + [len(parameter)]


def extract_reversed_x_positions(parameter):
    """
    find the positions within a string which are X and return as list reversed, including length of list

    :params and return: see extract_x_positions above
    """
    result = extract_x_positions(parameter)
    result.reverse()
    return result


def find_stem(stratified_string, joining_string="X"):
    """
    find the stem of the compartment name as the text leading up to the first occurrence of the joining string
        (usually "X")
    should run slightly faster than using find_name_components

    :param stratified_string: str
        the stratified string for the compartment or parameter name
    :param joining_string: str
        the string of interest whose character position will be used to truncate the stratified string
    :return: int
        the point at which the first occurrence of the joining string occurs
    """
    index = stratified_string.find(joining_string) if joining_string in stratified_string else len(stratified_string)
    return stratified_string[:index]


def find_name_components(compartment):
    """
    extract all the components of a stratified compartment or parameter name, including the stem

    :param compartment: str
        name of the compartment or parameter to be interrogated
    :return: list
        the extracted compartment components
    """

    # add -1 at the start, which becomes zero to represent the start when one is added
    x_positions = [-1] + extract_x_positions(compartment)

    # add one to the first index to go past the joining character, but not at the end
    return [compartment[x_positions[x_pos] + 1: x_positions[x_pos + 1]] for x_pos in range(len(x_positions) - 1)]


def find_stratum_index_from_string(compartment, stratification, remove_stratification_name=True):
    """
    finds the stratum which the compartment (or parameter) name falls in when provided with the compartment name and the
        name of the stratification of interest
    for example, if the compartment name was infectiousXhiv_positiveXdiabetes_none and the stratification of interest
        provided through the stratification argument was hiv, then return positive

    :param compartment: str
        name of the compartment (or parameter) to be interrogated
    :param stratification: str
        the stratification of interest
    :param remove_stratification_name: bool
        whether to remove the stratification name and its trailing _ from the string to return
    :return: str
        the name of the stratum within which the compartment falls
    """
    stratum_name = [name for n_name, name in enumerate(find_name_components(compartment)) if stratification in name][0]
    return stratum_name[stratum_name.find("_") + 1:] if remove_stratification_name else stratum_name


"""
data manipulation functions
"""


def increment_list_by_index(list_to_increment, index_to_increment, increment_value):
    """
    very simple but general method to increment the odes by a specified value

    :param list_to_increment: list
        the list to be incremented, expected to be the list of ODEs
    :param index_to_increment: int
        the index of the list that needs to be incremented
    :param increment_value: float
        the value to increment the list by
    :return: list
        the list after incrementing
    """
    list_to_increment[index_to_increment] += increment_value
    return list_to_increment


def normalise_dict(value_dict):
    """
    normalise the values from a list using the total of all values, i.e. returning dictionary whose keys sum to one with
        same ratios between all the values in the dictionary after the function has been applied

    :param value_dict: dict
        dictionary whose values will be adjusted
    :return: dict
        same dictionary after values have been normalised to the total of the original values
    """
    return {key: value_dict[key] / sum(value_dict.values()) for key in value_dict}


def order_dict_by_keys(input_dict):
    """
    sort the input dictionary keys and return two separate lists with keys and values as lists with corresponding
        elements

    :param input_dict: dict
        dictionary to be sorted
    :return:
        :dict_keys: list
            sorted list of what were the dictionary keys
        : list
            values applicable to the sorted list of dictionary keys
    """
    dict_keys = list(input_dict.keys())
    dict_keys.sort()
    return dict_keys, [input_dict[key] for key in dict_keys]


def element_list_multiplication(list_1, list_2):
    """
    multiply elements of two lists to return another list with the same dimensions

    :param list_1: list
        first list of numeric values to be multiplied
    :param list_2: list
        second list of numeric values to be multiplied
    :return: list
        resulting list populated with the multiplied values
    """
    return [a * b for a, b in zip(list_1, list_2)]


def element_list_division(list_1, list_2):
    """
    divide elements of two lists to return another list with the same dimensions

    :param list_1: list
        first list of numeric values for numerators of division
    :param list_2: list
        second list of numeric values for denominators of division
    :return: list
        list populated with quotient values
    """
    return [a / b for a, b in zip(list_1, list_2)]


def convert_boolean_list_to_indices(list_of_booleans):
    """
    take a list of boolean values and return the indices of the elements containing True

    :param list_of_booleans: list
        sequence of booleans
    :return: list with integer values
        list of the values that were True in the input list_of_booleans
    """
    return [n_element for n_element, element in enumerate(list_of_booleans) if element]


"""
functions needed for dealing with age stratification
"""


def add_zero_to_age_breakpoints(breakpoints):
    """
    append a zero on to a list if there isn't one already present, for the purposes of age stratification

    :param breakpoints: list
        integers for the age breakpoints requested
    :return: list
        age breakpoints with the zero value included
    """
    return [0] + breakpoints if 0 not in breakpoints else breakpoints


def split_age_parameter(age_breakpoints, parameter):
    """
    creates a dictionary to request splitting of a parameter according to age breakpoints, but using values of 1 for
        each age stratum
    allows that later parameters that might be age-specific can be modified for some age strata

    :param age_breakpoints: list
        list of the age breakpoints to be requested, with breakpoints as string
    :param parameter: str
        name of parameter that will need to be split
    :return: dict
        dictionary with age groups as string as keys and ones for all the values
    """
    age_breakpoints = ["0"] + age_breakpoints if "0" not in age_breakpoints else age_breakpoints
    return {parameter: {str(age_group): 1.0 for age_group in age_breakpoints}}


"""
functions of functions for use in stratified models
"""


def create_multiplicative_function(multiplier):
    """
    return multiplication by a fixed value as a function

    :param multiplier: float
        value that the returned function multiplies by
    :return: function
        function that can multiply by the multiplier parameter when called
    """
    return lambda input_value, time: multiplier * input_value


def create_time_variant_multiplicative_function(time_variant_function):
    """
    similar to create_multiplicative_function, except that the value to multiply by can be a function of time, rather
        than a single value

    :param time_variant_function: function
        a function with the independent variable of time that returns the value that the input should be multiplied by
    :return: function
        function that will multiply the input value by the output value of the time_variant_function
    """
    return lambda input_value, time: time_variant_function(time) * input_value


def create_additive_function(increment):
    """
    return the addition of a fixed value as a function

    :param increment: float
        value that the returned function increments by
    :return: function
        function that can increment by the value parameter when called
    """
    return lambda value: value + increment


def create_sloping_step_function(start_x, start_y, end_x, end_y):
    """
    create sloping step function, returning start y-value for input values below starting x-value, ending y-value
        for input values above ending x-value and connecting slope through the middle

    :param start_x: float
        starting x-value
    :param start_y: float
        starting y-value
    :param end_x: float
        ending x-value
    :param end_y: float
        ending y-value
    :return: function
        sloping function as described above
    """
    gradient = (end_y - start_y) / (end_x - start_x)

    def step_function(age):
        if age < start_x:
            return start_y
        elif start_x <= age < end_x:
            return gradient * age + start_y - gradient * start_x
        elif end_x <= age:
            return end_y

    return step_function


def create_function_of_function(outer_function, inner_function):
    """
    function that can itself return a function that sequentially apply two functions and so can be used recursively to
        create a series of functions

    :param outer_function: function
        last function to be called
    :param inner_function: function
        first function to be called
    :return: function
        composite function that applies the inner and then the outer function, allowing the time parameter to be passed
            through if necessary
    """
    return lambda time: outer_function(inner_function(time), time)


"""
flow diagram creation
"""


def create_flowchart(model_object, strata=None, name="flow_chart"):
    """
    use graphviz module to create flow diagram of compartments and inter-compartmental flows

    :param model_object: summer object
        model whose inter-compartmental flows need to be graphed
    :param strata: int
        number of stratifications that have been implemented at the point that diagram creation requested
    :param name: str
        filename for the image to be put out as
    """

    # find the stratification level of interest, with the fully stratified model being the default
    if strata is None:
        strata = len(model_object.all_stratifications)

    # set styles for graph
    styles = {"graph": {"label": "",
                        "fontsize": "16", },
              "nodes": {"fontname": "Helvetica",
                        "style": "filled",
                        "fillcolor": "#CCDDFF", },
              "edges": {"style": "dotted",
                        "arrowhead": "open",
                        "fontname": "Courier",
                        "fontsize": "10", }}

    # colour dictionary for different nodes indicating different stages of infection
    default_colour_dict = {"susceptible": "#F0FFFF", "early_latent": "#A64942", "late_latent": "#A64942",
                           "infectious": "#FE5F55", "recovered": "#FFF1C1"}

    def apply_styles(graph, _styles):
        graph.graph_attr.update(("graph" in _styles and _styles["graph"]) or {})
        graph.node_attr.update(("nodes" in _styles and _styles["nodes"]) or {})
        graph.edge_attr.update(("edges" in _styles and _styles["edges"]) or {})
        return graph

    # find input nodes and edges
    type_of_flow = model_object.transition_flows[model_object.transition_flows.implement == strata]

    # find compartment names to be used, from all compartments listed as origins or destinations in transition flows
    new_labels = list(set().union(type_of_flow["origin"].values, type_of_flow["to"].values))

    # start building graph
    model_object.flow_diagram = Digraph(format="png")

    # inputs are sectioned according to the stem value so colours can be added to each type
    for label in new_labels:
        comp_name = find_name_components(label)[0]
        node_color = default_colour_dict[comp_name] if comp_name in default_colour_dict.keys() else "#F0FFFF"
        model_object.flow_diagram.node(label, fillcolor=node_color)

    # build the graph edges
    for row in type_of_flow.iterrows():
        model_object.flow_diagram.edge(row[1]["origin"], row[1]["to"], row[1]["parameter"])
    model_object.flow_diagram = apply_styles(model_object.flow_diagram, styles)
    model_object.flow_diagram.render(name)


def store_database(outputs, table_name="outputs",  database_name="../databases/outputs.db", append=True):
    """
    store outputs from the model in sql database for use in producing outputs later
    """
    engine = create_engine("sqlite:///" + database_name, echo=False)
    outputs.to_sql(table_name, con=engine, if_exists="append",  index=False)


class EpiModel:
    """
    general epidemiological model for constructing compartment-based models, typically of infectious disease
        transmission
    see README.md file for full description of purpose and approach of this class

    :attribute all_stratifications: dict
        will remain an empty dictionary in the unstratified version
        needed for stratified model only
    :attribute birth_approach: str
        approach to allowing entry flows into the model
        currently must be add_crude_birth_rate, replace_deaths or no_births
    :attribute compartment_names: list
        list of the strings representing the model compartments
    :attribute compartment_types: list
        copy of compartment_names, which is kept unchanged when stratification begins to increase the number of
            compartments
    :attribute compartment_values: list
        list of floats for the working values of the compartment sizes
    :attribute customised_flow_functions: dict
        user-defined functions that calculate specific model quantities that are needed to determine the rate of
            specific flows - for example, transitions that need to be implemented as absolute rates regardless of the
            size of the origin compartment
    :attribute death_flows: pandas data frame
        grid containing the information for the compartment-specific death flows to be implemented
        columns are type, parameter, origin, implement
    :attribute death_indices_to_implement: list
        indices of the death indices to be implemented because applicable to the final level of stratification
    :attribute derived_output_functions: dict
        functions that can be used during the process of integration to calculate quantities emerging from the model
            that may not be as simple as that specified in output_connections
    :attribute derived_outputs: dict
        quantities whose values are to be recorded throughout integration, i.e. tracked_quantities
        collated as lists for each time step with descriptive keys equivalent to those for tracked_quantities
    :attribute entry_compartment: str
        name of the compartment that births come in to
    :attribute equilibrium_stopping_tolerance: float
        value at which relative changes in compartment size trigger stopping when equilibrium reached
    :attribute infectious_compartment: tuple
        name(s) of the infectious compartment for calculation of inter-compartmental infection flows
    :attribute infectious_denominators: float64
        in the unstratified version, the size of the total population, relevant for frequency-dependent transmission
    :attribute infectious_indices: list
        elements are booleans to indicate whether that compartment index represents an infectious compartment
    :attribute infectious_populations: float64
        working size of the infectious population
    :attribute initial_conditions: dict
        keys are compartment types, values are starting population values for each compartment
        note that not all compartment_types must be included as keys in requests
    :attribute initial_conditions_to_total: bool
        whether to sum the initial conditions up to a certain total if this value hasn't yet been reached through the
            initial_conditions argument
    :attribute integration_type: str
        integration approach for numeric solution to odes
        currently must be odeint or solveivp, but will likely be extended as this module is developed
    :attribute output_connections: dict
        keys are the names of the quantities to be tracked
        value is dict containing the origin and the destination ("to") compartments on which to base these calculations
            as well as conditional strings that should appear in the origin and destination compartment names
    :attribute outputs: numpy array
        array containing all the evaluated compartment sizes
    :attribute parameters: dict
        string keys for each parameter, with values either string to refer to a time-variant function or float
        becomes more complicated in the stratified version below
    :attribute reporting_sigfigs: int
        number of significant figures to output to when reporting progress
    :attribute requested_flows: list
        list with each element being a model flow that contains fixed key names according to the type of flow requested
        flows are dicts in standard format
    :attribute starting_compartment: str
        optional name of the compartment to add population recruitment to
    :attribute starting_population: numeric (int or float)
        value for the total starting population to be supplemented to if initial_conditions_to_total requested
    :attribute time_variants: dict
        keys parameter names, values functions with independent variable being time and returning parameter value
    :attribute times: list
        time steps at which outputs are to be evaluated
    :attribute tracked_quantities: dict
        keys tracked quantities, which are also the keys of output_connections and derived_outputs
        values are the current working value for this quantity during integration
    :attribute transition_flows: pandas data frame
        grid containing the information for the inter-compartmental transition flows to be implemented
        columns are type, parameter, origin, to, implement, strain
    :attribute transition_indices_to_implement: list
        indices of the transition indices to be implemented because applicable to the final level of stratification
    :attribute unstratified_flows:
    :attribute verbose: bool
        whether to output progress in model construction as this process proceeds
    """

    """
    general method
    """

    def output_to_user(self, comment):
        """
        output some information
        short function just to save the if statement in every call

        :param: comment: string for the comment to be displayed to the user
        """
        if self.verbose:
            print(comment)

    """
    model construction methods
    """

    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment=("infectious",), birth_approach="no_birth",
                 verbose=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 starting_compartment="", equilibrium_stopping_tolerance=1e-6, integration_type="odeint",
                 output_connections={}, derived_output_functions={}):
        """
        construction method to create a basic compartmental model
        includes checking that the arguments have been provided correctly by the user
        at this stage the model remains unstratified, but has characteristics required to support the stratification
            process

        :params: all arguments essentially become object attributes and are described in the first main docstring to
            this object class above
        """

        # set flow attributes as pandas data frames with fixed column names
        self.transition_flows = \
            pd.DataFrame(columns=("type", "parameter", "origin", "to", "implement", "strain", "force_index"))
        self.death_flows = pd.DataFrame(columns=("type", "parameter", "origin", "implement"))

        # attributes with specific format that are independent of user inputs
        self.tracked_quantities, self.time_variants, self.all_stratifications, self.customised_flow_functions,\
            self.derived_outputs_shadow = ({} for _ in range(5))
        self.derived_outputs = {"times": []}
        self.compartment_values, self.compartment_names, self.infectious_indices = ([] for _ in range(3))

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
            infectious_compartment, birth_approach, verbose, reporting_sigfigs, entry_compartment,
            starting_population, starting_compartment, equilibrium_stopping_tolerance, integration_type,
            output_connections, derived_output_functions)

        # stop ide complaining about attributes being defined outside __init__, even though they aren't
        self.times, self.compartment_types, self.initial_conditions, self.parameters, self.requested_flows, \
            self.initial_conditions_to_total, self.infectious_compartment, self.birth_approach, self.verbose, \
            self.reporting_sigfigs, self.entry_compartment, self.starting_population, \
            self.starting_compartment, self.equilibrium_stopping_tolerance, self.outputs, self.integration_type, \
            self.output_connections, self.infectious_populations, self.infectious_denominators, \
            self.derived_output_functions, self.transition_indices_to_implement, self.death_indices_to_implement = \
            (None for _ in range(22))

        # for storing derived output in db
        self.step = 0

        # convert input arguments to model attributes
        for attribute in \
                ("times", "compartment_types", "initial_conditions", "parameters", "initial_conditions_to_total",
                 "infectious_compartment", "birth_approach", "verbose", "reporting_sigfigs", "entry_compartment",
                 "starting_population", "starting_compartment", "infectious_compartment",
                 "equilibrium_stopping_tolerance", "integration_type", "output_connections",
                 "derived_output_functions"):
            setattr(self, attribute, eval(attribute))

        # keep copy of the compartment types in case the compartment names are stratified later
        self.compartment_names = copy.copy(self.compartment_types)

        # set initial conditions and implement flows
        self.set_initial_conditions(initial_conditions_to_total)

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.initialise_default_quantities()

    def check_and_report_attributes(
            self, _times, _compartment_types, _initial_conditions, _parameters, _requested_flows,
            _initial_conditions_to_total, _infectious_compartment, _birth_approach, _verbose, _reporting_sigfigs,
            _entry_compartment, _starting_population, _starting_compartment, _equilibrium_stopping_tolerance,
            _integration_type, _output_connections, _derived_output_functions):
        """
        check all input data have been requested correctly

        :parameters: all parameters have come directly from the construction (__init__) method unchanged and have been
            renamed with a preceding _ character to indicate different scope
        """

        # check variables are of the expected type
        for expected_numeric_variable in ("_reporting_sigfigs", "_starting_population"):
            if not isinstance(eval(expected_numeric_variable), int):
                raise TypeError("expected integer for %s" % expected_numeric_variable)
        for expected_float_variable in ("_equilibrium_stopping_tolerance",):
            if not isinstance(eval(expected_float_variable), float):
                raise TypeError("expected float for %s" % expected_float_variable)
        for expected_list in ("_times", "_compartment_types", "_requested_flows"):
            if not isinstance(eval(expected_list), list):
                raise TypeError("expected list for %s" % expected_list)
        for expected_tuple in ("_infectious_compartment",):
            if not isinstance(eval(expected_tuple), tuple):
                raise TypeError("expected tuple for %s" % expected_tuple)
        for expected_string in ("_birth_approach", "_entry_compartment", "_starting_compartment", "_integration_type",):
            if not isinstance(eval(expected_string), str):
                raise TypeError("expected string for %s" % expected_string)
        for expected_boolean in ("_initial_conditions_to_total", "_verbose"):
            if not isinstance(eval(expected_boolean), bool):
                raise TypeError("expected boolean for %s" % expected_boolean)
        for expected_dict in ("_derived_output_functions",):
            if not isinstance(eval(expected_dict), dict):
                raise TypeError("expected dictionary for %s" % expected_dict)

        # check some specific requirements
        if any(_infectious_compartment) not in _compartment_types:
            ValueError("infectious compartment name is not one of the listed compartment types")
        if _birth_approach not in ("add_crude_birth_rate", "replace_deaths", "no_births"):
            ValueError("requested birth approach unavailable")
        if sorted(_times) != _times:
            self.output_to_user("requested integration times are not sorted, now sorting")
            self.times = sorted(self.times)
        for output in _output_connections:
            if any(item not in ("origin", "to", "origin_condition", "to_condition") for item in _output_connections[output]):
                raise ValueError("output connections incorrect specified, need an 'origin' and possibly 'to'," 
                                 "'origin_condition and 'to_condition' keys")

        # report on characteristics of inputs
        if _verbose:
            print("integration times are from %s to %s (time units are always arbitrary)"
                  % (round(_times[0], _reporting_sigfigs), round(_times[-1], _reporting_sigfigs)))
            print("unstratified, unprocessed initial conditions requested are:")
            for compartment in _initial_conditions:
                print("\t%s: %s" % (compartment, _initial_conditions[compartment]))
            print("infectious compartment(s) are '%s'" % _infectious_compartment)
            print("birth approach is %s" % _birth_approach)

    def set_initial_conditions(self, _initial_conditions_to_total):
        """
        set starting compartment values according to user request

        :param _initial_conditions_to_total: bool
            unchanged from argument to __init__
        """

        # first set all compartments to zero
        self.compartment_values = [0.0] * len(self.compartment_names)

        # set starting values of (unstratified) compartments to requested value
        for comp in self.initial_conditions:
            if comp in self.compartment_types:
                self.compartment_values[self.compartment_names.index(comp)] = self.initial_conditions[comp]
            else:
                raise ValueError("compartment %s requested in initial conditions not found in model compartment types")

        # sum to a total value if requested
        if _initial_conditions_to_total:
            self.sum_initial_compartments_to_total()

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain user-requested value
        """
        remainder_compartment = self.find_remainder_compartment()
        remaining_population = self.find_remainder_population(remainder_compartment)
        self.output_to_user("requested that total population sum to %s" % self.starting_population)
        self.output_to_user("remaining population of %s allocated to %s compartment"
                            % (remaining_population, remainder_compartment))
        self.compartment_values[self.compartment_names.index(remainder_compartment)] = remaining_population

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total

        :return: str
            name of the compartment to assign the remaining population size to
        """
        if len(self.starting_compartment) == 0:
            self.output_to_user("no default starting compartment requested for unallocated population, " +
                                "so will be allocated to entry compartment %s" % self.entry_compartment)
            return self.entry_compartment
        elif self.starting_compartment not in self.compartment_types:
            raise ValueError("starting compartment to populate with initial values not found in available compartments")
        else:
            return self.starting_compartment

    def find_remainder_population(self, _remainder_compartment):
        """
        start by calculating the unassigned starting population as the difference between the specified compartment
            values and the requested starting population
        if the remainder compartment has been specified by the user as needing some initial values, add these back on,
            because they would have contributed to the calculation of the remainder, which they shouldn't do

        :param _remainder_compartment: str
            name of the compartment to assign the remaining population size to
        :return: remaining_population: float
            value of the population size to assign to the initial conditions remainder population compartment
        """
        remainder = self.starting_population - sum(self.compartment_values)
        if _remainder_compartment in self.initial_conditions:
            remainder += self.initial_conditions[_remainder_compartment]
        if remainder < 0.0:
            raise ValueError("total of requested compartment values is greater than the requested starting population")
        return remainder

    def implement_flows(self, _requested_flows):
        """
        add all flows to create data frames from input lists

        :param _requested_flows: dict
            unchanged from argument to __init__
        """
        for flow in _requested_flows:

            # check flow requested correctly
            if flow["parameter"] not in self.parameters and flow["parameter"] not in self.time_variants:
                raise ValueError("flow parameter not found in parameter list")
            if flow["origin"] not in self.compartment_types:
                raise ValueError("from compartment name not found in compartment types")
            if "to" in flow and flow["to"] not in self.compartment_types:
                raise ValueError("to compartment name not found in compartment types")

            # add flow to appropriate data frame
            if flow["type"] == "compartment_death":
                self.add_death_flow(flow)
            else:
                self.add_transition_flow(flow)

    def initialise_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed during integration
        """

        # universal death rate
        if "universal_death_rate" not in self.parameters:
            self.parameters["universal_death_rate"] = 0.0

        # birth approach-specific parameters
        if self.birth_approach == "add_crude_birth_rate" and "crude_birth_rate" not in self.parameters:
            self.parameters["crude_birth_rate"] = 0.0
        elif self.birth_approach == "replace_deaths":
            self.tracked_quantities["total_deaths"] = 0.0

        # for each derived quantity to be recorded, initialise derived outputs, and a tracked quantity if needed
        for output in self.output_connections:
            self.tracked_quantities[output] = 0.0
            self.derived_outputs[output] = []
            # add empty string as conditions if no condition provided
            for comp in ["origin", "to"]:
                if comp + "_condition" not in self.output_connections[output].keys():
                    self.output_connections[output][comp + "_condition"] = ""
        for output in self.derived_output_functions:
            self.derived_outputs[output] = []



        # parameters essential for later stratification, if requested
        self.parameters["entry_fractions"] = 1.0

    def add_transition_flow(self, _flow):
        """
        add a flow (row) to the data frame storing the transition flows

        :param _flow: dict
            user-submitted flow with keys that must match existing format
        """

        # implement value starts at zero for unstratified that can be progressively incremented during stratification
        _flow["implement"] = len(self.all_stratifications) if "implement" not in _flow else _flow["implement"]
        self.transition_flows = self.transition_flows.append(
            {key: value for key, value in _flow.items() if key != "function"}, ignore_index=True)

        # record the associated function if the flow being considered is a customised flow
        if _flow["type"] == "customised_flows":
            if "function" not in _flow.keys():
                raise ValueError("a customised flow requires a function to be specified in user request dictionary")
            elif not callable(_flow["function"]):
                raise ValueError("value of 'function' key must be a function")
            self.customised_flow_functions[self.transition_flows.shape[0] - 1] = _flow["function"]

    def add_death_flow(self, _flow):
        """
        same as previous method, but for compartment-specific death flows

        :param _flow: as for previous method
        """
        _flow["implement"] = len(self.all_stratifications) if "implement" not in _flow else _flow["implement"]
        self.death_flows = self.death_flows.append(_flow, ignore_index=True)

    """
    pre-integration methods
    """

    def prepare_to_run(self):
        """
        primarily for use in the stratified version when over-written
        here just find all of the compartments that are infectious and prepare some list indices to speed integration
        """
        self.infectious_indices = self.find_all_infectious_indices()
        self.transition_indices_to_implement = self.find_transition_indices_to_implement()
        self.death_indices_to_implement = self.find_death_indices_to_implement()

    def find_all_infectious_indices(self):
        """
        find all the compartment names that begin with one of the requested infectious compartments

        :return: list
            booleans for whether each compartment is infectious or not
        """
        return [find_stem(comp) in self.infectious_compartment for comp in self.compartment_names]

    def find_transition_indices_to_implement(self):
        """
        primarily for use in the stratified version when over-written
        here just returns the indices of all the transition flows, as they all need to be implemented

        :return: list
            integers for all the rows of the transition matrix
        """
        return list(range(len(self.transition_flows)))

    def find_death_indices_to_implement(self):
        """
        for over-writing in stratified version, here just returns the indices of all the transition flows, as they all
            need to be implemented

        :return: list
            integers for all the rows of the death matrix
        """
        return list(range(len(self.death_flows)))

    """
    model running methods
    """

    def run_model(self):
        """
        main function to integrate model odes, called externally in the master running script
        """
        self.output_to_user("\n-----\nnow integrating")
        self.prepare_to_run()

        # basic default integration method
        if self.integration_type == "odeint":
            def make_model_function(compartment_values, time):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes([0.0] * len(self.compartment_names), compartment_values, time)

            self.outputs = odeint(make_model_function, self.compartment_values, self.times, atol=1.e-3, rtol=1.e-3)

        # alternative integration method
        elif self.integration_type == "solve_ivp":

            # solve_ivp requires arguments to model function in the reverse order
            def make_model_function(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes([0.0] * len(self.compartment_names), compartment_values, time)

            # add a stopping condition, which was the original purpose of using this integration approach
            def set_stopping_conditions(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                return max(list(map(abs, self.apply_all_flow_types_to_odes(
                    [0.0] * len(self.compartment_names), compartment_values, time)))) - \
                       self.equilibrium_stopping_tolerance
            set_stopping_conditions.terminal = True

            # solve_ivp returns more detailed structure, with (transposed) outputs (called "y") being just one component
            self.outputs = solve_ivp(
                make_model_function,
                (self.times[0], self.times[-1]), self.compartment_values, t_eval=self.times,
                events=set_stopping_conditions)["y"].transpose()

        else:
            raise ValueError("integration approach requested not available")

        # check that all compartment values are >= 0
        if numpy.any(self.outputs < 0.):
            print("warning, compartment or compartments with negative values")

        self.output_to_user("integration complete")
        self.calculate_outputs_post_integration()

    def store_derived_outputs_to_db(self):
        """
        currently inactive code to store outputs to sql database as they are created
        """
        derived_output_df = pd.DataFrame.from_dict(self.derived_outputs)
        derived_output_df["step"] = self.step
        store_database(derived_output_df, table_name="derived_outputs")
        self.step += 1

    def apply_all_flow_types_to_odes(self, _ode_equations, _compartment_values, _time):
        """
        apply all flow types sequentially to a vector of zeros
        note that deaths must come before births in case births replace deaths

        :param _ode_equations: list
            comes in as a list of zeros with length equal to that of the number of compartments for integration
        :param _compartment_values: numpy.ndarray
            working values of the compartment sizes
        :param _time: float
            current integration time
        :return: ode equations as list
            updated ode equations in same format but with all flows implemented
        """
        _ode_equations = self.apply_transition_flows(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_compartment_death_flows(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_universal_death_flow(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_birth_rate(_ode_equations, _compartment_values, _time)
        return _ode_equations

    def apply_transition_flows(self, _ode_equations, _compartment_values, _time):
        """
        apply fixed or infection-related inter-compartmental transition flows to odes

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_flow in self.transition_indices_to_implement:
            net_flow = self.find_net_transition_flow(n_flow, _time, _compartment_values)

            # update equations
            _ode_equations = increment_list_by_index(
                _ode_equations, self.compartment_names.index(self.transition_flows.origin[n_flow]), -net_flow)
            _ode_equations = increment_list_by_index(
                _ode_equations, self.compartment_names.index(self.transition_flows.to[n_flow]), net_flow)

            # track any quantities dependent on flow rates
            if _time not in self.derived_outputs["times"]:
                self.track_derived_outputs(n_flow, net_flow)

        # add another element to the derived outputs vector
        if _time not in self.derived_outputs["times"]:
            self.extend_derived_outputs(_time)

        # return flow rates
        return _ode_equations

    def find_net_transition_flow(self, n_flow, _time, _compartment_values):
        """
        common code to finding transition flows during and after integration packaged into single function

        :param n_flow: int
            row of interest in transition flow dataframe
        :param _time: float
            time step, which may be time during integration or post-integration time point of interest
        :param _compartment_values: list
            list of current compartment sizes
        :return: float
            net transition between the two compartments being considered
        """

        # find adjusted parameter value
        parameter_value = self.get_parameter_value(self.transition_flows.parameter[n_flow], _time)

        # find from compartment and the "infectious population" (which equals one for non-infection-related flows)
        infectious_population = self.find_infectious_multiplier(n_flow)

        # find the index of the origin or from compartment
        from_compartment = self.compartment_names.index(self.transition_flows.origin[n_flow])

        # implement flows according to whether customised or standard/infection-related
        return parameter_value * self.customised_flow_functions[n_flow](self, n_flow) if \
            self.transition_flows.type[n_flow] == "customised_flows" else \
            parameter_value * _compartment_values[from_compartment] * infectious_population

    def track_derived_outputs(self, _n_flow, _net_flow):
        """
        calculate derived quantities to be tracked, which are stored as the self.derived_outputs dictionary for the
            current working time step

        :param _n_flow: int
            row number of the flow being considered in the preceding method
        :param _net_flow: float
            previously calculated magnitude of the transition flow
        """
        for output_type in self.output_connections:
            if self.output_connections[output_type]["origin"] in self.transition_flows.origin[_n_flow] \
                    and self.output_connections[output_type]["to"] in self.transition_flows.to[_n_flow]\
                    and self.output_connections[output_type]["origin_condition"] in self.transition_flows.origin[_n_flow]\
                    and self.output_connections[output_type]["to_condition"] in self.transition_flows.to[_n_flow]:
                self.tracked_quantities[output_type] += _net_flow

    def extend_derived_outputs(self, _time):
        """
        add the derived quantities being tracked to the end of the tracking vector, taking the self.derived_outputs
        dictionary for a single time point and updating the derived outputs dictionary of lists for all time points

        :param _time: float
            current time in integration process
        """
        self.derived_outputs["times"].append(_time)
        for output_type in self.output_connections:
            self.derived_outputs[output_type].append(self.tracked_quantities[output_type])

        for output_type in self.derived_output_functions:
            self.derived_outputs[output_type].append(self.derived_output_functions[output_type](self))

    def apply_compartment_death_flows(self, _ode_equations, _compartment_values, _time):
        """
        equivalent method to for transition flows above, but for deaths

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_flow in self.death_indices_to_implement:
            parameter_value = self.get_parameter_value(self.death_flows.parameter[n_flow], _time)
            from_compartment = self.compartment_names.index(self.death_flows.origin[n_flow])
            net_flow = parameter_value * _compartment_values[from_compartment]
            _ode_equations = increment_list_by_index(_ode_equations, from_compartment, -net_flow)
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return _ode_equations

    def apply_universal_death_flow(self, _ode_equations, _compartment_values, _time):
        """
        apply the population-wide death rate to all compartments

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_comp, compartment in enumerate(self.compartment_names):
            net_flow = self.get_compartment_death_rate(compartment, _time) * _compartment_values[n_comp]
            _ode_equations = increment_list_by_index(_ode_equations, n_comp, -net_flow)

            # track deaths in case births need to replace deaths
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return _ode_equations

    def get_compartment_death_rate(self, _compartment, _time):
        """
        calculate rate of non-disease-related deaths from the compartment being considered
        needs to be separated out from the previous method for later stratification

        :param _compartment: str
            name of the compartment being considered, which is ignored here
        :param _time: float
            time in model integration
        :return: float
            rate of death from the compartment of interest
        """
        return self.get_parameter_value("universal_death_rate", _time)

    def apply_birth_rate(self, _ode_equations, _compartment_values, _time):
        """
        apply a birth rate to the entry compartment

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        return increment_list_by_index(_ode_equations, self.compartment_names.index(self.entry_compartment),
                                       self.find_total_births(_compartment_values, _time))

    def find_total_births(self, _compartment_values, _time):
        """
        calculate the total births to apply dependent on the approach requested

        :param _compartment_values:
            as for preceding methods
        :param _time:
            as for preceding methods
        :return: float
            total rate of births to be implemented in the model
        """
        if self.birth_approach == "add_crude_birth_rate":
            return self.get_single_parameter_component("crude_birth_rate", _time) * sum(_compartment_values)
        elif self.birth_approach == "replace_deaths":
            return self.tracked_quantities["total_deaths"]
        else:
            return 0.0

    def find_infectious_multiplier(self, n_flow):
        """
        find the multiplier to account for the infectious population in dynamic flows
        using plural for infectious_populations because may be stratified later

        :param n_flow: int
            index for the row of the transition_flows data frame
        :return:
            the total infectious quantity, whether that be the number or proportion of infectious persons
            needs to return as one for flows that are not transmission dynamic infectiousness flows
        """
        if self.transition_flows.type[n_flow] == "infection_density":
            return self.infectious_populations
        elif self.transition_flows.type[n_flow] == "infection_frequency":
            return self.infectious_populations / self.infectious_denominators
        else:
            return 1.0

    def update_tracked_quantities(self, _compartment_values):
        """
        update quantities that emerge during model running (not pre-defined functions of time)

        :param _compartment_values:
            as for preceding methods
        """
        self.find_infectious_population(_compartment_values)
        for quantity in self.tracked_quantities:
            self.tracked_quantities[quantity] = 0.0

    def find_infectious_population(self, _compartment_values):
        """
        calculations to find the effective infectious population

        :param _compartment_values:
            as for preceding methods
        """
        self.infectious_populations = sum([_compartment_values[compartment]
                                           for compartment in convert_boolean_list_to_indices(self.infectious_indices)])
        self.infectious_denominators = sum(_compartment_values)

    def get_parameter_value(self, _parameter, _time):
        """
        primarily for use in the stratified version when over-written

        :param _parameter: str
            parameter name
        :param _time: float
            current integration time
        :return: float
            parameter value
        """
        return self.get_single_parameter_component(_parameter, _time)

    def get_single_parameter_component(self, _parameter, time):
        """
        find the value of a parameter with time-variant values trumping constant ones

        :param _parameter: str
            string for the name of the parameter of interest
        :param time: float
            model integration time (if needed)
        :return: float
            parameter value, whether constant or time variant
        """
        return self.time_variants[_parameter](time) if _parameter in self.time_variants \
            else self.parameters[_parameter]

    """
    post-integration collation of user-requested output values
    """

    def calculate_outputs_post_integration(self):
        """
        find outputs for each requested time point, rather than at the time points that the model integration steps
            occurred at, which are arbitrary and determined by the integration routine used
        """
        self.derived_outputs_shadow["times"] = self.times
        for output in self.output_connections:
            self.derived_outputs_shadow[output] = [0.0] * len(self.times)
            transition_indices = self.find_output_transition_indices(output)
            for n_time, time in enumerate(self.times):
                self.restore_past_state(time)
                for n_flow in transition_indices:
                    net_flow = self.find_net_transition_flow(n_flow, time, self.compartment_values)
                    self.derived_outputs_shadow[output][n_time] += net_flow
        for output in self.derived_output_functions:
            self.derived_outputs_shadow[output] = [0.0] * len(self.times)
            for n_time, time in enumerate(self.times):
                self.restore_past_state(time)
                self.derived_outputs_shadow[output][n_time] = self.derived_output_functions[output](self)

    def restore_past_state(self, time):
        """
        return compartment values and tracked quantities to the values current at a particular time during model
            integration from the returned outputs structure
        the model need not have been evaluated at this particular time point

        :param time: float
            time point to go back to
        """
        self.compartment_values = self.outputs[self.times.index(time)]
        self.update_tracked_quantities(self.compartment_values)

    def find_output_transition_indices(self, output):
        """
        find the transition indices that are relevant to a particular output evaluation request from the
            output_connections dictionary created from the user's request

        :param output: str
            name of the output of interest
        :return: list
            integers referencing the transition flows relevant to this output connection
        """
        return [row for row in range(len(self.transition_flows)) if
                self.transition_flows.implement[row] == len(self.all_stratifications) and
                find_stem(self.transition_flows.origin[row]) == self.output_connections[output]["origin"] and
                find_stem(self.transition_flows.to[row]) == self.output_connections[output]["to"] and
                self.output_connections[output]["origin_condition"] in self.transition_flows.origin[row] and
                self.output_connections[output]["to_condition"] in self.transition_flows.to[row]
                ]

    """
    simple output methods, although most outputs will be managed outside of this module
    """

    def get_total_compartment_size(self, compartment_tags):
        """
        find the total values of the compartments listed by the user

        :param compartment_tags: list
            list of string variables for the compartment stems of interest
        """
        indices_to_plot = [i_comp for i_comp in range(len(self.compartment_names)) if
                           find_stem(self.compartment_names[i_comp]) in compartment_tags]
        return self.outputs[:, indices_to_plot].sum(axis=1)

    def plot_compartment_size(self, compartment_tags, multiplier=1.):
        """
        plot the aggregate population of the compartments, the name of which contains all items of the list
            compartment_tags
        kept very simple for now, because output visualisation will generally be done outside of python

        :param compartment_tags: list
            string variables for the compartments to plot
        :param multiplier: float
            scalar value to multiply the compartment values by
        """
        matplotlib.pyplot.plot(self.times, multiplier * self.get_total_compartment_size(compartment_tags))
        matplotlib.pyplot.show()


class StratifiedModel(EpiModel):
    """
    stratified version of the epidemiological model that inherits from EpiModel above, which is a concrete class and
        could in theory run stratified models independently
    however, this class should make the stratification process more algorithmic, easier and more reliable

    :attribute adaptation_functions: dict
        one-stage functions for each parameter sub-component from which to build the final functions (which are in
            final_parameter_functions)
    :attribute all_stratifications: dictionary
        keys are all the stratification names implemented so far. values are the list of strata for each stratification
    :attribute available_death_rates: list
        single strata names for which population_wide mortality will be adjusted (or over-written)
    :attribute compartment_types_to_stratify: list
        the compartments that are being stratified at this round of model stratification
    :attribute final_parameter_functions: dict
        a function for every parameter to be implemented during integration, constructed recursively for stratification
    :attribute full_stratifications_list: list
        all the stratification names implemented so far that apply to all of the compartment types
    :attribute heterogeneous_mixing: bool
        whether any stratification has requested heterogeneous mixing, such that it will be implemented
    :attribute infectious_compartments: tuple
        all of the compartment stems that represent compartments with some degree of infectiousness
    :attribute infectious_indices: dict
        keys are strains being implemented with "all_strains" an additional standard key, such that models that are not
            stratified by strain will only have the key "all_strains"
        values are lists of the indices of the compartments that are infectious for that strain (or overall)
    :attribute infectious_denominators: float
        total size of the population, which effective infectious population will be divided through by in the case of
            frequency-dependent transmission
    :attribute infectious_populations: dict
        keys are strains
        values are lists with each list element representing a mixing category, so that this can be multiplied through
            by a row of the mixing matrix
    :attribute infectiousness_adjustments: dict
        user-submitted adjustments to infectiousness for the stratification currently being implemented
    :attribute infectiousness_levels: dict
        keys are any strata for any stratification for which infectiousness will be adjusted, which does not need to be
            exhaustive
        values are their relative multipliers
    :attribute infectiousness_multipliers: list
        multipliers for the relative infectiousness of each compartment attributable to stratification, regardless of
            whether they are actually infectious compartments or not and with arbitrary values which start from one and
            are then modified by the user requests
    :attribute mixing_categories: list
        the effective mixing categories, which consists of all the possible combinations of all the strata within the
            model's full stratifications that incorporate heterogeneous mixing
        contents are strings joined with the standard linking character
    :attribute mixing_denominator_indices: dict
        keys are te mixing categories
        values are lists of the indices that should be used to calculate the infectious population for that mixing
            category
    :attribute mixing_matrix: numpy array
        array formed by taking the kronecker product of all the mixing matrices provided for full stratifications for
            which heterogeneous mixing was requested
    :attribute mortality_components: dict
        keys for the name of each compartment, values the list of functions needed to recursively create the functions
            to calculate the mortality rates for each compartment
    :attribute overwrite_character: str
        standard string (usually single character and currently "W") to indicate that a stratum request is intended to
            over-write less stratified parameters
    :attribute overwrite_key: str
        standard string used by model to identify the dictionary element that represents the over-write parameters,
            rather than a request to a particular stratum
    :attribute overwrite_parameters: list
        parameters which will result in all the less stratified parameters closer to the stratification tree's trunk
            being ignored
    :attribute parameter_components: dict
        keys for the name of each transition parameter, values the list of functions needed to recursively create the
            functions to create these parameter values
    :attribute parameters: dict
        same format as for EpiModel, but described here again given the other parameter-related attributes
        unprocessed parameters, which may be either float values or strings pointing to the keys of adaptation functions
    :attribute removed_compartments: list
        all unstratified compartments that have been removed through the stratification process
    :attribute overwrite_parameters: list
        any parameters that are intended as absolute values to be applied to that stratum and not multipliers for the
            unstratified parameter further up the tree
    :attribute strain_mixing_elements: dict
        first tier of keys is strains
        second tier of keys is mixing categories
        content of lists at lowest/third tier is the indices of the compartments that are relevant to this strain and
            category
    :attribute strain_mixing_multipliers: dict
        first tier of keys is strains
        second tier of keys is mixing categories
        content of lists at lowest/third tier is the final infectiousness multiplier for the compartments for this
            strain and category
    :attribute strains: list
        the strata to the strains stratification with specific behaviour
    """

    """
    general methods
    """

    def add_compartment(self, new_compartment_name, new_compartment_value):
        """
        add a compartment by specifying its name and the starting value for it to take

        :param new_compartment_name: str
            name of the new compartment to be created
        :param new_compartment_value: float
            initial value to be assigned to the new compartment before integration
        """
        self.compartment_names.append(new_compartment_name)
        self.compartment_values.append(new_compartment_value)
        self.output_to_user("adding compartment: %s" % new_compartment_name)

    def remove_compartment(self, compartment_name):
        """
        remove a compartment by taking the element out of the compartment_names and compartment_values attributes
        store name of removed compartment in a separate attribute

        :param compartment_name: str
            name of compartment to be removed
        """
        self.removed_compartments.append(compartment_name)
        del self.compartment_values[self.compartment_names.index(compartment_name)]
        del self.compartment_names[self.compartment_names.index(compartment_name)]
        self.output_to_user("removing compartment: %s" % compartment_name)

    """
    model construction and methods
    """

    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment=("infectious",), birth_approach="no_birth",
                 verbose=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 starting_compartment="", equilibrium_stopping_tolerance=1e-6, integration_type="odeint",
                 output_connections={}, derived_output_functions={}):
        """
        constructor mostly inherits from parent class, with a few additional attributes that are required for the
        stratified version

        :parameters: all parameters coming in as arguments are those that are also attributes of the parent class
        """
        EpiModel.__init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                          initial_conditions_to_total=initial_conditions_to_total,
                          infectious_compartment=infectious_compartment, birth_approach=birth_approach,
                          verbose=verbose, reporting_sigfigs=reporting_sigfigs, entry_compartment=entry_compartment,
                          starting_population=starting_population, starting_compartment=starting_compartment,
                          equilibrium_stopping_tolerance=equilibrium_stopping_tolerance,
                          integration_type=integration_type, output_connections=output_connections,
                          derived_output_functions=derived_output_functions)

        self.full_stratification_list, self.removed_compartments, self.overwrite_parameters, \
            self.compartment_types_to_stratify, self.infectious_denominators, self.strains, self.mixing_categories = \
            ([] for _ in range(7))
        self.all_stratifications, self.infectiousness_adjustments, self.final_parameter_functions,\
            self.adaptation_functions, self.infectiousness_levels, self.infectious_indices, \
            self.infectious_compartments, self.infectiousness_multipliers, self.parameter_components, \
            self.mortality_components, self.infectious_populations, self.strain_mixing_elements, \
            self.strain_mixing_multipliers = ({} for _ in range(13))
        self.overwrite_character, self.overwrite_key = "W", "overwrite"
        self.heterogeneous_mixing, self.mixing_matrix, self.available_death_rates = False, None, [""]

    """
    stratification methods
    """

    def stratify(
            self, stratification_name, strata_request, compartment_types_to_stratify, requested_proportions,
            entry_proportions={}, adjustment_requests=(), infectiousness_adjustments={}, mixing_matrix=None,
            verbose=True):
        """
        calls to initial preparation, checks and methods that stratify the various aspects of the model

        :param stratification_name:
            see prepare_and_check_stratification
        :param strata_request:
            see find_strata_names_from_input
        :param compartment_types_to_stratify:
            see check_compartment_request
        :param adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :param requested_proportions:
            see prepare_starting_proportions
        :param entry_proportions:

        :param infectiousness_adjustments:

        :param mixing_matrix:
            see check_mixing
        :param verbose: bool
            whether to report on progress
            note that this can be changed at this stage from what was requested at the original unstratified model
                construction
        """

        # check inputs correctly specified
        strata_names, adjustment_requests = self.prepare_and_check_stratification(
            stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, verbose)

        # work out ageing flows - comes first, so that the compartment names remain in the unstratified form
        if stratification_name == "age":
            self.set_ageing_rates(strata_names)

        # stratify the compartments
        requested_proportions = self.prepare_starting_proportions(strata_names, requested_proportions)
        self.stratify_compartments(stratification_name, strata_names, requested_proportions)

        # stratify the flows
        self.stratify_transition_flows(stratification_name, strata_names, adjustment_requests)
        self.stratify_entry_flows(stratification_name, strata_names, entry_proportions, requested_proportions)
        if self.death_flows.shape[0] > 0:
            self.stratify_death_flows(stratification_name, strata_names, adjustment_requests)
        self.stratify_universal_death_rate(
            stratification_name, strata_names, adjustment_requests, compartment_types_to_stratify)

        # if stratifying by strain
        self.strains = strata_names if stratification_name == "strain" else self.strains

        # check submitted mixing matrix and combine with existing matrix, if any
        self.prepare_mixing_matrix(mixing_matrix, stratification_name, strata_names)

        # prepare infectiousness levels attribute
        self.prepare_infectiousness_levels(stratification_name, strata_names, infectiousness_adjustments)

    """
    stratification checking methods
    """

    def prepare_and_check_stratification(
            self, _stratification_name, _strata_names, _compartment_types_to_stratify, _adjustment_requests, _verbose):
        """
        initial preparation and checks of user-submitted arguments

        :param _stratification_name: str
            the name of the stratification - i.e. the reason for implementing this type of stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _compartment_types_to_stratify:
            see check_compartment_request
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :param _verbose:
            see stratify
        :return:
            _strata_names: list
                revised version of user request after adaptation to class requirements
            adjustment_requests:
                revised version of _adjustment_requests after adaptation to class requirements
        """

        # collate all the stratifications that have been implemented so far
        if not _compartment_types_to_stratify:
            self.full_stratification_list.append(_stratification_name)

        # report progress
        self.verbose = _verbose
        self.output_to_user("\n___________________\nimplementing stratification for: %s" % _stratification_name)

        # deal with stratifications that have specific behaviour
        if _stratification_name == "age":
            _strata_names = self.check_age_stratification(_strata_names, _compartment_types_to_stratify)
        elif _stratification_name == "strain":
            self.output_to_user("implementing strain stratification with specific behaviour")

        # make sure the stratification name is a string
        if not isinstance(_stratification_name, str):
            _stratification_name = str(_stratification_name)
            self.output_to_user("converting stratification name %s to string" % _stratification_name)

        # ensure requested stratification hasn't previously been implemented
        if _stratification_name in self.all_stratifications.keys():
            raise ValueError("requested stratification has already been implemented, please choose a different name")

        # record stratification as model attribute, find the names to apply strata and check requests
        _strata_names = self.find_strata_names_from_input(_strata_names)
        self.all_stratifications[_stratification_name] = _strata_names
        _adjustment_requests = self.incorporate_alternative_overwrite_approach(_adjustment_requests)
        self.check_compartment_request(_compartment_types_to_stratify)
        self.check_parameter_adjustment_requests(_adjustment_requests, _strata_names)
        return _strata_names, _adjustment_requests

    def check_age_stratification(self, _strata_names, _compartment_types_to_stratify):
        """
        check that the user request meets the requirements for stratification by age

        :parameters: all parameters have come directly from the stratification (stratify) method unchanged and have been
            renamed with a preceding _ character
        :return: _strata_names: list
            revised names of the strata tiers to be implemented
        """
        self.output_to_user("implementing age stratification with specific behaviour")
        if len(_compartment_types_to_stratify) > 0:
            raise ValueError("requested age stratification, but compartment request should be passed as empty vector " +
                             "in order to apply to all compartments")
        elif not all([isinstance(stratum, (int, float)) for stratum in _strata_names]):
            raise ValueError("inputs for age strata breakpoints are not numeric")
        if 0 not in _strata_names:
            _strata_names.append(0)
            self.output_to_user("adding age stratum called '0' because not requested, which represents those aged " +
                                "less than %s" % min(_strata_names))
        if _strata_names != sorted(_strata_names):
            _strata_names = sorted(_strata_names)
            self.output_to_user("requested age strata not ordered, so have been sorted to: %s" % _strata_names)
        return _strata_names

    def find_strata_names_from_input(self, _strata_names):
        """
        find the names of the strata to be implemented from a particular user request

        :parameters: list or alternative format to be adapted
            strata requested in the format provided by the user (except for age, which is dealth with in the preceding
            method)
        :return: strata_names: list
            modified list of strata to be implemented in model
        """
        if type(_strata_names) == int:
            _strata_names = numpy.arange(1, _strata_names + 1)
            self.output_to_user("single integer provided as strata labels for stratification, hence strata " +
                                "implemented will be integers from one to %s" % _strata_names)
        elif type(_strata_names) == float:
            raise ValueError("single value passed as request for strata labels, but not an integer greater than " +
                             "one, so unclear what to do - stratification failed")
        elif type(_strata_names) == list and len(_strata_names) > 0:
            pass
        else:
            raise ValueError("requested to stratify, but strata-level names not submitted in correct format")
        for name in range(len(_strata_names)):
            _strata_names[name] = str(_strata_names[name])
            self.output_to_user("adding stratum: %s" % _strata_names[name])
        return _strata_names

    def incorporate_alternative_overwrite_approach(self, _adjustment_requests):
        """
        alternative approach to working out which parameters to overwrite
        can put a capital W at the string's end to indicate that it is an overwrite parameter, as an alternative to
        submitting a separate dictionary key to represent the strata which need to be overwritten

        :param _adjustment_requests: dict
            user-submitted version of adjustment requests
        :return: revised_adjustments: dict
            modified version of _adjustment_requests after working out whether any parameters began with W
        """

        # has to be constructed as a separate dictionary to avoid change of size during iteration
        revised_adjustments = {}
        for parameter in _adjustment_requests:
            revised_adjustments[parameter] = {}

            # ignore overwrite if submitted with the standard approach
            for stratum in _adjustment_requests[parameter]:
                if stratum == "overwrite":
                    continue

                # if the parameter ends in W, interpret as an overwrite parameter and added to this key
                elif stratum[-1] == self.overwrite_character:
                    if "overwrite" not in revised_adjustments[parameter]:
                        revised_adjustments[parameter]["overwrite"] = []
                    revised_adjustments[parameter][stratum[: -1]] = _adjustment_requests[parameter][stratum]
                    revised_adjustments[parameter]["overwrite"].append(stratum[: -1])

                # otherwise just accept the parameter in its submitted form
                else:
                    revised_adjustments[parameter][stratum] = _adjustment_requests[parameter][stratum]
            if "overwrite" not in revised_adjustments:
                revised_adjustments["overwrite"] = []
        return revised_adjustments

    def check_compartment_request(self, _compartment_types_to_stratify):
        """
        check the requested compartments to be stratified has been requested correctly

        :param _compartment_types_to_stratify: list
            the names of the compartment types that the requested stratification is intended to apply to
        """

        # if list of length zero passed, stratify all the compartment types in the model
        if len(_compartment_types_to_stratify) == 0:
            self.compartment_types_to_stratify = self.compartment_types
            self.output_to_user("no compartment names specified for this stratification, " +
                                "so stratification applied to all model compartments")

        # otherwise check all the requested compartments are available and implement the user request
        elif any([compartment not in self.compartment_types for compartment in self.compartment_types_to_stratify]):
            raise ValueError("requested compartment or compartments to be stratified are not available in this model")
        else:
            self.compartment_types_to_stratify = _compartment_types_to_stratify

    def check_parameter_adjustment_requests(self, _adjustment_requests, _strata_names):
        """
        check parameter adjustments have been requested appropriately and add parameter for any strata not referred to

        :param _adjustment_requests: dict
            version of the submitted adjustment_requests modified by incorporate_alternative_overwrite_approach
        :param _strata_names:
            see find_strata_names_from_input
        """
        for parameter in _adjustment_requests:
            if any(requested_stratum not in _strata_names + [self.overwrite_key]
                   for requested_stratum in _adjustment_requests[parameter]):
                raise ValueError("a stratum was requested in adjustments that is not available in this stratification")

    """
    stratification preparation methods
    """

    def set_ageing_rates(self, _strata_names):
        """
        set inter-compartmental flows for ageing from one stratum to the next as the reciprocal of the width of the age
            bracket

        :param _strata_names:
            see find_strata_names_from_input
        """
        for stratum_number in range(len(_strata_names[: -1])):
            start_age = int(_strata_names[stratum_number])
            end_age = int(_strata_names[stratum_number + 1])
            ageing_parameter_name = "ageing%sto%s" % (start_age, end_age)
            ageing_rate = 1.0 / (end_age - start_age)
            self.output_to_user("ageing rate from age group %s to %s is %s"
                                % (start_age, end_age, round(ageing_rate, self.reporting_sigfigs)))
            self.parameters[ageing_parameter_name] = ageing_rate
            for compartment in self.compartment_names:
                self.transition_flows = self.transition_flows.append(
                    {"type": "standard_flows",
                     "parameter": ageing_parameter_name,
                     "origin": create_stratified_name(compartment, "age", start_age),
                     "to": create_stratified_name(compartment, "age", end_age),
                     "implement": len(self.all_stratifications)},
                    ignore_index=True)

    def prepare_starting_proportions(self, _strata_names, _requested_proportions):
        """
        prepare user inputs for starting proportions for the initial conditions to apply to the exact set of strata
            requested
        if one or more strata not specified, the proportion of the initial conditions allocated to that group will be
            the total unallocated population divided by the number of strata for which no request was specified

        :param _strata_names:
            see find_strata_names_from_input
        :param _requested_proportions: dict
            dictionary with keys for the stratum to assign starting population to and values the proportions to assign
        :return: dict
            revised dictionary of starting proportions after cleaning
        """
        self.output_to_user(
            "\n-----\ncalculating proportions of initial conditions to assign to each stratified starting compartment")
        if any(stratum not in _strata_names for stratum in _requested_proportions):
            raise ValueError("requested starting proportion for stratum that does not appear in requested strata")
        if sum(_requested_proportions.values()) > 1.0:
            raise ValueError("requested starting proportions sum to a value greater than one")

        # assuming an equal proportion of the unallocated population if no request specified
        unrequested_strata = [stratum for stratum in _strata_names if stratum not in _requested_proportions]
        unrequested_proportions = {}
        for stratum in unrequested_strata:
            starting_proportion = (1.0 - sum(_requested_proportions.values())) / len(unrequested_strata)
            unrequested_proportions[stratum] = starting_proportion
            self.output_to_user(
                "no starting proportion requested for %s stratum so provisionally allocated %s of total"
                % (stratum, round(starting_proportion, self.reporting_sigfigs)))

        # update specified proportions with inferred unspecified proportions
        _requested_proportions.update(unrequested_proportions)
        return _requested_proportions

    def stratify_compartments(self, _stratification_name, _strata_names, _requested_proportions):
        """
        stratify the model compartments, which affects the compartment_names and the compartment_values attributes

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _requested_proportions:
            see prepare_starting_proportions
        """

        # find the existing compartments that need stratification
        self.output_to_user("\n-----\ndetermining which compartments to add and which to remove")
        for compartment in \
                [comp for comp in self.compartment_names if find_stem(comp) in self.compartment_types_to_stratify]:

            # add and remove compartments
            for stratum in _strata_names:
                self.add_compartment(create_stratified_name(compartment, _stratification_name, stratum),
                                     self.compartment_values[self.compartment_names.index(compartment)] *
                                     _requested_proportions[stratum])
            self.remove_compartment(compartment)

    def stratify_transition_flows(self, _stratification_name, _strata_names, _adjustment_requests):
        """
        stratify flows depending on whether inflow, outflow or both need replication

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        """
        self.output_to_user("\n-----\nstratifying transition flows and calculating associated parameters")
        for n_flow in self.find_transition_indices_to_implement(back_one=1):
            self.add_stratified_flows(
                n_flow, _stratification_name, _strata_names,
                find_stem(self.transition_flows.origin[n_flow]) in self.compartment_types_to_stratify,
                find_stem(self.transition_flows.to[n_flow]) in self.compartment_types_to_stratify,
                _adjustment_requests)
        self.output_to_user("\n-----\nstratified transition flows matrix\n%s" % self.transition_flows)

    def add_stratified_flows(
            self, _n_flow, _stratification_name, _strata_names, stratify_from, stratify_to, _adjustment_requests):
        """
        add additional stratified flow to the transition flow data frame attribute of the class

        :param _n_flow: int
            row location of the unstratified flow within the transition flow attribute
        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param stratify_from: bool
            whether to stratify the from/origin compartment
        :param stratify_to:
            whether to stratify the to/destination compartment
        :param _adjustment_requests:
            see incorporate _alternative_overwrite_approach and check_parameter_adjustment_requests
        """
        if stratify_from or stratify_to:
            self.output_to_user(
                "for flow from %s to %s in stratification %s"
                % (self.transition_flows.origin[_n_flow], self.transition_flows.to[_n_flow], _stratification_name))

            # loop over each stratum in the requested stratification structure
            for stratum in _strata_names:

                # find parameter name
                parameter_name = self.add_adjusted_parameter(
                    self.transition_flows.parameter[_n_flow], _stratification_name, stratum, _adjustment_requests)
                if not parameter_name:
                    parameter_name = self.sort_absent_transition_parameter(
                        _stratification_name, _strata_names, stratum, stratify_from, stratify_to,
                        self.transition_flows.parameter[_n_flow])
                self.output_to_user("\t\tadding parameter %s" % parameter_name)

                # determine whether to and/or from compartments are stratified
                from_compartment = \
                    create_stratified_name(self.transition_flows.origin[_n_flow], _stratification_name, stratum) if \
                    stratify_from else self.transition_flows.origin[_n_flow]
                to_compartment = \
                    create_stratified_name(self.transition_flows.to[_n_flow], _stratification_name, stratum) if \
                    stratify_to else self.transition_flows.to[_n_flow]

                # add the new flow
                strain = stratum if _stratification_name == "strain" else self.transition_flows.strain[_n_flow]
                self.transition_flows = self.transition_flows.append(
                    {"type": self.transition_flows.type[_n_flow],
                     "parameter": parameter_name,
                     "origin": from_compartment,
                     "to": to_compartment,
                     "implement": len(self.all_stratifications),
                     "strain": strain},
                    ignore_index=True)

                # update the customised flow function storage dictionary
                if self.transition_flows.type[_n_flow] == 'customised_flows':
                    self.update_customised_flow_function_dict(_n_flow)

        # if flow applies to a transition not involved in the stratification, still increment to ensure implemented
        else:
            new_flow = self.transition_flows.loc[_n_flow, :].to_dict()
            new_flow["implement"] += 1
            self.transition_flows = self.transition_flows.append(new_flow, ignore_index=True)
            if self.transition_flows.type[_n_flow] == 'customised_flows':
                self.update_customised_flow_function_dict(_n_flow)

    def add_adjusted_parameter(self, _unadjusted_parameter, _stratification_name, _stratum, _adjustment_requests):
        """
        find the adjustment request that is relevant to a particular unadjusted parameter and stratum and add the
            parameter value (str for function or float) to the parameters dictionary attribute
        otherwise allow return of None

        :param _unadjusted_parameter:
            name of the unadjusted parameter value
        :param _stratification_name:
            see prepare_and_check_stratification
        :param _stratum:
            stratum being considered by the method calling this method
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :return: parameter_adjustment_name: str or None
            if returned as None, assumption will be that the original, unstratified parameter should be used
            otherwise create a new parameter name and value and store away in the appropriate model structure
        """
        parameter_adjustment_name = None
        relevant_adjustment_request = self.find_relevant_adjustment_request(_adjustment_requests, _unadjusted_parameter)
        if relevant_adjustment_request is not None:
            parameter_adjustment_name = \
                create_stratified_name(_unadjusted_parameter, _stratification_name, _stratum) if \
                _stratum in _adjustment_requests[relevant_adjustment_request] else _unadjusted_parameter
            self.output_to_user("\t parameter for %s stratum of %s stratification is called %s"
                                % (_stratum, _stratification_name, parameter_adjustment_name))
            if _stratum in _adjustment_requests[relevant_adjustment_request]:
                self.parameters[parameter_adjustment_name] = _adjustment_requests[relevant_adjustment_request][_stratum]

            # record the parameters that over-write the less stratified parameters closer to the trunk of the tree
            if self.overwrite_key in _adjustment_requests[relevant_adjustment_request] and \
                    _stratum in _adjustment_requests[relevant_adjustment_request][self.overwrite_key]:
                self.overwrite_parameters.append(parameter_adjustment_name)
        return parameter_adjustment_name

    def find_relevant_adjustment_request(self, _adjustment_requests, _unadjusted_parameter):
        """
        find the adjustment requests that are extensions of the base parameter type being considered
        expected behaviour is as follows:
        * if there are no submitted requests (keys to the adjustment requests) that are extensions of the unadjusted
            parameter, will return None
        * if there is one submitted request that is an extension of the unadjusted parameter, will return that parameter
        * if there are multiple submitted requests that are extensions to the unadjusted parameter and one is more
            stratified than any of the others (i.e. more instances of the "X" string), will return this most stratified
            parameter
        * if there are multiple submitted requests that are extensions to the unadjusted parameter and several of them
            are equal in having the greatest extent of stratification, will return the first one with the greatest
            length in the order of looping through the keys of the request dictionary

        :param _unadjusted_parameter:
            see add_adjusted_parameter
        :param _adjustment_requests:
            see prepare_and_check_stratification
        :return: str or None
            the key of the adjustment request that is applicable to the parameter of interest if any, otherwise None
        """

        # find all the requests that start with the parameter of interest and their level of stratification
        applicable_params = [param for param in _adjustment_requests if _unadjusted_parameter.startswith(param)]
        applicable_param_lengths = [len(find_name_components(param)) for param in applicable_params]

        # find the first most stratified parameter
        return applicable_params[applicable_param_lengths.index(max(applicable_param_lengths))] \
            if applicable_param_lengths else None

    def sort_absent_transition_parameter(
            self, _stratification_name, _strata_names, _stratum, _stratify_from, _stratify_to, unstratified_name):
        """
        work out what to do if a specific transition parameter adjustment has not been requested

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _stratum:
        :param _stratify_from:
            see add_stratified_flows
        :param _stratify_to:
            see add_stratified_flows
        :param unstratified_name: str
            the name of the parameter before the stratification is implemented
        :return: str
            parameter name for revised parameter than wasn't provided
        """

        # default behaviour if not specified is to split the parameter into equal parts if to compartment is split
        if not _stratify_from and _stratify_to:
            self.output_to_user("\t splitting existing parameter value %s into %s equal parts"
                                % (unstratified_name, len(_strata_names)))
            parameter_name = create_stratified_name(unstratified_name, _stratification_name, _stratum)
            self.parameters[parameter_name] = 1.0 / len(_strata_names)
            self.adaptation_functions[parameter_name] = create_multiplicative_function(1.0 / len(_strata_names))
            return parameter_name

        # otherwise if no request, retain the existing parameter
        else:
            self.output_to_user("\tretaining existing parameter value %s" % unstratified_name)
            return unstratified_name

    def update_customised_flow_function_dict(self, _n_flow):
        """
        when a stratified flow is being created and if the original unstratified flow was customised, we need to update
            the dictionary listing the functions associated with the customised flows with a new key to refer to the
            same function

        :param _n_flow: int
            the index of the unstratified flow
        """
        self.customised_flow_functions[self.transition_flows.shape[0] - 1] = self.customised_flow_functions[_n_flow]

    def stratify_entry_flows(self, _stratification_name, _strata_names, _entry_proportions, _requested_proportions):
        """
        stratify entry/recruitment/birth flows according to requested entry proportion adjustments
        again, may need to revise behaviour for what is done if some strata are requested but not others

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _entry_proportions: dict
            user requested proportions to enter to each stratum
        :param _requested_proportions:
            see prepare_starting_proportions
        :return:
            normalised dictionary of the compartments that the new entry flows should come in to
        """
        if self.entry_compartment in self.compartment_types_to_stratify:
            self.output_to_user(
                "\n-----\ncalculating proportions of births/recruitment to assign to each stratified entry compartment")
            for stratum in _strata_names:
                entry_fraction_name = create_stratified_name("entry_fraction", _stratification_name, stratum)

                # specific behaviour for age stratification
                if _stratification_name == "age" and str(stratum) == "0":
                    self.parameters[entry_fraction_name] = 1.0
                    continue
                elif _stratification_name == "age":
                    self.parameters[entry_fraction_name] = 0.0
                    continue

                # where a request for splitting entry rates has been submitted
                elif stratum in _entry_proportions and type(_entry_proportions[stratum]) == float:
                    self.parameters[entry_fraction_name] = _entry_proportions[stratum]
                    self.output_to_user("assigning requested proportion %s of births/recruitment to %s stratum"
                                        % (_entry_proportions[stratum], stratum))

                # if an incorrect string has been submitted by the user
                elif stratum in _entry_proportions and type(_entry_proportions[stratum]) == str and \
                        _entry_proportions[stratum] not in self.time_variants:
                    raise ValueError("requested entry fraction function for %s stratum not available in time variants")

                # otherwise it must already be a defined function that can be called during integration
                elif stratum in _entry_proportions and type(_entry_proportions[stratum]) == str:
                    self.time_variants[entry_fraction_name] = self.time_variants[_entry_proportions[stratum]]
                    self.output_to_user("function %s submitted for proportion of births assigned to %s"
                                        % (_entry_proportions[stratum], stratum))
                    continue

                # otherwise if no request made
                else:
                    self.parameters[entry_fraction_name] = 1.0 / len(_strata_names)

    def stratify_death_flows(self, _stratification_name, _strata_names, _adjustment_requests):
        """
        add compartment-specific death flows to death_flows data frame attribute

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
             see find_strata_names_from_input
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        """
        for n_flow in self.find_death_indices_to_implement(back_one=1):

            # if the compartment with an additional death flow is being stratified
            if find_stem(self.death_flows.origin[n_flow]) in self.compartment_types_to_stratify:
                for stratum in _strata_names:

                    # get stratified parameter name if requested to stratify, otherwise use the unstratified one
                    parameter_name = self.add_adjusted_parameter(
                        self.death_flows.parameter[n_flow], _stratification_name, stratum, _adjustment_requests)
                    if not parameter_name:
                        parameter_name = self.death_flows.parameter[n_flow]

                    # add the stratified flow to the death flows data frame
                    self.death_flows = self.death_flows.append(
                        {"type": self.death_flows.type[n_flow],
                         "parameter": parameter_name,
                         "origin": create_stratified_name(self.death_flows.origin[n_flow], _stratification_name, stratum),
                         "implement": len(self.all_stratifications)},
                        ignore_index=True)

            # otherwise if not part of the stratification, accept the existing flow and increment the implement value
            else:
                new_flow = self.death_flows.loc[n_flow, :].to_dict()
                new_flow["implement"] += 1
                self.death_flows = self.death_flows.append(new_flow, ignore_index=True)

    def stratify_universal_death_rate(
            self, _stratification_name, _strata_names, _adjustment_requests, _compartment_types_to_stratify):
        """
        stratify the approach to universal, population-wide deaths (which can be made to vary by stratum)
        adjust every parameter that refers to the universal death rate, according to user request if submitted and
            otherwise populated with a value of one by default

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :param _compartment_types_to_stratify:
            see above
        """
        if _stratification_name not in self.full_stratification_list and "universal_death_rate" in _adjustment_requests:
            raise ValueError("universal death rate can only be stratified when applied to all compartment types")
        elif _stratification_name not in self.full_stratification_list:
            self.output_to_user("universal death rate not adjusted as stratification not applied to all compartments")
            return

        # ensure baseline function available for modification in universal death rates
        self.adaptation_functions["universal_death_rateX"] = self.time_variants["universal_death_rate"] \
            if "universal_death_rate" in self.time_variants else lambda time: self.parameters["universal_death_rate"]

        # if stratification applied to all compartment types
        for stratum in _strata_names:
            if "universal_death_rate" in _adjustment_requests and \
                    stratum in _adjustment_requests["universal_death_rate"]:
                stratum_name = create_stratum_name(_stratification_name, stratum, joining_string="")
                self.available_death_rates.append(stratum_name)

                # use existing function or create new one from constant as needed
                if type(_adjustment_requests["universal_death_rate"][stratum]) == str:
                    self.adaptation_functions["universal_death_rateX" + stratum_name] = \
                        self.time_variants[_adjustment_requests["universal_death_rate"][stratum]]
                elif isinstance(_adjustment_requests["universal_death_rate"][stratum], (int, float)):
                    self.adaptation_functions["universal_death_rateX" + stratum_name] = \
                        create_multiplicative_function(
                            self.time_variants[_adjustment_requests["universal_death_rate"][stratum]])

                # record the parameters that over-write the less stratified parameters closer to the trunk of the tree
                if self.overwrite_key in _adjustment_requests["universal_death_rate"] and \
                        stratum in _adjustment_requests["universal_death_rate"][self.overwrite_key]:
                    self.overwrite_parameters.append(
                        create_stratified_name("universal_death_rate", _stratification_name, stratum))

    def prepare_mixing_matrix(self, _mixing_matrix, _stratification_name, _strata_names):
        """
        check that the mixing matrix has been correctly specified and call the other relevant functions

        :param _mixing_matrix: numpy array
            must be square
            represents the mixing of the strata within this stratification
        :param _stratification_name: str
            the name of the stratification - i.e. the reason for implementing this type of stratification
        :param _strata_names: list
            see find_strata_names_from_input
        """
        if _mixing_matrix is None:
            return
        elif type(_mixing_matrix) != numpy.ndarray:
            raise ValueError("submitted mixing matrix is wrong data type")
        elif len(_mixing_matrix.shape) != 2:
            raise ValueError("submitted mixing matrix is not two-dimensional")
        elif _mixing_matrix.shape[0] != _mixing_matrix.shape[1]:
            raise ValueError("submitted mixing is not square")
        elif _mixing_matrix.shape[0] != len(_strata_names):
            raise ValueError("mixing matrix does not sized to number of strata being implemented")
        self.combine_new_mixing_matrix_with_existing(_mixing_matrix, _stratification_name, _strata_names)

    def combine_new_mixing_matrix_with_existing(self, _mixing_matrix, _stratification_name, _strata_names):
        """
        master mixing matrix function to take in a new mixing matrix and combine with the existing ones

        :param _mixing_matrix: numpy array
            array, which must be square representing the mixing of the strata within this stratification
        :param _stratification_name: str
            the name of the stratification - i.e. the reason for implementing this type of stratification
        :param _strata_names: list
            see find_strata_names_from_input
        """

        # if no mixing matrix yet, just convert the existing one to a dataframe
        if self.mixing_matrix is None:
            self.mixing_categories = [_stratification_name + "_" + i for i in _strata_names]
            self.mixing_matrix = _mixing_matrix

        # otherwise take the kronecker product to get the new mixing matrix
        else:
            self.mixing_categories = \
                [old_strata + "X" + _stratification_name + "_" + new_strata
                 for old_strata, new_strata in itertools.product(self.mixing_categories, _strata_names)]
            self.mixing_matrix = numpy.kron(self.mixing_matrix, _mixing_matrix)

    def prepare_infectiousness_levels(self, _stratification_name, _strata_names, _infectiousness_adjustments):
        """
        store infectiousness adjustments as dictionary attribute to the model object, with first tier of keys the
            stratification and second tier the strata to be modified

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
             see find_strata_names_from_input
        :param _infectiousness_adjustments: dict
            requested adjustments to infectiousness for this stratification
        """
        if type(_infectiousness_adjustments) != dict:
            raise ValueError("infectiousness adjustments not submitted as dictionary")
        elif not all(key in _strata_names for key in _infectiousness_adjustments.keys()):
            raise ValueError("infectiousness adjustment key not in strata being implemented")
        else:
            for stratum in _infectiousness_adjustments:
                self.infectiousness_levels[create_stratum_name(_stratification_name, stratum, joining_string="")] = \
                    _infectiousness_adjustments[stratum]

    """
    pre-integration methods
    """

    def prepare_to_run(self):
        """
        methods that can be run prior to integration to save various function calls being made at every time step
        """
        self.prepare_stratified_parameter_calculations()
        self.prepare_infectiousness_calculations()
        self.transition_indices_to_implement = self.find_transition_indices_to_implement()
        self.death_indices_to_implement = self.find_death_indices_to_implement()

        # ensure there is a universal death rate available even if the model hasn't been stratified at all
        if len(self.all_stratifications) == 0 and isinstance(self.parameters["universal_death_rate"], (float, int)):
            self.final_parameter_functions["universal_death_rate"] = \
                lambda time: self.parameters["universal_death_rate"]
        elif len(self.all_stratifications) == 0 and type(self.parameters["universal_death_rate"]) == str:
            self.final_parameter_functions["universal_death_rate"] = self.adaptation_functions["universal_death_rate"]

    def prepare_stratified_parameter_calculations(self):
        """
        prior to integration commencing, work out what the components are of each parameter being implemented
        populates self.parameter_components even though it is not needed elsewhere, to allow that the components that
            were used to create each given parameter can be determined later
        """

        # create list of all the parameters that we need to find the set of adjustment functions for
        parameters_to_adjust = []
        for n_flow in range(self.transition_flows.shape[0]):
            if self.transition_flows.implement[n_flow] == len(self.all_stratifications) and \
                    self.transition_flows.parameter[n_flow] not in parameters_to_adjust:
                parameters_to_adjust.append(self.transition_flows.parameter[n_flow])
        for n_flow in range(self.death_flows.shape[0]):
            if self.death_flows.implement[n_flow] == len(self.all_stratifications) and \
                    self.death_flows.parameter[n_flow] not in parameters_to_adjust:
                parameters_to_adjust.append(self.death_flows.parameter[n_flow])

        # and adjust
        for parameter in parameters_to_adjust:
            self.parameter_components[parameter] = self.find_transition_components(parameter)
            self.create_transition_functions(parameter, self.parameter_components[parameter])

        # similarly for all model compartments
        for compartment in self.compartment_names:
            self.mortality_components[compartment] = self.find_mortality_components(compartment)
            if len(self.all_stratifications) > 0:
                self.create_mortality_functions(compartment, self.mortality_components[compartment])

    def find_mortality_components(self, _compartment):
        """
        find the sub-parameters for population-wide natural mortality that are relevant to a particular compartment
        used in prepare_stratified_parameter_calculations for creating functions to find the mortality rate for each
            compartment
        similar to find_transition_components, except being applied by compartment rather than parameter

        :param _compartment: str
            name of the compartment of interest
        :return: all_sub_parameters: list
            list of all the mortality-related sub-parameters for the compartment of interest
        """
        all_sub_parameters = []
        compartments_strata = find_name_components(_compartment)[1:]
        compartments_strata.reverse()
        compartments_strata.append("")

        # loop through each stratification of the parameter and adapt if the parameter is available
        for stratum in compartments_strata:
            if stratum in self.available_death_rates:
                all_sub_parameters.append("universal_death_rateX" + stratum)
            if "universal_death_rateX" + stratum in self.overwrite_parameters:
                break
        all_sub_parameters.reverse()
        return all_sub_parameters

    def create_mortality_functions(self, _compartment, _sub_parameters):
        """
        loop through all the components to the population-wide mortality and create the recursive functions

        :param _compartment: str
            name of the compartment of interest
        :param _sub_parameters: list
            the names of the functions that need to update the upstream parameters
        :return:
        """
        self.final_parameter_functions["universal_death_rateX" + _compartment] = \
            self.adaptation_functions[_sub_parameters[0]]
        for component in _sub_parameters[1:]:

            # get the new function to act on the less stratified function (closer to the "tree-trunk")
            if component not in self.parameters:
                raise ValueError("parameter component %s not found in parameters attribute" % component)
            elif type(self.parameters[component]) == float:
                update_function = create_multiplicative_function(self.parameters[component])
            elif type(self.parameters[component]) == str:
                update_function = create_time_variant_multiplicative_function(self.adaptation_functions[component])
            else:
                raise ValueError("parameter component %s not appropriate format" % component)

            # create the composite function
            self.final_parameter_functions["universal_death_rateX" + _compartment] = create_function_of_function(
                update_function, self.final_parameter_functions["universal_death_rateX" + _compartment])

    def find_transition_components(self, _parameter):
        """
        finds each of the strings for the functions acting on the next function in the sequence

        :param _parameter: str
            full name of the parameter of interest
        """
        sub_parameters = []

        # work backwards to allow stopping for overwriting requests, then reverse in preparation for function creation
        for x_instance in extract_reversed_x_positions(_parameter):
            component = _parameter[: x_instance]
            sub_parameters.append(component)
            if component in self.overwrite_parameters:
                break
        sub_parameters.reverse()
        return sub_parameters

    def create_transition_functions(self, _parameter, _sub_parameters):
        """
        builds up each parameter to be implemented as a function, recursively creating an outer function that calls the
            inner function

        :param _parameter: str
            full name of the parameter of interest
        :param _sub_parameters: list
            list of the strings representing the sub-parameters, including the base parameter as the stem and with all
                of the relevant strata in the stratification sequence following
        """

        # start from base value as a function of time, even if the time argument is ignored
        if isinstance(self.parameters[_sub_parameters[0]], (float, int)):
            self.final_parameter_functions[_parameter] = lambda time: self.parameters[_sub_parameters[0]]
        elif type(self.parameters[_sub_parameters[0]]) == str:
            self.final_parameter_functions[_parameter] = self.adaptation_functions[_sub_parameters[0]]

        # then cycle through other applicable components and extend function recursively, only if component available
        for component in _sub_parameters[1:]:

            # get the new function to act on the less stratified function (closer to the "tree-trunk")
            if component not in self.parameters:
                raise ValueError("parameter component %s not found in parameters attribute" % component)
            elif isinstance(self.parameters[component], float):
                update_function = create_multiplicative_function(self.parameters[component])
            elif type(self.parameters[component]) == str:
                update_function = create_time_variant_multiplicative_function(self.adaptation_functions[component])
            else:
                raise ValueError("parameter component %s not appropriate format" % component)

            # create the composite function
            self.final_parameter_functions[_parameter] = create_function_of_function(
                update_function, self.final_parameter_functions[_parameter])

    def prepare_infectiousness_calculations(self):
        """
        master method to run all the code concerned with preparation for force of infection calculations
        """

        # infectiousness preparations
        self.prepare_all_infectiousness_multipliers()
        self.find_infectious_indices()

        # mixing preparations
        if self.mixing_matrix is not None:
            self.add_force_indices_to_transitions()
        mixing_indices = self.find_mixing_denominators()

        # reconciling the strains and the mixing attributes together into one structure
        self.find_strain_mixing_multipliers(mixing_indices)

    def prepare_all_infectiousness_multipliers(self):
        """
        find the infectiousness multipliers for each compartment being implemented in the model
        """

        # start from assumption that each compartment is fully and equally infectious
        self.infectiousness_multipliers = [1.] * len(self.compartment_names)

        # if infectiousness modification requested for the compartment type, multiply through by the current value
        for n_comp, compartment in enumerate(self.compartment_names):
            for modifier in self.infectiousness_levels:
                if modifier in find_name_components(compartment):
                    self.infectiousness_multipliers[n_comp] *= self.infectiousness_levels[modifier]

    def find_infectious_indices(self):
        """
        find the infectious indices by strain and overall, as opposed to just overall in EpiModel
        note that this changes the structure by one hierarchical level compared to EpiModel - in that previously we had
            self.infectious_indices a list of infectious indices and now it is has a dictionary structure at the highest
            level, followed by keys for each strain with values being lists that are equivalent to the
            self.infectious_indices list for the unstratified version
        """

        # find the indices for the compartments that are infectious across all strains
        self.infectious_indices["all_strains"] = self.find_all_infectious_indices()

        # then find the infectious compartment for each strain separately
        for strain in self.strains:
            self.infectious_indices[strain] = \
                convert_boolean_list_to_indices(
                    [create_stratum_name("strain", strain, joining_string="")
                     in find_name_components(comp) and self.infectious_indices["all_strains"][i_comp]
                     for i_comp, comp in enumerate(self.compartment_names)])

    def add_force_indices_to_transitions(self):
        """
        find the indices from the force of infection vector to be applied for each infection flow and populate to the
            force_index column of the flows frame
        """

        # identify the indices of all the infection-related flows to be implemented
        infection_flow_indices = \
            [n_flow for n_flow, flow in enumerate(self.transition_flows.type)
             if "infection" in flow and self.transition_flows.implement[n_flow] == len(self.all_stratifications)]

        # loop through and find the index of the mixing matrix applicable to the flow, of which there should be only one
        for n_flow in infection_flow_indices:
            found = False
            for n_group, force_group in enumerate(self.mixing_categories):
                if all(stratum in find_name_components(self.transition_flows.origin[n_flow])
                       for stratum in find_name_components(force_group)):
                    self.transition_flows.force_index[n_flow] = n_group
                    if found:
                        raise ValueError("mixing group found twice for transition flow number %s" % n_flow)
                    found = True
                    continue
            if not found:
                raise ValueError("mixing group not found for transition flow number %s" % n_flow)

    def find_mixing_denominators(self):
        """
        for each mixing category, create a list of the compartment numbers that are relevant

        :return mixing_indices: list
            indices of the compartments that are applicable to a particular mixing category
        """
        if self.mixing_matrix is None:
            return {"all_population": range(len(self.compartment_names))}
        else:
            mixing_indices = {}
            for category in self.mixing_categories:
                mixing_indices[category] = \
                    [i_comp for i_comp, compartment in enumerate(self.compartment_names)
                     if all([component in find_name_components(compartment)
                             for component in find_name_components(category)])]
            return mixing_indices

    def find_strain_mixing_multipliers(self, mixing_indices):
        """
        find the relevant indices to be used to calculate the force of infection contribution to each strain from each
            mixing category as a list of indices - and separately find multipliers as a list of the same length for
            their relative infectiousness extracted from self.infectiousness_multipliers
        """
        for strain in self.strains + ["all_strains"]:
            self.strain_mixing_elements[strain], self.strain_mixing_multipliers[strain] = {}, {}
            for category in ["all_population"] if self.mixing_matrix is None else self.mixing_categories:
                self.strain_mixing_elements[strain][category] = \
                    [index for index in mixing_indices[category] if index in self.infectious_indices[strain]]
                self.strain_mixing_multipliers[strain][category] = \
                    [self.infectiousness_multipliers[i_comp]
                     for i_comp in self.strain_mixing_elements[strain][category]]

    def find_transition_indices_to_implement(self, back_one=0):
        """
        find all the indices of the transition flows that need to be stratified
        separated out as very short method in order that it can over-ride the version in the unstratified EpiModel

        :param back_one: int
            number to subtract from self.all_stratification, which will be one if this method is being called after the
                stratification has been added
        :return: list
            list of indices of the flows that need to be stratified
        """
        return self.transition_flows[self.transition_flows.implement == len(self.all_stratifications) - back_one].index

    def find_death_indices_to_implement(self, back_one=0):
        """
        find all the indices of the death flows that need to be stratified
        separated out as very short method in order that it can over-ride the version in the unstratified EpiModel

        :param back_one: int
            number to subtract from self.all_stratification, which will be one if this method is being called after the
                stratification has been added
        :return: list
            list of indices of the flows that need to be stratified
        """
        return self.death_flows[self.death_flows.implement == len(self.all_stratifications) - back_one].index

    """
    methods to be called during the process of model running
    """

    def get_parameter_value(self, _parameter, _time):
        """
        returns a parameter value by calling the function represented by its string within the parameter_functions
            attribute

        :param _parameter: str
            name of the parameter to be called (key to the parameter_functions dictionary)
        :param _time: float
            current time of model integration
        :return: float
            the parameter value needed
        """
        return self.final_parameter_functions[_parameter](_time)

    def find_infectious_population(self, _compartment_values):
        """
        find vectors for the total infectious populations and the total population that is needed in the case of
            frequency-dependent transmission

        :param _compartment_values: numpy array
            current values for the compartment sizes
        """
        mixing_categories = ["all_population"] if self.mixing_matrix is None else self.mixing_categories
        for strain in self.strains if self.strains else ["all_strains"]:
            self.infectious_populations[strain] = []
            for category in mixing_categories:
                self.infectious_populations[strain].append(sum(
                    element_list_multiplication(
                        [_compartment_values[i_comp] for i_comp in self.strain_mixing_elements[strain][category]],
                        self.strain_mixing_multipliers[strain][category])))
        self.infectious_denominators = sum(_compartment_values)

    def find_infectious_multiplier(self, n_flow):
        """
        find the multiplier to account for the infectious population in dynamic flows

        :param n_flow: int
            index for the row of the transition_flows data frame
        :return:
            the total infectious quantity, whether that is the number or proportion of infectious persons
            needs to return as one for flows that are not transmission dynamic infectiousness flows
        """
        if "infection" not in self.transition_flows.type[n_flow]:
            return 1.0
        strain = "all_strains" if not self.strains else self.transition_flows.strain[n_flow]
        mixing_elements = [1.0] if self.mixing_matrix is None else \
            list(self.mixing_matrix[self.transition_flows.force_index[n_flow], :])
        denominator = 1.0 if "_density" in self.transition_flows.type[n_flow] else self.infectious_denominators
        return sum(element_list_multiplication(self.infectious_populations[strain], mixing_elements)) / denominator

    def get_compartment_death_rate(self, _compartment, _time):
        """
        find the universal or population-wide death rate for a particular compartment

        :param _compartment: str
            name of the compartment
        :param _time: float
            current integration time
        :return: float
            death rate
        """
        return self.get_parameter_value("universal_death_rateX" + _compartment, _time) if \
            len(self.all_stratifications) > 0 else self.get_parameter_value("universal_death_rate", _time)

    def apply_birth_rate(self, _ode_equations, _compartment_values, _time):
        """
        apply a population-wide death rate to all compartments
        all the entry_fraction proportions should be present in either parameters or time_variants given how they are
            created in the process of implementing stratification

        :parameters: all parameters have come directly from the apply_all_flow_types_to_odes method unchanged
        """

        # find the total number of births entering the system at the current time point
        total_births = self.find_total_births(_compartment_values, _time)

        # split the total births across entry compartments
        for compartment in [comp for comp in self.compartment_names if find_stem(comp) == self.entry_compartment]:

            # calculate adjustment to original stem entry rate
            entry_fraction = 1.0
            for stratum in find_name_components(compartment)[1:]:
                entry_fraction *= self.get_single_parameter_component("entry_fractionX%s" % stratum, _time)

            # apply to that compartment
            _ode_equations = increment_list_by_index(
                _ode_equations, self.compartment_names.index(compartment), total_births * entry_fraction)
        return _ode_equations


if __name__ == "__main__":

    def get_total_popsize(model):
        return sum(model.compartment_values)

    sir_model = StratifiedModel(
        numpy.linspace(0, 60 / 365, 61).tolist(),
        ["susceptible", "infectious", "recovered"],
        {"infectious": 0.001},
        {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
        [{"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
         {"type": "infection_density", "parameter": "beta", "origin": "susceptible", "to": "infectious"},
         {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}],
        output_connections={
            "incidence": {"origin": "susceptible", "to": "infectious"},
            "incidence_hiv_positive": {"origin": "susceptible", "to": "infectious",
                                       "origin_condition": "hiv_positive", "to_condition": "hiv_positive"}},
        verbose=False, integration_type="odeint", derived_output_functions={'population': get_total_popsize}
    )
    # sir_model.adaptation_functions["increment_by_one"] = create_additive_function(1.)

    hiv_mixing = numpy.ones(4).reshape(2, 2)
    # hiv_mixing = None

    sir_model.stratify("hiv", ["negative", "positive"], [], {"negative": 0.6},
                       {"recovery": {"negative": "increment_by_one", "positive": 0.5},
                        "infect_death": {"negative": 0.5},
                        "entry_fraction": {"negative": 0.6, "positive": 0.4}},
                       adjustment_requests={"recovery": {"negative": 0.7}},
                       infectiousness_adjustments={"positive": 0.5},
                       mixing_matrix=hiv_mixing,
                       verbose=False)

    sir_model.stratify("strain", ["sensitive", "resistant"], ["infectious"],
                       adjustment_requests={"recoveryXhiv_negative": {"sensitive": 0.9},
                                            "recovery": {"sensitive": 0.8}},
                       requested_proportions={}, verbose=False)

    age_mixing = None
    sir_model.stratify("age", [1, 10, 3], [], {}, {"recovery": {"1": 0.5, "10": 0.8}},
                       infectiousness_adjustments={"1": 0.8},
                       mixing_matrix=age_mixing, verbose=False)

    sir_model.run_model()

    create_flowchart(sir_model, name="sir_model_diagram")

    # create_flowchart(sir_model)
    #
    sir_model.plot_compartment_size(['infectious', 'hiv_positive'])

