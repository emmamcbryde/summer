import theano.tensor as tt
from python_source_code.tb_model import *
import python_source_code.post_processing as post_proc
from itertools import chain

import pymc3 as pm
import theano
import numpy as np
import logging
logger = logging.getLogger("pymc3")
logger.setLevel(logging.DEBUG)

_logger = logging.getLogger("theano.gof.compilelock")
_logger.setLevel(logging.DEBUG)

theano.config.optimizer = 'None'


def build_model_for_calibration(start_time=1800.):
    input_database = InputDB()

    integration_times = numpy.linspace(start_time, 2020.0, 201).tolist()

    # set basic parameters, flows and times, then functionally add latency
    case_fatality_rate = 0.4
    untreated_disease_duration = 3.0
    parameters = \
        {"contact_rate": 20.,
         "recovery": case_fatality_rate / untreated_disease_duration,
         "infect_death": (1.0 - case_fatality_rate) / untreated_disease_duration,
         "universal_death_rate": 1.0 / 50.0,
         "case_detection": 0.0,
         "crude_birth_rate": 20.0 / 1e3}
    parameters.update(change_parameter_unit(provide_aggregated_latency_parameters(), 365.251))

    # sequentially add groups of flows
    flows = add_standard_infection_flows([])
    flows = add_standard_latency_flows(flows)
    flows = add_standard_natural_history_flows(flows)

    # compartments
    compartments = ["susceptible", "early_latent", "late_latent", "infectious", "recovered"]

    # define model
    _tb_model = StratifiedModel(
        integration_times, compartments, {"infectious": 1e-3}, parameters, flows, birth_approach="replace_deaths")

     # add case detection process to basic model
    _tb_model.add_transition_flow(
        {"type": "standard_flows", "parameter": "case_detection", "origin": "infectious", "to": "recovered"})

    # age stratification
    age_breakpoints = [5, 15]
    age_infectiousness = get_parameter_dict_from_function(logistic_scaling_function(15.0), age_breakpoints)
    age_params = get_adapted_age_parameters(age_breakpoints)
    age_params.update(split_age_parameter(age_breakpoints, "contact_rate"))

    _tb_model.stratify("age", copy.deepcopy(age_breakpoints), [], {}, adjustment_requests=age_params,
                       infectiousness_adjustments=age_infectiousness, verbose=False)


     # get bcg coverage function
    _tb_model = get_bcg_functions(_tb_model, input_database, 'MNG')

    # stratify by vaccination status
    bcg_wane = create_sloping_step_function(15.0, 0.7, 30.0, 0.0)
    age_bcg_efficacy_dict = get_parameter_dict_from_function(lambda value: bcg_wane(value), age_breakpoints)
    bcg_efficacy = substratify_parameter("contact_rate", "vaccinated", age_bcg_efficacy_dict, age_breakpoints)
    # _tb_model.stratify("bcg", ["vaccinated", "unvaccinated"], ["susceptible"],
    #                    requested_proportions={"vaccinated": 0.0},
    #                    entry_proportions={"vaccinated": "bcg_coverage",
    #                                       "unvaccinated": "bcg_coverage_complement"},
    #                    adjustment_requests=bcg_efficacy,
    #                    verbose=False)

    # loading time-variant case detection rate
    input_database = InputDB()
    res = input_database.db_query("gtb_2015", column="c_cdr", is_filter="country", value="Mongolia")

    # add scaling case detection rate
    cdr_adjustment_factor = 1.
    cdr_mongolia = res["c_cdr"].values / 1e2 * cdr_adjustment_factor
    cdr_mongolia = numpy.concatenate(([0.0], cdr_mongolia))
    res = input_database.db_query("gtb_2015", column="year", is_filter="country", value="Mongolia")
    cdr_mongolia_year = res["year"].values
    cdr_mongolia_year = numpy.concatenate(([1950.], cdr_mongolia_year))
    cdr_scaleup = scale_up_function(cdr_mongolia_year, cdr_mongolia, smoothness=0.2, method=5)
    prop_to_rate = convert_competing_proportion_to_rate(1.0 / untreated_disease_duration)
    detect_rate = return_function_of_function(cdr_scaleup, prop_to_rate)
    _tb_model.time_variants["case_detection"] = detect_rate

    _tb_model.stratify("strain", ["ds", "mdr"], ["early_latent", "late_latent", "infectious"], {}, verbose=False)

    _tb_model.stratify("smear", ["smearpos", "smearneg", "extrapul"], ["infectious"],
                       adjustment_requests={}, verbose=False, requested_proportions={})

    return _tb_model


class Calibration:
    """
    this class handles model calibration using an MCMC algorithm if sampling from the posterior distribution is
    required, or using maximum likelihood estimation if only one calibrated parameter set is required.
    """
    def __init__(self, model_builder, priors, targeted_outputs):
        self.model_builder = model_builder  # a function that builds a new model without running it
        self.base_model = model_builder()  # a built model that has not been run
        self.running_model = None  # a model that will be run during calibration
        self.post_processing = None  # a PostProcessing object containing the required outputs of a model that has been run
        self.priors = priors  # a list of dictionaries. Each dictionary describes the prior distribution for a parameter
        self.param_list = [self.priors[i]['param_name'] for i in range(len(self.priors))]
        self.targeted_outputs = targeted_outputs  # a list of dictionaries. Each dictionary describes a target
        self.data_as_array = None  # will contain all targeted data points in a single array

        self.loglike = None  # will store theano object

        self.format_data_as_array()
        self.workout_unspecified_sds()
        self.create_loglike_object()

        self.mcmc_trace = None  # will store the results of the MCMC model calibration
        self.mle_estimates = {}  # will store the results of the maximum-likelihood calibration

    def update_post_processing(self):
        """
        updates self.post_processing attribute based on the newly run model
        :return:
        """
        if self.post_processing is None:  # we need to initialise a PostProcessing object
            requested_outputs = [self.targeted_outputs[i]['output_key'] for i in range(len(self.targeted_outputs))]
            requested_times = {}

            for output in self.targeted_outputs:
                requested_times[output['output_key']] = output['years']

            self.post_processing = post_proc.PostProcessing(self.running_model, requested_outputs, requested_times)
        else:  # we just need to update the post_processing attribute and produce new outputs
            self.post_processing.model = self.running_model
            self.post_processing.generated_outputs = {}

            self.post_processing.generate_outputs()

    def run_model_with_params(self, params):
        """
        run the model with a set of params.
        :param params: a dictionary containing the parameters to be updated
        """
        if 'start_time' in self.param_list:  # we need to re-build a model
            this_param_index = self.param_list.index('start_time')
            self.running_model = self.model_builder(start_time=params[this_param_index])
        else:  # we cjust need to copy the existing base_model
            self.running_model = copy.deepcopy(self.base_model)  # reset running model

        # update parameter values
        for i in range(len(params)):
            param_name = self.priors[i]['param_name']
            if param_name != 'start_time':
                value = params[i]
                self.running_model.parameters[param_name] = value

        # run the model
        self.running_model.run_model()

        # perform post-processing
        self.update_post_processing()

    def loglikelihood(self, params):
        """
        defines the loglikelihood
        :param params: model parameters
        :return: the loglikelihood
        """
        # run the model
        self.run_model_with_params(params)

        ll = 0
        for target in self.targeted_outputs:
            key = target['output_key']
            data = np.array(target['values'])
            model_output = np.array(self.post_processing.generated_outputs[key])

            print("############")
            print(data)
            print(model_output)

            ll += -(0.5/target['sd']**2)*np.sum((data - model_output)**2)

        return ll

    def format_data_as_array(self):
        """
        create a list of data values based on the target outputs
        """
        data = []
        for target in self.targeted_outputs:
            data.append(target['values'])

        data = list(chain(*data))  # create a simple list from a list of lists
        self.data_as_array = np.asarray(data)

    def workout_unspecified_sds(self):
        """
        If the sd parameter of the targeted output is not specified, it will automatically be calculated such that the
        95% CI of the associated normal distribution covers 50% of the mean value of the target.
        :return:
        """
        for i, target in enumerate(self.targeted_outputs):
            if 'sd' not in target.keys():
                self.targeted_outputs[i]['sd'] = 0.5 / 4. * np.mean(target['values'])
                print(self.targeted_outputs[i]['sd'])

    def create_loglike_object(self):
        """
        create a 'theano-type' object to compute the likelihood
        """
        self.loglike = LogLike(self.loglikelihood)

    def run_fitting_algorithm(self, run_mode='mle', mcmc_method='Metropolis', n_iterations=100, n_burned=10, n_chains=1,
                              parallel=True):
        """
        master method to run model calibration.

        :param run_mode: either 'mcmc' (for sampling from the posterior) or 'mle' (maximum likelihood estimation)
        :param mcmc_method: if run_mode == 'mcmc' , either 'Metropolis' or 'DEMetropolis'
        :param n_iterations: number of iterations requested for sampling (excluding burn-in phase)
        :param n_burned: number of burned iterations before effective sampling
        :param n_chains: number of chains to be run
        :param parallel: boolean to trigger parallel computing
        """
        basic_model = pm.Model()
        with basic_model:

            fitted_params = []
            for prior in self.priors:
                if prior['distribution'] == 'uniform':
                    value = pm.Uniform(prior['param_name'], lower=prior['distri_params'][0],
                                       upper=prior['distri_params'][1])
                    fitted_params.append(value)

            theta = tt.as_tensor_variable(fitted_params)

            pm.DensityDist('likelihood', lambda v: self.loglike(v), observed={'v': theta})

            if run_mode == 'mle':
                self.mle_estimates = pm.find_MAP()
            elif run_mode == 'mcmc':  # full MCMC requested
                if mcmc_method == 'Metropolis':
                    mcmc_step = pm.Metropolis()
                elif mcmc_method == 'DEMetropolis':
                    mcmc_step = pm.DEMetropolis()
                else:
                    ValueError("requested mcmc mode is not supported. Must be one of ['Metropolis', 'DEMetropolis']")

                self.mcmc_trace = pm.sample(draws=n_iterations, step=mcmc_step, tune=n_burned, chains=n_chains,
                                            progressbar=False, parallelize=parallel)
            else:
                ValueError("requested run mode is not supported. Must be one of ['mcmc', 'lme']")


class LogLike(tt.Op):
    """
    Define a theano Op for our likelihood function.
    Specify what type of object will be passed and returned to the Op when it is
    called. In our case we will be passing it a vector of values (the parameters
    that define our model) and returning a single "scalar" value (the
    log-likelihood)
    """
    itypes = [tt.dvector]  # expects a vector of parameter values when called
    otypes = [tt.dscalar]  # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike):
        """
        Initialise the Op with various things that our log-likelihood function
        requires. Below are the things that are needed in this particular
        example.

        :param loglike: The log-likelihood function
        """

        # add inputs as class attributes
        self.likelihood = loglike

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        theta, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta)

        outputs[0][0] = np.array(logl)  # output the log-likelihood


if __name__ == "__main__":

    par_priors = [{'param_name': 'contact_rate', 'distribution': 'uniform', 'distri_params': [2., 100.]},
                  {'param_name': 'late_progression', 'distribution': 'uniform', 'distri_params': [.001, 0.003]},
                  {'param_name': 'start_time', 'distribution': 'uniform', 'distri_params': [1800., 1850.]}
                  ]
    target_outputs = [{'output_key': 'prevXinfectiousXamongXage_15', 'years': [2015, 2016], 'values': [0.005, 0.004],
                       'sd': 0.0005},
                      {'output_key': 'prevXlatentXamongXage_5', 'years': [2014], 'values': [0.096], 'sd': 0.012}
                     ]
    calib = Calibration(build_model_for_calibration, par_priors, target_outputs)

    calib.run_fitting_algorithm(run_mode='mle')  # for maximum-likelihood estimation

    print(calib.mle_estimates)
    #
    # calib.run_fitting_algorithm(run_mode='mcmc', mcmc_method='DEMetropolis', n_iterations=100, n_burned=10,
    #                             n_chains=4, parallel=True)  # for mcmc

