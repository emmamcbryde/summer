# SUMMER Handbook
Scalable, Universal Mathematical Modelling in R and Python.

We encourage wide use of this repository with acknowledgement. Although it is anticipated that the code base will be
developed by our team with time, the general functionality of the base code should not change markedly.

# Introduction

## Philosophical approach
The purpose of this repository is to provide a codebase to facilitate the construction of mathematical models producing
epidemiological simulations of infectious disease transmission.
Many aspects of models of infectious disease transmission are built using high-level packages. These include numeric
integrators (e.g. ODE solvers), packages to manipulate data structures, Bayesian calibration tools and graphing packages
to visualise outputs.

However, it is our group's experience that pre-built packages are less frequently used for defining the model structures
themselves (although this is not to say that this is never done and we are aware of some other packages that exist).
Nevertheless, it remains common modelling practice to define compartmental models of infectious disease transmission by
hand-coding a series of ordinary differential equations representing the rates of change of each modelled population
compartment.
In our view, this process can be facilitated with a modular code base to allow greater focus of modellers
on the epidemiological results of their simulations.

## Benefits
It is hoped that this approach to model construction will have the following advantages:
* Construction of more complicated models without the need for a large number of hand-coded ODEs or code loops
* Avoidance of errors that may occur with manually written code (e.g. transitions between compartments that do not have
an equivalent inflow to the destination compartment to balance the outflow from the origin compartment, etc.)
* Easy visualisation of the process of model construction through flow-diagrams of inter-compartmental flows
* Accessibility of modelling code and model construction to epidemiologists, policy-makers, clinicians, etc. who do not
have an extensive background in mathematical modelling
* Future extensibility to a broader range of capabilities as our group implements further elaboration of the modelling
platform

## Code
The code base is provided in R and Python as two open-source platforms that are among the most commonly used in
infectious disease modelling applications.
Python is more naturally suited to object-oriented programming (OOP), such that the R code requires the R6 package as a
dependency through much of its structure.
The code in the two languages is intended to be as equivalent as possible.
This document is intended to explain the SUMMER structure using language and terminology familiar to readers with a
background in OOP and familiarity with these two languages.
Therefore, although it has been developed with infectious disease epidemiologists in mind, it is not intended to be
readily readable by epidemiologists without also reading some background on these subjects.
Nevertheless, a very brief introduction to OOP follows.
In general Python lists are implemented as R vectors and Python dictionaries are implemented as R lists.

## Object-oriented structure
The OOP used in SUMMER is a commonly used paradigm in computer programming.
It involves creating objects, which have attributes (collections of data) and methods (which govern the object's
behaviour).
Objects are instances of classes, with SUMMER defining a number of classes that we have developed.
Therefore, a class defines the basic structure of an object and can be thought of as a template for objects.
An important feature of OOP is inheritance.
Inheritance allows a new class to be developed from an existing class, and to inherit its structure and features.

SUMMER considers epidemiological models as objects with attributes and methods relevant to the general construction of
epidemiological models typically used to simulate infectious disease transmission dynamics. When a model object is
instantiated, it is created with a set of features that allow for the construction of a compartmental model in a
standard form.

The current code structure allows for the creation of stratified or unstratified models, with the StratifiedModel class
inheriting from the unstratified basic epidemiological model class called EpiModel.
In this way, EpiModel allows for the manual construction of models of any degree of complexity, while the additional
features provided in the StratifiedModel class allow for rapid and reliable stratification of the model compartments
created in EpiModel without the need for the repetitive code or the use of loops.

The general approach to using the model objects is described below, but the details of how individual arguments should
be formatted, etc. is provided in the docstrings to each method of the model object(s).
The attributes of the classes are described in single docstrings at the start of these class definitions.

# Workflow of model creation
## Base model construction
As mentioned, the base EpiModel class allows for the construction of standard compartmental epidemiological models
implemented in ODEs. The object constructor must be provided with arguments that include:
1. The times at which the model outputs are to be evaluated
    * Note that time units are arbitrary and the user should recall what time unit is being used and ensure parameter
    value units and requested time are requested with the same time unit in mind (e.g. day, year, etc.)
2. The names of the compartments (if unstratified) or types of compartments (if stratified) to be used by the model
3. The distribution of the initial conditions of the population
4. Model parameters for the calculation of flows
5. Inter-compartmental flows
6. A range of optional arguments pertaining to modelling and reporting features

The first five arguments are required, but inter-compartmental flows can be provided as an empty list and populated
with specific flows later.
Note also that not all compartments must have an initial value provided at the construction stage, as default
assumptions can control this distribution later.

Get more information on the process of model construction by setting verbose to True in object instantiation.

### Compartment names and user-submitted strings
Compartment names can be any strings, although in general any user-submitted strings provided to SUMMER should:
* Avoid the use of the characters "X", "W" and "_", which trigger specific behaviours by SUMMER
* Be lower case (as future development may use further upper case letters in addition to X and W)
* Not exactly duplicate strings that were previously used
* Strings that are as descriptive as possible are encouraged
Therefore, typical compartment names might include "susceptible", "infectious", "recovered".

### Initial conditions
The process of setting the initial values of the model runs as follows:
1. All compartments requested in construction are set to values of zero by default.
Alternative initial values should be requested for specific compartments in the initial conditions, but if no request is
made, each compartment will start with a value of zero.
2. Each requested compartment in the keys/names of the initial conditions request are looped through and checked to
ensure that they are available compartments. If they are, the value is set to the request.
3. If the user has requested to sum initial conditions to a total value, the remainder of the total population that has
not yet been allocated will be assigned to the requested starting compartment. If this is a negative value (i.e. if more
population has been allocated than the already requested starting population), an error will result. The starting
compartment is assigned as a user input, but if no user input is requested, then the compartment to which population
births accrue will be used.

### Flow requests
Inter-compartmental or death flows are requested as a dict/list. The keys/names of the request must be:
* type: specifies the type of the flow being requested as standard_flows, infection_density, infection_frequency,
compartment_death or customised_flows
* parameter: provides the parameter name that the model should use to calculate the transition rate during integration
* origin: the name of the compartment from which the flow arises
* to: the name of the compartment that the transition goes towards, not applicable to death flows
standard_flows will calculate the flow rate as the parameter value by the size of the origin compartment.

infection_density flows will calculate the flow rate as the parameter value multiplied by the size of the infectious
population multiplied by the size of the origin compartment.

infection_frequency flows will calculate the flow rate as the parameter value multiplied by the size of the infectious
population divided by the total population size multiplied by the size of the origin compartment.

customised_flows use a user-defined function to define a rate of transition, which may not necessarily be multiplied by
the size of the origin compartment, such that negative compartment values could result. We consider this a more advanced
feature that would typically be used by experienced programmers when developing more complicated models and is not a
necessary feature of most standard compartmental models.

### Model running
Model running is called through the run_model method to the model object once constructed
In Python, odeint and solve_ivp have been implemented, where solve_ivp can be used to stop model integration once
equilibrium has been reached (if requested).
Both integration functions are available from the scipy.integrate module.

### Infectious compartment(s)
Differing from earlier releases, this is now specified as a list of all the compartment types that are infectious, in
order to allow for multiple infectious compartment types.
For example, patients may be infectious when active undiagnosed, as well as after diagnosis.
This can be requested by specifying the infectious_compartments with multiple list elements
(e.g. infectious_compartment=["infectious", "detected"] as the Python argument in construction).
The default argument to infectious_compartment is ["infectious"] (so the model will not run if no alternative argument
is submitted and "infectious" is not included in the compartment types requested).

### Tracking outputs
For outputs such as incidence and total mortality, a dictionary/list can be passed (output_connections) to request
that specific model quantities emerging from the model are tracked during integration.
The key/name of the quantity is the name of the indicator, while the value is a dictionary/list specifying the origin
and destination ("to") compartments as well as conditional strings that should appear in the origin and destination 
compartment names. For example, in an SIR model stratified by HIV status in Python, this would be specified as: 
{"incidence": {"origin": "susceptible", "to": "infectious"},
 "incidence_hiv_positive": {"origin": "susceptible", "to": "infectious", "origin_condition": "hiv_positive", 
							"to_condition": "hiv_positive"}
}.

Alternatively, user-defined functions can be tracked with any function used to track quantities within the model. These
are submitted as dictionaries with keys to name the quantities and values to define the functions. 

### Birth rates
There are currently three options for calculating the rate of birth implemented in SUMMER (although further approaches
may be introduced in the future).
The currently implemented options are:
* no_birth: No births/recruitment is implemented.
This will result in a declining population over time if either infection-related deaths or population-wide mortality
are also implemented.
* add_crude_birth_rate: The parameter "crude_birth_rate" is multiplied by the total population size to determine the
birth/recruitment rate.
This will result in the population increasing over time if the birth rate is higher than all death rates (including
infection-related deaths and population-wide mortality).
* replace_deaths: Under this approach the total rate of all deaths to be tracked throughout model integration (which is
implemented automatically if this approach is requested) and these deaths re-enter the model as births.
This will result in a stable population over time, or a "closed" model in epidemiological terms.
Whichever approach is adopted, births then accrue to the compartment nominated by the user as the entry_compartment,
which is assumed to be "susceptible" by default.

### Death rates
In addition to the disease/infection-specific death rates that are implemented individually as mentioned above, users
may request a universal, population-wide death rate to be added to all compartments.
This can be achieved by specifying the parameter "universal_death_rate" in the parameters submitted to the model, which
may be specified as a constant or time-variant parameter.
If this parameter name is not provided, it will be taken to be zero by default.
The universal death rate can also be altered for population groups (eg. age groups) in stratified versions of the model.

### Outputs
Once integration has completed, the class will have an attribute called outputs to represent the compartment sizes
estimated at each time step, as well as any additional outputs collated into the derived_output attribute.
Visualisation and other functions for outputting results have deliberately been kept minimal, as it is anticipated that
other software platforms/modules would be used for these purposes.
However, simple functions have been provided to calculate the total size of several compartments, which may be
stratified and to plot the total value of this group of compartments. 

## Stratification
If a stratified model is required, this can be built through the class StratifiedModel.
This provides additional methods for automatically stratifying the original model built using the parent (EpiModel)
class's methods only (although the model object should be created as an instance of StratifiedModel throughout to ensure
all methods are available).
Stratification could also be achieved manually by specifying compartments and flows through loops, etc. and using the
EpiModel class only.
However, this would require longer and more complicated code, with the associated potential for errors.
Enabling a simpler and faster approach to stratification was a major part of the reason for the development of SUMMER.

Model stratification is requested through the stratify method to StratifiedModel, which is the master method that calls
all other processes that need to be applied to the base model structure to construct the stratified mode.
After construction, the model is then run using the run_model method as it would be for the unstratified EpiModel
version, with some additional methods run in preparation for integration at this stage.
These methods do not need to be run at every integration step and include processes such as determining which
compartment indices fall within which strata, pre-assigning dictionary structures to calculate forces of infection for
multiple strains and for multiple population groups under heterogeneous mixing, etc.

Get more information on the process of stratification as it proceeds by requesting verbose to be True when calling the
stratify method.
(Note that this argument differs from the verbose argument to the model construction process, such that verbose can be
turned on for model construction but off for stratification, or vice versa). 

### Stratification naming
The name of the stratification to be applied to the model is requested as the first argument to the stratify method.
This should be provided as a string and have the features described above in "Compartment names and user-submitted
strings". In addition, the strings "age" and "strain" have specific behaviours as described below, and so should not be
used unless this behaviour is specifically desired.

### Naming of strata
The strata levels to be implemented would typically be submitted as a list of string/character variables for the names
of these levels as the second argument to the stratify method.
Non-string stratum names may be submitted and will be converted to strings/characters during the stratification process
if this is possible.
Alternatively, if an integer value is submitted, the levels will be automatically considered to be all the integer
values from one to this value (e.g. 1, 2, 3 if 3 is submitted).
Otherwise, if an alternative single value (e.g float, boolean) is submitted, stratification will fail. 

### Request for compartments to which stratification should apply
The third argument to the stratify method is the last required argument and is a list of the base compartments to which
the stratification should apply.
For example, for features of the active infection state, only the name of this state would be applied.
If an empty list/vector is supplied, the assumption is made that the stratification should apply to all compartments,
which would be the typical assumption for conditions that affect or are features of the entire population (e.g.
comorbidities, geographical regions, etc.).
Accordingly, for the case of age, the compartments for stratification must be provided as an empty list or an error will
be raised.

### Age stratification
As noted above, the stratification name should be "age" when age stratification is desired, as several features have
been built in to the stratification process when this stratification name is used.

The names of the strata levels within this stratification must be submitted as numeric values (integer or float).
These values represent "age breakpoints", such that the age classes are separated at each of these points.
Each age stratum is named by the lower limit of the age range.
Therefore, if 0 is not included as an age breakpoint, it will be automatically added to represent those aged less than
the lowest age breakpoint.
Requests submitted un-ordered will be sorted for implementation.
For example, if the user submits the values 15 and 5 in this (or any) order, the model will be stratified into the age
ranges: 0 to 4, 5 to 14 and 15+, with these three stratum names referred to as 0, 5 and 15 in the model code.

Ageing flows (that is transitions from one age category to the next, older age category) are implemented automatically
for every age stratum except for the oldest one.
The rate of the flow is equal to the inverse of the width of the age range.
For example, the rate of transition from the 0 to 4 bracket to the 5 to 14 bracket in the example above would be 0.2.

### Starting proportions for initial conditions
The population assigned to each compartment under the user request for initial conditions (as described above) must be
split between the strata of the stratification being requested - if the stratification applies to that compartment.
This can be set through the user request, using a dictionary/list with keys being the strata names and values being the
proportion of the population to be assigned to each stratum.
The starting value provided for each stratum must be associated with a key representing one of the strata names, and the
total of the values assigned must be less than or equal to one.

For any stratum not provided in this argument, an equal proportion of the unallocated starting population will be
assigned to that stratum (i.e. the reciprocal of the number of unassigned strata).
For example, if the strata are "positive" and "negative" and the user requests 0.6 for "positive", 0.4 will be assigned
to "negative".
If the strata are "a", "b" and "c" and the user requests 0.4 for "a", then "0.3" will be assigned to "b" and "c".

### Entry flows
The distribution of new births (or "entries/recruitment") into strata may also be specified by the user, which is done
through a separate argument to that for initial conditions (earlier SUMMER releases differ in this behaviour).
This consideration only applies if the entry compartment is stratified.
The proportion of births assigned to each stratum can be submitted as a constant or time-variant parameter, as for
parameters pertaining to other model processes, such as inter-compartmental transitions.

If no request is submitted for any particular stratum, the proportion of births assigned to that stratum is calculated
as the reciprocal of the number of strata being applied in the current stratification process.
Note that this will result in correct calculations if no request is submitted for any stratum.
Alternatively, if requests are submitted for some or all strata, the user must ensure these sum to one (although we may
implement code in future releases to ensure this problem cannot occur).
The reason for this is that the proportion of births allocated to the unspecified proportion will still take on the
value of the reciprocal of the number of strata being implemented.
This may or may not result in the total proportions summing to one.
It is difficult to ensure that the proportions must always sum to one because user-submitted time-variant functions may
be used for some of the strata, as this cannot be done with a simple normalisation function.

### Parameter stratification
Adjustments to existing parameters can be applied through the process of stratification through the optional
adjustment_requests argument to the stratify method.

This approach can be applied to any existing parameter in the model, including those that have previously been
stratified.
(Request verbose at the previous level of stratification to see the names of the parameters after stratification has
been applied.)
Additionally, this can be applied to the universal_death_rate parameter, which is a specifically named parameter that
depletes all model compartments.

By default, any user request for a particular parameter and stratum is interpreted as a multiplier.
There are two ways to avoid this behaviour and instead supply a new parameter that will ignore the value of the
unstratified parameter and instead overwrite the stratified parameter value for the stratum requested with a new
parameter value:
1. Add the string "W" to the end of the stratum string (e.g. "negativeW") within the parameter being requested
2. Provide an additional component in the request for the parameter which provides all the strata for which overwriting
is desired, which should be named "overwrite"

When stratified model parameters are implemented, implementation can be thought of conceptually as working through the
smallest branches of the tree of stratification towards the trunk.
At each branching point, the user may have specified an adjustment to the parameter at the level below (closer to the
trunk) in the hierarchy.
If this is the case, the lower parameter is modified by the function specified in the call to stratification, with the
function being a simple multiplication if the parameter was submitted as a float.
If no value was specified, the parameter remains unadjusted and if an overwrite parameter was specified the process
stops at that point in the branching process with the parameter specified at that point of stratification.

If a string is provided rather than a numeric/float value, this will be interpreted as representing a function, that
will then be looked for as the key to the parameter_function attribute of the model.
Therefore, a function must be provided before integration as the corresponding value to the parameter_function
attribute.

### Heterogeneous mixing
Heterogeneous mixing between specific compartment types can be requested by submitting a square numpy array with values
representing the group-specific contact rates between the strata being requested, such that this array must also have
both dimensions equal to the number of strata being implemented.
Note that this feature can only be applied if the stratification is being applied to all compartments.

If this heterogeneous mixing is applied at multiple stages of stratification, the new mixing matrix will be square with
dimensions equal to n x m, where n is the number of strata in the first heterogeneous mixing stratification and m is the
number of strata in the second heterogeneous mixing stratification.
In this case, the composite mixing matrix is calculated as the kronecker product of the two mixing matrices for the two
stratifications.
The contents of the matrices may not strictly be per capita contact rates, because the rates of contact are still
calculated by multiplying through by the parameter associated with that transition (although if this were set to one, 
then this would be the case).

### Heterogeneous infectiousness
Heterogeneous infectiousness of different strata can be specified through the infectiousness_adjustments argument to the
stratify method.
These are specified as a dictionary of relative infectiousness values with keys the strata names and values the relative
infectiousness of that strata.

### Multi-strain models
Multi-strain models can be created by specifying "strain" as the stratification name when calling stratify.
As with "age", using this stratification name will result in specific model behaviours.
When this is done, the force of infection will be calculated separately for each stratum, representing a strain of the
organism of interest.
Only compartments representing infection with the organism should be included in the compartments being stratified.

# Model structure visualisation

## Model attributes
These are described in detail in the docstrings to the EpiModel class and the StratifiedModel class.
However, note that the two attributes that are stored as data frames can be easily examined.
For example, in Python, calling the to_csv method to the transition and death data frames allows easy visualisation of
these structures in Excel.
This has the following columns:
* type: the type of flow being implemented, specifying the behaviour of the flow
* parameter: the name of the parameter being used for the transition calculation
* origin: the compartment that the flow comes from or depletes
* to: the compartment that the flow goes to or increments - note that this applies to transition flows only and not to
death flows
* implement: an integer value specifying which stratification this applies to - note that only flows for which the value
of this field is equal to the number of stratifications implemented are actually applied, the others are recorded for
use in flow diagram plotting, etc.

## Flow diagram
A flow diagram of the model compartments and inter-compartmental flows can be created with the create_flowchart function
which takes the model object as an argument.



