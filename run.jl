### Low regret strategies for the use of biomass in the European energy transition ###
### Master's Thesis Oliver Keller 2024 ###

# EuSysMod modified to contain the different biomass conversion pathways
# uncertainty is adressed with an additional file where the uncertainties can be specified
# this script creates a latin hyper cube samplie and varies the uncertain inputs
# the optimization of the model is executed for every input combination created with the LHS

## the functions are defined in the document "support_functions.jl"
## the output is written in the results folder. There a subfolder is created with date/time.
## The output contains the file "outputsOfInterest.csv" which contains the values how much biomass (TWh/a) is used for which purpose. 
## Additionaly for each run a table containing the different biomass carriers and the outputs of interests is saved.

### The file "outputsOfInterest.csv" is used as input for the decision tree creation with the DECIDE tool written in python (the adapted version can be found under https://github.com/oliver-keller/Decide)
### using the output of DECIDE the regret for the different strategies can be calculated




using AnyMOD, Gurobi, CSV, Statistics, Dates
include("support_functions.jl")

# b = "C:/Git/EuSysMod/"
# par_df = CSV.read(b * "settings.csv", DataFrame)

# if isempty(ARGS)
#     id_int = 1
#     t_int = 4
# else
#     id_int = parse(Int,ARGS[1])
#     t_int = parse(Int,ARGS[2]) # number of threads
# end


t_int = 8 # Number of threads (4)

# define input and output directories
inputMod_arr = ["_basis","timeSeries/96hours_2008"]
resultDir_str = "./results/" * Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
if !isdir(resultDir_str)
    mkdir(resultDir_str)
end

# name of the modelrun
obj_str = "Sampling_1000"


# read uncertain parameters
uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])
# foreach(x -> uncertain_parameters[!,x] = string.(uncertain_parameters[!,x]), [:parameter, :carrier, :region, :technology, :timestep])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0) # creates a new Latin Hypercube Sample and saves as csv
lhs = CSV.read("lhs.csv", DataFrame, header=false) # uses an already created lhc sample

# set iteration boundaries
start_iteration = 1
end_iteration = 100

# initialize model (deepcopy of initialized model might lead to a a terminal crash. In this case initialize the model separately for every iteration.)
anyM0 = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true)
outputsOfInterest = nothing


for iteration in range(start_iteration, end_iteration)
    println("[info] iteration $iteration/$end_iteration")

    anyM = deepcopy(anyM0) # deepcopy of initialized model might lead to a a terminal crash => initialize the model for every iteration
    # anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true) 

    # modify uncertain parameters
    modify_parameters(anyM, uncertain_parameters)

    # create optimization model
    createOptModel!(anyM)
    setObjective!(:cost, anyM)

    # solve model
    set_optimizer(anyM.optModel, Gurobi.Optimizer)  # select a solver
    set_optimizer_attribute(anyM.optModel, "Method", 2);
    set_optimizer_attribute(anyM.optModel, "Crossover", 0);
    set_optimizer_attribute(anyM.optModel, "Threads",t_int);
    set_optimizer_attribute(anyM.optModel, "BarConvTol", 1e-3);  # 1e-5
    optimize!(anyM.optModel) # solve the model

    global outputsOfInterest
    outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
    
    # Save outputsOfInterest as a CSV file
    CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
end
