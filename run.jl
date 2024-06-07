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
resultDir_str = "./results/" * Dates.format(now(), "yyyy-mm-dd_HH-MM")
if !isdir(resultDir_str)
    mkdir(resultDir_str)
end

# name of the modelrun
obj_str = "Sampling_300"


# read uncertain parameters
uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])
# foreach(x -> uncertain_parameters[!,x] = string.(uncertain_parameters[!,x]), [:parameter, :carrier, :region, :technology, :timestep])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0) # creates a new Latin Hypercube Sample and saves as csv
lhs = CSV.read("lhs.csv", DataFrame, header=false) # uses an already created lhc sample

# set iteration boundaries
start_iteration = 1
end_iteration = 3

# initialize model (deepcopy of initialized model might lead to a a terminal crash. In this case initialize the model separately for every iteration.)
anyM0 = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true)

outputsOfInterest = nothing # initialize parameter to store the outputs of interest
objective = DataFrame(iteration = Int[], objectiveValue = Float64[]) # initialize dataframe to store the objective values of the optimization runs

# set additional constraint for regret calculations -> Set the minimum use of biomass for the following categories (biomass usage in GWh/a)
new_constraint = Dict(
    "bioConversionOil" => 0,
    "bioConversionSyngas" => 0,
    "bioConversionBiogas" => 0,
    "bioConversionChp" => 0,
    "bioConversionHvc" => 0,
    "bioConversionCoal" => 0,
    "networkHeat" => 0,
    "spaceHeat" => 0,
    "proHeat" => 0
)

for row in eachrow(anyM0.parts.lim.par[:useLow].data)
    row.val = get(new_constraint, findTechnology(anyM0, row.Te), row.val)
end


for iteration in range(start_iteration, end_iteration)
    println("[info] iteration $iteration/$end_iteration")

    anyM = deepcopy(anyM0) # deepcopy of initialized model might lead to a a terminal crash => initialize the model for every iteration
    # anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true) 

    # modify uncertain parameters
    modify_parameters(anyM, uncertain_parameters, lhs, iteration)

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

    # update the outputs of interest and the objective value dataframes
    global outputsOfInterest
    global objective
    outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
    push!(objective, [iteration, objective_value(anyM.optModel)])
    
    # Save outputsOfInterest and total_cost as a CSV file
    CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
    CSV.write(resultDir_str * "/total_cost.csv", DataFrame(objective))
end
