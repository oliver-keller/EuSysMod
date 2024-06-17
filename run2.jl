


using AnyMOD, Gurobi, CSV, Statistics, Dates
include("support_functions.jl")



t_int = 8 # Number of threads (4)

# define input and output directories
inputMod_arr = ["_basis","timeSeries/96hours_2008"]
# resultDir_str = "./results/151-300/" 
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
start_iteration = 100
end_iteration = 200

# initialize model (deepcopy of initialized model might lead to a a terminal crash. In this case initialize the model separately for every iteration.)
anyM0 = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true)

outputsOfInterest = nothing # initialize parameter to store the outputs of interest
objective = DataFrame(iteration = Int[], objectiveValue = Float64[]) # initialize dataframe to store the objective values of the optimization runs

# set additional constraint for regret calculations -> Set the minimum use of biomass for the following categories (biomass usage in GWh/a)
set_constraint = true # set to true if additional constraints should be set # make sure that the file "df_input_with_final_cluster.csv" is available in EuSysMod. The file is created by the python script in Decide.

if set_constraint
    cluster_index = 1
    new_constraint = Dict(
        "bioConversionOil" => 0, # 253.5 TWh/a
        "bioConversionSyngas" =>  94050, # 94 TWh/a
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

    #determine if a scenario is already within the respective cluster and should be skipped in the optimization
    cluster_of_scenarios = CSV.read("df_input_with_final_cluster.csv", DataFrame)[!, "cluster_final"] # make sure that the file "df_input_with_final_cluster.csv" is available in EuSysMod. The file is created by the python script in Decide.
end


for iteration in range(start_iteration, end_iteration)
    println("[info] iteration $iteration/$end_iteration")

    if set_constraint && cluster_of_scenarios[iteration] == cluster_index
        continue
    end

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
    set_optimizer_attribute(anyM.optModel, "BarConvTol", 1e-5);  # 1e-5
    optimize!(anyM.optModel) # solve the model

    # update the outputs of interest and the objective value dataframes
    global outputsOfInterest
    global objective
    outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
    push!(objective, [iteration, objective_value(anyM.optModel)])
    
    # Save outputsOfInterest and total_cost as a CSV file
    CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
    set_constraint ? CSV.write(resultDir_str * "/total_cost_cluster_$cluster_index.csv", DataFrame(objective)) : CSV.write(resultDir_str * "/total_cost.csv", DataFrame(objective))
end
