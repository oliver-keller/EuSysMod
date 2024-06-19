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

MODELRUN = 1

using AnyMOD, Gurobi, CSV, Statistics, Dates
include("support_functions.jl")

par_df = CSV.read("settings.csv", DataFrame)

if isempty(ARGS)
    id_int = MODELRUN # id of the modelrun
    t_int = 8 # number of threads
else
    id_int = parse(Int,ARGS[1]) # id of the modelrun
    t_int = parse(Int,ARGS[2]) # number of threads
end

## Define the input parameter ##
START_ITERATION = par_df[id_int,:START_ITERATION]               # set iteration boundaries
END_ITERATION = par_df[id_int,:END_ITERATION]                   # set iteration boundaries
OBJ_STR_INPUT = par_df[id_int,:OBJ_STR]                         # name of the modelrun
MIN_BIOMASS_FOR_OIL = par_df[id_int,:MIN_BIOMASS_FOR_OIL]       # constraint for using certain amount of biomass for oil produiction (TWh/a)
MIN_BIOMASS_FOR_HVC = par_df[id_int,:MIN_BIOMASS_FOR_HVC]       # constraint for using certain amount of biomass for HVC produiction (TWh/a)
MIN_BIOMASS_FOR_SYNGAS = par_df[id_int,:MIN_BIOMASS_FOR_SYNGAS] # constraint for using certain amount of biomass for syngas produiction (TWh/a)
CLUSTER_INDEX = par_df[id_int,:CLUSTER_INDEX]                   # cluster index for the additional constraint (if no constraint parameter is not used)
BARRIER_CONV_TOL = par_df[id_int,:BARRIER_CONV_TOL]             # barrier convergence tolerance
NUMERIC_FOCUS = par_df[id_int,:NUMERIC_FOCUS]                   # numeric focus for the solver)

# define the name of the modelrun
OBJ_STR = OBJ_STR_INPUT * "_iteraiton" * string(START_ITERATION) * "-" * string(END_ITERATION) * "_oil" * string(MIN_BIOMASS_FOR_OIL) * "_hvc" * string(MIN_BIOMASS_FOR_HVC) * "_syngas" * string(MIN_BIOMASS_FOR_SYNGAS) * "_cluster" * string(CLUSTER_INDEX)

# define input and output directories
inputMod_arr = ["./_basis","./timeSeries/96hours_2008"]
resultDir_str = "./results/" * OBJ_STR * Dates.format(now(), "_yyyy-mm-dd_HH-MM")
mkdir(resultDir_str)

# read uncertain parameters
uncertain_parameters = CSV.read("./_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0) # creates a new Latin Hypercube Sample and saves as csv
lhs = CSV.read("./lhs.csv", DataFrame, header=false) # uses an already created lhc sample

# initialize model (deepcopy of initialized model might lead to a a terminal crash. In this case initialize the model separately for every iteration.)
anyM0 = anyModel(inputMod_arr, resultDir_str, objName = OBJ_STR, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true)

outputsOfInterest = nothing # initialize parameter to store the outputs of interest
objective = DataFrame(iteration = Int[], objectiveValue = Float64[]) # initialize dataframe to store the objective values of the optimization runs

# set additional constraint for regret calculations -> Set the minimum use of biomass for the following categories (biomass usage in GWh/a)
# make sure that the file "df_input_with_final_cluster.csv" is available in EuSysMod. The file is created by the python script in Decide.
set_constraint = (MIN_BIOMASS_FOR_HVC > 0) || (MIN_BIOMASS_FOR_OIL > 0) || (MIN_BIOMASS_FOR_SYNGAS > 0) # set to true if additional constraints exists

if set_constraint
    new_constraint = Dict(
        "bioConversionOil" => MIN_BIOMASS_FOR_OIL*1000,
        "bioConversionSyngas" => MIN_BIOMASS_FOR_SYNGAS*1000, 
        "bioConversionBiogas" => 0, 
        "bioConversionChp" => 0,
        "bioConversionHvc" => MIN_BIOMASS_FOR_HVC*1000,
        "bioConversionCoal" => 0,
        "networkHeat" => 0,
        "spaceHeat" => 0,
        "proHeat" => 0
    )

    for row in eachrow(anyM0.parts.lim.par[:useLow].data)
        row.val = get(new_constraint, findTechnology(anyM0, row.Te), row.val)
    end

    #determine if a scenario is already within the respective cluster and should be skipped in the optimization
    cluster_of_scenarios = CSV.read("df_filtered_with_cluster.csv", DataFrame)[!, "cluster_final"] # make sure that the file "df_filtered_with_cluster.csv" is available in EuSysMod. The file is created by the python script in Decide.
end


for iteration in range(START_ITERATION, END_ITERATION)
    start_time_iteration = @elapsed begin
        println("[info] iteration $iteration/$END_ITERATION")

        if set_constraint && cluster_of_scenarios[iteration] == CLUSTER_INDEX
            continue
        end

        anyM = deepcopy(anyM0) # deepcopy of initialized model might lead to a a terminal crash => initialize the model for every iteration
        # anyM = anyModel(inputMod_arr, resultDir_str, objName = OBJ_STR, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true) 

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
        set_optimizer_attribute(anyM.optModel, "BarConvTol", BARRIER_CONV_TOL); 
        set_optimizer_attribute(anyM.optModel, "NumericFocus", NUMERIC_FOCUS); 
        optimize!(anyM.optModel) # solve the model
        
    
        # update the outputs of interest and the objective value dataframes
        global outputsOfInterest
        global objective
        outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
        push!(objective, [iteration, objective_value(anyM.optModel)])
        
        # Save outputsOfInterest and total_cost as a CSV file
        CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
        set_constraint ? CSV.write(resultDir_str * "/total_cost_cluster_$CLUSTER_INDEX.csv", DataFrame(objective)) : CSV.write(resultDir_str * "/total_cost.csv", DataFrame(objective))
    end
    println("Iteration $iteration took $start_time_iteration seconds")
end
