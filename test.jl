############### Test Version ####################


## Define the input parameter ##
START_ITERATION = 1               # set iteration boundaries
END_ITERATION = 1                 # set iteration boundaries
OBJ_STR_INPUT = "test"                # name of the modelrun
MIN_BIOMASS_FOR_OIL = 0         # constraint for using certain amount of biomass for oil produiction (TWh/a)
MIN_BIOMASS_FOR_HVC = 0         # constraint for using certain amount of biomass for HVC produiction (TWh/a)
MIN_BIOMASS_FOR_SYNGAS = 0      # constraint for using certain amount of biomass for syngas produiction (TWh/a)
CLUSTER_INDEX = 0               # cluster index for the additional constraint (if no constraint parameter is not used) 
BARRIER_CONV_TOL = 1e-3         # barrier convergence tolerance
NUMERIC_FOCUS = 0               # numeric focus for the solver


using AnyMOD, Gurobi, CSV, Statistics, Dates
include("support_functions.jl")

t_int = 8 # number of threads

# define the name of the modelrun
OBJ_STR = OBJ_STR_INPUT * "_iteraiton" * string(START_ITERATION) * "-" * string(END_ITERATION) * "_oil" * string(MIN_BIOMASS_FOR_OIL) * "_hvc" * string(MIN_BIOMASS_FOR_HVC) * "_syngas" * string(MIN_BIOMASS_FOR_SYNGAS) * "_cluster" * string(CLUSTER_INDEX)

# define input and output directories
inputMod_arr = ["./_basis","./timeSeries/96hours_2008"]
resultDir_str = "./results/" * OBJ_STR * Dates.format(now(), "_yyyy-mm-dd_HH-MM")
mkdir(resultDir_str)


# read uncertain parameters
uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0) # creates a new Latin Hypercube Sample and saves as csv
lhs = CSV.read("lhs.csv", DataFrame, header=false) # uses an already created lhc sample


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


anyM = anyModel(inputMod_arr, resultDir_str, objName = OBJ_STR, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true) 

# modify uncertain parameters
iteration = 1
modify_parameters(anyM, uncertain_parameters, lhs, iteration)

# set exchange costs to zero to see if biomass is exchanged between regions
for rows in eachrow(anyM.parts.cost.par[:costVarExc].data)
    rows.val = 0
end

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



##### get output of the model #####


# update the outputs of interest and the objective value dataframes
outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
push!(objective, [iteration, objective_value(anyM.optModel)])

# Save outputsOfInterest and total_cost as a CSV file
CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
set_constraint ? CSV.write(resultDir_str * "/total_cost_cluster_$CLUSTER_INDEX.csv", DataFrame(objective)) : CSV.write(resultDir_str * "/total_cost.csv", DataFrame(objective))


total_cost = objective_value(anyM.optModel) # objective value (in mio Euro)


reportResults(:summary,anyM, addRep = (:capaConvOut,), addObjName = true)
reportResults(:exchange,anyM, addObjName = true)
reportResults(:cost,anyM, addObjName = true)
reportTimeSeries(:electricity,anyM)


plotSankeyDiagram(anyM, dropDown = (:timestep,)) # sankey for the whole europe
# plotSankeyDiagram(anyM) # sankey with dropdown for the regions and contires
plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", dropDown = (:timestep, ), fontSize=16, name="biomass_europe") 
# plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", fontSize=20, name="biomass") 
plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), fontSize=16, name="all_europe") 
# plotSankeyDiagram(anyM; ymlFilter = "all.yml", fontSize=20, name="all") 
# plotTree(:region, anyM)
# plotTree(:carrier, anyM)
# plotTree(:technology, anyM)
# plotTree(:timestep, anyM)
# plotNetworkGraph(anyM) # flow diagramm

