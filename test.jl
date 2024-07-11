############### Test Version ####################


## Define the input parameter ##
iteration = 800 # 164, 338, 800              # set iteration
OBJ_STR_INPUT = "test"                # name of the modelrun
MIN_BIOMASS_FOR_OIL = 0      # constraint for using certain amount of biomass for oil produiction (TWh/a)
MAX_BIOMASS_FOR_OIL = 9999       # constraint for using certain amount of biomass for oil produiction (TWh/a)
MIN_BIOMASS_FOR_HVC = 0       # constraint for using certain amount of biomass for HVC produiction (TWh/a)
MAX_BIOMASS_FOR_HVC = 9999      # constraint for using certain amount of biomass for HVC produiction (TWh/a)
CLUSTER_INDEX = 0               # cluster index for the additional constraint (if no constraint parameter is not used) 
BARRIER_CONV_TOL = 1e-6        # barrier convergence tolerance


using Gurobi, AnyMOD, CSV, Statistics, Dates
include("support_functions.jl")

t_int = 8 # number of threads


# define the name of the modelrun
OBJ_STR = OBJ_STR_INPUT * "_iteraiton" * string(iteration) * "_oil" * string(MIN_BIOMASS_FOR_OIL) * "_hvc" * string(MIN_BIOMASS_FOR_HVC) * "_cluster" * string(CLUSTER_INDEX)

# define input and output directories
inputMod_arr = ["./_basis","./timeSeries/96hours_2008"]
resultDir_str = "C:/Users/olive/Desktop/Projects/EuSysMod/results/V1.2_median/" * string(iteration)
#resultDir_str = "./results/"
#mkdir(resultDir_str)


# read uncertain parameters
uncertain_parameters = CSV.read("./_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0) # creates a new Latin Hypercube Sample and saves as csv
lhs = CSV.read("./lhs.csv", DataFrame, header=false) # uses an already created lhc sample


outputsOfInterest = nothing # initialize parameter to store the outputs of interest
objective = DataFrame(iteration = Int[], objectiveValue = Float64[]) # initialize dataframe to store the objective values of the optimization runs


obj_str_model = string(iteration) * "_" * OBJ_STR

# anyM = deepcopy(anyM0) # deepcopy of initialized model might lead to a a terminal crash => initialize the model for every iteration
anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str_model, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true) 


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
set_optimizer_attribute(anyM.optModel, "BarHomogeneous", 1)
set_optimizer_attribute(anyM.optModel, "PreDual", -1)
set_optimizer_attribute(anyM.optModel, "Presolve", 2) 
set_optimizer_attribute(anyM.optModel, "BarOrder", 0) 


numFoc_int = 0
while true
    set_optimizer_attribute(anyM.optModel, "NumericFocus", numFoc_int); 
    println("[info] NumericFocus: $numFoc_int")
    optimize!(anyM.optModel) # solve the model
    println("[info] Termination status: ", termination_status(anyM.optModel))
    if termination_status(anyM.optModel) == MOI.OPTIMAL || numFoc_int == 3
        break
    else
        numFoc_int = numFoc_int + 1
    end
end

# update the outputs of interest and the objective value dataframes
outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
push!(objective, [iteration, objective_value(anyM.optModel)])

# Save outputsOfInterest and total_cost as a CSV file
CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))





##### get output of the model #####


# update the outputs of interest and the objective value dataframes
# outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
# push!(objective, [iteration, objective_value(anyM.optModel)])

# Save outputsOfInterest and total_cost as a CSV file
# CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
# set_constraint ? CSV.write(resultDir_str * "/total_cost_cluster_$CLUSTER_INDEX.csv", DataFrame(objective)) : CSV.write(resultDir_str * "/total_cost.csv", DataFrame(objective))


total_cost = objective_value(anyM.optModel) # objective value (in mio Euro)


reportResults(:summary,anyM, addRep = (:capaConvOut,), addObjName = true)
# reportResults(:exchange,anyM, addObjName = true)
reportResults(:cost,anyM, addObjName = true)
# reportTimeSeries(:electricity,anyM)


plotSankeyDiagram(anyM, dropDown = (:timestep,), minVal = 1.) # sankey for the whole europe
# plotSankeyDiagram(anyM) # sankey with dropdown for the regions and contires
plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", dropDown = (:timestep, ), fontSize=16, name="biomass_europe", minVal = 1.) 
# plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", fontSize=20, name="biomass") 
plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), fontSize=16, name="all_europe", minVal = 1. ) 
# plotSankeyDiagram(anyM; ymlFilter = "all.yml", fontSize=20, name="all") 
# plotTree(:region, anyM)
# plotTree(:carrier, anyM)
# plotTree(:technology, anyM)
# plotTree(:timestep, anyM)
# plotNetworkGraph(anyM) # flow diagramm

 





df = reportResults(:summary, anyM, rtnOpt = (:csvDf,))
emission_variables = filter(row -> row.variable == :emission, df)
totEmissions = sum(emission_variables.value)
println("Total emissions: $totEmissions tCO2/a")
 