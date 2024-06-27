############### Test Version ####################


## Define the input parameter ##
START_ITERATION = 1               # set iteration boundaries
END_ITERATION = 1                 # set iteration boundaries
OBJ_STR_INPUT = "test"                # name of the modelrun
MIN_BIOMASS_FOR_OIL = 0         # constraint for using certain amount of biomass for oil produiction (TWh/a)
MIN_BIOMASS_FOR_HVC = 0         # constraint for using certain amount of biomass for HVC produiction (TWh/a)
MIN_BIOMASS_FOR_SYNGAS = 0      # constraint for using certain amount of biomass for syngas produiction (TWh/a)
CLUSTER_INDEX = 0               # cluster index for the additional constraint (if no constraint parameter is not used) 
BARRIER_CONV_TOL = 1e-8         # barrier convergence tolerance
NUMERIC_FOCUS = 2               # numeric focus for the solver


using AnyMOD, Gurobi, CSV, Statistics, Dates
include("support_functions.jl")

t_int = 8 # number of threads

# define the name of the modelrun
OBJ_STR = OBJ_STR_INPUT * "_iteraiton" * string(START_ITERATION) * "-" * string(END_ITERATION) * "_oil" * string(MIN_BIOMASS_FOR_OIL) * "_hvc" * string(MIN_BIOMASS_FOR_HVC) * "_syngas" * string(MIN_BIOMASS_FOR_SYNGAS) * "_cluster" * string(CLUSTER_INDEX)

# define input and output directories
inputMod_arr = ["./_basis","./timeSeries/96hours_2008"]
# resultDir_str = "./results/" * OBJ_STR * Dates.format(now(), "_yyyy-mm-dd_HH-MM")
# mkdir(resultDir_str)
resultDir_str = "./results/"

# read uncertain parameters
uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0) # creates a new Latin Hypercube Sample and saves as csv
lhs = CSV.read("lhs.csv", DataFrame, header=false) # uses an already created lhc sample


outputsOfInterest = nothing # initialize parameter to store the outputs of interest
objective = DataFrame(iteration = Int[], objectiveValue = Float64[]) # initialize dataframe to store the objective values of the optimization runs

anyM = anyModel(inputMod_arr, resultDir_str, objName = OBJ_STR, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true) 


# set exchange costs to zero to see if biomass is exchanged between regions
# for rows in eachrow(anyM.parts.cost.par[:costVarExc].data)
#     rows.val = 0
# end

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
# reportResults(:exchange,anyM, addObjName = true)
reportResults(:cost,anyM, addObjName = true)
# reportTimeSeries(:electricity,anyM)


plotSankeyDiagram(anyM, dropDown = (:timestep,), minVal = 1.) # sankey for the whole europe
# plotSankeyDiagram(anyM) # sankey with dropdown for the regions and contires
plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", dropDown = (:timestep, ), fontSize=16, name="biomass_europe", minVal = 1.) 
# plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", fontSize=20, name="biomass") 
plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), fontSize=16, name="all_europe", minVal = 1.) 
# plotSankeyDiagram(anyM; ymlFilter = "all.yml", fontSize=20, name="all") 
# plotTree(:region, anyM)
# plotTree(:carrier, anyM)
# plotTree(:technology, anyM)
# plotTree(:timestep, anyM)
# plotNetworkGraph(anyM) # flow diagramm

 





df = reportResults(:summary, model, rtnOpt = (:csvDf,))
use_variables = filter(row -> row.variable == :use, df)


# List of all carriers and technologies
includeWaste ? carriers = ["wood", "greenWaste", "manure", "sludge", "waste", "digestate"] : carriers = ["wood", "greenWaste", "manure", "sludge", "digestate"]
technologies = ["pyrolysisOil", "biomassToHvc", "chp", "boilerDh", "boilerSpace", "boilerProLow", "boilerProMed", "boilerProHigh", "biochemicalWoodOil", "liquefaction", "digestion", "gasification", "carbonization", "pyrolysisCoal"]

# Initialize the DataFrame with all zeros
technology_input = DataFrame([:carrier => carriers; Symbol.(technologies) .=> 0.0])

# Iterate over use_variables and update technology_input
for row in eachrow(use_variables)
    for tech in technologies
        if occursin(tech[2:end], row.technology)
            # carrier_index = findfirst(==(row.carrier), technology_input.carrier)
            carrier_index = findfirst(x -> occursin(x[2:end], row.carrier), technology_input.carrier)
            if !isnothing(carrier_index)
                technology_input[carrier_index, tech] += row.value
            end
        end
    end
end

# Save technology_input to CSV
csv_file = "biomassUsage_iteration_$iteration.csv"
CSV.write(joinpath(resultDir, csv_file), technology_input)


if outputsOfInterest === nothing
    outputsOfInterest = DataFrame(
        iteration = [iteration],
        crudeOil = [sum(technology_input.pyrolysisOil)+sum(technology_input.biochemicalWoodOil)+sum(technology_input.liquefaction)],
        HVC = [sum(technology_input.biomassToHvc)],
        CHP = [sum(technology_input.chp)],
        LTH = [sum(technology_input.boilerProLow)],
        MTH = [sum(technology_input.boilerProMed)],
        HTH = [sum(technology_input.boilerProHigh)],
        DH = [sum(technology_input.boilerDh)],
        SpaceHeating = [sum(technology_input.boilerSpace)],
        rawBiogas = [sum(technology_input.digestion)-sum(technology_input[technology_input.carrier .== "digestate", :][1, 2:end])],
        syngas = [sum(technology_input.gasification)],
        coal = [sum(technology_input.carbonization)+ sum(technology_input.pyrolysisCoal)]
    )
else
    push!(outputsOfInterest, [iteration, sum(technology_input.pyrolysisOil)+sum(technology_input.biochemicalWoodOil)+sum(technology_input.liquefaction), sum(technology_input.biomassToHvc), sum(technology_input.chp), sum(technology_input.boilerProLow), sum(technology_input.boilerProMed), sum(technology_input.boilerProHigh), sum(technology_input.boilerDh), sum(technology_input.boilerSpace), sum(technology_input.digestion)-sum(technology_input[technology_input.carrier .== "digestate", :][1, 2:end]), sum(technology_input.gasification), sum(technology_input.carbonization)+ sum(technology_input.pyrolysisCoal)])
end

