############### Test Version ####################

using AnyMOD, Gurobi, CSV, Statistics, Dates
include("support_functions.jl")

t_int = 8 # Number of threads (4)

inputMod_arr = ["_basis","timeSeries/96hours_2008"]
resultDir_str = "./results/" * Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
if !isdir(resultDir_str)
    mkdir(resultDir_str)
end

obj_str = "Test"

# read uncertain parameters
uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])
# uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame)
# foreach(x -> uncertain_parameters[!,x] = string.(uncertain_parameters[!,x]), [:parameter, :carrier, :region, :technology, :timestep])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0)
lhs = CSV.read("lhs.csv", DataFrame, header=false)

anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true);

iteration = 3

# modify uncertain parameters
modify_parameters(anyM, uncertain_parameters)

# set additional constraint for regret calculations
anyM.parts.lim.par






# solve model
createOptModel!(anyM)
setObjective!(:cost, anyM)
set_optimizer(anyM.optModel, Gurobi.Optimizer) # select a solver
set_optimizer_attribute(anyM.optModel, "Method", 2);
set_optimizer_attribute(anyM.optModel, "Crossover", 0);
set_optimizer_attribute(anyM.optModel, "Threads",t_int);
set_optimizer_attribute(anyM.optModel, "BarConvTol", 1e-3);  # 1e-5
optimize!(anyM.optModel) # solve the model


##### get output of the model #####

outputsOfInterest = nothing
outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false)
outputsOfInterest


reportResults(:summary, anyM)
plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), fontSize=16, name="all_europe") 
plotSankeyDiagram(anyM, dropDown = (:timestep,), name="test")
anyM.parts.tech[:biomassToHvc].cns





# reportResults(:summary,anyM, addRep = (:capaConvOut,), addObjName = true, rtnOpt = (:csvDf,))
# reportResults(:exchange,anyM, addObjName = true)
# reportResults(:cost,anyM, addObjName = true)
# reportTimeSeries(:electricity,anyM)


# plotSankeyDiagram(anyM, dropDown = (:timestep,)) # sankey for the whole europe
# plotSankeyDiagram(anyM) # sankey with dropdown for the regions and contires
# plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", dropDown = (:timestep, ), fontSize=16, name="biomass_europe") 
# plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", fontSize=20, name="biomass") 
# plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), fontSize=16, name="all_europe") 
# plotSankeyDiagram(anyM; ymlFilter = "all.yml", fontSize=20, name="all") 
# plotTree(:region, anyM)
# plotTree(:carrier, anyM)
# plotTree(:technology, anyM)
# plotTree(:timestep, anyM)
# plotNetworkGraph(anyM) # flow diagramm


