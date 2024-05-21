using AnyMOD, Gurobi, CSV, Statistics
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

t_int = 8 # 4

inputMod_arr = ["_basis","timeSeries/96hours_2008"]
resultDir_str = "results"

#region # * create and solve main model
objName = "biomass_included"
anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true);

createOptModel!(anyM)
setObjective!(:cost, anyM)

# solve model
set_optimizer(anyM.optModel, Gurobi.Optimizer)  # select a solver
set_optimizer_attribute(anyM.optModel, "Method", 2);
set_optimizer_attribute(anyM.optModel, "Crossover", 0);
set_optimizer_attribute(anyM.optModel, "Threads",t_int);
set_optimizer_attribute(anyM.optModel, "BarConvTol", 1e-3);  # 1e-5

optimize!(anyM.optModel) # solve the model

#endregion

#region # * write results

reportResults(:summary,anyM, addRep = (:capaConvOut,), addObjName = true)
# reportResults(:exchange,anyM, addObjName = true)
reportResults(:cost,anyM, addObjName = true)
# reportTimeSeries(:electricity,anyM)

# create plots
# plotSankeyDiagram(anyM, dropDown = (:timestep,)) # sankey for the whole europe
# plotSankeyDiagram(anyM) # sankey with dropdown for the regions and contires
plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", dropDown = (:timestep, ), digVal=0, fontSize=16, name="biomass_europe") 
plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", digVal=0, fontSize=20, name="biomass") 
plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), digVal=0, fontSize=16, name="all_europe") 
plotSankeyDiagram(anyM; ymlFilter = "all.yml", digVal=0, fontSize=20, name="all") 
# plotTree(:region, anyM)
# plotTree(:carrier, anyM)
# plotTree(:technology, anyM)
# plotTree(:timestep, anyM)
# plotNetworkGraph(anyM) # flow diagramm

#endregion

#anyM.parts.tech[:bevPsngRoadPrvt].var