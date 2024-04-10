using AnyMOD, Gurobi, CSV, Statistics

# par_df = CSV.read("settings.csv",DataFrame)

if isempty(ARGS)
    id_int = 1
    t_int = 4
else
    id_int = parse(Int,ARGS[1])
    t_int = parse(Int,ARGS[2]) # number of threads
end

# h = string(par_df[id_int,:h]) # resolution of time-series for actual solve, can be 96, 1752, 4392, or 8760
# h_heu = string(par_df[id_int,:h_heu]) # resolution of time-series for pre-screening, can be 96, 1752, 4392, or 8760
# grid = string(par_df[id_int,:grid]) # scenario 
h = "96"
h_heu = "96"
grid = "_gridExp"


obj_str = h * "hours_" * h_heu * "hoursHeu" * grid * "_updated"
inputMod_arr = ["_basis",grid,"timeSeries/" * h * "hours_2008"]

resultDir_str = "results"

#region # * create and solve main model


anyM = anyModel(inputMod_arr,resultDir_str, objName = obj_str, supTsLvl = 2, shortExp = 5, redStep = 1.0, emissionLoss = false, holdFixed = true)

createOptModel!(anyM)
setObjective!(:cost,anyM)

set_optimizer(anyM.optModel, Gurobi.Optimizer)  # select a solver
set_optimizer_attribute(anyM.optModel, "Method", 2);
set_optimizer_attribute(anyM.optModel, "Crossover", 0);
set_optimizer_attribute(anyM.optModel, "Threads",t_int);
set_optimizer_attribute(anyM.optModel, "BarConvTol", 1e-5);

optimize!(anyM.optModel) # solve the model

#endregion

#region # * write results

reportResults(:summary,anyM, addRep = (:capaConvOut,), addObjName = true)
reportResults(:exchange,anyM, addObjName = true)
reportResults(:cost,anyM, addObjName = true)

reportTimeSeries(:electricity,anyM)

# create plots
plotSankeyDiagram(anyM) # sankey
plotTree(:region, anyM)
plotTree(:carrier, anyM)
plotTree(:technology, anyM)
plotTree(:timestep, anyM)
#plotNetworkGraph(anyM) # flow diagramm

#endregion