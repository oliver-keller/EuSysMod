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


t_int = 8 # 4

inputMod_arr = ["_basis","timeSeries/96hours_2008"]
resultDir_str = "./results/" * Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
if !isdir(resultDir_str)
    mkdir(resultDir_str)
end
obj_str = "Sampling_1000"


#region # * uncertain parameter definition and create and solve model

# read uncertain parameters
uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])
# uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame)
# foreach(x -> uncertain_parameters[!,x] = string.(uncertain_parameters[!,x]), [:parameter, :carrier, :region, :technology, :timestep])

# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0)
lhs = CSV.read("lhs.csv", DataFrame, header=false)

# set iteration counter
iterations = 100

# initialize model
# anyM0 = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true);
outputsOfInterest = nothing

for iteration in range(1, iterations)
    println("[info] iteration $iteration/$iterations")

    anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true);
    # anyM = deepcopy(anyM0)

    # modify uncertain parameters
    for (index, row) in enumerate(eachrow(uncertain_parameters))
        ismissing(row[:minRel]) || ismissing(row[:maxRel]) ? factor = nothing : factor = row[:minRel]  + (row[:maxRel] - row[:minRel]) * lhs[iteration, index] 
        ismissing(row[:minAbs]) || ismissing(row[:maxAbs]) ? newValue = nothing : newValue = row[:minAbs] + (row[:maxAbs] - row[:minAbs]) * lhs[iteration, index]

        if row[:parameter] == "trdBuyUp"
            scaleBiomassPotential(anyM, factor=factor, newValue=newValue, carrier=row[:carrier], region=row[:region])
        elseif row[:parameter] == "capaConvUp"
            scaleRenewablePotential(anyM, factor=factor, newValue=newValue, technology=row[:technology], region=row[:region])
        elseif row[:parameter] == "costExpConv"
            scaleInvestmentCost(anyM, factor=factor, newValue=newValue, technology=row[:technology])
        end
    end

    createOptModel!(anyM)
    setObjective!(:cost, anyM)

    # solve model
    set_optimizer(anyM.optModel, Gurobi.Optimizer)  # select a solver
    set_optimizer_attribute(anyM.optModel, "Method", 2);
    set_optimizer_attribute(anyM.optModel, "Crossover", 0);
    set_optimizer_attribute(anyM.optModel, "Threads",t_int);
    set_optimizer_attribute(anyM.optModel, "BarConvTol", 1e-3);  # 1e-5

    optimize!(anyM.optModel) # solve the model
    # plotSankeyDiagram(anyM, dropDown = (:timestep,), name="iteration_" * string(iteration))
    if iteration == 1
        outputsOfInterest = calculateOutputs(anyM, iteration=iteration, includeWaste=false, resultDir = resultDir_str)
    else
        outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
    end
    
    # Save outputsOfInterest as a CSV file
    CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
end


#endregion

#region # * write results

# reportResults(:summary,anyM, addRep = (:capaConvOut,), addObjName = true, rtnOpt = (:csvDf,))
# reportResults(:exchange,anyM, addObjName = true)
# reportResults(:cost,anyM, addObjName = true)
# reportTimeSeries(:electricity,anyM)

# create plots
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

#endregion

#anyM.parts.tech[:bevPsngRoadPrvt].var







############### Test Version ####################

anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true);

anyM1 = deepcopy(anyM0)
iteration = 3

# modify uncertain parameters
for (index, row) in enumerate(eachrow(uncertain_parameters))
    ismissing(row[:minRel]) || ismissing(row[:maxRel]) ? factor = nothing : factor = row[:minRel]  + (row[:maxRel] - row[:minRel]) * lhs[iteration, index] 
    ismissing(row[:minAbs]) || ismissing(row[:maxAbs]) ? newValue = nothing : newValue = row[:minAbs] + (row[:maxAbs] - row[:minAbs]) * lhs[iteration, index]

    if row[:parameter] == "trdBuyUp"
        scaleBiomassPotential(anyM, factor=factor, newValue=newValue, carrier=row[:carrier], region=row[:region])
    elseif row[:parameter] == "capaConvUp"
        scaleRenewablePotential(anyM, factor=factor, newValue=newValue, technology=row[:technology], region=row[:region])
    elseif row[:parameter] == "costExpConv"
        scaleInvestmentCost(anyM, factor=factor, newValue=newValue, technology=row[:technology])
    end
end

createOptModel!(anyM)
setObjective!(:cost, anyM)


# solve model
set_optimizer(anyM.optModel, Gurobi.Optimizer) # select a solver
set_optimizer_attribute(anyM.optModel, "Method", 2);
set_optimizer_attribute(anyM.optModel, "Crossover", 0);
set_optimizer_attribute(anyM.optModel, "Threads",t_int);
set_optimizer_attribute(anyM.optModel, "BarConvTol", 1e-3);  # 1e-5

optimize!(anyM.optModel) # solve the model


##### get output of the model #####

# getAllVariables(:use, anyM, x -> x.C in (4,5,6,3)) # I didnt find the parameter values in the output of this function ...


if iteration == 1
    outputsOfInterest = calculateOutputs(anyM, iteration=iteration, includeWaste=false)
else
    outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false)
end

outputsOfInterest

reportResults(:summary, anyM)
plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), fontSize=16, name="all_europe") 
plotSankeyDiagram(anyM, dropDown = (:timestep,), name="test")
anyM.parts.tech[:biomassToHvc].cns
