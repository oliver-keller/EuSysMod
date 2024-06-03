using AnyMOD, Gurobi, CSV, Statistics, Dates
include("support_functions.jl")

t_int = 8 # 4

inputMod_arr = ["_basis","timeSeries/96hours_2008"]
resultDir_str = "./results/" * Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
if !isdir(resultDir_str)
    mkdir(resultDir_str)
end
obj_str = "Sampling_1000"



# read uncertain parameters
uncertain_parameters = CSV.read("_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])

lhs = CSV.read("lhs.csv", DataFrame, header=false)

# set iteration counter
iterations = 300
outputsOfInterest = nothing

for iteration in range(200, iterations)
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
    global outputsOfInterest
    outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
        
    # Save outputsOfInterest as a CSV file
    CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
end

