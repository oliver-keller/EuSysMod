using Random
using DelimitedFiles

function create_lhs(d::Int, N::Int, seed::Int)
    """
    Create a Latin Hypercube Sampling (LHS) of parameter trajectories.

    Parameters:
    - d (int): Number of dimensions.
    - N (int): Number of samples.
    - seed (int): Seed for random number generation.

    Returns:
    - parameter_trajectories (Matrix{Float64}): Array of parameter trajectories.
    """

    # Set the seed for reproducibility
    Random.seed!(seed)

    # Generate Latin Hypercube with N samples and d dimensions
    lhs = Matrix{Float64}(undef, N, d)
    for j in 1:d
        x = collect(1:N)
        shuffle!(x)
        for i in 1:N
            lhs[i, j] = (x[i] - rand()) / N
        end
    end

    # Save LHS to file (optional)
    writedlm("lhs.csv", lhs, ',')

    return lhs
end

"""
# Example usage:
d = 3
N = 10
seed = 1
parameter_trajectories = create_lhs(d, N, seed)
"""


function findRegion(model, region)
    """
    Find the index of a region or returns the name of a region given its index.
    
    Parameters:
    - region: Either a String representing the name of the region or an Int representing the index of the region.
    
    Returns:
    - If 'region' is a String, returns the index of the region.
    - If 'region' is an Int, returns the value of the region given its index.
    - If 'region' is neither a String nor an Int, returns "Invalid input".
    """
    if typeof(region) == String
        for i in 1:length(model.sets[:R].nodes)
            if model.sets[:R].nodes[i].val == region
                return i
            end
        end
    elseif typeof(region) == Int
        return model.sets[:R].nodes[region].val
    else
        return "Invalid input"
    end
end


function findCarrier(model, carrier)
    """
    Find the index of a carrier or returns the name of a carrier given its index.
    
    Parameters:
    - carrier: Either a String representing the name of the carrier or an Int representing the index of the carrier.
    
    Returns:
    - If 'carrier' is a String, returns the index of the carrier.
    - If 'carrier' is an Int, returns the value of the carrier given its index.
    - If 'carrier' is neither a String nor an Int, returns "Invalid input".
    """
    if typeof(carrier) == String
        for i in 1:length(model.sets[:C].nodes)
            if model.sets[:C].nodes[i].val == carrier
                return i
            end
        end
    elseif typeof(carrier) == Int
        return model.sets[:C].nodes[carrier].val
    else
        return "Invalid input"
    end
end


function findTechnology(model, technology)
    """
    Find the index of a technology or returns the name of a technology given its index.
    
    Parameters:
    - technology: Either a String representing the name of the technology or an Int representing the index of the technology.
    
    Returns:
    - If 'technology' is a String, returns the index of the technology.
    - If 'technology' is an Int, returns the value of the technology given its index.
    - If 'technology' is neither a String nor an Int, returns "Invalid input".
    """
    if typeof(technology) == String
        for i in 1:length(model.sets[:Te].nodes)
            if model.sets[:Te].nodes[i].val == technology
                return i
            end
        end
    elseif typeof(technology) == Int
        return model.sets[:Te].nodes[technology].val
    else
        return "Invalid input"
    end
end


function findTimestep(model, timestep)
    """
    Find the index of a timestep or returns the value of a timestep given its index.
    
    Parameters:
    - timestep: Either a String representing the name of the timestep or an Int representing the index of the timestep.
    
    Returns:
    - If 'timestep' is a String, returns the index of the timestep.
    - If 'timestep' is an Int, returns the value of the timestep given its index.
    - If 'timestep' is neither a String nor an Int, returns "Invalid input".
    """
    if typeof(timestep) == String
        for i in 1:length(model.sets[:Ts].nodes)
            if model.sets[:Ts].nodes[i].val == timestep
                return i
            end
        end
    elseif typeof(timestep) == Int
        return model.sets[:Ts].nodes[timestep].val
    else
        return "Invalid input"
    end
end


function children(model, type="technology", node=1)
    """
    Returns the children of a given node in the model tree.
    
    Parameters:
    - `model`: The model object.
    - `type`: The type of the tree to search for children. Default is 'technology'. Possible values are 'technology', 'region', 'carrier', and 'timestep'.
    - `node`: The node to find the children of. Default is 1.
    """
    if type == "technology"
        category = :Te
        if typeof(node) == String node = findTechnology(model, node) end
    elseif type == "region"
        category = :R
        if typeof(node) == String node = findRegion(model, node) end
    elseif type == "carrier"
        category = :C
        if typeof(node) == String node = findCarrier(model, node) end
    elseif type == "timestep"
        category = :Ts
        if typeof(node) == String node = findTimestep(model, node) end
    else
        return []
    end

    children = model.sets[category].nodes[node].down
    for _ in range(1, model.sets[category].height)
        for item in children children = [children; model.sets[category].nodes[item].down] end
    end

    return unique(children)
end



function scaleBiomassPotential(model; factor=nothing, newValue=nothing, carrier="all", region = "all")
    """
    Scale the biomass potential values in the given model by a factor or sets a new value.
    
    Parameters:
    - `model`: The model object.
    - `factor`: The scaling factor to apply to the biomass potential value(s). Default is `nothing`.
    - `newValue`: The new value to set for the biomass potential value(s). Default is `nothing`.
    - `carrier`: The carrier to scale the biomass potential values for. Default is "all".
    - `region`: The region to scale the biomass potential values for. Default is "all".
    
    If `region` is set to "all", the biomass potential values for all regions will be scaled.
    Otherwise, only the biomass potential values for the specified region will be scaled.
    """


    if typeof(region) == String # returns the index of the region if entered as sting
        if region == "all"
            region = 0
        else
            region = findRegion(model, region) 
        end
    end
    if typeof(carrier) == String # returns the index of the region if entered as sting
        if carrier == "all"
            carrier = findCarrier(model, "biomass")
        else
            carrier = findCarrier(model, carrier) 
        end
    end

    for row in eachrow(model.parts.lim.par[:trdBuyUp].data)
        if (row.R_dis in [region; children(model, "region", region)] || region == 0) && (row.C in [carrier; children(model, "carrier", carrier)] || carrier == 0)
            factor != nothing ? row.val *= factor : row.val = row.val
            newValue != nothing ? row.val = newValue : row.val = row.val
        end
    end
end


function scaleInvestmentCost(model; factor=nothing, newValue=nothing, technology="all")
    """
    Scale the investment cost of a given technology in the model.

    Arguments:
    - `model`: The model object.
    - `factor`: The scaling factor to multiply the investment cost by. Default is `nothing`.
    - `newValue`: The new value to set the investment cost to. Default is `nothing`.
    - `technology`: The technology to scale the investment cost for. Default is `"all"`.

    Returns:
    - Nothing.

    Examples:
    ```julia
    scaleInvestmentCost(model, factor=2.0, technology="heatpumpProLow")
    scaleInvestmentCost(model, newValue=1000.0, technology="chpGreenWasteProLow")
    ```
    """
    if typeof(technology) == String && technology != "all" technology = findTechnology(model, technology) end
    
    if typeof(technology) == String # returns the index of the technology if entered as sting
        if technology == "all"
            technology = 0
        else
            technology = findTechnology(model, technology, factor) 
        end
    end    
    for row in eachrow(model.parts.cost.par[:costExpConv].data)
        if row.Te in [technology; children(model, "technology", technology)]|| technology == 0
            factor != nothing ? row.val *= factor : row.val = row.val
            newValue != nothing ? row.val = newValue : row.val = row.val
        end
    end
end



function scaleRenewablePotential(model; factor=nothing, newValue=nothing, technology="all", region = "all")
    """
    Scale the renewable potential values in the given model by a factor or sets a new value.
    
    Parameters:
    - `model`: The model object.
    - `factor`: The scaling factor to apply to the biomass potential value(s). Default is `nothing`.
    - `newValue`: The new value to set for the biomass potential value(s). Default is `nothing`.
    - `technology`: The technology to scale the potential values for. Possible values are 'wind', 'solar' or 'all'. Default is "all".
    - `region`: The region to scale the renewable potential values for. Default is "all".
    
    If `region` is set to "all", the renewable potential values for all regions will be scaled.
    Otherwise, only the biomass potential values for the specified region will be scaled.
    Analogue for 'technology'.

    Caution: This function assumes that the parameter 'trdBuyUp' is only used for biomass potential.
    """
    if typeof(region) == String # returns the index of the region if entered as sting
        if region == "all"
            region = 0
        else
            region = findRegion(model, region) 
        end
    end
    if typeof(technology) == String # returns the index of the technology if entered as sting
        if technology == "all"
            technology = [findTechnology(model, "wind"), findTechnology(model, "solar")]
        else
            technology = findTechnology(model, technology) 
        end
    end
    for tec in technology
        for row in eachrow(model.parts.lim.par[:capaConvUp].data)
            if (row.R_exp in [region; children(model, "region", region)] || region == 0) && (row.Te in [tec; children(model, "technology", tec)])
                factor != nothing ? row.val *= factor : row.val = row.val
                newValue != nothing ? row.val = newValue : row.val = row.val
               
            end
        end
    end
end


function calculateOutputs(model;  iteration=1, outputsOfInterest=nothing, includeWaste=false, resultDir)
    """
    This function calculates the outputs of interest based on the given model and updates the `outputsOfInterest` DataFrame.

    ## Arguments
    - `model::AnyModel`: The model used for calculations.
    - `iteration::Int`: The iteration number. Default is 1.
    - `outputsOfInterest::DataFrame`: The DataFrame to store the calculated outputs. If not provided, a new DataFrame will be created.
    - `includeWaste::Bool`: Whether to include waste carriers in the calculations. Default is false.

    ## Returns
    - `outputsOfInterest::DataFrame`: The updated DataFrame with calculated outputs.
    """
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

    #println(technology_input)
    # Save technology_input to CSV
    csv_file = "biomassUsage_iteration_$iteration.csv"
    CSV.write(joinpath(resultDir, csv_file), technology_input)


    if outputsOfInterest === nothing #!@isdefined(outputsOfInterest)
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
    return outputsOfInterest
end


