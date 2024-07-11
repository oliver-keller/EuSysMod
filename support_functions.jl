### File containing all the support functions called from the main script called "run.jl" ###

using Random, DelimitedFiles

function create_lhs(d::Int, N::Int, seed::Int=0) ::Matrix{Float64}
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

    # Save LHS to file
    writedlm("lhs.csv", lhs, ',')

    return lhs
end


function findRegion(model::anyModel, region::Union{String, Int}) ::Union{Int, String}
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


function findCarrier(model::anyModel, carrier::Union{String, Int}) ::Union{Int, String}
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


function findTechnology(model::anyModel, technology::Union{Int, String}) ::Union{Int, String}
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


function findTimestep(model::anyModel, timestep::Union{Int, String}) ::Union{Int, String}
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


function children(model::anyModel, type::String="technology", node::Union{String, Int}=1) ::Array{Int, 1}
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



function scaleBiomassPotential(model::anyModel; factor::Union{Float64, Nothing}=nothing, newValue::Union{Float64, Int64, Nothing}=nothing, carrier::Union{String, Int}="all", region::Union{String, Int}="all") ::Nothing 
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





function scaleRenewablePotential(model::anyModel; factor::Union{Float64, Nothing}=nothing, newValue::Union{Float64, Int64, Nothing}=nothing, technology::Union{String, Int}="all", region::Union{String, Int}="all") ::Nothing
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


function calculateOutputs(model::anyModel;  iteration::Int64=1, outputsOfInterest::Union{DataFrame, Nothing}=nothing, includeWaste::Bool=false, resultDir::Union{String, Nothing}=nothing) ::DataFrame
    """
    This function calculates the outputs of interest based on the given model and updates the `outputsOfInterest` DataFrame.
    Furthermore, it saves the calculated outputs to a CSV file. 

    ## Arguments
    - `model::AnyModel`: The model used for calculations.
    - `iteration::Int`: The iteration number. Default is 1.
    - `outputsOfInterest::DataFrame`: The DataFrame to store the calculated outputs. If not provided, a new DataFrame will be created.
    - `includeWaste::Bool`: Whether to include waste carriers in the calculations. Default is false.
    - `resultDir::Union{String, Nothing}`: The directory to save the results. If not provided, the results will be saved in "./results/".

    ## Returns
    - `outputsOfInterest::DataFrame`: The updated DataFrame with calculated outputs.
    """
    if resultDir === nothing
        resultDir = "./results/"
        if !isdir(resultDir)
            mkdir(resultDir)
        end
    end

    df = reportResults(:summary, model, rtnOpt = (:csvDf,))
    use_variables = filter(row -> row.variable == :use, df)


    # List of all carriers and technologies
    includeWaste ? carriers = ["wood", "greenWaste", "manure", "sludge", "waste", "digestate"] : carriers = ["wood", "greenWaste", "manure", "sludge", "digestate"]
    technologies = ["pyrolysisOil", "biomassToHvc", "chp", "boilerDh", "boilerSpace", "boilerProLow", "boilerProMed", "boilerProHigh", "liquefaction", "digestion", "gasification", "carbonization", "pyrolysisCoal", "CCC"] # the comparison is made from the second to the last character. The first is omitted due to capitalization missmatch
    # CCC means Carbon Capture (second, to last char) ans is not an additional technology. This column is created that the fraction of biomass used for BECCS can be quantified

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
            crudeOil = [sum(technology_input.pyrolysisOil)+sum(technology_input.liquefaction)],
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
        push!(outputsOfInterest, [iteration, sum(technology_input.pyrolysisOil)+sum(technology_input.liquefaction), sum(technology_input.biomassToHvc), sum(technology_input.chp), sum(technology_input.boilerProLow), sum(technology_input.boilerProMed), sum(technology_input.boilerProHigh), sum(technology_input.boilerDh), sum(technology_input.boilerSpace), sum(technology_input.digestion)-sum(technology_input[technology_input.carrier .== "digestate", :][1, 2:end]), sum(technology_input.gasification), sum(technology_input.carbonization)+ sum(technology_input.pyrolysisCoal)])
    end


    ### same to find the usage of the created oil/syngas ###
    # List of all carriers and technologies
    carriers = ["crudeOil", "refinedOil", "syngas", "methane", "h2"]
    technologies = ["refinery", "oilPlant", "dieselEnginePeak", "dieselFrtRail", "dieselFrtRoadHeavy", "dieselFrtRoadLight", "dieselPsngRail", "dieselPsngRoadPub", "ottoPsngRoadPrvt", "avaAndNaviTech",
        "oilBoilerDh", "oilBoilerSpace", "oilBoilerProHigh", "oilBoilerProLow", "oilBoilerProMed",
        "Methanation", "fischerTropsch", "waterGasShift", "chpSyngasDh", "chpSyngasProLow", "syngasEngineDh", "syngasEngineProLow", "syngasEngineProMed",
        "pyrMethaneH2", "ccgtGasDh", "ccgtGasProMed", "ocgtGas", "ocgtGasChp", "methaneEngineDh", "methaneEngineProLow", "methaneEngineProMed", "methaneEnginePeak", "cngPsngRoadPrvt", "sofc", "gasToHvc", "stemMethaneReforming", "methaneBoilerDh", "methaneBoilerSpace", "methaneBoilerProHigh", "methaneBoilerProLow", "methaneBoilerProMed",
        "biogasUpgrade", "ccgtH2Dh", "ccgtH2ProMed", "ocgtH2", "fcevPsngRoadPrvt", "fcFrtRail", "fcFrtRoadHeavy", "fcFrtRoadLight", "fcPsngRail", "fcPsngRoadPub", "h2ToCrudeOil", "h2ToMethane", "pefc", "haberBoschH2AmmoniaDh", "haberBoschH2AmmoniaProLow", "h2BoilerDh", "h2BoilerProHigh", "h2BoilerProLow", "h2BoilerProMed"]
        


    # Initialize the DataFrame with all zeros
    technology_input = DataFrame([:carrier => carriers; Symbol.(technologies) .=> 0.0])

    # Iterate over use_variables and update technology_input
    for row in eachrow(use_variables)
        for tech in technologies
            if occursin(tech[2:end], row.technology)
                carrier_index = findfirst(x -> occursin(x[2:end], row.carrier), technology_input.carrier)
                if !isnothing(carrier_index)
                    technology_input[carrier_index, tech] += row.value
                end
            end
        end
    end

    # Save technology_input to CSV
    csv_file = "oilSyngas_usage_iteration_$iteration.csv"
    CSV.write(joinpath(resultDir, csv_file), technology_input)

    return outputsOfInterest
end




function modify_parameters(model::anyModel, uncertain_parameters::DataFrame, lhs::DataFrame, iteration::Int) ::Nothing
    """
    Modify the parameters of a model based on uncertain parameters.

    # Arguments
    - `model::anyModel`: The model to modify.
    - `uncertain_parameters::DataFrame`: A table of uncertain parameters.
    - `lhs::DataFrame`: A Latin Hypercube Sample.
    - `iteration::Int`: The iteration number.
    """
    for (index, row) in enumerate(eachrow(uncertain_parameters))
        ismissing(row[:minRel]) || ismissing(row[:maxRel]) ? factor = nothing : factor = row[:minRel]  + (row[:maxRel] - row[:minRel]) * lhs[iteration, index] 
        ismissing(row[:minAbs]) || ismissing(row[:maxAbs]) ? newValue = nothing : newValue = row[:minAbs] + (row[:maxAbs] - row[:minAbs]) * lhs[iteration, index]

        if row[:parameter] == "trdBuyUp"
            scaleBiomassPotential(model, factor=factor, newValue=newValue, carrier=row[:carrier], region=row[:region])
        elseif row[:parameter] == "capaConvUp"
            scaleRenewablePotential(model, factor=factor, newValue=newValue, technology=row[:technology], region=row[:region])
        elseif row[:parameter] == "costExpConv" || row[:parameter] == "costOprConv"|| row[:parameter] == "costOprStIn" || row[:parameter] == "costExpStSize"
            scaleCost(model, parameter=row[:parameter] , factor=factor, newValue=newValue, technology=row[:technology])
        elseif row[:parameter] == "trdBuyPrc"
            scalePrice(model, factor=factor, newValue=newValue, carrier=row[:carrier])
        end
    end
end


function scaleCost(model::anyModel;  parameter::String ,factor::Union{Float64, Nothing}=nothing, newValue::Union{Float64, Int64, Nothing}=nothing, technology::Union{String, Int}="all") ::Nothing
    """
    Scale the cost of a given technology in the model.

    Parameters:
    - `model`: The model object.
    - `parameter`: The parameter to scale the cost for. Possible values are 'costExpConv', 'costOprConv', "costOprStIn" and "costExpStSize"
    - `factor`: The scaling factor to multiply the cost by. Default is `nothing`.
    - `newValue`: The new value to set the cost to. Default is `nothing`.
    - `technology`: The technology to scale the investment cost for. Default is `"all"`.
    """ 
    if typeof(technology) == String # returns the index of the technology if entered as sting
        if technology == "all"
            technology = 0
        else
            technology = findTechnology(model, technology)
        end
    end    
    for row in eachrow(model.parts.cost.par[Symbol(parameter)].data)
        if row.Te in [technology; children(model, "technology", technology)]|| technology == 0
            factor != nothing ? row.val *= factor : row.val = row.val
            newValue != nothing ? row.val = newValue : row.val = row.val
        end
    end
end


function scalePrice(model::anyModel; factor::Union{Float64, Nothing}=nothing, newValue::Union{Float64, Int64, Nothing}=nothing, carrier::Union{String, Int}="all") ::Nothing
    """
    Scale the price of a given carrier in the model.

    Parameters:
    - `model`: The model object.
    - `factor`: The scaling factor to multiply the investment price by. Default is `nothing`.
    - `newValue`: The new value to set the investment price to. Default is `nothing`.
    - `carrier`: The carrier to scale the investment price for. Default is `"all"`.
    """ 
    if typeof(carrier) == String # returns the index of the carrier if entered as sting
        if carrier == "all"
            carrier = 0
        else
            carrier = findCarrier(model, carrier)
        end
    end    
    for row in eachrow(model.parts.bal.par[:trdBuyPrc].data)
        if row.C in [carrier; children(model, "carrier", carrier)]|| carrier == 0
            factor != nothing ? row.val *= factor : row.val = row.val
            newValue != nothing ? row.val = newValue : row.val = row.val
        end
    end
end
