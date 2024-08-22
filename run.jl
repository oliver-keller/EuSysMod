### Low regret strategies for the use of biomass in the European energy transition ###
### Master's Thesis Oliver Keller 2024 ###

# EuSysMod modified to contain the different biomass conversion pathways
# uncertainty is adressed with an additional file "uncertain_parameters" where the uncertainties can be specified
# this script creates a latin hyper cube samplie and varies the uncertain inputs
# the optimization of the model is executed for every input combination created with the LHS


## the functions are defined in the document "support_functions.jl"
## the output is written in the results folder. There a subfolder is created with date/time.
## The output contains the file "outputsOfInterest.csv" which contains the values how much biomass (TWh/a) is used for which purpose. 
## Additionaly for each run a table containing the different biomass carriers and the outputs of interests is saved.

### The file "outputsOfInterest.csv" is used as input for the decision tree creation with the DECIDE tool written in python (the adapted version can be found under https://github.com/oliver-keller/Decide)
### using the output of DECIDE the regret for the different strategies can be calculated. To do so the thresholds for the different strategies have to be defined in the settings.csv file. Furthermore the file "df_filtered_with_clusters.csv", which is an output of the DECIDE code, is needed to only run the scenarios that are not optimal under the respective strategy.


MODELRUN = 1 # id of the modelrun - scecifies the row in the settings.csv file - only needed if the script is run without arguments -> On the server the arguments are passed directly to the script from the runEuler.sh file.

using Gurobi, AnyMOD, CSV, Statistics, Dates
include("support_functions.jl")

par_df = CSV.read("settings.csv", DataFrame)

if isempty(ARGS)
    id_int = MODELRUN # id of the modelrun
    t_int = 4 # number of threads
else
    id_int = parse(Int,ARGS[1]) # id of the modelrun
    t_int = parse(Int,ARGS[2]) # number of threads
end

## Define the input parameter ##
START_ITERATION = par_df[id_int,:START_ITERATION]                               # set iteration boundaries
END_ITERATION = par_df[id_int,:END_ITERATION]                                   # set iteration boundaries
OBJ_STR_INPUT = par_df[id_int,:OBJ_STR]                                         # name of the modelrun
MIN_BIOMASS_FOR_OIL = par_df[id_int,:MIN_BIOMASS_FOR_OIL]                       # constraint for using certain amount of biomass for oil produiction (TWh/a)
MAX_BIOMASS_FOR_OIL = par_df[id_int,:MAX_BIOMASS_FOR_OIL]                       # constraint for using certain amount of biomass for oil produiction (TWh/a)
MIN_BIOMASS_FOR_HVC = par_df[id_int,:MIN_BIOMASS_FOR_HVC]                       # constraint for using certain amount of biomass for HVC produiction (TWh/a)
MAX_BIOMASS_FOR_HVC = par_df[id_int,:MAX_BIOMASS_FOR_HVC]                       # constraint for using certain amount of biomass for HVC produiction (TWh/a)
CLUSTER_INDEX = par_df[id_int,:CLUSTER_INDEX]                                   # cluster index for the additional constraint (if no constraint parameter is not used)
BARRIER_CONV_TOL = par_df[id_int,:BARRIER_CONV_TOL]                             # barrier convergence tolerance
BIOMASS_ALLOWED = par_df[id_int,:BIOMASS_ALLOWED]                               # biomass carrier allowed or disabled
ALL_BIOMASS_FOR_LTH = par_df[id_int,:ALL_BIOMASS_FOR_LTH]                     # biomass is only used for district heating
ONLY_BIOMASS_FOR_HVC_AND_OIL = par_df[id_int,:ONLY_BIOMASS_FOR_HVC_AND_OIL]     # biomass is only used for for hvc and oil


# define the name of the modelrun
OBJ_STR = OBJ_STR_INPUT * "_iteraiton" * string(START_ITERATION) * "-" * string(END_ITERATION)

# define input and output directories
inputMod_arr = ["./_basis","./timeSeries/96hours_2008"]
resultDir_str = "./results/" * OBJ_STR * Dates.format(now(), "_yyyy-mm-dd_HH-MM")
mkdir(resultDir_str)

# read uncertain parameters
uncertain_parameters = CSV.read("./_basis/uncertain_parameters.csv", DataFrame, types=[String, String, String, String, String, Float64, Float64, Float64, Float64, String])

# create or read Latin Hypercube Sample. If you change the number of uncertain parameters you have to create a new LHS, otherwise it is recommended to always use the same lhs which is stored as a csv file.
# lhs = create_lhs(size(uncertain_parameters)[1], 1000, 0) # creates a new Latin Hypercube Sample and saves as csv
lhs = CSV.read("./lhs.csv", DataFrame, header=false) # uses an already created lhc sample

# initialize model (deepcopy of initialized model might lead to a a terminal crash. In this case initialize the model separately for every iteration.)
anyM0 = anyModel(inputMod_arr, resultDir_str, objName = OBJ_STR, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true)

outputsOfInterest = nothing # initialize parameter to store the outputs of interest
objective = DataFrame(iteration = Int[], objectiveValue = Float64[]) # initialize dataframe to store the objective values of the optimization runs


# loop over the iterations to create and optimize the model
for iteration in range(START_ITERATION, END_ITERATION)
    start_time_iteration = @elapsed begin # measure the time for each iteration
        println("[info] iteration $iteration/$END_ITERATION")

        obj_str_model = string(iteration) * "_" * OBJ_STR

        anyM = deepcopy(anyM0) # deepcopy of initialized model might lead to a a terminal crash => initialize the model for every iteration
        # anyM = anyModel(inputMod_arr, resultDir_str, objName = obj_str_model, supTsLvl = 2, repTsLvl = 3, shortExp = 5, emissionLoss = false, holdFixed = true) 

        # modify uncertain parameters
        modify_parameters(anyM, uncertain_parameters, lhs, iteration)

        # set biomass potential to zero if biomass is not allowed
        if !BIOMASS_ALLOWED
            scaleBiomassPotential(anyM, newValue=0, carrier="all", region="all")
        end


        # use biomass only for HVC and oil technologies -> set all other technologies to zero
        if ONLY_BIOMASS_FOR_HVC_AND_OIL
            for row in eachrow(anyM.parts.lim.par[:useUp].data)
                if findTechnology(anyM, row.Te) != "bioConversionOil" && findTechnology(anyM, row.Te) != "bioConversionHvc"
                    row.val = 0
                end
            end
        end

        # use all biomass for district heating technologies (greenWaste, manure and sludge are not used because there is no district heating technology for these carriers) -> all wood is used for woodBoilerDh
        if ALL_BIOMASS_FOR_LTH
            scaleBiomassPotential(anyM, newValue=0, carrier="sludge", region="all")
            scaleBiomassPotential(anyM, newValue=0, carrier="manure", region="all")
            scaleBiomassPotential(anyM, newValue=0, carrier="digestate", region="all")
            scaleBiomassPotential(anyM, newValue=0, carrier="greenWaste", region="all")

            wood_potential = sum(anyM.parts.lim.par[:trdBuyUp].data[anyM.parts.lim.par[:trdBuyUp].data.C .== findCarrier(anyM, "wood"), :].val)

            for row in eachrow(anyM.parts.lim.par[:useLow].data)
                if findTechnology(anyM, row.Te) == "woodBoilerDh"
                    row.val = wood_potential
                end
            end
        end


        # set additional constraint for regret calculations -> Set the minimum use of biomass for the following categories (biomass usage in GWh/a)
        # make sure that the file "df_filtered_with_clusters.csv" is available in EuSysMod. The file is created by the python script in Decide.
        set_constraint = (MIN_BIOMASS_FOR_HVC > 0) || (MIN_BIOMASS_FOR_OIL > 0) || (MAX_BIOMASS_FOR_HVC<9999) || (MAX_BIOMASS_FOR_OIL<9999) # set to true if additional constraints exists

        if set_constraint
            # define the new constraints for the minimum use of biomass for HVC and oil
            new_min_constraint = Dict(
                "bioConversionOil" => MIN_BIOMASS_FOR_OIL*1000, # *1000 to convert TWh to GWh
                "bioConversionHvc" => MIN_BIOMASS_FOR_HVC*1000
            )

            # replace the values in the useLow table with the new values
            for row in eachrow(anyM.parts.lim.par[:useLow].data)
                row.val = get(new_min_constraint, findTechnology(anyM, row.Te), row.val)
            end

            # define the new constraints for the maximum use of biomass for HVC and
            new_max_constraint = Dict(
                "bioConversionOil" => MAX_BIOMASS_FOR_OIL*1000,
                "bioConversionHvc" => MAX_BIOMASS_FOR_HVC*1000,
            )

            # replace the values in the useUp table with the new values
            for row in eachrow(anyM.parts.lim.par[:useUp].data)
                row.val = get(new_max_constraint, findTechnology(anyM, row.Te), row.val)
            end

            #determine if a scenario is already within the respective cluster and should be skipped in the optimization
            cluster_of_scenarios = CSV.read("df_filtered_with_cluster.csv", DataFrame)[!, "cluster_final"] # make sure that the file "df_filtered_with_cluster.csv" is available in EuSysMod. The file is created by the python script in Decide.
        end

        # skip the optimization if the scenario is already within the respective cluster
        if set_constraint && cluster_of_scenarios[iteration] == CLUSTER_INDEX
            continue
        end

        # create optimization model
        createOptModel!(anyM)
        setObjective!(:cost, anyM)

        ### Solve Model ### 
        # The gurobi-parameters are described in the word document "Gurobi_parameters_for_LP" or in the gurobi documentation https://www.gurobi.com/documentation/9.1/refman/parameters.html
        set_optimizer(anyM.optModel, Gurobi.Optimizer)  # select a solver
        set_optimizer_attribute(anyM.optModel, "Method", 2); # Barrier (= interior point) method
        set_optimizer_attribute(anyM.optModel, "Crossover", 0); # no crossover -> only needed to have results that lie on the hull of the feasible region. For our purposes it is not needed.
        set_optimizer_attribute(anyM.optModel, "Threads",t_int); # number of threads 4 or 8 are reasonable values
        set_optimizer_attribute(anyM.optModel, "BarConvTol", BARRIER_CONV_TOL); # barrier convergence tolerance 10-6 needed to have accuracy to calculate the regret
        set_optimizer_attribute(anyM.optModel, "BarHomogeneous", 1) # homogeneous barrier -> faster convergence in my experience
        set_optimizer_attribute(anyM.optModel, "Presolve", 2) # aggressive presolve -> takes time to presolve but faster optimization
        set_optimizer_attribute(anyM.optModel, "BarOrder", 0) # AMD -> Appriximate Minimum Degree ordering -> faster in my experience

        # numeric focus should be low to have fast optimization but if the gets stuck in a suboptimal solution increase the numeric focus and resolve the model
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
        global outputsOfInterest
        global objective
        outputsOfInterest = calculateOutputs(anyM, iteration=iteration, outputsOfInterest=outputsOfInterest, includeWaste=false, resultDir = resultDir_str)
        push!(objective, [iteration, objective_value(anyM.optModel)])
        
        # Save outputsOfInterest and total_cost as a CSV file
        CSV.write(resultDir_str * "/outputsOfInterest.csv", DataFrame(outputsOfInterest))
        typeof(CLUSTER_INDEX) == Int64 ? CSV.write(resultDir_str * "/total_cost_cluster_$CLUSTER_INDEX.csv", DataFrame(objective)) : CSV.write(resultDir_str * "/total_cost.csv", DataFrame(objective))


        ## Possible outputs of the model - used mostly for single model runs and for debugging
        # reportResults(:summary,anyM, addRep = (:capaConvOut,), addObjName = true)
        # reportResults(:exchange,anyM, addObjName = true)
        # reportResults(:cost,anyM, addObjName = true)
        # plotSankeyDiagram(anyM, dropDown = (:timestep,), minVal = 1.) # sankey for the whole europe
        # plotSankeyDiagram(anyM) # sankey with dropdown for the regions and contires
        # plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", dropDown = (:timestep, ), fontSize=16, name="biomass_europe", minVal = 1.) 
        # plotSankeyDiagram(anyM; ymlFilter = "biomass.yml", fontSize=20, name="biomass") 
        # plotSankeyDiagram(anyM; ymlFilter = "all.yml", dropDown = (:timestep, ), fontSize=16, name="all_europe", minVal = 1. ) 
        # plotSankeyDiagram(anyM; ymlFilter = "all.yml", fontSize=20, name="all") 
    end
    println("Iteration $iteration took $start_time_iteration seconds") 
end
