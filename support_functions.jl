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
