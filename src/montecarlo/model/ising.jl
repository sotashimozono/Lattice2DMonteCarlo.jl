@kwdef struct IsingModel <: AbstractModel{Int}
    J::Float64 = 1.0
    h::Float64 = 0.0
end
export IsingModel

function propose(
    rng::AbstractRNG,
    ::SpinFlip,
    model::IsingModel,
    lat::Lattice,
    grids::AbstractVector{Int},
    site::Int,
)
    current = grids[site]
    proposed = -current
    return (LocalChange(site, proposed, current),)
end

function propose(
    rng::AbstractRNG,
    ::SpinExchange,
    model::IsingModel,
    lat::Lattice,
    grids::AbstractVector{Int},
    site1::Int,
)
    site2 = rand(rng, lat.nearest_neighbors[site1])

    v1 = grids[site1]
    v2 = grids[site2]
    if v1 == v2
        return ()
    end

    return (LocalChange(site1, v2, v1), LocalChange(site2, v1, v2))
end

function local_hamiltonian(
    grids::AbstractVector{Int},
    lat::Lattice,
    site::Int,
    model::IsingModel;
    val::Int=grids[site],
)
    energy = 0.0
    for neighbor in lat.nearest_neighbors[site]
        energy -= model.J * val * grids[neighbor]
    end
    energy -= model.h * val
    return energy
end