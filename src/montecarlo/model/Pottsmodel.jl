@kwdef struct PottsModel <: AbstractModel{Int}
    q::Int = 3
    J::Float64 = 1.0
end
export PottsModel

function propose_state(rng, ::SpinFlip, model::PottsModel, current::Int)
    shift = rand(rng, 1:(model.q - 1))
    return mod1(current + shift, model.q)
end

function local_hamiltonian(
    grids::AbstractVector{Int},
    lat::Lattice,
    site::Int,
    model::PottsModel;
    val::Int=grids[site],
)
    energy = 0.0
    for neighbor in lat.nearest_neighbors[site]
        if val == grids[neighbor]
            energy -= model.J
        end
    end
    return energy
end