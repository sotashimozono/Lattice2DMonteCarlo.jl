@kwdef struct PottsModel <: AbstractModel{Int}
    q::Int = 3
    J::Float64 = 1.0
end
export PottsModel

function local_hamiltonian(
    grids::AbstractVector{Int},
    lat::Lattice,
    model::PottsModel,
    site::Int;
    val::Int=grids[site],
    kwargs...,
)
    energy = 0.0
    for neighbor in lat.nearest_neighbors[site]
        if val == grids[neighbor]
            energy -= model.J
        end
    end
    return energy
end

function propose(
    rng::AbstractRNG,
    ::SpinFlip,
    grids::AbstractVector{Int},
    lat::Lattice,
    model::PottsModel,
    site::Int;
    kwargs...,
)
    current = grids[site]
    shift = rand(rng, 1:(model.q - 1))
    new_val = mod1(current + shift, model.q)
    return (LocalChange(site, new_val, current),)
end
get_binder_coeff(::PottsModel) = 3.0

function measure_magnetization(grids, lat, model::PottsModel)
    counts = zeros(Int, model.q)
    for s in grids
        counts[s] += 1
    end
    m_frac = maximum(counts) / length(grids)
    return (model.q * m_frac - 1.0) / (model.q - 1.0)
end