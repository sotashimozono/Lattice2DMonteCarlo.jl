@kwdef struct XYModel <: AbstractModel{Float64}
    J::Float64 = 1.0
end
export XYModel

function local_hamiltonian(
    grids::AbstractVector{Float64},
    lat::Lattice,
    model::XYModel,
    site::Int;
    val::Float64=grids[site],
    kwargs...,
)
    energy = 0.0
    for neighbor in lat.nearest_neighbors[site]
        energy -= model.J * cos(val - grids[neighbor])
    end
    return energy
end

function propose(
    rng::AbstractRNG,
    alg::UniformShift,
    grids::AbstractVector{Float64},
    lat::Lattice,
    model::XYModel,
    site::Int;
    kwargs...,
)
    current = grids[site]
    delta = (rand(rng) - 0.5) * alg.width
    new_val = current + delta
    return (LocalChange(site, new_val, current),)
end
get_binder_coeff(::XYModel) = 2.0

function measure_magnetization(grids, lat, ::XYModel)
    M_x = sum(cos, grids)
    M_y = sum(sin, grids)
    return sqrt(M_x^2 + M_y^2) / length(grids)
end
