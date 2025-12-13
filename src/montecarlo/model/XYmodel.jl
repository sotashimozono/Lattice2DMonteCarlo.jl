@kwdef struct XYModel <: AbstractModel{Float64}
    J::Float64 = 1.0
end
export XYModel
function propose_state(rng, p::UniformShift, ::XYModel, current::Float64)
    # random shift from -width/2 to +width/2
    # we don't need to wrap around [0, 2π) because cos() is periodic
    theta = current + (rand(rng) - 0.5) * p.width
    return theta
end

function local_hamiltonian(
    grids::AbstractVector{Float64},
    lat::Lattice,
    site::Int,
    model::XYModel;
    val::Float64=grids[site],
)
    energy = 0.0
    for neighbor in lat.nearest_neighbors[site]
        # J = -J * cos(θ_i - θ_j)
        energy -= model.J * cos(val - grids[neighbor])
    end

    return energy
end