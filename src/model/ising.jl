@kwdef struct IsingModel <: AbstractModel{Int}
    J::Float64 = 1.0
    h::Float64 = 0.0
end
export IsingModel

function local_hamiltonian(
    grids::AbstractVector{Int},
    lat::Lattice,
    model::IsingModel,
    site::Int;
    val::Int=grids[site],
    kwargs...,
)
    energy = 0.0
    for neighbor in lat.nearest_neighbors[site]
        energy -= model.J * val * grids[neighbor]
    end
    energy -= model.h * val
    return energy
end
function Lattice2DMonteCarlo.total_energy(
    grids::AbstractVector{Int}, lat::Lattice, model::IsingModel
)
    E_bond_sum = 0.0
    E_field_sum = 0.0

    for site in 1:(lat.N)
        s = grids[site]
        for n in lat.nearest_neighbors[site]
            E_bond_sum -= model.J * s * grids[n]
        end
        E_field_sum -= model.h * s
    end

    return (E_bond_sum / 2.0) + E_field_sum
end
function propose(
    rng::AbstractRNG,
    ::SpinFlip,
    grids::AbstractVector{Int},
    lat::Lattice,
    model::IsingModel,
    site::Int;
    kwargs...,
)
    current = grids[site]
    proposed = -current
    return (LocalChange(site, proposed, current),)
end

function propose(
    rng::AbstractRNG,
    ::SpinExchange,
    grids::AbstractVector{Int},
    lat::Lattice,
    model::IsingModel,
    site1::Int;
    kwargs...,
)
    site2 = rand(rng, lat.nearest_neighbors[site1])

    v1 = grids[site1]
    v2 = grids[site2]

    if v1 == v2
        return ()
    end

    return (LocalChange(site1, v2, v1), LocalChange(site2, v1, v2))
end

function Lattice2DMonteCarlo.calculate_diff_energy(
    grids::AbstractVector{Int},
    lat::Lattice,
    model::IsingModel,
    changes::Tuple{LocalChange{Int},LocalChange{Int}};
    kwargs...,
)
    c1 = changes[1]
    c2 = changes[2]

    is_neighbor = c2.index in lat.nearest_neighbors[c1.index]

    dE1 =
        local_hamiltonian(grids, lat, model, c1.index; val=c1.new_val) -
        local_hamiltonian(grids, lat, model, c1.index; val=c1.old_val)

    dE2 =
        local_hamiltonian(grids, lat, model, c2.index; val=c2.new_val) -
        local_hamiltonian(grids, lat, model, c2.index; val=c2.old_val)
    total_dE = dE1 + dE2
    if is_neighbor
        term1 = -model.J * (c1.new_val - c1.old_val) * c2.old_val
        term2 = -model.J * (c2.new_val - c2.old_val) * c1.old_val
        total_dE -= (term1 + term2)
    end

    return total_dE
end
get_binder_coeff(::IsingModel) = 3.0

measure_magnetization(grids, lat, ::IsingModel) = abs(sum(grids)) / length(grids)