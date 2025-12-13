# =========================================================
# common interfaces
# =========================================================
function run!(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
    observers::AbstractVector{<:AbstractObserver},
    nsteps::Int,
) where {T}
    # initial observation
    for obs in observers
        observe!(obs, grids, lat, kbT, model, 0)
    end

    for step in 1:nsteps
        update_step!(rng, grids, lat, kbT, model, alg)

        # observation
        for obs in observers
            if step % obs.interval == 0
                observe!(obs, grids, lat, kbT, model, step)
            end
        end
    end
end
function run(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
    observers::AbstractVector{<:AbstractObserver},
    nsteps::Int,
) where {T}
    grids_copy = copy(grids)
    run!(rng, grids_copy, lat, kbT, model, alg, observers, nsteps)
    return grids_copy
end
export run!, run

function propose(
    rng::AbstractRNG,
    alg::ProposalMethod,
    model::AbstractModel{T},
    lat::Lattice,
    grids::AbstractVector{T},
    site::Int,
) where {T}
    return error("propose is not implemented for $(typeof(alg)) and $(typeof(model))")
end
export propose

function local_hamiltonian(
    grids::AbstractVector{T},
    lat::Lattice,
    site::Int,
    model::AbstractModel{T};
    val::T=grids[site],
) where {T}
    return error("local_hamiltonian not implemented for $(typeof(model))")
end
export local_hamiltonian

function total_energy(
    grids::AbstractVector{T}, lat::Lattice, model::AbstractModel{T}
) where {T}
    energy = 0.0
    for site in 1:(lat.N)
        energy += local_hamiltonian(grids, lat, site, model)
    end
    return energy / 2.0
end
export total_energy

function calculate_diff_energy(
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    changes::Tuple{LocalChange{T},LocalChange{T},Vararg{LocalChange{T}}},
) where {T}
    return error(
        "Generic calculate_diff_energy is not implemented for multi-site updates on $(typeof(model)). Please implement a specific method in models/$(typeof(model)).jl.",
    )
end
export calculate_diff_energy

function update_step!(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
) where {T}
    return error("update_step! not implemented for algorithm $(typeof(alg))")
end
export update_step!

function process_site_selection!(
    rng::AbstractRNG,
    S::SiteSelectionMethod,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
) where {T}
    return error("Unsupported site selection method $(S) for algorithm.")
end
export process_site_selection!

function check_acceptance(rng::AbstractRNG, rule::AcceptanceRule, dE::Float64, kbT::Float64)
    return error("check_acceptance not implemented for $(typeof(rule))")
end
export check_acceptance

function observe!(
    obs::AbstractObserver,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    step::Int,
) where {T}
    return error("observe! not implemented for $(typeof(obs))")
end
export observe!