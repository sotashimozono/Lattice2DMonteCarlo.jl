# Order of arguments:
# rng, (dispatch target), grids, lat, model, alg, observers
# others will be passed with kwargs... = kbT, nsteps
# =========================================================
# common interfaces
# =========================================================

function run!(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
    observers::AbstractVector{<:AbstractObserver};
    kbT::Float64=1.0,
    nsteps::Int=1000,
    kwargs...,
) where {T}
    for obs in observers
        observe!(obs, grids, lat, model, 0; kwargs..., kbT=kbT)
    end

    for step in 1:nsteps
        update_step!(rng, grids, lat, model, alg; kwargs..., kbT=kbT)

        # observation
        for obs in observers
            if step % obs.interval == 0
                observe!(obs, grids, lat, model, step; kwargs..., kbT=kbT)
            end
        end
    end
end

function run(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
    observers::AbstractVector{<:AbstractObserver};
    kwargs...,
) where {T}
    grids_copy = copy(grids)
    run!(rng, grids_copy, lat, model, alg, observers; kwargs...)
    return grids_copy
end
export run!, run

function propose(
    rng::AbstractRNG,
    alg::ProposalMethod,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    site::Int,
) where {T}
    return error("propose is not implemented for $(typeof(alg)) and $(typeof(model))")
end
export propose

function local_hamiltonian(
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    site::Int;
    val::T=grids[site],
) where {T}
    return error("local_hamiltonian not implemented for $(typeof(model))")
end
export local_hamiltonian

# total_energy
function total_energy(
    grids::AbstractVector{T}, lat::Lattice, model::AbstractModel{T}
) where {T}
    energy = 0.0
    for site in 1:(lat.N)
        energy += local_hamiltonian(grids, lat, model, site)
    end
    return energy / 2.0
end
export total_energy

function calculate_diff_energy(
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    changes::Tuple{LocalChange{T},LocalChange{T}, Vararg{LocalChange{T}}},
) where {T}
    return error(
        "Generic calculate_diff_energy is not implemented for multi-site updates on $(typeof(model)).",
    )
end
export calculate_diff_energy

function update_step!(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::UpdateAlgorithm;
    kbT::Float64=1.0,
) where {T}
    return error("update_step! not implemented for algorithm $(typeof(alg))")
end
export update_step!

function process_site_selection!(
    rng::AbstractRNG,
    S::SiteSelectionMethod, # dispatch target
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::UpdateAlgorithm;
    kbT::Float64=1.0,
) where {T}
    return error("Unsupported site selection method $(S) for algorithm.")
end
export process_site_selection!

function check_acceptance(rng::AbstractRNG, rule::AcceptanceRule, dE::Float64; kbT::Float64)
    return error("check_acceptance not implemented for $(typeof(rule))")
end
export check_acceptance

function observe!(
    obs::AbstractObserver,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    step::Int;
    kwargs...,
) where {T}
    return error("observe! not implemented for $(typeof(obs))")
end
export observe!