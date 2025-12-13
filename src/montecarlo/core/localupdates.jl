# TODO: local update
function process_site_selection!(
    rng::AbstractRNG,
    ::SequentialSweep,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::LocalUpdate,
) where {T}
    for site in 1:(lat.N)
        update_single_site!(rng, site, grids, lat, kbT, model, alg)
    end
end
function process_site_selection!(
    rng::AbstractRNG,
    S::SiteSelectionMethod,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
) where {T}
    for _ in 1:(lat.N)
        site = rand(rng, 1:(lat.N))
        update_single_site!(rng, site, grids, lat, kbT, model, alg)
    end
end

# TODO: energy difference calculation
function calculate_diff_energy(
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    changes::Tuple{LocalChange{T}},
) where {T}
    c = changes[1]
    E_old = local_hamiltonian(grids, lat, c.index, model; val=c.old_val)
    E_new = local_hamiltonian(grids, lat, c.index, model; val=c.new_val)
    return E_new - E_old
end

function calculate_diff_energy(
    grids,
    lat,
    model::AbstractModel{T},
    changes::Tuple{Vararg{LocalChange{T}}}, # 要素数2以上
) where {T}
    return error(
        "Generic calculate_diff_energy is not implemented for multi-site updates on $(typeof(model)). Please implement a specific method in models/$(typeof(model)).jl.",
    )
end
# TODO: accept criteria
function check_acceptance(rng::AbstractRNG, ::Metropolis, dE::Float64, kbT::Float64)
    if dE <= 0
        return true
    elseif kbT <= 0
        return false
    else
        return rand(rng) < exp(-dE / kbT)
    end
end
function check_acceptance(rng::AbstractRNG, ::Glauber, dE::Float64, kbT::Float64)
    if kbT <= 0
        return dE < 0
    else
        return rand(rng) < 1.0 / (1.0 + exp(dE / kbT))
    end
end

# TODO: update process
function update_step!(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::LocalUpdate,
) where {T}
    return process_site_selection!(rng, alg.selection, grids, lat, kbT, model, alg)
end
function update_single_site!(
    rng::AbstractRNG,
    site::Int,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::LocalUpdate,
) where {T}
    changes = propose(rng, alg.proposal, model, lat, grids, site)
    if isempty(changes)
        return nothing
    end
    dE = calculate_diff_energy(grids, lat, model, changes)
    if check_acceptance(rng, alg.rule, dE, kbT)
        for c in changes
            grids[c.index] = c.new_val
        end
    end
end
