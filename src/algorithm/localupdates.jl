function process_site_selection!(
    rng::AbstractRNG,
    ::SequentialSweep, # dispatch target
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::UpdateAlgorithm;
    kwargs...
) where {T}
    for site in 1:(lat.N)
        update_single_site!(rng, site, grids, lat, model, alg; kwargs...)
    end
end

function process_site_selection!(
    rng::AbstractRNG,
    ::RandomSiteSelection,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::UpdateAlgorithm;
    kwargs...
) where {T}
    for _ in 1:(lat.N)
        site = rand(rng, 1:(lat.N))
        update_single_site!(rng, site, grids, lat, model, alg; kwargs...)
    end
end


function calculate_diff_energy(
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    changes::Tuple{LocalChange{T}};
    kwargs...
) where {T}
    c = changes[1]
    E_old = local_hamiltonian(grids, lat, model, c.index; val=c.old_val, kwargs...)
    E_new = local_hamiltonian(grids, lat, model, c.index; val=c.new_val, kwargs...)
    return E_new - E_old
end

function calculate_diff_energy(
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    changes::Tuple{Vararg{LocalChange{T}}};
    kwargs...
) where {T}
    return error(
        "Generic calculate_diff_energy is not implemented for multi-site updates on $(typeof(model)). Please implement a specific method in models/$(typeof(model)).jl.",
    )
end


function check_acceptance(rng::AbstractRNG, ::Metropolis, dE::Float64, kbT::Float64; kwargs...)
    if dE <= 0
        return true
    elseif kbT <= 0
        return false
    else
        return rand(rng) < exp(-dE / kbT)
    end
end

function check_acceptance(rng::AbstractRNG, ::Glauber, dE::Float64, kbT::Float64; kwargs...)
    if kbT <= 0
        return dE < 0
    else
        return rand(rng) < 1.0 / (1.0 + exp(dE / kbT))
    end
end

function update_step!(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::LocalUpdate;
    kwargs...
) where {T}
    return process_site_selection!(rng, alg.selection, grids, lat, model, alg; kwargs...)
end

function update_single_site!(
    rng::AbstractRNG,
    site::Int,
    grids::AbstractVector{T},
    lat::Lattice,
    model::AbstractModel{T},
    alg::LocalUpdate;
    kbT::Float64=1.0,
    kwargs...
) where {T}
    changes = propose(rng, alg.proposal, grids, lat, model, site)

    if isempty(changes)
        return nothing
    end

    dE = calculate_diff_energy(grids, lat, model, changes; kwargs...)

    if check_acceptance(rng, alg.rule, dE, kbT; kwargs...)
        for c in changes
            grids[c.index] = c.new_val
        end
    end
end