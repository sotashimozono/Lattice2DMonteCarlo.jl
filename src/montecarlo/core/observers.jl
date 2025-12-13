# =========================================================
# Observer Interface
# =========================================================

# =========================================================
# Concrete Observers
# =========================================================
@kwdef mutable struct MagnetizationObserver <: AbstractObserver
    interval::Int = 100
    steps::Vector{Int} = Int[]
    history::Vector{Float64} = Float64[]
end
export MagnetizationObserver

function observe!(
    obs::MagnetizationObserver,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    step::Int,
) where {T}
    M = abs(sum(grids)) / length(grids)
    push!(obs.steps, step)
    return push!(obs.history, M)
end

@kwdef mutable struct EnergyObserver <: AbstractObserver
    interval::Int = 100
    steps::Vector{Int} = Int[]
    history::Vector{Float64} = Float64[]
end
export EnergyObserver

function observe!(
    obs::EnergyObserver,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    step::Int,
) where {T}
    E_total = total_energy(grids, lat, model)
    push!(obs.steps, step)
    return push!(obs.history, E_total / lat.N)
end

@kwdef mutable struct ConfigurationObserver{T} <: AbstractObserver
    interval::Int = 1000
    steps::Vector{Int} = Int[]
    history::Vector{Vector{T}} = Vector{T}[]
end
export ConfigurationObserver
function observe!(
    obs::ConfigurationObserver{T},
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    step::Int,
) where {T}
    snapshot = copy(grids)
    push!(obs.steps, step)
    return push!(obs.history, snapshot)
end