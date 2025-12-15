measure_energy(grids, lat, model) = total_energy(grids, lat, model)
measure_magnetization(grids, lat, model) = error("measure_magnetization not implemented for $(typeof(model))")

mutable struct FunctionObserver <: AbstractObserver
    name::String
    measure_func::Function
    interval::Int
    steps::Vector{Int}
    history::Vector{Float64}
end
export FunctionObserver, measure_magnetization, measure_energy

function FunctionObserver(name::String, func::Function; interval::Int=100)
    return FunctionObserver(name, func, interval, Int[], Float64[])
end

function observe!(
    obs::FunctionObserver,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    step::Int,
) where {T}
    val = obs.measure_func(grids, lat, model)
    push!(obs.steps, step)
    push!(obs.history, val)
    return nothing
end

function to_dataframe(obs::FunctionObserver)
    return DataFrame(:Step => obs.steps, Symbol(obs.name) => obs.history)
end

mutable struct ThermodynamicObserver <: AbstractObserver
    interval::Int
    n_samples::Int
    sum_E::Float64
    sum_E2::Float64
    sum_M::Float64
    sum_M2::Float64
    sum_M4::Float64
end

function ThermodynamicObserver(; interval::Int=10)
    return ThermodynamicObserver(interval, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
end
function observe!(
    obs::ThermodynamicObserver,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    step::Int,
) where {T}
    E = measure_energy(grids, lat, model)
    M = measure_magnetization(grids, lat, model)

    obs.n_samples += 1
    obs.sum_E += E
    obs.sum_E2 += E^2
    obs.sum_M += M

    M2 = M^2
    obs.sum_M2 += M2
    return obs.sum_M4 += M2^2 # M^4
end
function get_thermodynamics(
    obs::ThermodynamicObserver, kbT::Float64, N::Int, model::AbstractModel
)
    if obs.n_samples == 0
        return Dict{String,Float64}()
    end
    mean_E = obs.sum_E / obs.n_samples
    mean_E2 = obs.sum_E2 / obs.n_samples
    mean_M = obs.sum_M / obs.n_samples
    mean_M2 = obs.sum_M2 / obs.n_samples
    mean_M4 = obs.sum_M4 / obs.n_samples

    var_E = mean_E2 - mean_E^2
    Cv = var_E / (kbT^2) / N

    var_M = mean_M2 - mean_M^2
    Chi = var_M / kbT * N

    bindercoeff = get_binder_coeff(model)
    denominator = bindercoeff * (mean_M2^2)
    U4 = (denominator == 0.0) ? 0.0 : (1.0 - mean_M4 / denominator)

    return Dict(
        "Energy" => mean_E / N,
        "Magnetization" => mean_M,
        "SpecificHeat" => Cv,
        "Susceptibility" => Chi,
        "BinderParam" => U4,
        "Samples" => Float64(obs.n_samples),
    )
end
export ThermodynamicObserver, observe!, get_thermodynamics, to_dataframe
