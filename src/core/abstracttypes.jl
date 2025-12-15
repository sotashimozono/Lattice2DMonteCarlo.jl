# =========================================================
# definition of Abstract Types
# =========================================================
abstract type AbstractMonteCarlo end

abstract type AbstractModel{T} <: AbstractMonteCarlo end
# T (Int, Float64, Vector etc.) represents the state type of the model
abstract type AbstractObserver <: AbstractMonteCarlo end
export AbstractObserver

abstract type UpdateAlgorithm <: AbstractMonteCarlo end
abstract type LocalUpdateAlgorithm <: UpdateAlgorithm end
abstract type ClusterUpdateAlgorithm <: UpdateAlgorithm end

abstract type AcceptanceRule <: AbstractMonteCarlo end
abstract type SiteSelectionMethod <: AbstractMonteCarlo end
abstract type ProposalMethod <: AbstractMonteCarlo end
export AcceptanceRule, SiteSelectionMethod, ProposalMethod

export AbstractModel, UpdateAlgorithm, LocalChange
# --- Local Changes ---
# Represents a local change at a single site, used in energy difference calculations
struct LocalChange{T}
    index::Int
    new_val::T
    old_val::T
end
# --- Proposals ---
struct SpinFlip <: ProposalMethod end
export SpinFlip

struct SpinExchange <: ProposalMethod end
export SpinExchange

@kwdef struct UniformShift <: ProposalMethod
    width::Float64 = 0.1
end
export UniformShift

# --- Rules ---
struct Metropolis <: AcceptanceRule end
export Metropolis

struct Glauber <: AcceptanceRule end
export Glauber

# --- Selections ---
struct RandomSiteSelection <: SiteSelectionMethod end
export RandomSiteSelection

struct SequentialSweep <: SiteSelectionMethod end
export SequentialSweep

# --- Algorithms ---
@kwdef struct LocalUpdate{R<:AcceptanceRule,S<:SiteSelectionMethod,P<:ProposalMethod} <:
              LocalUpdateAlgorithm
    rule::R = Metropolis()
    selection::S = RandomSiteSelection()
    proposal::P = SpinFlip()
end
export LocalUpdate