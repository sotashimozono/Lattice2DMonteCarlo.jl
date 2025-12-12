using Lattice2D
using Random, Statistics, Plots

abstract type MonteCarloMethod end
abstract type MetropolisMethod <: MonteCarloMethod end
struct Sweep <: MetropolisMethod end
struct RandomChoice <: MetropolisMethod end
const sweep = Sweep()
const random = RandomChoice()

function local_hamiltonian(grids, lat, site; J=1.0, H=0.0)
    local_h = H * grids[site]
    for neighbor in lat.nearest_neighbors[site]
        local_h -= J * grids[site] * grids[neighbor]
    end
    return local_h
end
function total_hamiltonian(grids, lat; J=1.0, H=0.0)
    total_h = 0.0
    for site in 1:(lat.N)
        total_h += local_hamiltonian(grids, lat, site; J=J, H=H)
    end
    return total_h / 2
end
function metropolis_step!(grids, lat, kbT; J=1.0, H=0.0, method=sweep)
    return metropolis_step!(method, grids, lat, kbT; J=J, H=H)
end
function metropolis_step!(::Sweep, grids, lat, kbT; J=1.0, H=0.0)
    for site in 1:(lat.N)
        dE = -2 * local_hamiltonian(grids, lat, site; J=J, H=H)
        if dE <= 0 || rand() < exp(-dE / kbT)
            grids[site] *= -1
        end
    end
end
function metropolis_step!(::RandomChoice, grids, lat, kbT; J=1.0, H=0.0)
    for _ in 1:(lat.N)
        site = rand(1:(lat.N))
        dE = -2 * local_hamiltonian(grids, lat, site; J=J, H=H)
        if dE <= 0 || rand() < exp(-dE / kbT)
            grids[site] *= -1
        end
    end
end
function Metropolis(
    grids, lat, kbT; J=1.0, H=0.0, method=sweep, n_thermal=1000, n_steps=1000
)
    energies = zeros(Float64, n_steps)
    magnetizations = zeros(Float64, n_steps)
    for step in 1:n_thermal
        metropolis_step!(grids, lat, kbT; J=J, H=H, method=method)
    end
    for step in 1:n_steps
        metropolis_step!(grids, lat, kbT; J=J, H=H, method=method)
        E = total_hamiltonian(grids, lat; J=J, H=H) / lat.N
        M = mean(grids)
        energies[step] = E
        magnetizations[step] = M
    end
    return grids, energies, magnetizations
end
