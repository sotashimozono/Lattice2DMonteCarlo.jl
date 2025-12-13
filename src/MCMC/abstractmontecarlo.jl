
# =========================================================
# 3. 実行エンジン (Execution Engine)
# =========================================================

# --- Top Level: run! ---
function run!(
    rng::AbstractRNG,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::UpdateAlgorithm,
    observers::AbstractVector{<:AbstractObserver}, # 型を緩める
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
            if step % obs.interval == 0 # intervalはObserverが持つ
                observe!(obs, grids, lat, kbT, model, step)
            end
        end
    end
    return grids
end

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

# --- Layer 3: Site Selection Loop ---
# sequencial
function process_site_selection!(rng, ::SequentialSweep, grids, lat, kbT, model, alg)
    for site in 1:(lat.N)
        update_single_site!(rng, site, grids, lat, kbT, model, alg)
    end
end

# random
function process_site_selection!(rng, ::RandomSiteSelection, grids, lat, kbT, model, alg)
    for _ in 1:(lat.N)
        site = rand(rng, 1:(lat.N))
        update_single_site!(rng, site, grids, lat, kbT, model, alg)
    end
end

# --- Layer 4: Single Site Update (Core Logic) ---
function update_single_site!(
    rng::AbstractRNG,
    site::Int,
    grids::AbstractVector{T},
    lat::Lattice,
    kbT::Float64,
    model::AbstractModel{T},
    alg::LocalUpdate,
) where {T}
    current_state = grids[site]
    proposed_state = propose_state(rng, alg.proposal, model, current_state)

    dE = calculate_local_dE(grids, lat, site, proposed_state, current_state, model)

    if check_acceptance(rng, alg.rule, dE, kbT)
        grids[site] = proposed_state
    end
end
# Isingの場合: 単純反転

# Pottsの場合: 別の色を選ぶ

# XYの場合: 回転させる
function calculate_local_dE(grids, lat, site, new_val, old_val, model::IsingModel)
    # Isingの高速化: 周囲との相互作用だけ計算
    # E_new - E_old
    # ここでは概念的に書きますが、実際は local_hamiltonian を使えばOK
    E_old = local_hamiltonian(grids, lat, site, model) # ここは現在値を使う

    # グリッドを一時的に書き換えて計算する実装と、
    # 引数で new_val を受け取れる local_hamiltonian を作る実装がある。
    # 計算コスト的には引数渡し推奨。

    # 簡易実装:
    grids[site] = new_val
    E_new = local_hamiltonian(grids, lat, site, model)
    grids[site] = old_val # 戻す

    return E_new - E_old
end