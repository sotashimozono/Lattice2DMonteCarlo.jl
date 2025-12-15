using ColorSchemes
const p_init = plot(;
    aspect_ratio=:equal,
    grid=false,
    axis=false,
    ticks=false,
    legend=:bottomright,
)

function plot_setup(lat::Lattice; title="")
    ms = find_marker_size(lat)
    p = plot(;
        aspect_ratio=:equal, grid=false, axis=false, ticks=false, legend=false, title=title
    )
    visualize_bonds(p, lat)
    return p, ms
end
function find_marker_size(lat::Lattice; ms_scale=80.0)
    min_dist = 0.0
    if isempty(lat.bonds)
        min_dist = 1.0
    else
        min_dist = minimum([
            norm(lat.positions[b.src] - lat.positions[b.dst]) for b in lat.bonds
        ])
    end
    xs = [p[1] for p in lat.positions]
    ys = [p[2] for p in lat.positions]
    area = [(maximum(xs) - minimum(xs)), (maximum(ys) - minimum(ys))]
    scaling = min_dist / norm(area)
    marker_size = ms_scale * scaling
    return marker_size
end
function visualize_bonds(p, lat::Lattice)
    threhold = 1.5 * maximum(norm.(lat.unit_cell.basis))
    seg_x, seg_y = Float64[], Float64[]
    for bond in lat.bonds
        src_pos = lat.positions[bond.src]
        dst_pos = lat.positions[bond.dst]
        if norm(dst_pos - src_pos) < threhold
            push!(seg_x, src_pos[1], dst_pos[1], NaN)
            push!(seg_y, src_pos[2], dst_pos[2], NaN)
        end
    end
    return plot!(p, seg_x, seg_y; color=:black, lw=1.0, label="")
end


function plot_state!(
    p::Plots.Plot,
    lat::Lattice,
    grids::AbstractVector,
    model::AbstractModel;
    marker_size=10.0,
)
    xs = [pos[1] for pos in lat.positions]
    ys = [pos[2] for pos in lat.positions]
    colors = get_state_colors(model, grids)
    return scatter!(p, xs, ys; ms=marker_size, mc=colors, markerstrokewidth=0, label="")
end

function get_state_colors(model::IsingModel, grids::AbstractVector)
    return [s > 0 ? :red : :blue for s in grids]
end
function get_state_colors(model::PottsModel, grids::AbstractVector)
    q = model.q
    palette = cgrad(:tab10, q; categorical=true)

    return [palette[(s - 1) / (q - 1)] for s in grids]
end
# XY model : visualize as arrows
function plot_state!(
    p::Plots.Plot,
    lat::Lattice,
    grids::AbstractVector,
    model::XYModel;
    marker_size=10.0,
)
    xs = [pos[1] for pos in lat.positions]
    ys = [pos[2] for pos in lat.positions]
    u = cos.(grids) .* (marker_size * 0.01)
    v = sin.(grids) .* (marker_size * 0.01)
    return quiver!(p, xs, ys; quiver=(u, v), color=:black)
end

function visualize_snapshot(grids, lat, model)
    p, ms = plot_setup(lat; title="$(typeof(model)) Step")
    plot_state!(p, lat, grids, model; marker_size=ms)
    return p
end

export plot_setup,
    find_marker_size, visualize_bonds, plot_state!, visualize_snapshot, get_state_colors