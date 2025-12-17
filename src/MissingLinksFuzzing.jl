module MissingLinksFuzzing
import ArchGDAL as AG
import MissingLinks: graph_from_gdal
import StableRNGs: StableRNG
import Distributions: Normal, Uniform, truncated
import LinearAlgebra: norm2
import Graphs: Graph, complete_graph, prim_mst, dijkstra_shortest_paths, ne, rem_edge!, edges, nv, connected_components,
    neighbors
import MetaGraphsNext: label_for
import StatsBase: sample
import Logging: @debug, @warn, NullLogger, with_logger
import DataFrames: DataFrame, metadata!, nrow
import GLMakie
import Makie: Figure, lines!, scatter!, Observable, Figure, Axis, @lift, GridLayout, Button, on, Label,
    Slider, SliderGrid, colsize!, rowsize!, Fixed, Auto, Relative, events, Mouse, mouseposition, Point2f, limits!,
    Textbox, text!, Checkbox
import ArgParse: ArgParseSettings, @add_arg_table!, parse_args
import Random: rand
import MissingLinks: graph_from_gdal, remove_tiny_islands, fill_distance_matrix!, identify_potential_missing_links, deduplicate_links,
    links_to_gis, get_xy, score_links, realize_graph, create_graph_weights
import Tables: eachrow
import LibGEOS
import ProgressMeter: ProgressUnknown, next!
import ProgressMeter

include("graph_construction.jl")
include("ui.jl")
include("visualizer.jl")
include("fuzz.jl")
include("progressheadless.jl")

export FuzzedGraphSettings, build_fuzzed_graph
end