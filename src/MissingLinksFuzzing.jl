module MissingLinksFuzzing
import ArchGDAL as AG
import MissingLinks: graph_from_gdal
import StableRNGs: StableRNG
import Distributions: Normal, Uniform, truncated
import LinearAlgebra: norm2
import Graphs: Graph, complete_graph, prim_mst, dijkstra_shortest_paths, ne, rem_edge!, edges, nv, connected_components
import StatsBase: sample
import Logging: @debug
import DataFrames: DataFrame, metadata!
import GLMakie
import Makie: Figure, lines!, scatter!, Observable, Figure, Axis, @lift, GridLayout, Button, on, Label,
    Slider, SliderGrid, colsize!, rowsize!, Fixed, Auto, Relative
import ArgParse: ArgParseSettings, @add_arg_table!, parse_args
import Random: rand
import MissingLinks: graph_from_gdal, remove_tiny_islands, fill_distance_matrix!, identify_potential_missing_links, deduplicate_links, links_to_gis

include("graph_construction.jl")
include("ui.jl")
include("visualizer.jl")

export FuzzedGraphSettings, build_fuzzed_graph
end