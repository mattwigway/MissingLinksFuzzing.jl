## Graph 8140876093691447795

This graph gives this error:

┌ Error: Error occurred with seed 8140876093691447795
│   exception =
│    AssertionError: (1091.8928786054205, 1047.781996754713) not found in weights
│    Stacktrace:
│      [1] get_weights(rng::StableRNGs.LehmerRNG, fuzzed::MissingLinksFuzzing.FuzzedGraph, G::MetaGraphsNext.MetaGraph{Int64, Graphs.SimpleGraphs.SimpleGraph{Int64}, MissingLinks.VertexID, Tuple{Float64, Float64}, @NamedTuple{length_m::Float64, link_type::Union{Missing, String}, geom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}}, Nothing, MissingLinks.var"#15#16", Float64}, poly::Bool)

which is coming from the assertion that all graph vertices are present in the original weights. The location flagged does not have a vertex in the original dataset, but two lines do cross there. I think what is happening is that these two geometries happen to get unioned by combine_geometries for that part of the test, and when you union two linestrings that cross, GDAL splits them at the point where they cross.