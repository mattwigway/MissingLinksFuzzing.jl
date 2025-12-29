"Randomly combine some geometries into multigeometries"
function combine_geometries(ps::Vector{AG.IGeometry{AG.wkbLineString}})
    Base.require_one_based_indexing(ps)

    ps = Vector{Union{AG.IGeometry{AG.wkbLineString}, AG.IGeometry{AG.wkbMultiLineString}}}(ps)

    for i in 1:(length(ps) ÷ 10 + 1)
        g1 = pop!(ps)
        g2idx = min(i * 5, length(ps))
        g2 = ps[g2idx]
        if !AG.intersects(g1, g2)
            # if they intersect, ArchGDAL will create a node where they cross,
            # which will affect at the very least weight calculation since there
            # will be nodes that are not expected. So if they cross, just don't
            # merge them
            ps[g2idx] = AG.union(g1, g2)
        else
            # just put things back how they were if they cross
            push!(ps, g1)
        end
    end
    
    return ps
end

function combine_geometries(df::DataFrame)
    newdf = DataFrame(geometry=combine_geometries(df.geometry))
    metadata!(newdf, "geometrycolumns", ("geometry",))
    newdf
end