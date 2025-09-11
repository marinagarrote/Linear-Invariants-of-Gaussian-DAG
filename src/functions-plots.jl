using Graphs, Colors, GraphPlot
using Gadfly

function generate_n_colors(n)
    return [HSV(i * 360.0 / n, 0.8, 0.9) for i in 0:(n - 1)]
end

function convert_to_Oscar(::Type{T}, G::DiGraph{Int64}) where {T <:Union{Directed, Undirected}}
    edges = []
    for e in Graphs.edges(G)
        push!(edges, Oscar.Edge(e.src, e.dst))
    end

    return Oscar.graph_from_edges(T, Vector{Edge}(edges))
end

convert_to_Oscar(G::DiGraph{Int64}) = convert_to_Oscar(Directed, G)

function covert_to_DiGraph(G::Oscar.Graph)
     A =  matrix(ZZ, Oscar.adjacency_matrix(G))
    return Graphs.DiGraph(Matrix{Int}(A))
end

function create_fixed_layout(G::Oscar.Graph)
    graph = covert_to_DiGraph(G)
    positions = spring_layout(graph; C=20)
    fixed_layout = (g) -> positions
    return fixed_layout
end

function plot_colored_graph(G::Oscar.Graph, edges_color, nodes_color, 
                           layout=layout=(args...)->spring_layout(args...; C=20), title="")

    n = Oscar.nv(G)
    G = covert_to_DiGraph(G)
    
    # Node labels
    # nodelabel = collect(values(nodes_color)) #[nodes_color[k] for k in 1:n]
    # nodelabel_unique = unique(nodelabel)
    n_node_colors = maximum(unique(collect(values(nodes_color))))
    
    # Edge labels
    edge_list = collect(Graphs.edges(G))  # Collect into a vector
    edge_matrix = zeros(Int64, length(edge_list), 3)
    for (i, edge) in enumerate(edge_list)
        edge_matrix[i, 1] = edge.src
        edge_matrix[i, 2] = edge.dst
        edge_matrix[i, 3] = edges_color[Oscar.Edge(edge.src, edge.dst)]
    end
    edgelabel_unique = sort(unique(edge_matrix[:, 3]))

    # Colors
    n_colors = n_node_colors + maximum(edgelabel_unique)
    # print(n_colors)
    colors = generate_n_colors(n_colors)

    nodefillc = [nodes_color[i] for i in 1:n] #vcat([findall(x -> x ==i, nodelabel_unique) for i in nodelabel]...)
    nodefillc = colors[nodefillc]

    #edgefillc = vcat([findall(x -> x == i, edgelabel_unique) for i in edge_matrix[:, 3]]...)
    #edgefillc = vcat([findall(x -> x == i, edgelabel_unique) for i in edge_matrix[:, 3]]...)
    edgefillc = edge_matrix[:, 3]
    edgefillc = colors[edgefillc .+ n_node_colors]

    #layout=(args...)->spring_layout(args...; C=20) 
    p = gplot(G, layout=layout, nodefillc=nodefillc, nodelabel=1:n, 
    edgestrokec=edgefillc, EDGELINEWIDTH=0.5,
    title=title, title_color="white", title_size=6,
    pad=1cm,
    plot_size=(30cm, 15cm))
    #, 
    #edgelabelc=edgefillc,
    #edgelabel=edge_matrix[:,3],  edgelabeldistx=0.8, edgelabeldisty=0.8)
    return(p)
end

function hierarchical_layout(G::Oscar.Graph, poset::Bool = false)

    # Ensure DAG
    if is_cyclic(G)
        error("Graph must be a DAG for hierarchical layout.")
    end

    A =  matrix(ZZ, Oscar.adjacency_matrix(G))

    # Find sources (no incoming edges) and sinks (no outgoing edges)
    indegrees = sum(A[i,:] for i in 1:size(A)[1])
    outdegree = sum(A[:,i] for i in 1:size(A)[1])
    sources =  findall((indegrees .== 0) .& (outdegree .> 0))
    sinks   = findall(outdegree .== 0)

    # Compute longest distance from any source (layer assignment)
    layers = Dict(v => 0 for v in Oscar.vertices(G))

   
    for v in Oscar.vertices(G)
        ds = maximum(distance_to_sinks(G, v))
        dh = maximum(height(G,v))
        layers[v] = ifelse(poset, dh, ds)
    end
   
    

    # Group vertices by layer
    layer_groups = Dict{Int, Vector{Int}}()
    for (v, l) in layers
        push!(get!(layer_groups, l, Int[]), v)
    end

    # # Assign coordinates
    # positions = Dict{Int, Tuple{Float64, Float64}}()
    # for (l, group) in sort(collect(layer_groups), by=x->x[1])
    #     n = length(group)
    #     for (i, v) in enumerate(group)
    #         # Spread nodes horizontally, y = -layer for downward layout
    #         positions[v] = (i - (n+1)/2, -l)
    #     end
    # end

    #  # --- Return a callable layout function ---
    # fixed_layout = (g) -> positions
    # return fixed_layout


     # Assign coordinates
    positions = Dict{Int, Tuple{Float64, Float64}}()
    max_width = maximum(length(g) for g in values(layer_groups))

    for (l, group) in sort(collect(layer_groups), by=x->x[1])
        n = length(group)
        for (i, v) in enumerate(sort(group)) # Sort group for consistent ordering
            # Spread nodes horizontally, y = -layer for downward layout
            x_pos = (i - 1) * (max_width / (n > 1 ? n - 1 : 1)) - (max_width / 2)
            if n == 1
                x_pos = 0.0
            end
            positions[v] = ifelse(poset, (x_pos, -l), (x_pos, -l))
             
        end
    end

    # --- This is the corrected part ---
    # Convert the Dict of positions to a Vector ordered by vertex index
    n_v = Oscar.nv(G)
    pos_vector_x = [positions[i][1] for i in 1:n_v]
    pos_vector_y = [positions[i][2] for i in 1:n_v]

    # Return a callable layout function that provides the ordered vector
    pos_vector = (pos_vector_x, pos_vector_y)
    fixed_layout = (g) -> pos_vector
    return fixed_layout
end


#########################
#
#        POSETS
#
#########################


function plot_poset(P::labeledPoset{T, K, V}, G::Oscar.Graph, type::Symbol, title::String) where {T, K, V<: Vector{Int}}
    G_posets = cover_digraph(poset(P))
    lbs = p_labels(P)

    node_lbs = []
    for i in 1:length(Graphs.vertices(G_posets))
        w = width(G, lbs[i], type)
        push!(node_lbs, string(lbs[i])*" w="*string(w))
    end

    lay = hierarchical_layout(convert_to_Oscar(G_posets), true)

    p = gplot(G_posets, layout=lay, 
    nodelabel=node_lbs, nodelabeldist=4, nodelabelsize=5,
    EDGELINEWIDTH=0.3, nodelabelangleoffset = π/4.7,
    nodelabelc = "white", NODESIZE = 0.05,
    title=uppercasefirst(string(type))*" "*title*"\n", title_color="white", title_size=5,
    pad=1.2cm, # Add 10% padding on all sides
    plot_size=(15cm, 20cm)
    )
    return p

end

