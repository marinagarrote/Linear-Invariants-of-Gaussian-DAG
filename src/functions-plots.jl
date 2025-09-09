using Graphs, Colors, GraphPlot

function generate_n_colors(n)
    return [HSV(i * 360.0 / n, 0.8, 0.9) for i in 0:(n - 1)]
end

function create_fixed_layout(G)
    A =  matrix(ZZ, Oscar.adjacency_matrix(G))
    graph = Graphs.DiGraph(Matrix{Int}(A))
    positions = spring_layout(graph; C=20)
    fixed_layout = (g) -> positions
    return fixed_layout
end

function plot_colored_grap(G, edges_color, nodes_color, 
                           layout=layout=(args...)->spring_layout(args...; C=20), title="")

    A =  matrix(ZZ, Oscar.adjacency_matrix(G))
    G = Graphs.DiGraph(Matrix{Int}(A))
    n = size(A,1)

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
    edgestrokec=edgefillc, EDGELINEWIDTH=0.7,
    title=title, title_color="white", title_size=6)#, 
    #edgelabelc=edgefillc,
    #edgelabel=edge_matrix[:,3],  edgelabeldistx=0.8, edgelabeldisty=0.8)
    return(p)
end