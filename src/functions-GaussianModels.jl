using Oscar

##################################
#
#   COLORING
#
##################################

## COLOR DICTIONARIES

function edge_random_coloring(G::Oscar.Graph, n_colors::Int64)
    edgs = collect(Oscar.edges(G))

    # Initatte Dictionary with edge-colors and node-colors
    edges_color = Dict{Oscar.Edge, Int64}(e => 0 for e in edgs)

    for e in edgs
        edges_color[e] = rand(1:n_colors)
    end

    return edges_color
end

function node_random_coloring(G::Oscar.Graph, n_colors::Int64)
    n = n_vertices(G)

    nodes_color = Dict{Int64, Int64}(i => 0 for i in 1:n)

    for i in 1:length(nodes_color)
       nodes_color[i] = rand(1:n_colors)
    end

    return nodes_color
end

function edge_constant_coloring(G::Oscar.Graph)
    edgs = collect(Oscar.edges(G))

    edges_color = Dict{Oscar.Edge, Int64}(e => 0 for e in edgs)

    for e in edgs
        edges_color[e] = 1
    end

    return edges_color
end

function node_constant_coloring(G::Oscar.Graph)
    n = n_vertices(G)
   
    nodes_color = Dict{Int64, Int64}(i => 0 for i in 1:n)

    for i in 1:length(nodes_color)
       nodes_color[i] = 1
    end

    nodes_color
end

function edge_coloring(edges::Vector{Vector{Vector{Int64}}})
    n = n_vertices(G)

    edges_color = Dict{Oscar.Edge, Int64}()
    c = 1
    for color in edges
        for e in color
            edges_color[Edge(e[1],e[2])] = c 
        end
        c = c + 1
    end

    return(edges_color)
end

function edge_coloring(edges::Vector{Vector{Edge}})
    n = n_vertices(G)

    edges_color = Dict{Oscar.Edge, Int64}()
    c = 1
    for color in edges
        for e in color
            edges_color[e] = c 
        end
        c = c + 1
    end

    return(edges_color)
end

function node_coloring(vertices::Vector{Vector{Vector{Int}}})

    vertices_color = Dict{Int, Int}()

    c = 1
    for color in vertices
        for v in color
            vertices_color[v] = c 
        end
        c = c + 1
    end

    return(vertices_color)
end

## UPDATE PARAMETERS

function apply_edge_coloring(L, edges_color)

    colors = sort(unique(collect(values(edges_color))))

    for col in colors
        edges_col = sort([e for (e, c) in edges_color if c == col])
        if length(edges_col) <= 1
            continue
        end

        idx₀ = [Oscar.src(edges_col[1]), Oscar.dst(edges_col[1])]
        for i in 2:length(edges_col)
            idx = [Oscar.src(edges_col[i]), Oscar.dst(edges_col[i])]
            L[idx[1], idx[2]] = L[idx₀[1], idx₀[2]]
        end
    end
    return L
end

function apply_node_coloring(W, nodes_color)

    colors = sort(unique(collect(values(nodes_color))))

    for col in colors
        nodes_col = sort([n for (n, c) in nodes_color if c == col])
        if length(nodes_col) <= 1
            continue
        end

        idx₀ = nodes_col[1]
        for i in 2:length(nodes_col)
            idx = nodes_col[i]
            W[idx, idx] = W[idx₀, idx₀]
        end
    end
    return W
end

function directed_edges_matrix_coloring(M, edges_color)
    L = directed_edges_matrix(M)
    L = apply_edge_coloring(L, edges_color)
    return L
end

function error_covariance_matrix_coloring(M, nodes_color)
    W = error_covariance_matrix(M)
    W = apply_node_coloring(W, nodes_color)
    return W
end


##################################
#
#   PARAMETRIZATION
#
##################################

function parameterization(M::GraphicalModel{Oscar.Graph{Directed}, GaussianRing}, 
    L::MatElem=directed_edges_matrix(M), W::MatElem=error_covariance_matrix(M))

    S = ring(M)
    R = param_ring(M)

    #L = directed_edges_matrix(M)
    Id = identity_matrix(R, n_vertices(M.graph))
    #W = error_covariance_matrix(M)
    Sigma = transpose(inv(Id-L))*W*inv(Id-L)
    
    hom(S.ring, R, reduce(vcat, [[Sigma[i,j] for j in i:n_vertices(M.graph)] for i in 1:n_vertices(M.graph)]))
end

function parameterization(M::GraphicalModel{Oscar.Graph{Directed}, GaussianRing}, 
                          edges_color::Union{<: Dict, Nothing}, nodes_color::Union{<: Dict, Nothing}) 

    if isnothing(edges_color)
        L = directed_edges_matrix(M)
    else
        L = directed_edges_matrix_coloring(M, edges_color)
    end

    if isnothing(edges_color)
        W = error_covariance_matrix(M)
    else
        W = error_covariance_matrix_coloring(M, nodes_color)
    end

    pG = parameterization(M, L, W)

    return pG
end


##################################
#
#   OTHERS
#
##################################


function linear_constraints(I::Ideal)
    gR = grade(base_ring(I))[1]
    mI = minimal_generating_set(ideal(gR, gens(I)))

    linear_pols = findall(x -> x == 1, [total_degree(I[i]) for i in 1:length(mI)]);
    gens(I)[linear_pols]
end

