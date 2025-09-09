using Posets

################################
#
#  Intersect trek-polynomials
#
################################

function intersect_trek_polynomials(f::MPolyRingElem, g::MPolyRingElem)
    if f == g
        return f
    end

    mf = collect(monomials(f))
    mg = collect(monomials(g))

    monom = sort(unique(vcat(mf, mg)))

    cf = [coeff(f, m) for m in monom]
    cg = [coeff(g, m) for m in monom]
    
    c = [minimum([i,j]) for (i,j) in zip(cf, cg)]
        
    sum([c[i]*monom[i] for i in 1:length(c)])
end

function isless(f::QQMPolyRingElem, g::QQMPolyRingElem)
    if f == base_ring(f)(0)
        return true
    end

    mf = collect(monomials(f))
    mg = collect(monomials(g))

    if !all([m in mg for m in mf])
        return false
    end

    cf = [coeff(f, m) for m in mg]
    cg = [coeff(g, m) for m in mg]
    
    all([i <= j for (i,j) in zip(cf, cg)])
end

function new_intersections(p::Oscar.MPolyAnyMap)
    polys = p.img_gens
    resulting_intersections = []

    for i in 1:length(polys)
        for j in i:length(polys)
            q = intersect_trek_polynomials(polys[i], polys[j])
            
            if !(q in polys)
                push!(resulting_intersections, q)
            end
        end
    end
    return unique(resulting_intersections)
end  

all_intersections(p :: Oscar.MPolyAnyMap) = vcat(p.img_gens, new_intersections(p))

function is_pi_graph(M::GraphicalModel{Oscar.Graph{Directed}, GaussianRing}, 
                     edges_color::Union{<: Dict, Nothing}, nodes_color::Union{<: Dict, Nothing})
    p = parameterization(M, edges_color, nodes_color)
    int = new_intersections(p)
    return length(int[int .!= 0]) == 0
end

####################################
#
#  Build poset of trek-polynomials
#
####################################

struct BiDict{K, V}
    fwd::Dict{K, V}
    rev::Dict{V, K}

    function BiDict{K, V}() where {K, V}
        new(Dict{K, V}(), Dict{V, K}())
    end
end

# Add a method to set a key-value pair
function Base.setindex!(b::BiDict, k, v)
    b.fwd[k] = v
    b.rev[v] = k
end

# Add a method to get a value from a key
function Base.getindex(b::BiDict{K,V}, k::K) where {K, V}
    b.fwd[k]
end

function Base.getindex(b::BiDict{K,V}, v::V) where {K, V}
    b.rev[v]
end

function unique_value_keys(p::Oscar.MPolyAnyMap)
    
    seen_values = Set{MPolyRingElem}()
    unique_keys = Vector{MPolyRingElem}()

    keys = gens(p.domain)
    for k in keys
        if p(k) ∉ seen_values
            push!(seen_values, p(k))
            push!(unique_keys, k)
        end
    end

    return unique_keys
end

function poset(M, p::Oscar.MPolyAnyMap)
    s = ring(M).gens

    poset_labels = BiDict{Int, MPolyRingElem}()

    ks = unique_value_keys(p)

    for k in 1:length(ks)
        setindex!(poset_labels, k, ks[k])
    end

    inters = new_intersections(p)
    for k in 1:length(inters)
        setindex!(poset_labels, k+length(ks), inters[k])
    end

    n = length(ks) + length(inters)
    P = Poset(n)

    # Compare sij's vs sij's
    for i in 1:length(ks)
        for j in (i+1):length(ks)
            print(i, " ", j, "\n")
            s1 = ks[i]
            s2 = ks[j]
            if isless(p(s1), p(s2))
                add_relation!(P, poset_labels[s1], poset_labels[s2])
            elseif isless(p(s2), p(s1))
                add_relation!(P, poset_labels[s2], poset_labels[s1])
            end
        end
    end

    # Compare sij's vs intersections
    for i in 1:length(ks)
        for j in 1:length(inters)
            #print(i, " ", j, "\n")
            s1 = ks[i]
            if isless(p(s1), inters[j])
                add_relation!(P, poset_labels[s1], poset_labels[inters[j]])
            elseif isless(inters[j], p(s1))
                add_relation!(P, poset_labels[inters[j]], poset_labels[s1])
            end
        end
    end

    # Compare intersections vs intersections
    for i in 1:length(inters)
        for j in (i+1):length(inters)
            if isless(inters[i], inters[j])
                add_relation!(P, poset_labels[inters[i]], poset_labels[inters[j]])
            elseif isless(inters[j], inters[i])
                add_relation!(P, poset_labels[inters[j]], poset_labels[inters[i]])
            end
        end
    end

    return(P, poset_labels)

end

####################################
#
#  Graph aux functions
#
####################################

function is_cyclic(g::Oscar.Graph{Directed})
    n = Oscar.nv(g)
    if n == 0
    return false
    end

    # 0: unvisited (white set)
    # 1: visiting  (currently in the recursion stack)
    # 2: visited   (finished visiting this node and its descendants)
    visited_states = zeros(Int, n)

    # We need a recursive helper function for the DFS traversal.
    # It's defined inside the main function to capture `g`.
    function dfs_visit(u::Int)
        visited_states[u] = 1 # Mark as "visiting" (gray)

        for v in Oscar.outneighbors(g, u)
            if visited_states[v] == 1
                # Found a back edge to a node in the current recursion stack. This means there is a cycle.
                return true
            end

            if visited_states[v] == 0
                # If the neighbor is unvisited, recursively call DFS on it.
                if dfs_visit(v)
                    # If the recursive call found a cycle, propagate the result up.
                    return true
                end
            end
        end

        # Finished visiting all descendants of u. Mark it as "visited" (black).
        visited_states[u] = 2
        return false # No cycle found from this path
    end

    # Iterate through all vertices to handle disconnected graphs
    for i in 1:n
        if visited_states[i] == 0 # If the vertex is unvisited (white)
            if dfs_visit(i)
                return true # A cycle was found
            end
        end
    end

    return false # No cycles were found in the entire graph
end

function is_cyclic(g::Oscar.Graph{Undirected})
    n = Oscar.nv(g)
    # A graph with 0, 1, or 2 vertices cannot have a cycle.
    if n <= 2
    return false
    end

    visited = falses(n)

    # Recursive helper function for the DFS traversal.
    # `parent` is the vertex from which `u` was discovered.
    function dfs_visit(u::Int, parent::Int)
        visited[u] = true

        # For undirected graphs, use Oscar.neighbors to get all adjacent vertices.
        for v in Oscar.neighbors(g, u)
            if v == parent
                # This is the edge we just came from, so ignore it.
                continue
            end

            if visited[v]
                # If we find a visited neighbor that is not our parent, we've found a back edge.
                # This means a cycle exists.
                return true
            end

            # Recurse for an unvisited neighbor.
            if dfs_visit(v, u)
                # If a cycle was found in the recursive call, propagate the result.
                return true
            end
        end

        return false # No cycle found starting from this vertex `u`.
    end

    # Iterate through all vertices to handle disconnected graphs (forests).
    for i in 1:n
        if !visited[i]
            # Start a new DFS traversal. The parent of the starting node can be an invalid
            # vertex index like 0, since vertex indices are 1-based.
            if dfs_visit(i, 0)
                return true # Cycle found
            end
        end
    end

    return false # No cycles were found in any component of the graph.
end

function vertex_descendants(gr::Oscar.Graph{Directed}, v::Int, desc::Vector{Any} = [])
  outn = Oscar.outneighbors(gr, v)

  d = unique(append!(desc, outn))
  
  if length(outn) > 0
    for i in outn
      d = vertex_descendants(gr, i, d)
    end
    return d 
  end

  return d
end

function sink_descendants(gr::Oscar.Graph{Directed}, v::Int)
  lvs = findall(iszero, Oscar.outdegree(gr))
  v_desc = vertex_descendants(gr, v)
  
  intersect(lvs, v_desc)
end

function descendant_subgraph(G::Oscar.Graph{Directed}, v::Int)
    E = Oscar.edges(G)
    desc = vcat(v, vertex_descendants(G, v))

    E_desc = []
    for e in E
        i =  Oscar.src(e)
        j = Oscar.dst(e)
        if i in desc && j in desc
            push!(E_desc, e)
        end
    end

    return Oscar.graph_from_edges(Directed, Vector{Edge}(E_desc))
end

function source_nodes(G::Oscar.Graph{Directed})
    A = matrix(ZZ, Oscar.adjacency_matrix(G))
    indegrees = sum(A[i,:] for i in 1:size(A)[1])
    outdegree = sum(A[:,i] for i in 1:size(A)[1])

    findall((indegrees .== 0) .& (outdegree .> 0))
end

####################################
#
#  Graph - tableaux properties
#
####################################


function is_directed_tree(G::Oscar.Graph{Directed})
   
    roots = source_nodes(G)
    if length(roots) > 1
        return false
    end

    E = Oscar.edges(G)
    G_und = graph_from_edges(Undirected, E)
    
    if is_cyclic(G_und)
        return false
    end
    
    return true
end

function is_polytree(G::Oscar.Graph)
    roots = source_nodes(G)

    for r in roots
        subG = descendant_subgraph(G, r)
        if !is_directed_tree(subG)
            return false
        end
    end

    return true
end

function is_standard(G::Oscar.Graph)
    A = matrix(ZZ, Oscar.adjacency_matrix(G))
    indegrees = sum(A[i,:] for i in 1:size(A)[1])
    outdegree = sum(A[:,i] for i in 1:size(A)[1])

    v = findall((indegrees .> 1) .& (outdegree .> 0))

    return length(v) == 0
end

#TODO: Finish
function is_semistandard(G::Oscar.Graph)
    error("Not implemented yet, sorry!")
    if is_standard(G)
        return true
    end

end

function is_pi_graph_constant_coloring(G::Oscar.Graph)
    edges_color = edge_constant_coloring(G) 
    nodes_color = node_constant_coloring(G) 
    M = graphical_model(G, gaussian_ring(n_vertices(G)));
    is_pi_graph(M, edges_color, nodes_color)
end

function height(G::Oscar.Graph, v::Int)
    sn = source_nodes(G)

    if v in sn
        return 0
    end

    h = []
    for s in sn
        push!(h, length(shortest_path_dijkstra(G, s, v)) -1)
    end

    h = h[h .> 0]

    return(sort(h))
end

function height(G::Oscar.Graph)
    lvs = findall(iszero, Oscar.outdegree(G))
    h = 0
    for v in lvs
        hv = height(G, v)
        if length(hv) > 0
            h = maximum([h, maximum(hv)])
        end
    end

    return h
end

function width(G::Oscar.Graph, v::Int)
    if !is_polytree(G) 
        error("$G is not a polytree. The method is only implemented for polytrees.")
    end

    # if Oscar.outdegree(G,v) < 2
    #     error("Node $v is not the top of a simple trek")
    # end
    d = Oscar.outdegree(G,v)

    sink_dsc = sink_descendants(G, v)
    if length(sink_dsc) == 0
    return 0
    end

    w = []
    for i in sink_dsc
        l = length(shortest_path_dijkstra(G, v, i))
        push!(w, ifelse(l==0, l, l-1))
    end

    w = sort(w)
    return sum(w[1:d])
end

function vertices_with_height(G::Oscar.Graph, h::Vector{Int})
    V = Oscar.vertices(G)

    v_height = []
    for v in V
        hv = height(G, v)
        if h == hv
            push!(v_height, v)
        end
    end

    return v_height

end

function vertices_with_content(G::Oscar.Graph, c::Vector{Int})
    V = Oscar.vertices(G)

    v_content = []
    for v in V
        cv = content(G, v)
        if all(c .== cv)
            push!(v_content, v)
        end
    end

    return v_content

end

function width_height(G::Oscar.Graph, h::Vector{Int})
    vh = vertices_with_height(G, h)

    w = 0
    for v in vh
        w = max(width(G,v), w)
    end

    return w
end

function width_content(G::Oscar.Graph, c::Vector{Int})
    vh = vertices_with_content(G, c)

    w = 0
    for v in vh
        w = max(width(G,v), w)
    end

    return w
end

function width(G::Oscar.Graph)
    V = Oscar.vertices(G)
    w = 0
    for v in V
        wv = width(G, v)
        w = max(w, wv)
    end
    return w
end

function meet(v1::Vector{T}, v2::Vector{T}) where {T}
    L = max(length(v1), length(v2))
    v1 = ifelse(length(v1) < L, vcat(v1, zeros(T, L-length(v1))), v1)
    v2 = ifelse(length(v2) < L, vcat(v2, zeros(T, L-length(v2))), v2)

    m = zeros(T, L)
    for i in 1:L
        m[i] = min(v1[i], v2[i])
    end
    return m
end

function content(G::Oscar.Graph, v::Int, c::Vector{Any} = [], dist::Int = 1)
    
    in_vert = Oscar.inneighbors(G, v)
    if length(in_vert) == 0
        return c_invariants
    end

    if length(c) >= dist
        c[dist] = c[dist] + length(in_vert)
    else 
        push!(c, length(in_vert))
    end
    #c[dist] = c[dist] + length(in_vert)
    
    dist = dist + 1

    for i in in_vert
        content(G, i, c, dist)
    end
    return c
end

function content(G::Oscar.Graph)

end

function ancestral_polynomial(G::Oscar.Graph, i::Int)
end

