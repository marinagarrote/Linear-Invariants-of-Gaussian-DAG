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


##########################################
#
#  Intersect vectors of heights or content
#
##########################################

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

function isless(v1::Vector{T}, v2::Vector{T}) where {T}

    L = max(length(v1), length(v2))
    v1 = ifelse(length(v1) < L, vcat(v1, zeros(T, L-length(v1))), v1)
    v2 = ifelse(length(v2) < L, vcat(v2, zeros(T, L-length(v2))), v2)

    
    all(v1 .<= v2)
end

####################################
#
#  Structs Posets
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

struct labeledPoset{T, K, V}
    poset::Poset{T}
    labels::BiDict{K, V}

    function labeledPoset{T, K, V}() where {T, K, V}
        new(Poset{T}, BiDict{K,V})
    end

    function labeledPoset{T, K, V}(p::Poset{T}, l::BiDict{K,V}) where {T, K, V}
        new{T, K, V}(p, l)
    end
end

function labeledPoset(P::Poset{T}, labels::BiDict{K, V}) where {T, K, V}
    # This calls the inner constructor labeledPoset{T, K, V}(P, labels)
    return labeledPoset{T, K, V}(P, labels)
end

poset(P::labeledPoset) = P.poset
p_labels(P::labeledPoset) = P.labels.fwd



####################################
#
#  Trek-polynomial Posets
#
####################################

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

    return labeledPoset(P, poset_labels)
end

function plot_poset(P::labeledPoset{T, K, V}) where {T, K, V<: MPolyRingElem}
    G = convert_to_Oscar(cover_digraph(P))
    lay = hierarchical_layout(G)

    n = Oscar.nv(G)
    G = covert_to_DiGraph(G)
    
    # Node labels
    n_node_colors = 2
    
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
    colors = generate_n_colors(2)

    nodefillc = [nodes_color[i] for i in 1:n] #vcat([findall(x -> x ==i, nodelabel_unique) for i in nodelabel]...)
    nodefillc = colors[nodefillc]

    #edgefillc = vcat([findall(x -> x == i, edgelabel_unique) for i in edge_matrix[:, 3]]...)
    #edgefillc = vcat([findall(x -> x == i, edgelabel_unique) for i in edge_matrix[:, 3]]...)
    edgefillc = edge_matrix[:, 3]
    edgefillc = colors[edgefillc .+ n_node_colors]

    #layout=(args...)->spring_layout(args...; C=20) 
    p = gplot(G, layout=lay, nodefillc=nodefillc, nodelabel=1:n, 
    edgestrokec=edgefillc, EDGELINEWIDTH=0.7,
    title=title, title_color="white", title_size=6)#, 
    #edgelabelc=edgefillc,
    #edgelabel=edge_matrix[:,3],  edgelabeldistx=0.8, edgelabeldisty=0.8)
    return(p)
end



####################################
#
#  Height/Content Posets
#
####################################

# Define an abstract type for all our plot styles
abstract type PlotPosetStyle end

# Define a concrete type for each symbol value
struct HeightLayout <: PlotPosetStyle end
struct ContentLayout <: PlotPosetStyle end


function poset(G::Oscar.Graph, type::Symbol)
    if type == :height
        return _poset(G, HeightLayout())
    elseif type == :content
        return _poset(G, ContentLayout())
    else
        error("Unsupported poset type: ", type)
    end
end

function _poset(G::Oscar.Graph, layout::HeightLayout)
    V = Oscar.vertices(G)

    h = []
    
    for v in V
        push!(h, height(G,v))
    end
    h = unique(h)
    sort!(h, by = v -> (length(v), v))

    poset_labels = BiDict{Int, Vector{Int}}()

    for k in 1:length(h)
        setindex!(poset_labels, k, h[k])
    end

    P = Poset(length(h))

    # Compare sij's vs sij's
    for i in 1:length(h)
        for j in (i+1):length(h)
            h1 = h[i]
            h2 = h[j]
            if isless(h1, h2)
                add_relation!(P, poset_labels[h1], poset_labels[h2])
            elseif isless(h2,h1)
                add_relation!(P, poset_labels[h2], poset_labels[h1])
            end
        end
    end

    return labeledPoset(P, poset_labels)

end

function _poset(G::Oscar.Graph, layout::ContentLayout)
    
   V = Oscar.vertices(G)

    c = []
    
    for v in V
        push!(c, content(G,v))
    end
    c = unique(c)
    sort!(c, by = v -> (length(v), v))

    poset_labels = BiDict{Int, Vector{Int}}()

    for k in 1:length(c)
        setindex!(poset_labels, k, c[k])
    end

    P = Poset(length(c))

    for i in 1:length(c)
        for j in (i+1):length(c)
            c1 = c[i]
            c2 = c[j]
            if isless(c1, c2)
                add_relation!(P, poset_labels[c1], poset_labels[c2])
            elseif isless(c2,c1)
                add_relation!(P, poset_labels[c2], poset_labels[c1])
            end
        end
    end

    return labeledPoset(P, poset_labels)


end

