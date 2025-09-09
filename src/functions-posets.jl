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



