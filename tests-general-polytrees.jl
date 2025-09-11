include("src/functions-GaussianModels.jl")
include("src/functions-graphs.jl")
include("src/functions-posets.jl")
include("src/functions-plots.jl")



### Find randomly polytrees that are pi-graphs and not standard

G_list = []
n_iterations = 20
for n in 10:15
    print(n)
    for i in 1:n_iterations
        # Random polytree with Prim's algorithm
        # G = create_random_polytree(n) # 
        
        # Modified algorithm to get polytrees with longer height (i.e. longer content vectors)
        G = create_random_polytree(n; backbone_fraction=0.6, edge_reversal_fraction=0.3) 


        if is_pi_graph_constant_coloring(G) && !is_standard(G)
            push!(G_list, G)
        end
    end
end

### Plot those polytrees, together with the poset of heights and contents
i = rand(1:length(G_list))
G = G_list[i]

is_pi_graph_constant_coloring(G)
is_standard(G)

p = plot_colored_graph(G, edge_constant_coloring(G), 
                       node_constant_coloring(G), 
                       hierarchical_layout(G), "G$i"); 


Ph = poset(G, :height)
p_height = plot_poset(Ph, G, :height, "G$i"); #p_height

Pc = poset(G, :content)
p_content = plot_poset(Pc, G, :content, "G$i"); #p_content

P = vstack(p, hstack(p_height, p_content))







# g = [[1,3],[2,3],[3,4],[4,5],[4,9],[9,10],[1,6],[6,7],[6,8],[7,11]]

g = [[1,3],[2,3],[3,4],[2,5],[3,6]]
G =  Oscar.graph_from_edges(Directed, g)

p = plot_colored_graph(G, edge_constant_coloring(G), 
                       node_constant_coloring(G), 
                       hierarchical_layout(G), "G"); 

Ph = poset(G, :height)
p_height = plot_poset(Ph, G, :height, "G"); #p_height

Pc = poset(G, :content)
p_content = plot_poset(Pc, G, :content, "G"); #p_content

P = vstack(p, hstack(p_height, p_content))
