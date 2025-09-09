include("src/functions-GaussianModels.jl")
include("src/functions-plots.jl")
include("src/functions-posets.jl")


## Define a graph from a vector of edges
g = [[1,3],[2,3],[3,4],[2,5]]
G =  Oscar.graph_from_edges(Directed, g)

## Set edges and nodes colors
edges_color = edge_constant_coloring(G) # Same color in all edges
nodes_color = node_constant_coloring(G) # Same color in all nodes

## Plot the graph
p = plot_colored_graph(G, edges_color, nodes_color, hierarchical_layout(G), "G"); p

######################################
#
#   Properties from the graph
#
######################################

## Properties
is_directed_tree(G)
is_polytree(G)
is_standard(G)
is_pi_graph_constant_coloring(G)

## Statistics of a node in the graph:
v = 4
width(G, v) 
height(G, v)
content(G, v)

## General statistics of the graph
height(G)

## Given a vector representing a height: h = [h₁, …, hₙ] 
h = [1,1]
vertices_with_height(G, h)
width_height(G, h) # wₕ(G)

## Given a vector representing a content: c = [c₁, …, cₙ] 
c = [2]
vertices_with_content(G, c)
width_content(G, c)

## Given two vectors x = [x₁, …, xₙ], y = [y₁, …, yₘ]
## meet(x,y) return the minimum entrywise
meet([1,2,5,4], [2,1]) 


#########################################
#
#   Properties from the Graphical Model
#
#########################################

## Define Gaussian DAG Model
M = graphical_model(G, gaussian_ring(n_vertices(G)));

## Compute the parametrization of sᵢⱼ
p = parameterization(M, edges_color, nodes_color)

## Ideal of the model
I = kernel(p)

## Extract the linear constraints
L = linear_constraints(I)

## Intersections of the polynomials s that are not an s
## If the list is not empy, then (G,c) is not a π-graph
new_intersections(p)

## You can also use the function is_pi_graph given a 
## graphical model and the Dicts of colors
is_pi_graph(M, edges_color, nodes_color)


