# Linear Invariants of Colored Gaussian DAG Models

This repository contains a Julia package for studying **Colored Gaussian graphical models on Directed Acyclic Graphs (DAGs)** using the [Oscar.jl](https://oscar.computeralgebra.de/) ecosystem. It implements the methods and combinatorial objects introduced in the paper:

- 📄 **Linear Invariants of Gaussian DAG Models with Partial Homogeneity Constraints** by *Marina Garrote-López*, *Nataliia Kushnerchuk*, and *Liam Solus*.


### 🔬 Core Features

This package offers a suite of utilities to:

-   **Build Models:** Easily construct colored Gaussian graphical models from any DAG.
-   **Analyze Algebraically:** Compute the vanishing ideal of a model to find its linear invariants.
-   **Calculate Combinatorics:** Determine key combinatorial properties of the graph and its nodes, such as **width**, **height**, and **content**.
-  **Verify Properties:** Check if a graph is a **directed tree**, **polytree**, **standard graph**, or a **$\pi$-graph**.


### 🚀 Upcoming Features

Soon it will include tools for computing and plotting various associated posets:

-   Poset of Heights
-   Poset of Contents
-   Poset of Trek-Polynomials


---

## Installation

### 1. Install Julia

Download and install Julia from the official website:  
👉 [https://julialang.org/downloads/](https://julialang.org/downloads/)

Follow the instructions for your operating system (Windows, macOS, Linux).  
After installation, you can check your version by running:

```bash
julia --version
```

---

### 2. Install Oscar and other packages

Open Julia and enter the package manager (press `]` in the REPL). Then install `Oscar`:

```julia
using Pkg
Pkg.add("Oscar")
```

You should also install `Graphs`, `Colors`, `GraphPlot`, if you plan to visualize graphs:

```julia
Pkg.add("Graphs")
Pkg.add("Colors")
Pkg.add("GraphPlot")
```

To work with _Posets_ you need to install one more package:

```julia
Pkg.add("Posets")
```
---

### 3. Recommended: Set Up VS Code for Julia Development 🧑‍💻

While you can run all the code directly in the Julia terminal, using a code editor like **Visual Studio Code (VS Code)** is highly recommended for a much better development experience. It is not strictly necessary but will make coding much easier.

1.  **Install VS Code:** Download and install it for free from the official website:
    👉 [https://code.visualstudio.com/](https://code.visualstudio.com/)

2.  **Install the Julia Extension:**
    * Open VS Code.
    * Navigate to the **Extensions** view on the side panel (or press `Ctrl+Shift+X`).
    * Search for `Julia` in the marketplace.
    * Click **Install** on the extension provided by **julialang**. 

This setup provides powerful features like syntax highlighting, intelligent code completion, an integrated Julia REPL, plotting capabilities, and debugging tools directly within the editor.

---

## Usage

### Include helper functions

```julia
include("functions-GaussianModels.jl")
include("functions-plots.jl")
include("functions-posets.jl")
```

---

### Example: Define a Graph

```julia
g = [[1,3],[2,3],[3,4],[2,5]]
G = Oscar.graph_from_edges(Directed, g)
```

-----

### Set Edge and Node Colors

You can control how edges and nodes are colored in the graph `G`.

The functions `edge_constant_coloring` and `node_constant_coloring` assign a **single uniform color** to all edges or all nodes:

```julia
edges_color = edge_constant_coloring(G)  # all edges same color
nodes_color = node_constant_coloring(G)  # all nodes same color
```

---

#### Alternative colorings

- **Random coloring:** Given an integer `n`, the functions `edge_random_coloring` and `node_random_coloring` assign up to `n` different colors at random to edges or nodes.

```julia
edges_color = edge_random_coloring(G, 3)   # Up to 3 differnet random colors for edges
nodes_color = node_random_coloring(G, 3)   # Up to 3 differnetrandom colors for nodes
```

- **Grouped coloring:**  The functions `edge_coloring` and `node_coloring` let you manually assign colors by **grouping** elements. Each sub-list represents a group, and all edges or nodes in the same group receive the same color.

```julia
# Grouped edge coloring
edges_color = edge_coloring([[[1,3],[2,3]], [[3,4],[2,5]]])
# Edges (1→3, 2→3) share one color; edges (3→4, 2→5) share another.

# Grouped node coloring
nodes_color = node_coloring([[1, 2, 3], [4, 5]])
# Nodes 1, 2, 3 have one color; nodes 4 and 5 another.
```


---

### Plot the Graph

The function `plot_colored_grap` is used to visualize the graph with custom colors.  
It takes as input:
- the graph `G`,
- a coloring scheme for the edges,
- a coloring scheme for the nodes,
- a layout function (here `create_fixed_layout`), and
- a title for the plot.

```julia
p = plot_colored_grap(G, edges_color, nodes_color, create_fixed_layout(G), "G")
p
```

---

### Graph Properties

The following functions allow you to test whether the graph satisfies certain properties:

```julia
is_directed_tree(G)              # true if G is a directed tree
is_polytree(G)                   # true if G is a polytree (directed acyclic graph with no undirected cycles)
is_standard(G)                   # checks whether G is standard
is_pi_graph_constant_coloring(G) # tests whether G with constant coloring is a π-graph
```
---

### Node and Graph Statistics

These functions compute numerical characteristics of the graph or its nodes.

```julia
v = 4
width(G, v)    # width of the node v in G
height(G, v)   # height of node v in G
content(G, v)  # content vector associated with node v

# Global graph statistic
height(G)      # maximum height of the entire graph
```

You can also specify vectors to query properties:
```julia
h = [1,1]
vertices_with_height(G, h)  # returns vertices with given height vector
width_height(G, h)          # width associated with height vector

c = [2]
vertices_with_content(G, c) # returns vertices with given content vector
width_content(G, c)         # width associated with content vector
```

The **meet operation** combines two vectors entrywise, taking the minimum of each pair of entries:

```julia
meet([1,2,5,4], [2,1])  # returns [1,1,0,0]
```

---

### Gaussian DAG Models

Gaussian Directed Acyclic Graph (DAG) models can be constructed from a graph.  

```julia
M = graphical_model(G, gaussian_ring(n_vertices(G)))
```

Once the model is defined, you can compute and analyze its algebraic structure:

```julia
# Parametrization of the model
p = parameterization(M, edges_color, nodes_color)

# Ideal generated by the relations among parameters
I = kernel(p)

# Extract linear constraints from the ideal
L = linear_constraints(I)

# Check for new intersections in the trek-polynomials (if it's not empty, then indicates non-π-graph)
new_intersections(p)

# Direct check whether the model is a π-graph given the coloring
is_pi_graph(M, edges_color, nodes_color)
```

---


## Project Structure
- `functions-GaussianModels.jl` → Functions for Gaussian graphical models  
- `functions-plots.jl` → Plotting utilities  
- `functions-posets.jl` → Functions for graph and poset properties  