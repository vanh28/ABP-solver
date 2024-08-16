/*********************************************
 * OPL 22.1.1.0 Model
 * Author: hoang
 * Creation Date: Jul 31, 2024 at 5:17:11 PM
 *********************************************/
// abp_model.mod

using CP;
int num_vertices = ...; // Number of vertices
int num_edges = ...;    // Number of edges
int LB = ...; //Lower-bound
int UB = ...; //Upper-bound
range vertices = 1..num_vertices; // Set of vertices
tuple edge {
    int edge_from;
    int edge_to;
}
{edge} edges = ...; // Set of edges

// Decision variable: labels for each vertex
dvar int labels[vertices] in 1..num_vertices;
dvar int b;

// Degree of each vertex (for calculation only)
int vertex_degree[vertices];
int current_max_degree = -1;
// Calculate degrees and find the vertex with the maximum degree
int max_degree_vertex = -1;
execute {
    // Initialize vertex_degree array
    for (var v in vertices) {
        vertex_degree[v] = 0;
    }

    // Compute degree for each vertex
    for (var e in edges) {
        vertex_degree[e.edge_from] += 1;
        vertex_degree[e.edge_to] += 1;
    }

    // Find the vertex with the maximum degree
    for (var v in vertices) {
        if (vertex_degree[v] > current_max_degree) {
            current_max_degree = vertex_degree[v];
            max_degree_vertex = v;
        }
    }
    max_degree_vertex;
}

// Objective: maximize the minimum distance
maximize b;
// Constraints
subject to {
  	// b in range (LB, UB)
  	LB <= b <= UB;
  	
    // Constraint (C.2): b ≤ abs(li − li')
    forall(e in edges) {
        b <= abs(labels[e.edge_from] - labels[e.edge_to]);
    }
    // Constraint (C.3): allDifferent(labels)
    allDifferent(labels);
    // Symmetry breaking constraint: restrict the label of the vertex with maximum degree
    labels[max_degree_vertex] <= num_vertices div 2;
}

// Solution printing
execute {
    writeln("Objective value (b): ", b);
    writeln("Labels: ", labels);
    writeln("Vertex with maximum degree: ", max_degree_vertex);
}
