include "globals.mzn";
include "all_different.mzn";
int: num_vertices;
int: num_edges;
int: UB;
int: LB;
set of int: vertices = 1..num_vertices;
array[1..num_edges, 1..2] of int: edges;
array[vertices] of var 1..num_vertices: labels;
var int: b;

solve maximize b;

constraint
  LB <= b;

constraint
  b <= UB;

constraint
  forall(e in 1..num_edges) (
    b <= abs(labels[edges[e, 1]] - labels[edges[e, 2]])
  );

constraint all_different(labels);


output [
  "Objective value (b): ", show(b), "\n",
  "Labels: ", show(labels), "\n"
];