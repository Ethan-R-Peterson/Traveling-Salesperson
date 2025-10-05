# Traveling-Salesperson-MST-Solver
This project implements algorithms for solving graph problems related to the Traveling Salesperson Problem (TSP) and Minimum Spanning Tree (MST). It supports three modes of operation, each modeling different optimization strategies:

* MST – Uses Prim’s algorithm to compute a minimum spanning tree. Includes constraints where USA and Canada regions must connect through border points.

* FASTTSP – Generates an approximate TSP tour using the Furthest Insertion heuristic.

* OPTTSP – Finds the optimal TSP tour using a branch-and-bound search with MST-based lower bounds for pruning.

Features

* Command-line interface with --mode {MST|FASTTSP|OPTTSP}

* Region-aware MST construction (USA, Canada, Border constraints)

* Efficient graph data structures and distance matrix precomputation

* Heuristic (FASTTSP) and exact (OPTTSP) TSP solutions

* Formatted output of total cost and solution paths

Technical Highlights

* Implements Prim’s MST, Furthest Insertion heuristic, and branch-and-bound TSP solver

* Uses STL features (vector, algorithm, numeric, getopt) for performance and clarity

* Employs pruning strategies with MST-based bounds to make OPTTSP feasible on larger inputs

* Outputs reproducible results with fixed precision formatting
