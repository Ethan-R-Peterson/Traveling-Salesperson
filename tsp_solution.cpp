#include <iostream>
#include <getopt.h>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <numeric>

using namespace std;

// ---STRUCTURES---
enum class Region {
    USA,
    CANADA,
    BORDER
};

struct Point {
    int x, y;
    Region region;
    size_t index;
};

struct Node {
    bool inMST = false;
    double distance = numeric_limits<double>::infinity();
    int parent = -1;
};

// ---GLOBALS---
string mode;
vector<Point> points;
vector<vector<double>> optDistanceMatrix;  // Only for OPTTSP
bool hasUSA = false;
bool hasCanada = false;
bool hasBorder = false;

// ---FUNCTION DECLARATIONS---
void printHelp(const char* command);
void getOptions(int argc, char** argv);
void readInput();
double calculateDistanceMST(size_t a, size_t b);
double calculateDistanceTSP(size_t a, size_t b);
void precomputeOptDistanceMatrix();
void primMST();
vector<size_t> furthestInsertionTSP();
void fastTSP_output(vector<size_t> const & tour);

// ---OPTIONS IMPLEMENTATION---
void printHelp(const char* command) {
    cout << "Usage: " << command << " ./donut {OPTION}";
    cout << "Options:\n";
    cout << "  -h, --help           Show this help message\n";
    cout << "  -m, --mode           Precedes required argument of {MST|FASTTSP|OPTTSP}\n";
}

void getOptions(int argc, char** argv) {
    struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"mode", required_argument, nullptr, 'm'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "hm:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                printHelp(argv[0]);
                exit(0);
            case 'm':
                mode = optarg;
                if (mode != "MST" && mode != "FASTTSP" && mode != "OPTTSP") {
                    cerr << "Error: invalid mode " << mode << "\n";
                    exit(1);
                }
                break;
            default:
                cerr << "Invalid option. Use --help to see usage.\n";
                exit(1);
        }
    }

    if (mode.empty()) {
        cerr << "Error: no mode specified\n";
        exit(1);
    }
}

void readInput() {
    size_t n;
    cin >> n;
    points.resize(n);

    for (size_t i = 0; i < n; ++i) {
        cin >> points[i].x >> points[i].y;
        points[i].index = i;
        
        // Determine region
        if (points[i].x > 0 && points[i].y > 0) {
            points[i].region = Region::CANADA;
            hasCanada = true;
        } else if ((points[i].x == 0 && points[i].y == 0) || 
                  (points[i].x > 0 && points[i].y == 0) || 
                  (points[i].x == 0 && points[i].y > 0)) {
            points[i].region = Region::BORDER;
            hasBorder = true;
        } else {
            points[i].region = Region::USA;
            hasUSA = true;
        }
    }
}

double calculateDistanceMST(size_t a, size_t b) {
    // If both points are in different countries and neither is on the border
    if (points[a].region != points[b].region && 
        points[a].region != Region::BORDER && 
        points[b].region != Region::BORDER) {
        return numeric_limits<double>::infinity();
    }
    
    double dx = static_cast<double>(points[a].x - points[b].x);
    double dy = static_cast<double>(points[a].y - points[b].y);
    return sqrt(dx*dx + dy*dy);
}

double calculateDistanceTSP(size_t a, size_t b) { 
    double dx = static_cast<double>(points[a].x - points[b].x);
    double dy = static_cast<double>(points[a].y - points[b].y);
    return sqrt(dx*dx + dy*dy);
}

void precomputeOptDistanceMatrix() {
    size_t n = points.size();
    optDistanceMatrix.resize(n, vector<double>(n));
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double dx = static_cast<double>(points[i].x - points[j].x);
            double dy = static_cast<double>(points[i].y - points[j].y);
            optDistanceMatrix[i][j] = sqrt(dx*dx + dy*dy);
        }
    }
}

void primMST() {
    // Check if MST is possible
    if (hasUSA && hasCanada && !hasBorder) {
        cerr << "Cannot construct MST" << "\n";
        exit(1);
    }

    size_t num_vertices = points.size();
    vector<Node> nodes(num_vertices);
    vector<pair<size_t, size_t>> mstEdges;
    double totalWeight = 0.0;

    // Start with the first node (index 0)
    nodes[0].distance = 0.0;

    for (size_t count = 0; count < num_vertices; ++count) {
        // Find vertex with minimum distance not yet in MST
        size_t pni = 0;
        double min_dist = numeric_limits<double>::infinity();
        
        for (size_t i = 0; i < num_vertices; ++i) {
            if (!nodes[i].inMST && nodes[i].distance < min_dist) {
                min_dist = nodes[i].distance;
                pni = i;
            }
        }

        // If no valid vertex found, graph is disconnected
        if (!(min_dist < numeric_limits<double>::infinity())) {
            cerr << "Cannot construct MST" << "\n";
            exit(1);
        }

        // Add to MST
        nodes[pni].inMST = true;
        if (nodes[pni].parent != -1) {
            // Store edge with smaller index first
            if (nodes[pni].parent < static_cast<int>(pni)) {
                mstEdges.emplace_back(nodes[pni].parent, pni);
            } else {
                mstEdges.emplace_back(pni, nodes[pni].parent);
            }
            totalWeight += min_dist;
        }

        // Update adjacent vertices
        for (size_t v = 0; v < num_vertices; ++v) {
            if (!nodes[v].inMST) {
                double dist = calculateDistanceMST(pni, v);
                if (dist < nodes[v].distance) {
                    nodes[v].distance = dist;
                    nodes[v].parent = static_cast<int>(pni);
                }
            }
        }
    }

    // Output results
    cout << fixed << setprecision(2) << totalWeight << "\n";
    for (const auto& edge : mstEdges) {
        cout << edge.first << " " << edge.second << "\n";
    }
}

vector<size_t> furthestInsertionTSP() {
    size_t n = points.size();
    vector<size_t> tour;
    tour.reserve(n);
    
    // Start with just the first node
    tour.push_back(0);
    
    // Track which nodes are in the tour
    vector<bool> in_tour(n, false);
    in_tour[0] = true;
    
    // For nodes not in tour, track their closest distance to tour and nearest node
    vector<double> min_dist_to_tour(n);
    vector<size_t> nearest_in_tour(n, 0);
    
    // Initialize distances for the starting node
    for (size_t i = 1; i < n; ++i) {
        min_dist_to_tour[i] = calculateDistanceTSP(0, i);
    }
    
    while (tour.size() < n) {
        // Find the node not in tour with maximum min_dist_to_tour
        size_t furthest_node = 0;
        double max_dist = -1.0;
        
        for (size_t i = 0; i < n; ++i) {
            if (!in_tour[i] && min_dist_to_tour[i] > max_dist) {
                max_dist = min_dist_to_tour[i];
                furthest_node = i;
            }
        }
        
        // Find the best position to insert the furthest node
        double min_increase = numeric_limits<double>::infinity();
        size_t best_pos = 0;
        
        for (size_t i = 0; i < tour.size(); ++i) {
            size_t j = (i + 1) % tour.size();
            double original_dist = calculateDistanceTSP(tour[i], tour[j]);
            double new_dist1 = calculateDistanceTSP(tour[i], furthest_node);
            double new_dist2 = calculateDistanceTSP(furthest_node, tour[j]);
            double increase = new_dist1 + new_dist2 - original_dist;
            
            if (increase < min_increase) {
                min_increase = increase;
                best_pos = i + 1;
            }
        }
        
        // Insert the node at the best position 
        tour.insert(tour.begin() + static_cast<ptrdiff_t>(best_pos), furthest_node);
        in_tour[furthest_node] = true;
        
        // Update min_dist_to_tour for remaining nodes
        for (size_t i = 0; i < n; ++i) {
            if (!in_tour[i]) {
                double new_dist = calculateDistanceTSP(furthest_node, i);
                if (new_dist < min_dist_to_tour[i]) {
                    min_dist_to_tour[i] = new_dist;
                    nearest_in_tour[i] = furthest_node;
                }
            }
        }
    }
    return tour;
}

void fastTSP_output(vector<size_t> const & tour) {
    // Calculate total tour length
    double total_length = 0.0;
    for (size_t i = 0; i < tour.size(); ++i) {
        size_t j = (i + 1) % tour.size();
        total_length += calculateDistanceTSP(tour[i], tour[j]);
    }
    
    // Output results
    cout << fixed << setprecision(2) << total_length << "\n";
    for (size_t node : tour) {
        cout << node << " ";
    }
    cout << "\n";
}

class OPTTSP {
private:
    vector<Point>& points;
    vector<vector<double>>& distMatrix;
    vector<size_t> bestPath;
    double bestCost;
    vector<size_t> currentPath;
    double currentCost;
    
public:
    OPTTSP(vector<Point>& pts, vector<vector<double>>& matrix) 
        : points(pts), distMatrix(matrix), bestCost(numeric_limits<double>::max()) {
        // Initialize with FASTTSP solution
        bestPath = furthestInsertionTSP();
        bestCost = calculateTourCost(bestPath);
    }
    
    double calculateTourCost(const vector<size_t>& path) {
        double cost = 0.0;
        for (size_t i = 0; i < path.size(); ++i) {
            size_t j = (i + 1) % path.size();
            cost += distMatrix[path[i]][path[j]];
        }
        return cost;
    }
    
    double calculateMSTCost(const vector<size_t>& vertices) {
        if (vertices.size() <= 1) return 0.0;
        
        vector<bool> inMST(vertices.size(), false);
        vector<double> minDist(vertices.size(), numeric_limits<double>::max());
        minDist[0] = 0.0;
        double totalCost = 0.0;
        
        for (size_t count = 0; count < vertices.size(); ++count) {
            size_t u = vertices.size();
            double minVal = numeric_limits<double>::max();
            
            for (size_t i = 0; i < vertices.size(); ++i) {
                if (!inMST[i] && minDist[i] < minVal) {
                    minVal = minDist[i];
                    u = i;
                }
            }
            
            if (u == vertices.size()) break;
            inMST[u] = true;
            totalCost += minVal;
            
            for (size_t v = 0; v < vertices.size(); ++v) {
                if (!inMST[v]) {
                    double dist = distMatrix[vertices[u]][vertices[v]];
                    if (dist < minDist[v]) {
                        minDist[v] = dist;
                    }
                }
            }
        }
        
        return totalCost;
    }
    
    bool promising(size_t permLength) {
        if (permLength == currentPath.size()) return true;
            // Early exit if we've already found a good solution
        if (currentCost >= bestCost) return false;
        if (currentPath.size() - permLength <= 3) return true;
        
        size_t numUnvisited = currentPath.size() - permLength;
        
        // For small numbers of unvisited nodes, just continue
        if (numUnvisited <= 5) return true;
        
        // Get unvisited nodes
        vector<size_t> unvisited;
        unvisited.reserve(currentPath.size() - permLength);
        for (size_t i = permLength; i < currentPath.size(); ++i) {
            unvisited.push_back(currentPath[i]);
        }
        
        // Calculate MST of unvisited nodes
        double mstCost = calculateMSTCost(unvisited);
        
        // Calculate connections from first and last in partial path to unvisited nodes mst
        double minArm1 = numeric_limits<double>::max();
        double minArm2 = numeric_limits<double>::max();
        
        for (size_t node : unvisited) {
            double dist1 = distMatrix[currentPath[0]][node];
            if (dist1 < minArm1) minArm1 = dist1;
            
            double dist2 = distMatrix[currentPath[permLength-1]][node];
            if (dist2 < minArm2) minArm2 = dist2;
        }
        
        double totalEstimate = currentCost + minArm1 + minArm2 + mstCost;
        return totalEstimate < bestCost;
    }
    
    void genPerms(size_t permLength) {
        if (currentCost >= bestCost) return;
        if (permLength == currentPath.size()) {
            double finalCost = currentCost + distMatrix[currentPath.back()][currentPath[0]];
            if (finalCost < bestCost) {
                bestCost = finalCost;
                bestPath = currentPath;
            }
            return;
        }
        
        if (!promising(permLength)) {
            return;
        }
        
        for (size_t i = permLength; i < currentPath.size(); ++i) {
            swap(currentPath[permLength], currentPath[i]);
            
            // Update current cost
            if (permLength > 0) {
                currentCost += distMatrix[currentPath[permLength-1]][currentPath[permLength]];
            }
            
            genPerms(permLength + 1);
            
            // Revert current cost
            if (permLength > 0) {
                currentCost -= distMatrix[currentPath[permLength-1]][currentPath[permLength]];
            }

            swap(currentPath[permLength], currentPath[i]);
        }
    }
    
    void solve() {
        currentPath.resize(points.size());
        iota(currentPath.begin(), currentPath.end(), 0);
        currentCost = 0.0;
        
        // Start with node 0 fixed
        swap(currentPath[0], currentPath[0]);
        genPerms(1);
        
        // Output results
        cout << fixed << setprecision(2) << bestCost << "\n";
        for (size_t node : bestPath) {
            cout << node << " ";
        }
        cout << "\n";
    }
};

int main(int argc, char** argv) {
    std::ios_base::sync_with_stdio(false);
    getOptions(argc, argv);
    readInput();
    
    if (mode == "MST") {
        primMST();
    }
    else if (mode == "FASTTSP") {
        vector<size_t> fast_tour = furthestInsertionTSP();
        fastTSP_output(fast_tour);
    }
    else if (mode == "OPTTSP") {
        precomputeOptDistanceMatrix();
        OPTTSP solver(points, optDistanceMatrix);
        solver.solve();
    }
    
    return 0;
}
