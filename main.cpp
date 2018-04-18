#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "GD.h" 
#include "Params.h"

using namespace std;

int main(int argc, char** argv){
    // Read Traces
    Params params(argc, argv);
    params.dump();

    int D = atoi(argv[2]);
    vector<vector<Point>> traces = read_trace(argv[1], D);
    
    Graph graph(traces, params.method, params.ord);
    
    graph.compute_matchings(params.K, params.output_folder);
    
    //Gradient Descent
    cout << "learning rate=" << params.lambda << endl;
    cout << "threshold for correspondence=" << params.threshold << endl;
    string output_folder(params.output_folder);
    GDSolver gdsolver = GDSolver(params.lambda);
    gdsolver.solve(traces, 100000, output_folder);
}
