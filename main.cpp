#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "GD.h" 

using namespace std;

int main(int argc, char** argv){
    // Read Traces
    int period = 0;
    if (argc < 6){
        cerr << "./main <input_file> <D> <lambda> <threshold for LCS> <output_dir>" << endl;
        exit(1);
    }
    int D = atoi(argv[2]);
    vector<vector<Point>> traces = read_trace(argv[1], D);
    augment_trace(traces, period);
    //Gradient Descent
    Float lambda = atof(argv[3]);
    Float threshold = atof(argv[4]);
    cout << "learning rate=" << lambda << endl;
    cout << "threshold for LCS=" << threshold << endl;
    GDSolver gdsolver = GDSolver(lambda, 100000, period, D, threshold, string(argv[5]));
    gdsolver.solve(traces);
}
