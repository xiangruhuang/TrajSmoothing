#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "GD.h" 

using namespace std;

int main(int argc, char** argv){
    // Read Traces
    vector<vector<Point>> traces = read_trace(argv[1]);
    augment_trace(traces, 10);
    //Gradient Descent
    Float lambda = atof(argv[2]);
    cout << "learning rate=" << lambda << endl;
    GDSolver gdsolver = GDSolver(lambda, 100000);
    gdsolver.solve(traces);
}
