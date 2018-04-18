#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "GD.h" 
#include "AM.h"
#include "Trace.h"
#include "Graph.h"
#include "Params.h"
#include<cassert>

using namespace std;

int main(int argc, char** argv){
    // Read Traces
    Params params(argc, argv);
    params.dump();
    vector<string> filelist = read_filelist(params.filelist_path);
    vector<Trace*> traces;
    for (int i = 0; i < filelist.size(); i++){
        cerr << "reading " << filelist[i];
        if (params.D == 2){
            GPSTrace* trace = new GPSTrace(filelist[i]);
            cerr << ", #frames=" << trace->traces.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
            traces.push_back(trace);
        } else {
            HumanMotion* trace = new HumanMotion(filelist[i]);
            cerr << ", #frames=" << trace->traces.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
            traces.push_back(trace);
        }
    }
    cerr << "learning rate=" << params.learning_rate << endl;
    cerr << "threshold for correspondence=" << params.threshold << endl;
    
    Graph graph(traces, &params);
    string output_folder = params.output_folder;

    

    //Gradient Descent
    if (params.solver == "GD"){
        GDSolver gdsolver = GDSolver(params.learning_rate);
        for (int iter = 0; iter < 10; iter++){
            graph.construct();
            graph.compute_matchings(params.K, params.output_folder);
            gdsolver.solve(graph, &params);
        }
    } else {
        AMSolver amsolver = AMSolver(params.learning_rate);
        for (int iter = 0; iter < 1; iter++){
            graph.construct();
            graph.compute_matchings(params.K, params.output_folder);
            amsolver.solve(graph, &params);
        }
        graph.construct();
    }
}
