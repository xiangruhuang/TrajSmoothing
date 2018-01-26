#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "GD.h" 
#include "Motion.h"

using namespace std;

int main(int argc, char** argv){
    // Read Traces
    int period = 0;
    if (argc < 5){
        cerr << "./motion_smooth <folder> <lambda> <threhold> <output_folder>" << endl;
        exit(1);
    }
    string folder(argv[1]);
    Float lambda = atof(argv[2]);
    Float threshold = atof(argv[3]);
    vector<string> filelist = read_filelist(folder+"/motion_list.txt");
    vector<Motion> motions;
    vector<vector<Point>> traces;
    for (int i = 0; i < filelist.size(); i++){
        cerr << "reading " << folder + "/" + filelist[i];
        Motion motion(folder+"/"+filelist[i]);
        cerr << ", #frames=" << motion.rotations.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
        motions.push_back(motion);
        traces.push_back(motion.rotations);
    }

    //Gradient Descent
    cout << "learning rate=" << lambda << endl;
    cout << "threshold for correspondence=" << threshold << endl;
    string output_folder = argv[4];
    GDSolver gdsolver = GDSolver(lambda, 100000, period, 93, threshold, output_folder);
    gdsolver.solve_motion(traces, motions, filelist);
}
