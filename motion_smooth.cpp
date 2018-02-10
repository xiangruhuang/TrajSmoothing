#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "GD.h" 
#include "Motion.h"
#include "Graph.h"

using namespace std;

class Params{
    public:
        Params(int argc, char** argv){
            parse_cmd_line(argc, argv);
        }

        void exit_with_help(){
            cerr << "./motion_smooth (Options) [filelist path] (output_folder, default to 'motion/smooth')" << endl;
            cerr << "Options:" << endl;
            cerr << "-t: set threshold (200.0)" << endl;
            cerr << "-c: set correspondence method, candidates are {euclidean, direction, equal_dist}" << endl;
            cerr << "-l: set lambda (1e-3)" << endl;
            cerr << "-o: set order (1)" << endl;
            cerr << "-k: set K (10)" << endl;
            exit(1);
        }

        void parse_cmd_line(int argc, char** argv){
            int i;
            for(i=1;i<argc;i++){
                if( argv[i][0] != '-' )
                    break;
                if( ++i >= argc )
                    exit_with_help();

                switch(argv[i-1][1]){
                    case 't': threshold = atof(argv[i]);
                              break;
                    case 'l': lambda = atof(argv[i]);
                              break;
                    case 'o': ord = atoi(argv[i]);
                              break;
                    case 'k': K = atoi(argv[i]);
                              break;
                    case 'c': method = argv[i];
                              break;
                    //case 'u': param->do_subSolve = false; --i;
                    //          break;
                    //case 'o': param->heldout_period = atoi(argv[i]);
                    //          break;
                    default:
                              cerr << "unknown option: -" << argv[i-1][1] << endl;
                              exit(0);
                }
            }

            if(i>=argc)
                exit_with_help();

            filelist_path = argv[i];
            if (i + 1 < argc){
                output_folder = argv[i+1];
            }
        }

        void dump(){
            cerr << "Parameters: {" << endl;
            cerr << "\tthreshold=" << threshold << endl;
            cerr << "\tfilelist_path=" << filelist_path << endl;
            cerr << "\tord=" << ord << endl;
            cerr << "\tlambda=" << lambda << endl;
            cerr << "\tk=" << K << endl;
            cerr << "\toutput_folder=" << output_folder << endl;
            cerr << "\tcorrespondence method=" << method << endl;
            cerr << "}" << endl;
        }

        Float threshold = 100.0;
        int K = 10, ord = 1;
        Float lambda = 1e-3;
        string filelist_path, output_folder;
        string method;
};

int main(int argc, char** argv){
    // Read Traces
    Params params(argc, argv);
    params.dump();
    vector<string> filelist = read_filelist(params.filelist_path);
    vector<Motion> motions;
    for (int i = 0; i < filelist.size(); i++){
        cerr << "reading " << filelist[i];
        Motion motion(filelist[i]);
        cerr << ", #frames=" << motion.rotations.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
        motions.push_back(motion);
    }
    Graph graph(motions, params.method, params.ord);
    //Graph graph(motions, "euclidean");
    graph.compute_matchings(params.K, params.output_folder);

    //Gradient Descent
    cout << "learning rate=" << params.lambda << endl;
    cout << "threshold for correspondence=" << params.threshold << endl;
    string output_folder = "motion/smooth";
    GDSolver gdsolver = GDSolver(params.lambda);
    gdsolver.solve_motion(graph, 100000, output_folder);
}
