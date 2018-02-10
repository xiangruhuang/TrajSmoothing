#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "Motion.h"
#include "Graph.h"

using namespace std;

class Params{
    public:
        Params(int argc, char** argv){
            parse_cmd_line(argc, argv);
        }
        void exit_with_help(){
            cerr << "./build_graph (options) [filelist path]" << endl;
            cerr << "Options:" << endl;
            cerr << "-t: set threshold" << endl;
            cerr << "-m: set motion number" << endl;
            cerr << "-f: set frame number" << endl;
            cerr << "-k: set k for knn search" << endl;
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
                    case 'm': mnum = atoi(argv[i]);
                              break;
                    case 'f': fnum = atoi(argv[i]);
                              break;
                    case 'k': K = atoi(argv[i]);
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
        }
        void dump(){
            cerr << "threshold=" << threshold << ", filelist_path=" << filelist_path << endl;
        }
        Float threshold = 100.0;
        int mnum = 0, fnum = 0, K;
        string filelist_path;
};

inline vector<Motion> read_motions(vector<string> filelist, Params& params){
    vector<Motion> motions;
    for (int i = 0; i < filelist.size(); i++){
        //cerr << "reading " << filelist[i];
        Motion motion(filelist[i]);
        //cerr << ", #frames=" << motion.rotations.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
        motions.push_back(motion);
    }
    return motions;
}


int main(int argc, char** argv){
    Params params(argc, argv);
    vector<string> filelist = read_filelist(params.filelist_path);
    vector<Motion> motions = read_motions(filelist, params);
    Graph graph(motions);
    
    vector<pair<int, int>> ans;
    int K = params.K;
    int found = 0;
    do{
        found = 0;
        ans = graph.search(K, motions[params.mnum].rotations[params.fnum]);
        for (int i = 0; i < K; i++){
            pair<int, int>& pair = ans[i];
            if (pair.first != params.mnum){
                found++;
            }
        }
        if (found == 0){
            K = K * 2;
        };
    } while (found < params.K);
    
    int count = 0;
    cout << "Found " << found << " (" << params.K << ") " <<  " Out of " << K << endl;
    for (int i = 0; i < K; i++){
        pair<int, int>& pair = ans[i];
        if (pair.first == params.mnum){
            continue;
        }
        count++;
        if (count <= params.K){
            cout << filelist[pair.first] << ": Frame=" << pair.second << ", dist=" << distance(motions[params.mnum].rotations[params.fnum], motions[pair.first].rotations[pair.second]);
            cout << endl;
        }
    }
}
