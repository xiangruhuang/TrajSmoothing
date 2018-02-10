#include <iostream>
#include "util.h"
#include <fstream>
#include <cassert>
#include "Motion.h"
#include <omp.h>

using namespace std;

int main(int argc, char** argv){
    if (argc < 4){
        cerr << "./align <filelist path> <id> <threshold> (nickname)" << endl;
        exit(1);
    }
    string filelist_path(argv[1]);
    vector<string> filelist = read_filelist(filelist_path);
    //vector<string> tokens = split(filelist, "/");
    //string folder = tokens[tokens.size()-1];
    string src_file(argv[2]);
    int id = -1;
    Float threshold = atof(argv[3]);
    vector<Motion> motions;
    for (int i = 0; i < filelist.size(); i++){
        cerr << "reading " << filelist[i];
        int pos = filelist[i].find(src_file);
        if (pos >= 0){
            assert(id == -1);
            id = i;
        }
        Motion motion(filelist[i]);
        cerr << ", #frames=" << motion.rotations.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
        motions.push_back(motion);
    }
    

    // Align one with others
    string nickname = to_string(id);
    if (argc > 4){
        nickname = argv[4];
    }
    //string output_folder = folder+"/align_"+nickname;
    Motion& src = motions[id];
    src.write_bvh(src.filename + ".project");
    #pragma omp parallel for
    for (int i = 0; i < filelist.size(); i++){
        if (i == id){
            continue;
        }
        Motion& tgt = motions[i];
        tgt.skeleton = src.skeleton;
        vector<Point> positions;
        vector<Point> rotations;
        vector<int> aligned_ids;
        for (int j = 0; j < src.rotations.size(); j++){
            Point& r_j = src.rotations[j];
            Point proj = project(r_j, tgt.rotations);
            Float dist = distance(proj, r_j);
            if (dist < threshold){
                // aligned
                positions.push_back(src.positions[j]);
                aligned_ids.push_back(j);
                rotations.push_back(proj);
            }
        }
        tgt.positions = positions;
        tgt.rotations = rotations;
        string bvh_out = tgt.filename + ".project";
        tgt.write_bvh(bvh_out);
        ofstream fout(tgt.filename+".aligned", fstream::out);
        for (int i = 0; i < aligned_ids.size(); i++){
            fout << aligned_ids[i] << endl;
        }
        fout.close();
        cerr << filelist[i];
        cerr << ", #aligned=" << aligned_ids.size();
        cerr << ", ( " << i << " / " << filelist.size() << " ) ";
        cerr << endl;
    }

    return 0;
}

