#include <iostream>
#include "util.h"
#include <fstream>
#include <cassert>
#include "Motion.h"
#include <omp.h>

using namespace std;



int main(int argc, char** argv){
    if (argc < 4){
        cerr << "./align <folder> <id> <threshold> (nickname)" << endl;
        exit(1);
    }
    string folder(argv[1]);
    vector<string> filelist = read_filelist(folder+"/motion_list.txt");
    int id = atoi(argv[2]);
    Float threshold = atof(argv[3]);
    vector<Motion> motions;
    for (int i = 0; i < filelist.size(); i++){
        cerr << "reading" << folder + "/" + filelist[i];
        Motion motion(folder+"/"+filelist[i]);
        cerr << ", #frames=" << motion.rotations.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
        motions.push_back(motion);
    }
    
    // Align one with others
    string nickname = to_string(id);
    if (argc > 4){
        nickname = argv[4];
    }
    string output_folder = folder+"/align_"+nickname;
    Motion& src = motions[id];
    src.write_bvh(output_folder+"/"+filelist[id]);
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
        tgt.write_bvh(output_folder + "/" + filelist[i]);
        ofstream fout(output_folder + "/" + filelist[i]+".aligned", fstream::out);
        for (int i = 0; i < aligned_ids.size(); i++){
            fout << aligned_ids[i] << endl;
        }
        fout.close();
        cerr << filelist[i] << "... output to " << output_folder +"/" + filelist[i];
        cerr << ", #aligned=" << aligned_ids.size() << endl;
    }

    return 0;
}

