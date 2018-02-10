#ifndef MOTION_H
#define MOTION_H

#include <fstream>
#include <string>
#include "util.h"
#include <iomanip>
#include <cassert>

class Motion{
    public:
    Motion(string filename){
        this->filename = filename;
        read_bvh(filename); 
    }

    vector<Point> directions(int ord){
        vector<Point> dirs;
        for (vector<Point>::iterator it = rotations.begin() + ord; it < rotations.end() - ord; it++){
            Point dir;
            for (int offset = -ord; offset < ord; offset++){
                Point diff = *(it+offset+1) - *(it+offset);
                if (norm(diff) >= 1e-5){
                    diff = diff / norm(diff);
                }
                for (int d = 0; d < D; d++){
                    dir.push_back(diff[d]);
                }
            }
            dir = dir/ norm(dir);
            for (int d = 0; d < D; d++){
                dir.push_back((*it)[d]);
            }
            assert(dir.size() == D*(2*ord+1));
            dirs.push_back(dir);
        }
        assert(dirs.size() == rotations.size() - 2*ord);
        return dirs;
    }
    
    vector<pair<int, Point>> equal_dist_samples(int ord, Float dist){
        vector<pair<int, Point>> sample_vecs;
        for (vector<Point>::iterator it = rotations.begin(); it < rotations.end(); it++){
            //left
            vector<Point> left_samples;
            vector<Point>::iterator left_it = it;
            Float rest = 0.0;
            //cerr << "left begin" << endl;
            while (left_it != rotations.begin()){
                vector<Point>::iterator next_it = left_it - 1;
                Float d = distance(*next_it, *left_it);
                int i = 0;
                while (i + 1 <= (rest + d) / dist){
                    i++;
                    Float move = i * dist - rest;
                    Float ratio = move / d;
                    Point p = (*left_it) * ratio + (*next_it) * (1.0-ratio) - *(it);
                    p = p / norm(p);
                    left_samples.push_back(p);
                    if (left_samples.size() == ord){
                        break;
                    }
                }
                if (left_samples.size() == ord){
                    break;
                }
                rest = rest + d - dist * i;
                left_it = next_it;
            }
            //cerr << "left done" << endl;
            if (left_samples.size() < ord){
                continue;
            }
            //cerr << "right begin" << endl;
            //right 
            vector<Point> right_samples;
            vector<Point>::iterator right_it = it;
            rest = 0.0;
            while (right_it != rotations.end()-1){
                vector<Point>::iterator next_it = right_it + 1;
                Float d = norm(*next_it - *right_it);
                for (int i = 1; i <= (rest + d) / dist; i++){
                    Float move = i * dist - rest;
                    Float ratio = move / d;
                    Point p = (*right_it) * ratio + (*next_it) * (1.0-ratio) - *(it);
                    p = p / norm(p);
                    right_samples.push_back(p);
                    if (right_samples.size() == ord){
                        break;
                    }
                }
                if (right_samples.size() == ord){
                    break;
                }
                rest = rest + d - dist * ((rest + d) / dist);
                right_it = next_it;
            }
            //cerr << "right done" << endl;
            if (right_samples.size() < ord){
                continue;
            }
            
            Point sample_vec;
            for (int i = 0; i < ord; i++){
                Point& l_i = left_samples[i];
                Point& r_i = right_samples[i];
                assert(l_i.size() == D);
                assert(r_i.size() == D);
                for (int d = 0; d < D; d++){
                    sample_vec.push_back(l_i[d]);
                }
                for (int d = 0; d < D; d++){
                    sample_vec.push_back(r_i[d]);
                }
            }
            sample_vec = sample_vec/ norm(sample_vec);
            assert((*it).size() == D);
            for (int d = 0; d < D; d++){
                sample_vec.push_back((*it)[d]);
            }
            assert(sample_vec.size() == D*(2*ord+1));
            sample_vecs.push_back(make_pair(it-rotations.begin(), sample_vec));
        }
        return sample_vecs;
    }
    
    void read_bvh(string filename){
        ifstream fin(filename, fstream::in);
        char* line = new char[LINE_LEN];
        int used_joints[] = { 0, 2, 3, 4, 7, 8, 9, 12, 13, 15, 16, 18,
            19, 20, 25, 26, 27};
        bool* useful = new bool[31];
        for (int i = 0; i < 31; i++){
            useful[i] = false;
        }
        for (int i = 0; i < 17; i++){
            useful[used_joints[i]] = true;
        }
        while (!fin.eof()){
            fin.getline(line, LINE_LEN);
            string line_str(line);
            if (line_str.size() < 2 && fin.eof()){
                break;
            }
            if (line_str.find("HIERARCHY") == 0){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("ROOT") == 0){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("{") != -1){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("}") != -1){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("OFFSET") != -1){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("CHANNELS") != -1){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("JOINT") != -1){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("End") != -1){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("MOTION") != -1){
                skeleton.push_back(line_str);
                continue;
            }
            if (line_str.find("Frame") != -1){
                continue;
            }
            vector<string> tokens = split(line_str, " ");
            Point position;
            Point rotation;
            Point useless;
            Point all;
            for (int count = 0; count < tokens.size(); count++){
                if (tokens[count].size() < 2){
                    continue;
                }
                Float fval = stod(tokens[count]);
                if (count <= 2){
                    position.push_back(fval);
                } else {
                    assert(count / 3 - 1 < 31);
                    if (useful[count / 3 - 1]){
                        rotation.push_back(fval);
                    } else {
                        useless.push_back(fval);
                    }
                }
            }
            positions.push_back(position);
            rotations.push_back(rotation);
            unused.push_back(useless);
        }
        D = rotations[0].size();
        delete[] line;
        delete[] useful;
        fin.close();
    }

    void write_bvh(string filename){
        int used_joints[] = { 0, 2, 3, 4, 7, 8, 9, 12, 13, 15, 16, 18,
            19, 20, 25, 26, 27};
        bool* useful = new bool[31];
        for (int i = 0; i < 31; i++){
            useful[i] = false;
        }
        for (int i = 0; i < 17; i++){
            useful[used_joints[i]] = true;
        }

        string command = "mkdir -p " + dir(filename);
        system(command.c_str());
        ofstream fout(filename);
        for (vector<string>::iterator it = skeleton.begin(); it != skeleton.end(); it++){
            fout << *it << endl;
        }
        int num_frame = rotations.size();
        fout << "Frames: " << rotations.size() << endl;
        fout << "Frame Time: .0083333" << endl;
        for (int i = 0; i < num_frame; i++){
            Point& p = positions[i];
            fout << fixed << setprecision(6) << p[0] << " " << p[1] << " " << p[2];
            Point r = rotations[i];
            if (centered){
                r = (r * std) + mean_pose;
            }
            Point& u = unused[i];
            int count_r = 0, count_u = 0;
            for (int j = 0; j < 31; j++){
                Float f1, f2, f3;
                if (useful[j]){
                    f1 = r[count_r*3+0];
                    f2 = r[count_r*3+1];
                    f3 = r[count_r*3+2];
                    count_r++;
                } else {
                    f1 = u[count_u*3+0];
                    f2 = u[count_u*3+1];
                    f3 = u[count_u*3+2];
                    count_u++;
                }
                fout << fixed << setprecision(6) << " " << f1 << " " << f2 << " " << f3;
            }
            fout << endl;
        }
        fout.close();
    }

    int size(){
        return rotations.size();
    }

    void center(Point& mean_pose, Point& std){
        this->mean_pose = mean_pose;
        this->std = std;
        for (int i = 0; i < rotations.size(); i++){
            rotations[i] = (rotations[i] - mean_pose) / std;
        }
        centered = true;
    }
    
    int D;
    string filename;
    vector<string> skeleton;
    vector<Point> positions;
    vector<Point> rotations;
    vector<Point> unused;
    bool centered = false;
    Point mean_pose, std;
};

#endif
