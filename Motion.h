#ifndef MOTION_H
#define MOTION_H

#include <fstream>
#include <string>
#include "util.h"
#include <iomanip>

class Motion{
    public:
    Motion(string filename){
        this->filename = filename;
        read_bvh(filename);
        
    }
    
    void read_bvh(string filename){
        ifstream fin(filename, fstream::in);
        char* line = new char[LINE_LEN];
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
            for (int count = 0; count < tokens.size(); count++){
                if (tokens[count].size() < 2){
                    continue;
                }
                Float fval = stod(tokens[count]);
                if (count <= 2){
                    position.push_back(fval);
                } else {
                    rotation.push_back(fval);
                }
            }
            positions.push_back(position);
            rotations.push_back(rotation);
        }
        delete line;
        fin.close();
    }

    void write_bvh(string filename){
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
            Point& r = rotations[i];
            for (Point::iterator it = r.begin(); it != r.end(); it++){
                fout << " ";
                fout << fixed << setprecision(6) << *it;
            }
            fout << endl;
        }
        fout.close();
    }

    string filename;
    vector<string> skeleton;
    vector<Point> positions;
    vector<Point> rotations;
};

#endif
