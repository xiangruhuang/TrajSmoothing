#ifndef UTIL_H
#define UTIL_H

#include<fstream>
#include<vector>
#include<cassert>
#include<iostream>
#include<ANN/ANN.h>

using namespace std;

typedef float Float;
typedef vector<Float> Point;
typedef ANNcoord* ANNpoint;
typedef ANNpoint* ANNpointArray;
typedef ANNdist* ANNdistArray;
typedef ANNidx* ANNidxArray;

const int LINE_LEN = 100000000;

vector<Point> lifting(vector<Point> trace, int n){
    vector<Point> new_trace;
    for (int i = n; i < trace.size() - n; i++){
        Point new_point;
        for (int j = i - n; j <= i + n; j++){
            Point p_j = trace[j];
            for (int k = 0; k < p_j.size(); k++){
                new_point.push_back(p_j[k]);
            }
        }
        assert(new_point.size() == trace[i].size()*(2*n+1));
        new_trace.push_back(new_point);
    }
    return new_trace;
}

vector<string> split(string str, string pattern){
    vector<string> str_split;
    size_t i=0;
    size_t index=0;
    while( index != string::npos ){

        index = str.find(pattern,i);
        str_split.push_back(str.substr(i,index-i));

        i = index+1;
    }

    if( str_split.back()=="" )
        str_split.pop_back();

    return str_split;
};

vector<vector<Point>> read_trace(string filename){
    ifstream fin(filename, fstream::in);
    char* line = new char[LINE_LEN];
    vector<vector<Point>> traces;
    int line_count = 0;
    while (!fin.eof()){
        fin.getline(line, LINE_LEN);
        string line_str(line);
        if (line_str.length() < 2 && fin.eof()){
            break;
        }
        vector<string> tokens = split(line_str, " ");
        vector<Point> trace;
        int len = tokens.size() / 2;

        assert(tokens.size() % 2 == 0);
        for (int i = 0; i < len; i++){
            Float lat_i = stod(tokens[i*2 + 0]);
            Float long_i = stod(tokens[i*2 + 1]);
            Point new_point;
            new_point.push_back(lat_i);
            new_point.push_back(long_i);
            trace.push_back(new_point);
        }
        traces.push_back(trace);
    }
    delete[] line;

    return traces;
};


#endif
