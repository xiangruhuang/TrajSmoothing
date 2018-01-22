#ifndef UTIL_H
#define UTIL_H

#include<fstream>
#include<vector>
#include<cassert>
#include<iostream>
#include<cmath>

using namespace std;

typedef float Float;
typedef vector<Float> Point;

const int LINE_LEN = 100000000;

inline void print_point(Point const& a, string comment){
    cerr << comment;
    cerr << ":\t";
    for (int i = 0; i < a.size(); i++){
        if (i != 0){
            cerr << " ";
        }
        cerr << a[i];
    }
    cerr << endl;
}

inline Point zero_point(int d){
    Point zero_point;
    for (int i = 0; i < d; i++){
        zero_point.push_back(0.0);
    }
    return zero_point;
}

inline Point operator + (Point const &a, Point const &b){
    vector<Float> c;
    for (int i = 0; i < a.size(); i++){
        c.push_back(a[i] + b[i]);
    }
    return c;
}

inline Point operator + (Point const &a, Float const &c){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] + c);
    }
    return ans;
}

inline Point operator - (Point const &a, Point const &b){
    Point c;
    for (int i = 0; i < a.size(); i++){
        c.push_back(a[i] - b[i]);
    }
    return c;
}

inline Point operator * (Point const &a, Float const &c){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] * c);
    }
    return ans;
}

inline Point operator / (Point const &a, Float const &c){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] / c);
    }
    return ans;
}

inline Point operator / (Point const &a, Point const &b){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] / b[i]);
    }
    return ans;
}

inline Point operator * (Point const &a, Point const &b){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] * b[i]);
    }
    return ans;
}

inline Point sqrt_point(Point const &a){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(sqrt(a[i]));
    }
    return ans;
}

Float normsq(Point const &a){
    Float ans = 0.0;
    for (int i = 0; i < a.size(); i++){
        ans += a[i] * a[i];
    }
    return ans;
}

Float distance(vector<Float> a, vector<Float> b){
    Float dist = 0.0;
    for (int i = 0; i < a.size(); i++){
        Float diff = a[i] - b[i];
        dist += diff * diff;
    }
    return dist;
};

Float distance(pair<Float, Float> a, pair<Float, Float> b){
    Float diff1 = a.first - b.first, diff2 = a.second - b.second;
    Float dist = diff1 * diff1 + diff2 * diff2;
    return dist;
};

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

vector<string> readlines(ifstream& fin){
    vector<string> ans;
    char* line = new char[LINE_LEN];
    while (!fin.eof()){
        fin.getline(line, LINE_LEN);
        string line_str(line);
        ans.push_back(line_str);
    }
    return ans;
}

vector<vector<Point>> read_trace(string filename, int D){
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
        int len = tokens.size() / D;

        assert(tokens.size() % D == 0);
        for (int i = 0; i < len; i++){
            Point new_point;
            for (int d = 0; d < D; d++){
                new_point.push_back(stof(tokens[i*D + d]));
            }
            trace.push_back(new_point);
        }
        traces.push_back(trace);
    }
    delete[] line;

    return traces;
};

void augment_trace(vector<vector<Point>>& traces, int n){
    for (int i = 0; i < traces.size(); i++){
        vector<Point> trace = traces[i];
        vector<Point> new_trace;
        new_trace.push_back(trace[0]);
        for (int j = 1; j < trace.size(); j++){
            for (int t = 0; t < n; t++){
                Float frac = t * 1.0 / n;
                Point middle_t = trace[j-1] * (1.0-frac) + trace[j] * frac;
                new_trace.push_back(middle_t);
            }
            new_trace.push_back(trace[j]);
        }
        traces[i] = new_trace;
    }
}

#endif


