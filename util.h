#ifndef UTIL_H
#define UTIL_H

#include<fstream>
#include<vector>
#include<cassert>
#include<iostream>
#include<cmath>

using namespace std;

typedef double Float;
typedef vector<Float> Point;
class Matrix{
    public:
    Matrix(int n, int m){
        this->n = n;
        this->m = m;
        entry = new Float[n*m];
    }
    ~Matrix(){
        delete entry;
    }
    int n, m;
    Float* entry;
    
};

const int LINE_LEN = 100000000;

inline double* to_double(Point& point){
    int D = point.size();
    double* ans = new double[D];

    for (int d = 0; d < D; d++){
        ans[d] = point[d];
    }
    return ans;
}

inline void print_point(Point & a, string comment){
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

Float dot(Point& a, Point& b){
    Float ans = 0.0;
    for (int i = 0; i < a.size(); i++){
        ans += a[i] * b[i];
    }
    return ans;
}

inline Point to_vector(int D, double* point){
    Point ans;
    for (int d = 0; d < D; d++){
        ans.push_back(point[d]);
    }
    return ans;
}

inline double diffsq(double* a, double* b, int D){
    Float ans = 0.0;
    for (int i = 0; i < D; i++){
        ans += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return ans;
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

inline Point operator - (Point const &a, double* b){
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

inline Float normsq(Point const &a){
    Float ans = 0.0;
    for (int i = 0; i < a.size(); i++){
        ans += a[i] * a[i];
    }
    return ans;
}

inline Float norm(Point const &a){
    Float ans = 0.0;
    for (int i = 0; i < a.size(); i++){
        ans += a[i] * a[i];
    }
    return sqrt(ans);
}

inline Float distance(vector<Float> a, vector<Float> b){
    Float dist = 0.0;
    for (int i = 0; i < a.size(); i++){
        Float diff = a[i] - b[i];
        dist += diff * diff;
    }
    return dist;
};

inline Float distance(Float* a, Float* b, int D){
    Float dist = 0.0;
    for (int i = 0; i < D; i++){
        Float diff = a[i] - b[i];
        dist += diff * diff;
    }
    return dist;
};

inline Float distance(pair<Float, Float> a, pair<Float, Float> b){
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

inline Point project_onto_segment(Point& p, Point& a, Point& b){
    Point ap = p - a;
    Point ab = b - a;
    Point ab_hat = ab / norm(ab);
    Float dp = dot(ap, ab_hat);
    if (dp <= 0.0){
        return a;
    }
    if (dp >= norm(ab)){
        return b;
    }
    return a + ab_hat * dp;
}

inline Point project(Point& p, vector<Point>& traj){
    Float min_dist = 1e10;
    Point ans = traj[0];
    for (int i = 0; i < traj.size()-1; i++){
        Point proj = project_onto_segment(p, traj[i], traj[i+1]);
        Float d = distance(proj, p);
        if (d < min_dist){
            min_dist = d;
            ans = proj;
        }
    }
    return ans;
}

inline string dir(string fullname){
    bool absolute = (fullname.find("/") == 0);
    vector<string> tokens = split(fullname, "/");
    string ans = "";
    if (absolute){
        ans += "/";
    }
    for (int i = 0; i < tokens.size()-1; i++){
        if (i != 0){
            ans += "/";
        }
        ans += tokens[i];
    }
    return ans;
}

vector<string> read_filelist(string filename){
    ifstream fin(filename, fstream::in);
    vector<string> filelist;
    char* line = new char[LINE_LEN];
    while (!fin.eof()){
        fin.getline(line, LINE_LEN);
        string line_str(line);
        if (line_str.size() < 2 && fin.eof()){
            break;
        }
        filelist.push_back(line_str);
    }
    fin.close();
    delete line;
    return filelist;
}

string replace(string a, string before, string after){
    int m = a.find(before);
    if (m < 0){
        return a;
    }
    string ans = a.replace(m, before.size(), after);
    return ans;
}


#endif


