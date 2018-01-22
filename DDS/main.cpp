#include "util.h"
#include <cassert>
#include <ANN/ANN.h>

using namespace std;

const bool verbose = false;

void print_traces(vector<vector<Point>>& traces){
    for (int i = 0; i < traces.size(); i++){
        for (int j = 0; j < traces[i].size(); j++){
            cout << "(";
            for (int t = 0; t < traces[i][j].size(); t++){
                cout << traces[i][j][t] << ",";
            }
            cout << ") ";
        }
        cout << endl;
    }
    
}

vector<vector<Point>> smoothing(vector<vector<Point>>& traces, int n, int K){
    vector<vector<Point>> lifted_traces;
    int num_point = 0;
    int d = 2;
    for (int i = 0; i < traces.size(); i++){
        lifted_traces.push_back(lifting(traces[i], n));
        num_point += lifted_traces[i].size();
    }
    if (verbose){
        cerr << "lifting done" << endl;
    }
    //print_traces(lifted_traces);
    
    double** points = new double*[num_point];
    int D = d*(2*n+1);
    int count = 0;
    for (int i = 0; i < lifted_traces.size(); i++){
        vector<Point>& ti = lifted_traces[i];
        for (int j = 0; j < ti.size(); j++){
            points[count] = new double[D];
            for (int t = 0; t < D; t++){
                points[count][t] = ti[j][t];
            }
            count++;
        }
    }

    if (verbose){
        cerr << "creating knn tree...";
    }
    ANNkd_tree knn(points, num_point, D);
    if (verbose){
        cerr << "done" << endl;
    }
    int* nn_idx = new int[K];
    double* dists = new double[K];
    vector<vector<Point>> smoothed_traces;
    
    if (verbose){
        cerr << "computing smooth trace...";
    }
    double eps = 1e-6;
    count = 0;
    for (int i = 0; i < lifted_traces.size(); i++){
        vector<Point> smoothed_trace;
        for (int j = 0; j < lifted_traces[i].size(); j++){
            knn.annkSearch(points[count], K, nn_idx, dists, eps);
            Point smoothed_point;
            for (int t = 0; t < D; t++){
                smoothed_point.push_back(0.0);
            }
            for (int k = 0; k < K; k++){
                double* kth = points[nn_idx[k]];
                for (int t = 0; t < D; t++){
                    smoothed_point[t] += kth[t];
                }
            }
            for (int t = 0; t < D; t++){
                smoothed_point[t] /= K;
            }
            count++;
            smoothed_trace.push_back(smoothed_point);
        }
        smoothed_traces.push_back(smoothed_trace);
    }
    if (verbose){
        cout << "projecting back...";
    }
    
    //project back
    vector<vector<Point>> projected_traces;
    for (vector<vector<Point>>::iterator it = smoothed_traces.begin(); it != smoothed_traces.end(); it++){
        vector<Point> projected_trace;
        int org_size = it->size() + n*2;
        int* down = new int[it->size() + 2*n];
        for (int i = 0; i < it->size() + 2*n; i++){
            Point p_i;
            down[i] = 0;
            for (int j = 0; j < d; j++){
                p_i.push_back(0.0);
            }
            projected_trace.push_back(p_i);
            assert(p_i.size() == 2);
        }

        for (int i = 0; i < it->size(); i++){
            int p = 0;
            for (int j = 0; j < 2*n+1; j++){
                int org_idx = i+j;
                down[org_idx]++;
                for (int t = 0; t < d; t++){
                    projected_trace[org_idx][t] += it->at(i)[p];
                    p++;
                }
            }
        }
        for (int i = 0; i < it->size() + 2*n; i++){
            for (int t = 0; t < d; t++){
                projected_trace[i][t] /= down[i];
            }
        }
        delete down; 
        projected_traces.push_back(projected_trace);
    }
    if (verbose){
        cerr << "done" << endl;
    }
    delete dists;
    delete nn_idx;
    for (int i = 0; i < num_point; i++){
        delete points[i];
    }
    delete[] points;
    return projected_traces;
}

int main(int argc, char** argv){
    if (argc < 4){
        cerr << "./dds <trace_filename> <delay> <K> <output_file>" << endl;
        exit(1);
    }
    vector<vector<Point>> traces = read_trace(argv[1]);
    int delay = atoi(argv[2]);
    cerr << "delay=" << delay << endl;
    int K = atoi(argv[3]);
    int d = 2;
    for (int iter = 0; iter < 100; iter++){ 
        vector<vector<Point>> projected_traces = smoothing(traces, delay, K);
        string name(argv[4]);
        name = name + "_iter" + to_string(iter) + ".out";
        //cerr << "printing to " << name << "...";
        ofstream fout(name, fstream::out);
        for (int i = 0; i < projected_traces.size(); i++){
            for (int j = 0; j < projected_traces[i].size(); j++){
                if (j != 0){
                    fout << " ";
                }
                Point & p_ij = projected_traces[i][j];
                //cerr << "(" << p_ij.size() << ")"<< endl;
                for (int t = 0; t < d; t++){
                    if (t != 0){
                        fout << " ";
                    }
                    fout << projected_traces[i][j][t];
                }
            }
            fout << endl;
        }
        fout.close();
        if (verbose){
            cerr << "done" << endl;
        }
        traces = projected_traces;
    }
}
