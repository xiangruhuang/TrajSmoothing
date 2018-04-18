#ifndef GRAPH_H
#define GRAPH_H

#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "Trace.h"
#include <ANN/ANN.h>
#include <omp.h>
#include "Params.h"

using namespace std;

typedef pair<int, int> FID;
struct Match{
    FID a, b;
    Float weight;
};

class Graph{
    public:
        //void construct_via_euclidean(){
        //    num_point = 0;
        //    for (int i = 0; i < traces.size(); i++){
        //        offset.push_back(num_point);
        //        num_point += traces[i].size();
        //    }
        //    D = traces[0].D;
        //    points = new double*[num_point];
        //    int count = 0;
        //    for (int i = 0; i < traces.size(); i++){
        //        vector<Point>& m = traces[i].traces;
        //        for (int j = 0; j < m.size(); j++){
        //            points[count] = new double[D];
        //            map.push_back(make_pair(i, j));
        //            Point& p_j = m[j];
        //            double* p = points[count];
        //            for (int d = 0; d < D; d++){
        //                p[d] = p_j[d];
        //            }
        //            count++;
        //        }
        //    }
        //}
        //void construct_via_direction(){
        //    cerr << "constructing via directions" << ", ord=" << ord << endl;
        //    num_point = 0;
        //    D = traces[0].D * (2*ord+1);
        //    for (int i = 0; i < traces.size(); i++){
        //        if (traces[i].size() - 2 * ord <= 0){
        //            cerr << "Graph.h: HumanMotion " << i << " is too short, #frames=" << traces[i].size() << ", ord=" << ord << endl;
        //            exit(1);
        //        }
        //        offset.push_back(num_point);
        //        num_point += traces[i].size() - 2 * ord;
        //    }
        //    cerr << "generating points" << endl;
        //    points = new double*[num_point];
        //    int count = 0;
        //    for (int i = 0; i < traces.size(); i++){
        //        vector<Point> m = traces[i].directions(ord);
        //        for (int j = 0; j < m.size(); j++){
        //            points[count] = new double[D];
        //            map.push_back(make_pair(i, j+ord));
        //            Point& p_j = m[j];
        //            if (p_j.size() != D){
        //                cerr << p_j.size() << " " << D << " " << traces[0].D << endl;
        //            }
        //            assert(p_j.size() == D);
        //            double* p = points[count];
        //            for (int d = 0; d < D; d++){
        //                p[d] = p_j[d];
        //            }
        //            count++;
        //        }
        //    }
        //}
        //
        

        void construct_via_equal_dist(){
            cerr << "traces.size = " << traces.size() << endl;
            Float avg_norm = 0.0;
            Float norm_down = 0.0;
            for (int i = 0; i < traces.size(); i++){
                Trace* m = traces[i];
                for (int j = 1; j < m->traces.size(); j++){
                    Float dist = norm(m->traces[j] - m->traces[j-1]);
                    avg_norm += dist;
                }
                //cerr << m->traces.size() << endl;
                norm_down += m->traces.size() - 1.0;
            }
            cerr << "average distance (in l2 norm) = " << avg_norm / norm_down << endl;
            cerr << "constructing via equal dist" << ", ord=" << ord << ", dist=" << interval << endl;
            D = traces[0]->D * (2*ord+1);
            cerr << "generating points, D=" << D << endl;
            vector<vector<pair<int, Point>>> rec;
            rec.resize(traces.size());
            cerr << traces.size() << endl;
            #pragma omp parallel for
            for (int i = 0; i < traces.size(); i++){
                vector<pair<int, Point>> m = traces[i]->equal_dist_samples(ord, interval, 1);
                rec[i] = m;
            }
            num_point = 0;
            for (int i = 0; i < traces.size(); i++){
                num_point += rec[i].size();
                offset.push_back(num_point);
            }
            cerr << "got in total " << num_point << " points " << endl;
            int count = 0;
            points = new double*[num_point];
            for (int i = 0; i < traces.size(); i++){
                vector<pair<int, Point>>& m = rec[i];
                for (int j = 0; j < m.size(); j++){
                    points[count] = new double[D];
                    map.push_back(make_pair(i, m[j].first));
                    Point& p_j = m[j].second;
                    if (p_j.size() != D){
                        cerr << p_j.size() << " " << D << " " << traces[0]->D << endl;
                    }
                    assert(p_j.size() == D);
                    double* p = points[count];
                    for (int d = 0; d < D; d++){
                        p[d] = p_j[d];
                    }
                    count++;
                }
            }
            rec.clear();
        }

        void construct_via_interpolate(){
            int num_inter = 50;
            cerr << "traces.size = " << traces.size() << endl;
            Float avg_norm = 0.0;
            Float norm_down = 0.0;
            for (int i = 0; i < traces.size(); i++){
                Trace* m = traces[i];
                for (int j = 1; j < m->traces.size(); j++){
                    Float dist = norm(m->traces[j] - m->traces[j-1]);
                    avg_norm += dist;
                }
                //cerr << m->traces.size() << endl;
                norm_down += m->traces.size() - 1.0;
            }
            cerr << "average distance (in l2 norm) = " << avg_norm / norm_down << endl;
            Float interval = 0.1;
            cerr << "constructing via equal dist" << ", ord=" << ord << ", dist=" << interval << endl;
            D = traces[0]->D * (2*ord+1);
            cerr << "generating points, D=" << D << endl;
            vector<vector<pair<int, Point>>> rec;
            rec.resize(traces.size());
            cerr << traces.size() << endl;
            cerr << "interpolating" << endl;
            #pragma omp parallel for
            for (int i = 0; i < traces.size(); i++){
                cerr << i << "/" << traces.size() << endl;
                vector<pair<int, Point>> m = traces[i]->interpolate(num_inter, all_point->knn, ord, interval, "GPS/dump/"+to_string(i)+".txt");
                //vector<pair<int, Point>> m = traces[i]->equal_dist_samples(ord, interval);
                rec[i] = m;
            }
            num_point = 0;
            for (int i = 0; i < traces.size(); i++){
                num_point += rec[i].size();
                offset.push_back(num_point);
            }
            cerr << "got in total " << num_point << " points " << endl;
            int count = 0;
            points = new double*[num_point];
            for (int i = 0; i < traces.size(); i++){
                vector<pair<int, Point>>& m = rec[i];
                for (int j = 0; j < m.size(); j++){
                    points[count] = new double[D];
                    map.push_back(make_pair(i, m[j].first));
                    Point& p_j = m[j].second;
                    if (p_j.size() != D){
                        cerr << p_j.size() << " " << D << " " << traces[0]->D << endl;
                    }
                    assert(p_j.size() == D);
                    double* p = points[count];
                    for (int d = 0; d < D; d++){
                        p[d] = p_j[d];
                    }
                    count++;
                }
            }
            rec.clear();
        }
        void construct(){
            if (num_point != 0){
                for (int i = 0; i < num_point; i++){
                    delete[] points[i];
                }
                delete[] points;
            }
            bool flag = false;
            //if (method == "direction"){
            //    construct_via_direction();
            //    flag = true;
            //} 
            //if (method == "euclidean"){
            //    construct_via_euclidean();
            //    flag = true;
            //} 
            if (method == "equal_dist"){
                construct_via_equal_dist();
                flag = true;
            }
            if (method == "interpolate"){
                all_point = new Graph(traces, "equal_dist");
                all_point->construct();
                construct_via_interpolate();
                flag = true;
            }
            if (!flag){
                cerr << "unknown method: " << method << endl;
                cerr << "candidates are: " << endl;
                cerr << "\tdirections" << endl;
                cerr << "\teuclidean" << endl;
                cerr << "\tequal_dist" << endl;
                exit(1);
            }
            knn = new ANNkd_tree(points, num_point, D); 
        }
        Graph(vector<Trace*> _traces, Params* params){
            traces = _traces;
            ord = params->ord;
            method = params->method;
            interval = params->interval;
        }
        Graph(vector<Trace*> _traces, string _method){
            traces = _traces;
            ord = 0;
            method = _method;
        }
        ~Graph(){
            for (int i = 0; i < num_point; i++){
                delete[] points[i];
            }
            delete[] points;
            delete knn;
        }
       
        void center_traces(){
            assert(traces[0]->D != -1);
            mean_pose = zero_point(traces[0]->D);
            int num_pose = 0;
            for (int i = 0; i < traces.size(); i++){
                for (int j = 0; j < traces[i]->traces.size(); j++){
                    mean_pose = mean_pose + traces[i]->traces[j];
                }
                num_pose += traces[i]->traces.size();
            }
            mean_pose = mean_pose / num_pose;
            Point var = zero_point(traces[0]->D);
            for (int i = 0; i < traces.size(); i++){
                for (int j = 0; j < traces[i]->traces.size(); j++){
                    Point diff = traces[i]->traces[j]-mean_pose;
                    Point diffsq = diff * diff;
                    var = var + diffsq;
                }
            }
            var = var / num_pose;
            std = sqrt_point(var);
            for (int i = 0; i < traces.size(); i++){
                traces[i]->center(mean_pose, std); 
            }
        }

        

        vector<pair<int, int>> search(int K, Point point){
            int* nn_idx = new int[K];
            double* dists = new double[K];
            double* query = to_double(point);
            knn->annkSearch(points[0], K, nn_idx, dists, eps);
            vector<pair<int, int>> ans;
            for (int i = 0; i < K; i++){
                ans.push_back(map[nn_idx[i]]);
            }
            
            delete[] nn_idx;
            delete[] dists;
            delete[] query;

            return ans;
        }
        
        inline int id(FID fid){
            return offset[fid.first] + fid.second + ord;
        }

        void compute_matchings(int K, string output_folder){
            cerr << "Compute matchings, K=" << K << endl;
            matchings.clear();
            for (int i = 0; i < num_point; i++){
                vector<Match> match_i;
                matchings.push_back(match_i);
            }
            matchings_per_motion.clear();
            for (int i = 0; i < traces.size(); i++){
                vector<Match> match_i;
                matchings_per_motion.push_back(match_i);
            }
            //ifstream fin(output_folder+"/matchings");
            //if (fin.good()){
            //    cerr << "reading from " << output_folder + "/matchings" << endl;
            //    while (!fin.eof()){
            //        Match m;
            //        fin >> m.a.first >> m.a.second >> m.b.first >> m.b.second;
            //        fin >> m.weight;
            //        matchings_per_motion[m.a.first].push_back(m);
            //    }
            //    return;
            //}
            cerr << "doing knn search" << endl;
            #pragma omp parallel for
            for (int i = 0; i < num_point; i++){
                vector<Match>& match_i = matchings[i];
                int found;
                int tK = K;
                if (i % 1000 == 0){
                    cerr << i << " / " << num_point << endl;
                }
                do{
                    int* nn_idx = new int[tK];
                    double* dists = new double[tK];
                    found = 0;
                    knn->annkSearch(points[i], tK, nn_idx, dists, eps);
                    for (int k = 0; k < tK; k++){
                        Match m;
                        if (map[i].first == map[nn_idx[k]].first){
                            //if (abs(map[i].second-map[nn_idx[k]].second) <= 10){
                                continue;
                            //}
                        }
                        found++;
                        m.a = map[i];
                        m.b = map[nn_idx[k]];
                        m.weight = dists[k];
                        //m.dist = distance(points[i], points[nn_idx[k]], D);
                        match_i.push_back(m);
                        if (found == K){
                            break;
                        }
                    }
                    if (found < K){
                        match_i.clear();
                        tK = tK * 2;
                    }
                    delete[] nn_idx;
                    delete[] dists;
                } while (found < K);
            }
            cerr << "done." << endl; 
            if (normalize_factor < 0.0){
                normalize_factor = 0.0;
                Float down = 1.0;
                for (int i = 0; i < num_point; i++){
                    vector<Match>& match_i = matchings[i];
                    for (int k = 0; k < K; k++){
                        matchings_per_motion[match_i[k].a.first].push_back(match_i[k]);
                        Match m;
                        m.a = match_i[k].b;
                        m.b = match_i[k].a;
                        m.weight = match_i[k].weight;
                        matchings_per_motion[match_i[k].b.first].push_back(m);
                        normalize_factor += match_i[k].weight;
                        down += 1.0;
                    }
                }
                normalize_factor /= (down * 1.0);
            }
            cerr << "normalizing" << endl;
            // normalization
            for (int i = 0; i < num_point; i++){
                vector<Match>& match_i = matchings[i];
                for (int j = 0; j < match_i.size(); j++){
                    match_i[j].weight = exp(-match_i[j].weight / normalize_factor);
                }
            }
            for (int i = 0; i < traces.size(); i++){
                vector<Match>& match_i = matchings_per_motion[i];
                for (int j = 0; j < match_i.size(); j++){
                    match_i[j].weight = 1.0; //exp(-match_i[j].weight / normalize_factor);
                }
            }

            ofstream fout("GPS/"+output_folder + "/matchings");
            cerr << "dumping matchings to " << ("GPS/"+output_folder + "/matchings") << endl;

            for (int i = 0; i < num_point; i++){
                vector<Match>& match_i = matchings[i];
                for (vector<Match>::iterator it = match_i.begin(); it != match_i.end(); it++){

                    fout << it->a.first << " " << it->a.second << " " << it->b.first << " " << it->b.second << " " << it->weight << endl;
                }
            }

            fout.close();

        }

        vector<pair<int, int>> search(int K, double* point){
            return search(K, to_vector(D, point));
        }
       
        Float normalize_factor = -10000;
        vector<vector<Match>> matchings;
        vector<vector<Match>> matchings_per_motion;
        vector<Trace*> traces;
        vector<FID> map;
        vector<int> offset;
        double** points = NULL;
        int num_point = 0, D;
        ANNkd_tree* knn;
        double eps = 1e-3;
        int ord = 0;
        Point mean_pose;
        Point std;
        Point ANN;
        string method;
        Graph* all_point;
        Float interval;
};

#endif
