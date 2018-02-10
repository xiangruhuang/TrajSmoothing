#ifndef GRAPH_H
#define GRAPH_H

#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "Motion.h"
#include <ANN/ANN.h>
#include <omp.h>

using namespace std;

typedef pair<int, int> FID;
struct Match{
    FID a, b;
    Float weight;
};

class Graph{
    public:
        void construct_via_euclidean(){
            num_point = 0;
            for (int i = 0; i < motions.size(); i++){
                offset.push_back(num_point);
                num_point += motions[i].size();
            }
            D = motions[0].D;
            points = new double*[num_point];
            int count = 0;
            for (int i = 0; i < motions.size(); i++){
                vector<Point>& m = motions[i].rotations;
                for (int j = 0; j < m.size(); j++){
                    points[count] = new double[D];
                    map.push_back(make_pair(i, j));
                    Point& p_j = m[j];
                    double* p = points[count];
                    for (int d = 0; d < D; d++){
                        p[d] = p_j[d];
                    }
                    count++;
                }
            }
        }
        void construct_via_direction(){
            cerr << "constructing via directions" << ", ord=" << ord << endl;
            num_point = 0;
            D = motions[0].D * (2*ord+1);
            for (int i = 0; i < motions.size(); i++){
                if (motions[i].size() - 2 * ord <= 0){
                    cerr << "Graph.h: Motion " << i << " is too short, #frames=" << motions[i].size() << ", ord=" << ord << endl;
                    exit(1);
                }
                offset.push_back(num_point);
                num_point += motions[i].size() - 2 * ord;
            }
            cerr << "generating points" << endl;
            points = new double*[num_point];
            int count = 0;
            for (int i = 0; i < motions.size(); i++){
                vector<Point> m = motions[i].directions(ord);
                for (int j = 0; j < m.size(); j++){
                    points[count] = new double[D];
                    map.push_back(make_pair(i, j+ord));
                    Point& p_j = m[j];
                    if (p_j.size() != D){
                        cerr << p_j.size() << " " << D << " " << motions[0].D << endl;
                    }
                    assert(p_j.size() == D);
                    double* p = points[count];
                    for (int d = 0; d < D; d++){
                        p[d] = p_j[d];
                    }
                    count++;
                }
            }
        }
        void construct_via_equal_dist(){
            Float avg_norm = 0.0;
            Float norm_down = 0.0;
            for (int i = 0; i < motions.size(); i++){
                Motion& m = motions[i];
                for (int j = 1; j < m.rotations.size(); j++){
                    Float dist = norm(m.rotations[j] - m.rotations[j-1]);
                    avg_norm += dist;
                }
                norm_down += m.rotations.size() - 1.0;
            }
            cerr << "average distance (in l2 norm) = " << avg_norm / norm_down << endl;
            Float interval = 1.0;
            cerr << "constructing via circle" << ", ord=" << ord << ", dist=" << interval << endl;
            D = motions[0].D * (2*ord+1);
            cerr << "generating points" << endl;
            vector<vector<pair<int, Point>>> rec;
            rec.resize(motions.size());
            #pragma omp parallel for
            for (int i = 0; i < motions.size(); i++){
                vector<pair<int, Point>> m = motions[i].equal_dist_samples(ord, interval);
                rec[i] = m;
            }
            num_point = 0;
            for (int i = 0; i < motions.size(); i++){
                num_point += rec[i].size();
                offset.push_back(num_point);
            }
            cerr << "got in total " << num_point << " points " << endl;
            int count = 0;
            points = new double*[num_point];
            for (int i = 0; i < motions.size(); i++){
                vector<pair<int, Point>>& m = rec[i];
                for (int j = 0; j < m.size(); j++){
                    points[count] = new double[D];
                    map.push_back(make_pair(i, m[j].first));
                    Point& p_j = m[j].second;
                    if (p_j.size() != D){
                        cerr << p_j.size() << " " << D << " " << motions[0].D << endl;
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
        Graph(vector<Motion>& motions, string method="euclidean", int ord=1){
            this->motions = motions;
            this->ord = ord;
            center_motions();
            bool flag = false;
            if (method == "direction"){
                construct_via_direction();
                flag = true;
            } 
            if (method == "euclidean"){
                construct_via_euclidean();
                flag = true;
            } 
            if (method == "equal_dist"){
                construct_via_equal_dist();
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
        ~Graph(){
            for (int i = 0; i < num_point; i++){
                delete[] points[i];
            }
            delete[] points;
            delete knn;
        }
       
        void center_motions(){
            mean_pose = zero_point(motions[0].D);
            int num_pose = 0;
            for (int i = 0; i < motions.size(); i++){
                for (int j = 0; j < motions[i].rotations.size(); j++){
                    mean_pose = mean_pose + motions[i].rotations[j];
                }
                num_pose += motions[i].rotations.size();
            }
            mean_pose = mean_pose / num_pose;
            Point var = zero_point(motions[0].D);
            for (int i = 0; i < motions.size(); i++){
                for (int j = 0; j < motions[i].rotations.size(); j++){
                    Point diff = motions[i].rotations[j]-mean_pose;
                    Point diffsq = diff * diff;
                    var = var + diffsq;
                }
            }
            var = var / num_pose;
            std = sqrt_point(var);
            for (int i = 0; i < motions.size(); i++){
                motions[i].center(mean_pose, std);
            }
        }

        inline double* to_double(Point& point){
            int D = point.size();
            double* ans = new double[D];

            for (int d = 0; d < D; d++){
                ans[d] = point[d];
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
            for (int i = 0; i < motions.size(); i++){
                vector<Match> match_i;
                matchings_per_motion.push_back(match_i);
            }
            ifstream fin(output_folder+"/matchings");
            if (fin.good()){
                cerr << "reading from " << output_folder + "/matchings" << endl;
                while (!fin.eof()){
                    Match m;
                    fin >> m.a.first >> m.a.second >> m.b.first >> m.b.second;
                    fin >> m.weight;
                    matchings_per_motion[m.a.first].push_back(m);
                }
                return;
            }
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
            Float normalize_factor = 0.0;
            Float down = 1.0;
            for (int i = 0; i < num_point; i++){
                vector<Match>& match_i = matchings[i];
                for (int k = 0; k < K; k++){
                    matchings_per_motion[match_i[k].a.first].push_back(match_i[k]);
                    normalize_factor += match_i[k].weight;
                    down += 1.0;
                }
            }
            cerr << "normalizing" << endl;
            // normalization
            normalize_factor /= (down * 5);
            for (int i = 0; i < num_point; i++){
                vector<Match>& match_i = matchings[i];
                for (int j = 0; j < match_i.size(); j++){
                    match_i[j].weight = exp(-match_i[j].weight / normalize_factor);
                }
            }
            for (int i = 0; i < motions.size(); i++){
                vector<Match>& match_i = matchings_per_motion[i];
                for (int j = 0; j < match_i.size(); j++){
                    match_i[j].weight = exp(-match_i[j].weight / normalize_factor);
                }
            }

            ofstream fout(output_folder + "/matchings");

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
        
        vector<vector<Match>> matchings;
        vector<vector<Match>> matchings_per_motion;
        vector<Motion> motions;
        vector<FID> map;
        vector<int> offset;
        double** points;
        int num_point, D;
        ANNkd_tree* knn;
        double eps = 1e-3;
        int ord = 0;
        Point mean_pose, std;
};

#endif
