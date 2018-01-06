#ifndef GD_H
#define GD_H

#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>

using namespace std;
// Objective function L(f):
// smooth = 0.5 * sum_{i,j} \|(f_{i,j} - f_{i,j+1}) - (g_{i,j} - g_{i, j+1}) \|_2^2
// align = 0.5 * sum_{i,j,p} \|f_{i,p.first} - f_{j, p.second}\|_2^2 = \|PF\|_2^2
// vanilla = 0.5 * sum_{i,j} \|f_{i,j} - g_{i,j}\|_2^2
// Objective = C_smooth * smooth + C_align * align + 1.0 * vanilla

class GDSolver{
    public:

        GDSolver(Float _lambda, int _max_iter){
            this->lambda = _lambda;
            this->max_iter = _max_iter;
        }
        void compute_LCS(vector<vector<Point>>& traces){
            int count = 0;
            float avg_len = 0;
            for (int i = 0; i < N; i++){
                for (int j = i+1; j < N; j++){
                    avg_len += min(traces[i].size(), traces[j].size());
                }
            }
            avg_len /= (N * (N - 1) / 2);
            match.clear();
            float avg_match = 0;
            for (int i = 0; i < N; i++){
                for (int j = i+1; j < N; j++){
                    vector<pair<int, int>> match_ij = LCS(traces[i], traces[j], 2.1e-2, 10);
                    avg_match += match_ij.size();
                    count++;
                    match.insert(make_pair(i*N+j, match_ij));
                }
            }
            cout << avg_match/count << "/" << avg_len << endl;
        }
        void solve(vector<vector<Point>>& traces){
            N = traces.size();
            Float C_smooth = 1.0;
            Float C_align = 1.0;
            Float C_vanilla = 1.0;
            compute_LCS(traces);
            vector<vector<Point>> recover = traces;
            vector<vector<Point>> grads;
            for (int i = 0; i < N; i++){
                vector<Point> grad;
                for (int t = 0; t < traces[i].size(); t++){
                    grad.push_back(make_pair(0.0, 0.0));
                }
                grads.push_back(grad);
            }
            int iter = 0;
            vector<int> indices;
            for (int i = 0; i < N; i++){
                indices.push_back(i);
            }
            while (iter < max_iter){
                random_shuffle(indices.begin(), indices.end());
                Float grad_l2 = 0.0;
                for (int _i = 0; _i < N; _i++){
                    int i = indices[_i];
                    int T = traces[i].size();
                    vector<Point>& grad = grads[i];
                    for (int t = 0; t < grad.size(); t++){
                        grad[t] = make_pair(0.0, 0.0);
                    }
                    // smooth
                    for (int t = 1; t < traces[i].size(); t++){
                        Point diff = (recover[i][t] - recover[i][t-1]) - (traces[i][t] - traces[i][t-1]);
                        grad[t] = grad[t] + diff * C_smooth;
                        grad[t-1] = grad[t-1] - diff * C_smooth;
                    }
                    // align
                    //vector<int> count;
                    //for (int j = 0; j < T; j++){
                    //    count.push_back(0);
                    //}
                    for (int j = 0; j < N; j++){
                        if (i == j){
                            continue;
                        }
                        vector<pair<int, int>> match_ij;
                        if (i < j){
                            match_ij = match.find(i*N + j)->second;
                        } else {
                            match_ij = match.find(j*N + i)->second;
                        }
                        for (int p = 0; p < match_ij.size(); p++){
                            int a = match_ij[p].first, b = match_ij[p].second;
                            if (i > j){
                                int temp = a; a = b; b = temp;
                            }
                            grad[a] = grad[a] + (recover[i][a] - recover[j][b]) * C_align;
                            //count[a]++;
                        }
                    }
                    //int nonzero = 0;
                    //for (int j = 0; j < T; j++){
                    //    if (count[j] > 3){
                    //        nonzero++;
                    //    }
                    //}
                    //cout << ((nonzero * 1.0) / T) << endl;
                    // vanilla
                    for (int j = 0; j < traces[i].size(); j++){
                        if (j % (11) != 0){
                            continue;
                        }
                        grad[j] = grad[j] + (recover[i][j] - traces[i][j]) * C_vanilla;
                    }
                    //Update
                    for (int j = 0; j < recover[i].size(); j++){
                        grad_l2 += normsq(grad[j]);
                        recover[i][j] = recover[i][j] - grad[j] * lambda;
                    }
                }

                //Compute Energy
                //Vanilla Part
                Float energy_vanilla = 0.0;
                for (int i = 0; i < N; i++){
                    for (int j = 0; j < traces[i].size(); j++){
                        if (j % (11) != 0){
                            continue;
                        }
                        energy_vanilla += normsq(traces[i][j] - recover[i][j]) * C_vanilla * 0.5;
                    }
                }
                //Smooth Part
                Float energy_smooth = 0.0;
                for (int i = 0; i < N; i++){
                    for (int t = 1; t < traces[i].size(); t++){
                        Point diff = (recover[i][t] - recover[i][t-1]) - (traces[i][t] - traces[i][t-1]);
                        energy_smooth += normsq(diff) * C_smooth * 0.5;
                    }
                }
                //Align Part
                Float energy_align = 0.0;
                for (int i = 0; i < N; i++){
                    for (int j = i+1; j < N; j++){
                        vector<pair<int, int>> match_ij;
                        match_ij = (match.find(i*N + j))->second;
                        for (int p = 0; p < match_ij.size(); p++){
                            int a = match_ij[p].first, b = match_ij[p].second;
                            energy_align += normsq(recover[i][a] - recover[j][b]) * C_align * 0.5;
                        }
                    } 
                }
                cout << "iter=" << iter;
                cout << ", grad_l2=" << grad_l2;
                cout << ", vanilla=" << energy_vanilla;
                cout << ", smooth=" << energy_smooth;
                cout << ", align=" << energy_align;
                cout << ", energy=" << (energy_vanilla + energy_smooth + energy_align);
                cout << endl;
                iter++;
                if (iter % 1000 == 0){
                    ofstream fout("recover/"+to_string(iter)+".txt", fstream::out);
                    for (int i = 0; i < N; i++){
                        for (int p = 0; p < recover[i].size(); p++){
                            if (p != 0){
                                fout << " ";
                            }
                            fout << recover[i][p].first << " " << recover[i][p].second;
                        }
                        fout << endl;
                    }
                    fout.close();
                }
            }
        }

        int N;
        unordered_map<int, vector<pair<int,int>>> match;
        vector<vector<Point>> recover;
        Float lambda;
        int max_iter;
};

#endif
