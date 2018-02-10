#ifndef GD_H
#define GD_H

#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include<string>
#include "Motion.h"
#include "Graph.h"

using namespace std;
// Objective function L(f):
// smooth = 0.5 * sum_{i,j} \|(f_{i,j} - f_{i,j+1}) - (g_{i,j} - g_{i, j+1}) \|_2^2
// align = 0.5 * sum_{i,j,p} \|f_{i,p.first} - f_{j, p.second}\|_2^2 = \|PF\|_2^2
// vanilla = 0.5 * sum_{i,j} \|f_{i,j} - g_{i,j}\|_2^2
// Objective = C_smooth * smooth + C_align * align + 1.0 * vanilla

class GDSolver{
    public:

        GDSolver(Float _lambda){
            this->lambda = _lambda;
        }

        //void read_correspondences(vector<vector<Point>>& traces, ifstream& fin){
        //    match.clear();
        //    cerr << "reading matchings from " << workdir + "/matchings" << endl;
        //    while (!fin.eof()){
        //        int i, j, T;
        //        fin >> i >> j >> T;
        //        //cerr << i << " " << j << " " << T << " " << traces[i].size() << " " << traces[j].size() << endl;
        //        vector<pair<pair<int, int>, Float>> match_ij;
        //        for (int t = 0; t < T; t++){
        //            int a,b;
        //            Float dist;
        //            fin >> a >> b >> dist;
        //            //cout << b << " " << traces[j].size() << endl;
        //            //cerr << b << " " << traces[j].size() << endl;
        //            assert(a < traces[i].size());
        //            assert(b < traces[j].size());
        //            match_ij.push_back(make_pair(make_pair(a, b), dist));
        //        }
        //        match.insert(make_pair(i*N+j, match_ij));
        //    }
        //    cerr << "done" << endl;
        //    //vector<string> lines = readlines(fin);
        //    //cerr << "... done" << endl;
        //    //match.clear();
        //    //vector<pair<pair<int, int>, Float>> match_ij;
        //    //for (int i = 0; i < lines.size(); i++){
        //    //    string& line = lines[i];
        //    //    vector<string> tokens = split(line, " ");
        //    //    int a = stoi(tokens[0]);
        //    //    int b = stoi(tokens[1]);
        //    //    int T = stoi(tokens
        //    //    
        //    //    vector<pair<pair<int, int>, Float>> match_ij;

        //    //    for (int j = 0; j < tokens.size() / 3; j++){
        //    //        int a = stoi(tokens[j*3+0]);
        //    //        int b = stoi(tokens[j*3+1]);
        //    //        Float dist = stod(tokens[j*3+2]);
        //    //        match_ij.push_back(make_pair(make_pair(a, b), dist));
        //    //    }
        //    //    match.insert(make_pair(i*N+j, match_ij));

        //    //}
        //    fin.close();
        //    
        //}
        
        //void compute_correspondences(vector<vector<Point>>& traces, bool use_LCS){
        //    ifstream fin(workdir+"/matchings");
        //    if (fin.good()){
        //        read_correspondences(traces, fin);
        //        return;
        //    }
        //    int count = 0;
        //    float avg_len = 0;
        //    for (int i = 0; i < N; i++){
        //        for (int j = i+1; j < N; j++){
        //            avg_len += min(traces[i].size(), traces[j].size());
        //        }
        //    }
        //    avg_len /= (N * (N - 1) / 2);
        //    match.clear();
        //    float avg_match = 0;
        //    ofstream fout(workdir+"/matchings"/*+to_string(threshold)*/, fstream::out);
        //    for (int i = 0; i < N; i++){
        //        for (int j = i+1; j < N; j++){
        //            vector<pair<pair<int, int>, Float>> match_ij;
        //            if (use_LCS){
        //                match_ij = LCS(traces[i], traces[j], threshold, period);
        //            } else {
        //                match_ij = gaussian(traces[i], traces[j], threshold); 
        //            }
        //            fout << i << " " << j << " " << match_ij.size() << endl;
        //            
        //            for (int t = 0; t < match_ij.size(); t++){
        //                if (t != 0){
        //                    fout << " ";
        //                }
        //                pair<int, int> pair = match_ij[t].first;
        //                int a = pair.first;
        //                int b = pair.second;
        //                fout << a << " " << b << " " << distance(traces[i][a], traces[j][b]);
        //            }
        //            fout << endl;
        //            avg_match += match_ij.size();
        //            count++;
        //            match.insert(make_pair(i*N+j, match_ij));
        //            cout << avg_match/count << "/" << avg_len << endl;
        //        }
        //    }
        //    fout.close();
        //    //auto temp = match;
        //    //ifstream fin2(workdir+"/matchings");
        //    //read_correspondences(traces, fin2);
        //    //for (int i = 0; i < N; i++){
        //    //    for (int j = i+1; j < N; j++){
        //    //        cerr << i << " " << j << endl;
        //    //        auto match_ij = match.find(i*N+j)->second;
        //    //        auto temp_ij = temp.find(i*N+j)->second;
        //    //        cerr << match_ij.size() << " " << temp_ij.size() << endl;
        //    //        assert(match_ij.size() == temp_ij.size());
        //    //        for (int t = 0; t < match_ij.size(); t++){
        //    //            cerr << t << ":";
        //    //            cerr << match_ij[t].first.first << " " << temp_ij[t].first.first << endl;
        //    //            cerr << match_ij[t].first.second << " " << temp_ij[t].first.second << endl;
        //    //            cerr << match_ij[t].second << " " << temp_ij[t].second << endl;
        //    //            assert(match_ij[t].first.first == temp_ij[t].first.first);
        //    //            assert(match_ij[t].first.second == temp_ij[t].first.second);
        //    //            assert(fabs(match_ij[t].second-temp_ij[t].second) < 1e-2);
        //    //        }
        //    //    }
        //    //}
        //    
        //    cout << avg_match/count << "/" << avg_len << endl; 
        //}

        void solve_motion(Graph& graph, int max_iter, string output_folder){
            vector<Motion>& motions = graph.motions;
            int N = motions.size();
            vector<vector<Point>> traces;
            for (int i = 0; i < N; i++){
                traces.push_back(motions[i].rotations);
            }
            int D = traces[0][0].size();
            if (shift){
                mean = zero_point(D);
                Float count = 0;
                for (int i = 0; i < N; i++){
                    for (int t = 0; t < traces[i].size(); t++){
                        mean = mean + traces[i][t];
                        count += 1.0;
                    }
                }
                mean = mean / count;
                std = zero_point(D) + 1e-5;
                for (int i = 0; i < N; i++){
                    for (int t = 0; t < traces[i].size(); t++){
                        Point diff = traces[i][t] - mean;
                        Point diff_sq = diff * diff;
                        std = std + diff_sq;
                    }
                }
                std = std / count;
                std = sqrt_point(std);
                for (int i = 0; i < N; i++){
                    for (int t = 0; t < traces[i].size(); t++){
                        traces[i][t] = traces[i][t] / std;
                    }
                }
                print_point(mean, "mean");
                print_point(std, "std");
            }
            
            Float C_smooth = 100.0;
            Float C_align = 1.0;
            Float C_vanilla = 1.0;
            //compute_correspondences(traces, false);
            vector<vector<Point>> recover = traces;
            vector<vector<Point>> grads;
            for (int i = 0; i < N; i++){
                vector<Point> grad;
                for (int t = 0; t < traces[i].size(); t++){
                    grad.push_back(zero_point(D));
                }
                grads.push_back(grad);
            }
            int iter = 0;
            vector<int> indices;
            for (int i = 0; i < N; i++){
                indices.push_back(i);
            }
            vector<vector<Match>>& matchings = graph.matchings_per_motion;
            while (iter < max_iter){
                random_shuffle(indices.begin(), indices.end());
                Float grad_l2 = 0.0;
                #pragma omp parallel for
                for (int _i = 0; _i < N; _i++){
                    int i = indices[_i];
                    int T = traces[i].size();
                    vector<Point>& grad = grads[i];
                    for (int t = 0; t < grad.size(); t++){
                        grad[t] = zero_point(D);
                    }
                    // smooth
                    for (int t = 1; t < traces[i].size(); t++){
                        Point diff = (recover[i][t] - recover[i][t-1]) - (traces[i][t] - traces[i][t-1]);
                        grad[t] = grad[t] + diff * C_smooth;
                        grad[t-1] = grad[t-1] - diff * C_smooth;
                    }
                    // align
                    for (vector<Match>::iterator it = matchings[i].begin(); it != matchings[i].end(); it++){
                        assert(i == it->a.first);
                        int fi = it->a.second;
                        int j = it->b.first;
                        int fj = it->b.first;
                        grad[fi] = grad[fi] + (recover[i][fi] - recover[j][fj]) * C_align * it->weight;
                    }
                    // vanilla
                    for (int j = 0; j < traces[i].size(); j++){
                        grad[j] = grad[j] + (recover[i][j] - traces[i][j]) * C_vanilla;
                    }
                }

                //Update
                #pragma omp parallel for
                for (int i = 0; i < N; i++){
                    vector<Point>& grad = grads[i];
                    for (int j = 0; j < recover[i].size(); j++){
                        #pragma omp atomic
                        grad_l2 += normsq(grad[j]);
                        recover[i][j] = recover[i][j] - grad[j] * lambda;
                    }

                }

                //Compute Energy
                //Vanilla Part
                Float energy_vanilla = 0.0;
                #pragma omp parallel for
                for (int i = 0; i < N; i++){
                    int j;
                    for (j = 0; j < traces[i].size(); j++){
                        #pragma omp atomic
                        energy_vanilla += normsq(traces[i][j] - recover[i][j]) * C_vanilla * 0.5;
                    }
                }
                //Smooth Part
                Float energy_smooth = 0.0;
                #pragma omp parallel for
                for (int i = 0; i < N; i++){
                    int t;
                    for (t = 1; t < traces[i].size(); t++){
                        Point diff = (recover[i][t] - recover[i][t-1]) - (traces[i][t] - traces[i][t-1]);
                        #pragma omp atomic
                        energy_smooth += normsq(diff) * C_smooth * 0.5;
                    }
                }
                //Align Part
                Float energy_align = 0.0;
                #pragma omp parallel for
                for (int i = 0; i < N; i++){
                    vector<Match>::iterator it;
                    for (it = matchings[i].begin(); it != matchings[i].end(); it++){
                        assert(i == it->a.first);
                        int fi = it->a.second;
                        int j = it->b.first;
                        int fj = it->b.first;
                        #pragma omp atomic
                        energy_align += normsq(recover[i][fi] - recover[j][fj]) * C_align * it->weight * 0.5;
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
                if (iter % 100 == 0){
                    if (shift){
                        for (int i = 0; i < N; i++){
                            for (int t = 0; t < traces[i].size(); t++){
                                traces[i][t] = traces[i][t] * std;
                                recover[i][t] = recover[i][t] * std;
                            }
                        }
                    }
                    ofstream smooth_list("motion/smooth/smoothed_" + to_string(iter)+".txt");
                    for (int i = 0; i < traces.size(); i++){
                        motions[i].rotations = recover[i];
                        motions[i].write_bvh(replace(motions[i].filename, "data", "smooth/"+to_string(iter)));
                        smooth_list << replace(motions[i].filename, "data", "smooth/"+to_string(iter)) << endl;
                    }
                    smooth_list.close();
                    if (shift){
                        for (int i = 0; i < N; i++){
                            for (int t = 0; t < traces[i].size(); t++){
                                traces[i][t] = traces[i][t] / std;
                                recover[i][t] = recover[i][t] / std;
                            }
                        }
                    }
                }
            }
            if (shift){
                for (int i = 0; i < N; i++){
                    for (int t = 0; t < traces[i].size(); t++){
                        traces[i][t] = traces[i][t] * std;
                        recover[i][t] = recover[i][t] * std;
                    }
                }
            }
        }
        //void solve(vector<vector<Point>>& traces){
        //    N = traces.size();
        //    mean = zero_point(D);
        //    Float count = 0;
        //    for (int i = 0; i < N; i++){
        //        for (int t = 0; t < traces[i].size(); t++){
        //            mean = mean + traces[i][t];
        //            count += 1.0;
        //        }
        //    }
        //    mean = mean / count;
        //    std = zero_point(D) + 1e-5;
        //    for (int i = 0; i < N; i++){
        //        for (int t = 0; t < traces[i].size(); t++){
        //            Point diff = traces[i][t] - mean;
        //            Point diff_sq = diff * diff;
        //            std = std + diff_sq;
        //        }
        //    }
        //    std = std / count;
        //    std = sqrt_point(std);
        //    for (int i = 0; i < N; i++){
        //        for (int t = 0; t < traces[i].size(); t++){
        //            traces[i][t] = traces[i][t] / std;
        //        }
        //    }
        //    print_point(mean, "mean");
        //    print_point(std, "std");
        //    
        //    Float C_smooth = 1e-1;
        //    Float C_align = 10;
        //    Float C_vanilla = 1.0;
        //    compute_correspondences(traces, true);
        //    vector<vector<Point>> recover = traces;
        //    vector<vector<Point>> grads;
        //    for (int i = 0; i < N; i++){
        //        vector<Point> grad;
        //        for (int t = 0; t < traces[i].size(); t++){
        //            grad.push_back(zero_point(D));
        //        }
        //        grads.push_back(grad);
        //    }
        //    int iter = 0;
        //    vector<int> indices;
        //    for (int i = 0; i < N; i++){
        //        indices.push_back(i);
        //    }
        //    while (iter < max_iter){
        //        random_shuffle(indices.begin(), indices.end());
        //        Float grad_l2 = 0.0;
        //        for (int _i = 0; _i < N; _i++){
        //            int i = indices[_i];
        //            int T = traces[i].size();
        //            vector<Point>& grad = grads[i];
        //            for (int t = 0; t < grad.size(); t++){
        //                grad[t] = zero_point(D);
        //            }
        //            // smooth
        //            for (int t = 1; t < traces[i].size(); t++){
        //                Point diff = (recover[i][t] - recover[i][t-1]) - (traces[i][t] - traces[i][t-1]);
        //                grad[t] = grad[t] + diff * C_smooth;
        //                grad[t-1] = grad[t-1] - diff * C_smooth;
        //            }
        //            // align
        //            //vector<int> count;
        //            //for (int j = 0; j < T; j++){
        //            //    count.push_back(0);
        //            //}
        //            for (int j = 0; j < N; j++){
        //                if (i == j){
        //                    continue;
        //                }
        //                vector<pair<pair<int, int>, Float>> match_ij;
        //                if (i < j){
        //                    match_ij = match.find(i*N + j)->second;
        //                } else {
        //                    match_ij = match.find(j*N + i)->second;
        //                }
        //                for (int p = 0; p < match_ij.size(); p++){
        //                    pair<int, int>& pair = match_ij[p].first;
        //                    int a = pair.first, b = pair.second;
        //                    if (i > j){
        //                        int temp = a; a = b; b = temp;
        //                    }
        //                    grad[a] = grad[a] + (recover[i][a] - recover[j][b]) * C_align;
        //                    //count[a]++;
        //                }
        //            }
        //            //int nonzero = 0;
        //            //for (int j = 0; j < T; j++){
        //            //    if (count[j] > 3){
        //            //        nonzero++;
        //            //    }
        //            //}
        //            //cout << ((nonzero * 1.0) / T) << endl;
        //            // vanilla
        //            for (int j = 0; j < traces[i].size(); j++){
        //                if (j % (period + 1) != 0){
        //                    continue;
        //                }
        //                grad[j] = grad[j] + (recover[i][j] - traces[i][j]) * C_vanilla;
        //            }
        //            //Update
        //            for (int j = 0; j < recover[i].size(); j++){
        //                grad_l2 += normsq(grad[j]);
        //                recover[i][j] = recover[i][j] - grad[j] * lambda;
        //            }
        //        }

        //        //Compute Energy
        //        //Vanilla Part
        //        Float energy_vanilla = 0.0;
        //        for (int i = 0; i < N; i++){
        //            for (int j = 0; j < traces[i].size(); j++){
        //                if (j % (period + 1) != 0){
        //                    continue;
        //                }
        //                energy_vanilla += normsq(traces[i][j] - recover[i][j]) * C_vanilla * 0.5;
        //            }
        //        }
        //        //Smooth Part
        //        Float energy_smooth = 0.0;
        //        for (int i = 0; i < N; i++){
        //            for (int t = 1; t < traces[i].size(); t++){
        //                Point diff = (recover[i][t] - recover[i][t-1]) - (traces[i][t] - traces[i][t-1]);
        //                energy_smooth += normsq(diff) * C_smooth * 0.5;
        //            }
        //        }
        //        //Align Part
        //        Float energy_align = 0.0;
        //        for (int i = 0; i < N; i++){
        //            for (int j = i+1; j < N; j++){
        //                vector<pair<pair<int, int>, Float>> match_ij;
        //                match_ij = (match.find(i*N + j))->second;
        //                for (int p = 0; p < match_ij.size(); p++){
        //                    pair<int, int> pair = match_ij[p].first;
        //                    int a = pair.first, b = pair.second;
        //                    energy_align += normsq(recover[i][a] - recover[j][b]) * C_align * 0.5;
        //                }
        //            } 
        //        }
        //        cout << "iter=" << iter;
        //        cout << ", grad_l2=" << grad_l2;
        //        cout << ", vanilla=" << energy_vanilla;
        //        cout << ", smooth=" << energy_smooth;
        //        cout << ", align=" << energy_align;
        //        cout << ", energy=" << (energy_vanilla + energy_smooth + energy_align);
        //        cout << endl;
        //        iter++;
        //        if (iter % 100 == 0){
        //            ofstream fout(workdir+"/"+to_string(iter)+".txt", fstream::out);
        //            for (int i = 0; i < N; i++){
        //                for (int p = 0; p < recover[i].size(); p++){
        //                    if (p != 0){
        //                        fout << " ";
        //                    }
        //                    Point rip = recover[i][p] * std;
        //                    for (int tt = 0; tt < rip.size(); tt++){
        //                        if (tt != 0){
        //                            fout << " ";
        //                        }
        //                        fout << rip[tt] ;
        //                    }
        //                }
        //                fout << endl;
        //            }
        //            fout.close();
        //        }
        //    }
        //    for (int i = 0; i < N; i++){
        //        for (int t = 0; t < traces[i].size(); t++){
        //            traces[i][t] = traces[i][t] * std;
        //            recover[i][t] = recover[i][t] * std;
        //        }
        //    }
        //}
        bool shift=false;
        Point mean, std;
        vector<vector<Point>> recover;
        Float lambda;
        string workdir;
};

#endif
