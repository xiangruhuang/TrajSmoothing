#ifndef AM_H
#define AM_H

#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>

#include "Graph.h"
#include "Params.h"
#include <omp.h>

using namespace std;
// Objective function L(f):
// smooth = 0.5 * sum_{i,j} \|(f_{i,j} - f_{i,j+1}) - (g_{i,j} - g_{i, j+1}) \|_2^2
// align = 0.5 * sum_{i,j,p} \|f_{i,p.first} - f_{j, p.second}\|_2^2 = \|PF\|_2^2
// vanilla = 0.5 * sum_{i,j} \|f_{i,j} - g_{i,j}\|_2^2
// Objective = C_smooth * smooth + C_align * align + 1.0 * vanilla

class AMSolver{
    public:

        AMSolver(Float _lambda){
            this->learning_rate = _lambda;
            iter = 0;
            traces.clear();
        }

        void shift(){
            int D = traces[0][0].size();
            mean = zero_point(D);
            Float count = 0;
            int N = traces.size();
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
                    recover[i][t] = recover[i][t] / std;
                }
            }
            //print_point(mean, "mean");
            //print_point(std, "std"); 
        }

        void shift_back(){
            int N = traces.size();
            for (int i = 0; i < N; i++){
                for (int t = 0; t < traces[i].size(); t++){
                    traces[i][t] = traces[i][t] * std;
                    recover[i][t] = recover[i][t] * std;
                }
            }
        }

        void solve(Graph& graph, Params* params){
            int max_iter = params->max_iter;
            string output_folder = params->output_folder;
            Float C_smooth = params->C_smooth;
            Float C_align = params->C_align;
            Float C_reg = params->C_reg;
            mean = graph.mean_pose;
            std = graph.std;
            


            int N = graph.traces.size();
            recover.clear();
            vector<Trace*>& motions = graph.traces;
            for (int i = 0; i < N; i++){
                recover.push_back(motions[i]->traces);
            }
            if (traces.size() == 0){
                traces.clear();
                for (int i = 0; i < N; i++){
                    traces.push_back(motions[i]->traces);
                }
            }
            int D = traces[0][0].size();

            vector<int> indices;
            for (int i = 0; i < N; i++){
                indices.push_back(i);
            }
            vector<vector<Match>>& matchings = graph.matchings_per_motion;
            vector<vector<Float>> reg_weight;
            for (int i = 0; i < N; i++){
                vector<Float> reg_weight_i;
                for (int j = 0; j < recover[i].size(); j++){
                    reg_weight_i.push_back(0.0);
                }
                reg_weight.push_back(reg_weight_i);
            }

            for (int I = 0; I < N; I++){
                vector<Match> next;
                for (vector<Match>::iterator it = matchings[I].begin(); it != matchings[I].end(); it++){
                    assert(I == it->a.first);
                    int fi = it->a.second;
                    int J = it->b.first;
                    int fj = it->b.second;
                    reg_weight[I][fi] += it->weight;
                }
                
            }
            
            //vector<vector<int>> adj;
            //bool* visit = new bool[N];
            //for (int i = 0; i < N; i++){
            //    visit[i] = false;
            //}
            //for (int i = 0; i < N; i++){
            //    vector<int> adj_i;
            //    adj_i.push_back(i);
            //    visit[i] = true;
            //    for (vector<Match>::iterator it = matchings[i].begin(); it != matchings[i].end(); it++){
            //        assert(i == it->a.first);
            //        int j = it->b.first;
            //        if (!visit[j]){
            //            adj_i.push_back(j);
            //        }
            //        visit[j] = true;
            //    }
            //    adj.push_back(adj_i);
            //    for (vector<Match>::iterator it = matchings[i].begin(); it != matchings[i].end(); it++){
            //        assert(i == it->a.first);
            //        int j = it->b.first;
            //        visit[j] = false;
            //    }
            //    visit[i] = false;
            //}
            //delete[] visit;
            //omp_lock_t* locks = new omp_lock_t[N];
            //for (int i = 0; i < N; i++){
            //    omp_init_lock(&(locks[i]));
            //}
            //omp_lock_t global;
            //omp_init_lock(&global); 
            //
            //cerr << "locks initialized" << endl;
            Float max_dist = 1e10;
            Float last_energy = 1e100;
            while (iter < max_iter){
                random_shuffle(indices.begin(), indices.end());
                //#pragma omp parallel for
                for (int _i = 0; _i < N; _i++){
                    int I = indices[_i];
                    //int I = rand() % N;
                    //bool flag = false;
                    
                    //while (!flag){
                    //    omp_set_lock(&global);
                    //    int j;
                    //    flag = true;
                    //    for (j = 0; j < adj[I].size(); j++){
                    //        if (!omp_test_lock(&(locks[adj[I][j]]))){
                    //            flag = false;
                    //            break;
                    //        }
                    //    }

                    //    if (!flag){
                    //        for (j = j-1; j >= 0; j--){
                    //            omp_unset_lock(&(locks[adj[I][j]]));
                    //        }
                    //        omp_unset_lock(&global);
                    //        I = rand() % N;
                    //        continue;
                    //    }

                    //    omp_unset_lock(&global);
                    //}
                        
                    //cerr << "entered" << omp_get_thread_num() << endl;

                    int T = traces[I].size();
                    assert(T > 1);
                    vector<Point>& u = recover[I];
                    vector<Point>& trace = traces[I];
                    // Summarize align and regularization part into this format
                    // below
                    // \sum_{i=1}^T a_i \|u_i\|_2^2 - 2 <v_i, u_i>
                    
                    ////////////////////////DEBUG 
                    //Float energy_before = 0.0;
                    //for (int i = 0; i < T; i++){
                    //    energy_before += normsq(u[i] - trace[i]) * C_reg * 0.5;
                    //}

                    //for (vector<Match>::iterator it = matchings[I].begin(); it != matchings[I].end(); it++){
                    //    assert(I == it->a.first);
                    //    int fi = it->a.second;
                    //    int J = it->b.first;
                    //    int fj = it->b.second;
                    //    energy_before += normsq(u[fi] - recover[J][fj]) * ( C_align*it->weight * 0.5 );
                    //}
                    //
                    //for (int i = 1; i < T; i++){
                    //    energy_before += normsq(u[i] - u[i-1] - (trace[i] - trace[i-1])) * C_smooth * 0.5;
                    //}
                    //cerr << "energy_before:" << energy_before << endl;
                    /////////////////////////DEBUG_END

                    vector<Float> a(T, 0.0);
                    vector<Point> v;
                    for (int i = 0; i < T; i++){
                        v.push_back(Point(D, 0.0));
                    }
                   
                    #pragma omp parallel for
                    for (int i = 0; i < T; i++){
                        a[i] += C_reg * reg_weight[I][i];
                        v[i] = v[i] + trace[i] * C_reg * reg_weight[I][i];
                    }

                    for (vector<Match>::iterator it = matchings[I].begin(); it != matchings[I].end(); it++){
                        assert(I == it->a.first);
                        int fi = it->a.second;
                        int J = it->b.first;
                        int fj = it->b.second;
                        //if (iter == 4){
                        //    cerr << I << " " << fi << " --- " << J << " " << fj << endl;
                        //}
                        a[fi] += C_align * it->weight;
                        
                        v[fi] = v[fi] + recover[J][fj] * (C_align * it->weight);
                    }
                    
                    for (int i = 1; i < T; i++){
                        Point diff = (trace[i] - trace[i-1]);
                        v[i] = v[i] + diff * (C_smooth);
                        v[i-1] = v[i-1] - diff * (C_smooth);
                        a[i] += C_smooth;
                        a[i-1] += C_smooth;
                    }

                    #pragma omp parallel for
                    for (int i = 0; i < T; i++){
                        v[i] = v[i] / C_smooth;
                        a[i] = a[i] / C_smooth;
                    }

                    Float c = a[0];
                    Point z = v[0];
                    vector<Float> cs;
                    vector<Point> zs;
                    for (int i = 1; i < T; i++){
                        // Here c u[i-1] - u[i] == z should be hold 
                        // which is equivalent to 
                        // u[i-1] - (1.0/c) u[i] == (z/c)
                        // (i.e. gradient == 0)
                        // Also, - u[i-1] + a[i] u[i] - u[i+1] == v[i]. 
                        // Adding them together gives us
                        // (-1.0/c + a[i]) u[i] - u[i+1] == v[i] + (z/c)
                        // Therefore, new c is -1.0/c + a[i]
                        // new z is v[i] + z/c
                        cs.push_back(c);
                        zs.push_back(z);
                        z = z / c + v[i];
                        c = -1.0/c + a[i];
                    }
                    u[T-1] = z / c;
                    //#pragma omp parallel for
                    for (int i = T-1; i >= 1; i--){
                        u[i-1] = (zs[i-1] + u[i]) / cs[i-1];
                    }

                    ////////////////////////DEBUG 
                    //Float energy_after = 0.0;
                    //for (int i = 0; i < T; i++){
                    //    energy_after += normsq(u[i] - trace[i]) * C_reg * 0.5;
                    //}

                    //for (vector<Match>::iterator it = matchings[I].begin(); it != matchings[I].end(); it++){
                    //    assert(I == it->a.first);
                    //    int fi = it->a.second;
                    //    int J = it->b.first;
                    //    int fj = it->b.second;
                    //    energy_after += normsq(u[fi] - recover[J][fj]) * ( C_align*it->weight * 0.5 );
                    //}
                    //
                    //for (int i = 1; i < T; i++){
                    //    energy_after += normsq(u[i] - u[i-1] - (trace[i] - trace[i-1])) * C_smooth * 0.5;
                    //}
                    //cerr << "energy_after:" << energy_after << endl;
                    //assert(energy_after <= energy_before + 1e-6);
                    /////////////////////////DEBUG_END
                    //
                    //omp_set_lock(&global);

                    //for (int j = 0; j < adj[I].size(); j++){
                    //    omp_unset_lock(&(locks[adj[I][j]]));
                    //}

                    //omp_unset_lock(&global);
                }
                
                //Compute Energy
                //Vanilla Part
                Float energy_vanilla = 0.0;
                //#pragma omp parallel for
                for (int i = 0; i < N; i++){
                    int j;
                    for (j = 0; j < traces[i].size(); j++){
                        //#pragma omp atomic
                        energy_vanilla += normsq(traces[i][j] - recover[i][j]) * C_reg * reg_weight[i][j] * 0.5;
                    }
                }
                //Smooth Part
                Float energy_smooth = 0.0;
                //#pragma omp parallel for
                for (int i = 0; i < N; i++){
                    int t;
                    for (t = 1; t < traces[i].size(); t++){
                        Point diff = (recover[i][t] - recover[i][t-1]) - (traces[i][t] - traces[i][t-1]);
                        //#pragma omp atomic
                        energy_smooth += normsq(diff) * C_smooth * 0.5;
                    }
                }
                //Align Part
                Float energy_align = 0.0;
                //#pragma omp parallel for
                int num_edge = 0;
                for (int i = 0; i < N; i++){
                    vector<Match>::iterator it;
                    num_edge += matchings[i].size();
                    for (it = matchings[i].begin(); it != matchings[i].end(); it++){
                        assert(i == it->a.first);
                        int fi = it->a.second;
                        int j = it->b.first;
                        int fj = it->b.second;
                        //#pragma omp atomic
                        //cerr << i << " " << fi << " --- " << j << " " << fj << endl;
                        energy_align += normsq(recover[i][fi] - recover[j][fj]) * C_align * it->weight * 0.5 * 0.5;
                    }
                }
                cout << "iter=" << iter;
                cout << ", vanilla=" << energy_vanilla;
                cout << ", smooth=" << energy_smooth;
                cout << ", align=" << energy_align;
                cout << ", #edges=" << num_edge;
                cout << ", energy=" << (energy_vanilla + energy_smooth + energy_align);
                cout << endl;

                iter++;
                Float delta = (last_energy - (energy_vanilla + energy_smooth + energy_align))/last_energy;
                if ((iter % 1 == 0) || (delta < params->tol)){
                    if (D != 2){
                        string list_file_name;

                        list_file_name = "motion/smooth/smoothed_" + to_string(iter)+".txt";
                        ofstream smooth_list(list_file_name);

                        for (int i = 0; i < traces.size(); i++){
                            motions[i]->traces = recover[i];
                            string filename_i;
                            filename_i = replace(motions[i]->filename, "data", "smooth/"+to_string(iter));
                            motions[i]->write_data(filename_i);
                            smooth_list << filename_i << endl;
                            //replace(motions[i].filename, "data", "smooth/"+to_string(iter)) << endl;
                        }
                        smooth_list.close();
                    } else {
                        string file_name;
                        if (delta < params->tol){
                            file_name = "GPS/" + output_folder + "/GPS.txt";
                        } else {
                            file_name = "GPS/" + output_folder + "/GPS_"+to_string(iter)+".txt";
                        }
                        string directory = "/home/xiangru/Projects/Qixing/TrajSmoothing/" + dir(file_name);
                        cerr << "\tdumping to " << directory << endl;
                        string command = "mkdir -p " + directory;
                        int e = system(command.c_str());
                        cerr << " directory created" << endl;
                        //ofstream fout(file_name, fstream::out);
                        //fout.precision(10);
                        for (int i = 0; i < traces.size(); i++){
                            motions[i]->traces = recover[i];
                            motions[i]->write_data(file_name, true);
                            //for (int j = 0; j < motions[i]->traces.size(); j++){
                            //    if (j != 0){
                            //        fout << " ";
                            //    }
                            //    Point r = motions[i]->traces[j];
                            //    if (motions[i]->centered){
                            //        r = (r * motions[i]->std) + motions[i]->mean_pose;
                            //    }
                            //    for (int d = 0; d < D; d++){
                            //        if (d != 0){
                            //            fout << " ";
                            //        }
                            //        fout << r[d];
                            //    }
                            //}
                            //fout << endl;
                        }
                        //fout.close();
                    }
                }
                if (delta < params->tol){
                    break;
                }
                last_energy = (energy_vanilla + energy_smooth + energy_align);

            }// end of while
            
            //for (int i = 0; i < N; i++){
            //    omp_destroy_lock(&(locks[i]));
            //}
            //omp_destroy_lock(&global);
            for (int i = 0; i < traces.size(); i++){
                motions[i]->traces = recover[i];
            }
        }

        int N;
        unordered_map<int, vector<pair<int,int>>> match;
        vector<vector<Point>> recover;
        Point mean, std;
        vector<vector<Point>> traces;
        bool do_shift = false;
        Float learning_rate;
        int iter;
};

#endif
