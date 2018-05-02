#ifndef TRACE_H
#define TRACE_H

#include <fstream>
#include <string>
#include "util.h"
#include <iomanip>
#include <cassert>
#include <ANN/ANN.h>

class GPSTrace;
class HumanMotion;

static Float average_num_point = 0.0;
static int down = 0;

class Trace{
    public:
    Trace(){};
    Trace(string _filename): filename(_filename){};
    
    inline vector<int> search(int K, Point point, ANNkd_tree* tree){
        int* nn_idx = new int[K];
        double* dists = new double[K];
        double* query = to_double(point);
        tree->annkSearch(query, K, nn_idx, dists, eps);
        vector<int> ans;
        for (int i = 0; i < K; i++){
            ans.push_back(nn_idx[i]);
        }

        delete[] nn_idx;
        delete[] dists;
        delete[] query;

        return ans;
    }
    
    inline vector<int> search(int K, double* query, ANNkd_tree* tree){
        int* nn_idx = new int[K];
        double* dists = new double[K];
        tree->annkSearch(query, K, nn_idx, dists, eps);
        vector<int> ans;
        for (int i = 0; i < K; i++){
            ans.push_back(nn_idx[i]);
        }

        delete[] nn_idx;
        delete[] dists;

        return ans;
    }
    
    inline vector<pair<int, Point>> shortest_path(ANNkd_tree* tree, int ord, Float dist, string output_folder, vector<vector<pair<Float, int>>>& adj){

        // all points from the tree
        double** points = tree->thePoints();

        // every point is assumed to be inside this ANN kd tree, so first
        // retrieve their indices
        vector<int> indices;
        int* nn_idx = new int[1];
        double* dists = new double[1];
        for (int i = 0; i < traces.size(); i++){
            double* query = to_double(traces[i]);
            tree->annkSearch(query, 1, nn_idx, dists, 1e-10);
            indices.push_back(nn_idx[0]);
            assert(norm(traces[i]-points[nn_idx[0]]) <= 1e-10);
            delete[] query;
        }
        delete[] nn_idx;
        delete[] dists;
        
        vector<Point> path;
        path.push_back(traces[0]);
        vector<bool> inside;
        inside.push_back(true);

        for (int i = 1; i < traces.size(); i++){
            vector<int> shortest = dijkstra(indices[i-1], indices[i], tree, adj);
            for (int t = 1; t < shortest.size() - 1; t++){
                Point pt = to_vector(D, points[shortest[t]]);
                path.push_back(pt);
                inside.push_back(false);
            }
            path.push_back(traces[i]);
            inside.push_back(true);
        }
        string command = "mkdir -p " + dir(output_folder);
        int e = system(command.c_str());
        ofstream fout(output_folder, fstream::out);
        fout.precision(15);
        for (int i = 0; i < path.size(); i++){
            if (i != 0){
                fout << " ";
            }
            Point r = path[i];
            for (int d = 0; d < D; d++){
                if (d != 0){
                    fout << " ";
                }
                fout << r[d];
            }
        }
        fout << endl;
        fout.close();

        vector<Point> temp = traces;
        traces = path;
        vector<pair<int, Point>> ret = equal_dist_samples(ord, dist, inside);
        traces = temp;
        
        return ret;
    }

    inline vector<pair<int, Point>> interpolate(int num_inter, ANNkd_tree* tree, int ord, Float dist, string output_folder){
        vector<Point> new_traces;
        vector<bool> inside;
        new_traces.push_back(traces[0]);
        inside.push_back(true);
        int K = 100;
        double** points = tree->thePoints();
        int num_point = tree->nPoints();
        if (K > num_point){
            K = num_point;
        }
        Float eps=1e-6;
        double** f = new double*[num_inter];
        int** rec_prev = new int*[num_inter];
        for (int i = 0; i < num_inter; i++){
            f[i] = new double[K];
            rec_prev[i] = new int[K];
        }
        assert(D == 2);
        for (int I = 1; I < traces.size(); I++){
            // get interpolated points
            vector<vector<int>> c;
            vector<Point> candidates;
            for (int t = 0; t < num_inter; t++){
                double ratio = t / (num_inter+1.0);
                Point pt = traces[I-1] * ratio + traces[I] * (1.0-ratio);
                candidates.push_back(pt);
            }

            for (int T = 0; T < 10; T++){
                c.clear();
                for (int t = 0; t < num_inter; t++){
                    vector<int> ans = search(K, candidates[t], tree);
                    c.push_back(ans);
                }
                for (int i = 0; i < K; i++){
                    f[0][i] = normsq(traces[I-1] - points[c[0][i]]);
                }
                for (int j = 1; j < num_inter; j++){
                    for (int i = 0; i < K; i++){
                        f[j][i] = 1e100;
                        for (int prev = 0; prev < K; prev++){
                            double cost = f[j-1][prev] + diffsq(points[c[j][i]], points[c[j-1][prev]], D);
                            if (j == num_inter-1){
                                cost += normsq(traces[I] - points[c[j][i]]);
                            }
                            if (cost < f[j][i]){
                                f[j][i] = cost;
                                rec_prev[j][i] = prev;
                            }
                        }
                    }
                }
                vector<int> ans;
                for (int i = 0; i < num_inter; i++){
                    ans.push_back(0);
                }
                int best_prev = 0;
                for (int i = 1; i < K; i++){
                    if (f[num_inter-1][i] < f[num_inter-1][best_prev]){
                        best_prev = i;
                    }
                }
                ans[num_inter-1] = best_prev;
                for (int i = num_inter-2; i >= 0; i--){
                    best_prev = rec_prev[i+1][best_prev];
                    ans[i] = best_prev;
                }
                candidates.clear();
                for (int i = 0; i < num_inter; i++){
                    candidates.push_back(to_vector(D, points[c[i][ans[i]]]));
                }
            }
            for (int i = 0; i < num_inter; i++){
                new_traces.push_back(candidates[i]);
                inside.push_back(false);
            }
            new_traces.push_back(traces[I]);
            inside.push_back(true);
        }
        ofstream fout(output_folder, fstream::out);
        for (int i = 0; i < new_traces.size(); i++){
            if (i != 0){
                fout << " ";
                //cerr << " " ;
            }
            Point r = new_traces[i];
            for (int d = 0; d < D; d++){
                if (d != 0){
                    fout << " ";
                    //cerr << " ";
                }
                fout << r[d];
                //cerr << r[d];
            }
        }
        //cerr << endl;
        fout << endl;
        fout.close();
        vector<Point> temp = traces;
        traces = new_traces;
        vector<pair<int, Point>> ret = equal_dist_samples(ord, dist, inside);
        traces = temp;
        
        return ret;
    }


    inline vector<pair<int, Point>> equal_dist_samples(int ord, Float dist, vector<bool> inside){
        vector<pair<int, Point>> sample_vecs;
        Float average_num_point_local = 0.0;
        Float num_point_covered = 0.0;
        int down_local = 0;
        //for (int i = 0; i < traces.size(); i++){
        //    if (i != 0){
        //        cerr << " ";
        //    }
        //    Point r = traces[i];
        //    for (int d = 0; d < D; d++){
        //        if (d != 0){
        //            cerr << " ";
        //        }
        //        cerr << r[d];
        //    }
        //}
        //cerr << "traces.size()=" << traces.size() << ", skip=" << skip << ", ord=" << ord << ", dist=" << dist << endl;
        Float total_dist = 0.0;

        for (vector<Point>::iterator it = traces.begin()+1; it < traces.end(); it++){
            Float d = norm((*(it-1))- (*it));
            total_dist += d;
        }
        Float left_dist = 0.0;
        int trace_index = 0;
        int num_original_point = 0;
        for (vector<Point>::iterator it = traces.begin(); it < traces.end(); it++, trace_index++){
            if (!inside[trace_index]){
                continue;
            }
            num_original_point++;
            if (it != traces.begin()){
                Float dd = norm((*(it-1)) - (*it));
                left_dist += dd;
                total_dist -= dd;
            }
            //left
            vector<Point> left_samples;
            vector<Point>::iterator left_it = it;
            Float rest = 0.0;
            //cerr << "left begin" << endl;
            int temp_trace_index = trace_index;
            while ((left_samples.size() < ord) && (left_it != traces.begin())){
                vector<Point>::iterator next_it = left_it - 1;
                Float d = norm((*next_it) - (*left_it));
                left_dist += d;
                int i = 0;
                while (i + 1 <= (rest + d) / dist){
                    i++;
                    Float move = i * dist - rest;
                    Float ratio = move / d;
                    Point p = (*left_it) * ratio + (*next_it) * (1.0-ratio);
                    //p = p / norm(p);
                    left_samples.push_back(p);
                    if (left_samples.size() == ord){
                        break;
                    }
                }
                if (left_samples.size() == ord){
                    break;
                }
                rest = rest + d - dist * i;
                left_it = next_it;
                temp_trace_index--;
                if (inside[temp_trace_index]){
                    average_num_point_local++;
                }
            }
            
            //cerr << " left loc=" << (it-traces.begin()) << endl;
            //average_num_point_local += (it - left_it)/(skip * 1.0);
            down_local++;
            //cerr << "right begin" << endl;
            //right 
            vector<Point> right_samples;
            vector<Point>::iterator right_it = it;
            rest = 0.0;
            temp_trace_index = trace_index;
            while ((right_samples.size() < ord) && (right_it != traces.end()-1)){
                //cerr << right_it->at(0) << " " << right_it->at(1) << endl;
                vector<Point>::iterator next_it = right_it + 1;
                Float d = norm((*next_it) - (*right_it));
                int i = 0;
                while (i+1 <= (rest + d) / dist){
                    i++;
                    Float move = i * dist - rest;
                    Float ratio = move / d;
                    Point p = (*right_it) * ratio + (*next_it) * (1.0-ratio);
                    //p = p / norm(p);
                    right_samples.push_back(p);
                    if (right_samples.size() == ord){
                        break;
                    }
                }
                if (right_samples.size() == ord){
                    break;
                }
                rest = rest + d - dist * i;
                right_it = next_it;
                temp_trace_index++;
                if (inside[temp_trace_index]){
                    average_num_point_local++;
                }
            }
            if ((left_samples.size() < ord) && (right_samples.size() < ord)){
                continue;
            }
            num_point_covered++;
            //average_num_point_local += (right_it - it) / (skip * 1.0);
            down_local++;
            
            Point sample_vec;
            for (int i = 0; i < ord; i++){
                Point l_i;
                if (left_samples.size() < ord){
                    l_i = zero_point(D) + 1e5;
                } else {
                    l_i = left_samples[i];
                }
                assert(l_i.size() == D);
                for (int d = 0; d < D; d++){
                    sample_vec.push_back(l_i[d]);
                }
            }
            assert((*it).size() == D);
            for (int d = 0; d < D; d++){
                sample_vec.push_back((*it)[d]);
            }
            for (int i = 0; i < ord; i++){
                Point r_i;
                if (right_samples.size() < ord){
                    r_i = zero_point(D) + 1e5;
                } else {
                    r_i = right_samples[i];
                }
                assert(r_i.size() == D);
                for (int d = 0; d < D; d++){
                    sample_vec.push_back(r_i[d]);
                }
            }
            //sample_vec = sample_vec/(2*ord);
            assert(sample_vec.size() == D*(2*ord+1));
            sample_vecs.push_back(make_pair(num_original_point-1, sample_vec));
        }
        //assert(sample_vecs.size() != 0);
        num_point_covered /= num_original_point;
        average_num_point += average_num_point_local;
        down += down_local;
        if (rand() % 10000 == 0){
            //cerr << "average_num_point covered=" << average_num_point / down << ", down=" << down << endl;
        }
        //cerr << "percentage of point covered=" << num_point_covered << endl;
        return sample_vecs;
    }
    
    void center(Point& _mean_pose, Point& _std){
        this->mean_pose = _mean_pose;
        this->std = _std;
        for (int i = 0; i < traces.size(); i++){
            traces[i] = (traces[i] - mean_pose) / std;
        }
        centered = true;
    }

    virtual void read_data(string filename) = 0;

    virtual void write_data(string filename, bool append = false) = 0;

    void dump(){
        cerr << "D=" << D << endl;
        cerr << "filename=" << filename << endl;
        cerr << "traces.size=" << traces.size() << endl;
        for (int i = 0; i < traces.size(); i++){
            cerr << "i=" << i << ", traces[i].size=" << traces[i].size() << endl;
            assert(traces[i].size() == 2);
            for (int j = 0; j < traces[i].size(); j++){
                cerr << traces[i][j] << " ";
            }
            cerr << " ";
        }
        cerr << "centered=" << centered << endl;
    }
    int D = -1;
    string filename;
    vector<Point> traces;
    bool centered = false;
    Point mean_pose;
    Point std;
    double eps = 1e-6;
};

class HumanMotion : public Trace{
    public:
    HumanMotion(string filename){
        this->filename = filename;
        D = 51;
        read_data(filename);
    }

    //vector<Point> directions(int ord){
    //    vector<Point> dirs;
    //    for (vector<Point>::iterator it = traces.begin() + ord; it < traces.end() - ord; it++){
    //        Point dir;
    //        for (int offset = -ord; offset < ord; offset++){
    //            Point diff = *(it+offset+1) - *(it+offset);
    //            if (norm(diff) >= 1e-5){
    //                diff = diff / norm(diff);
    //            }
    //            for (int d = 0; d < D; d++){
    //                dir.push_back(diff[d]);
    //            }
    //        }
    //        dir = dir/ norm(dir);
    //        for (int d = 0; d < D; d++){
    //            dir.push_back((*it)[d]);
    //        }
    //        assert(dir.size() == D*(2*ord+1));
    //        dirs.push_back(dir);
    //    }
    //    assert(dirs.size() == traces.size() - 2*ord);
    //    return dirs;
    //}
    
    void read_data(string filename){
        ifstream fin(filename, fstream::in);
        char* line = new char[LINE_LEN];
        int used_joints[] = { 0, 2, 3, 4, 7, 8, 9, 12, 13, 15, 16, 18,
            19, 20, 25, 26, 27};
        bool* useful = new bool[31];
        for (int i = 0; i < 31; i++){
            useful[i] = false;
        }
        for (int i = 0; i < 17; i++){
            useful[used_joints[i]] = true;
        }
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
            Point trace;
            Point useless;
            Point all;
            for (int count = 0; count < tokens.size(); count++){
                if (tokens[count].size() < 2){
                    continue;
                }
                Float fval = stod(tokens[count]);
                if (count <= 2){
                    position.push_back(fval);
                } else {
                    assert(count / 3 - 1 < 31);
                    if (useful[count / 3 - 1]){
                        trace.push_back(fval);
                    } else {
                        useless.push_back(fval);
                    }
                }
            }
            positions.push_back(position);
            traces.push_back(trace);
            unused.push_back(useless);
        }
        D = traces[0].size();
        delete[] line;
        delete[] useful;
        fin.close();
    }

    inline void write_data(string filename, bool append=false){
        int used_joints[] = { 0, 2, 3, 4, 7, 8, 9, 12, 13, 15, 16, 18,
            19, 20, 25, 26, 27};
        bool* useful = new bool[31];
        for (int i = 0; i < 31; i++){
            useful[i] = false;
        }
        for (int i = 0; i < 17; i++){
            useful[used_joints[i]] = true;
        }

        string command = "mkdir -p " + dir(filename);
        int e = system(command.c_str());
        ofstream fout(filename);
        for (vector<string>::iterator it = skeleton.begin(); it != skeleton.end(); it++){
            fout << *it << endl;
        }
        int num_frame = traces.size();
        fout << "Frames: " << traces.size() << endl;
        fout << "Frame Time: .0083333" << endl;
        for (int i = 0; i < num_frame; i++){
            Point& p = positions[i];
            fout << fixed << setprecision(6) << p[0] << " " << p[1] << " " << p[2];
            Point r = traces[i];
            if (centered){
                r = (r * std) + mean_pose;
            }
            Point& u = unused[i];
            int count_r = 0, count_u = 0;
            for (int j = 0; j < 31; j++){
                Float f1, f2, f3;
                if (useful[j]){
                    f1 = r[count_r*3+0];
                    f2 = r[count_r*3+1];
                    f3 = r[count_r*3+2];
                    count_r++;
                } else {
                    f1 = u[count_u*3+0];
                    f2 = u[count_u*3+1];
                    f3 = u[count_u*3+2];
                    count_u++;
                }
                fout << fixed << setprecision(6) << " " << f1 << " " << f2 << " " << f3;
            }
            fout << endl;
        }
        fout.close();
    }

    int size(){
        return traces.size();
    }

    string filename;
    vector<string> skeleton;
    vector<Point> positions;
    vector<Point> traces; // rotations
    vector<Point> unused;
};

class GPSTrace: public Trace{
    public:
    GPSTrace(string _filename) : Trace(_filename) {
        filename=_filename;
        D = 2;
        this->read_data(filename);
    }

    void read_data(string filename){
        ifstream fin(filename, fstream::in);
        char* line = new char[LINE_LEN];

        fin.getline(line, LINE_LEN);
        string line_str(line);
        if (line_str.size() < 2 && fin.eof()){
            assert(false);
        }
        vector<string> tokens = split(line_str, " ");
        Point trace;
        traces.clear();
        int count = 0;
        assert(tokens.size() % D == 0);
        for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++){
            Float fval = stod(*it);
            trace.push_back(fval);
            count++;
            if (count % D == 0){
                Point temp = trace;
                traces.push_back(temp);
                trace.clear();
            }
        }
    }

    inline void write_data(string filename, bool append = false){
        string command = "mkdir -p " + dir(filename);
        int e = system(command.c_str());
        
        ofstream fout;
        if (!append){
            fout.open(filename, fstream::out);
        } else {
            fout.open(filename, ios_base::app);
        }
        fout.precision(10);
        for (int i = 0; i < traces.size(); i++){
            if (i != 0){
                fout << " ";
            }
            Point r = traces[i];
            if (centered){
                r = (r * std) + mean_pose;
            }
            for (int d = 0; d < D; d++){
                if (d != 0){
                    fout << " ";
                }
                fout << r[d];
            }
        }
        fout << endl;
        fout.close();
    }
};

#endif
