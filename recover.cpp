#include "util.h"
#include "algo.h"
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include "GD.h" 
#include "AM.h"
#include "Trace.h"
#include "Graph.h"
#include "Params.h"
#include<cassert>
#include <unordered_map>
#include<vector>

using namespace std;

class GraphStat{
    public:
    vector<vector<pair<Float, int>>> adj;
    vector<Point> points;
    Float** _points;
    ANNkd_tree* tree;
    int N, D;

    GraphStat(vector<Point> given_points, Params* params){
        points = given_points;
        N = points.size();
        assert(N > 0);
        D = points[0].size();
        tree = build_ann_tree(points);
        _points = tree->thePoints();
        int K = params->K;
        cerr << "building adjacency map with K=" << K << ", radius=" << params->radius * params->ord << "\t...";
        adj = adjacency(min(K, N), tree, params->radius*params->ord);
        
        cerr << "done." << endl;
    }
    ~GraphStat(){
        for (int i = 0; i < points.size(); i++){
            delete[] _points[i];
        }
        delete[] _points;
    }

    vector<int> dijkstra_idx(Point src, Point tgt){
        // project to get index for src and tgt
        vector<Point> st;
        st.push_back(src);
        st.push_back(tgt);
        vector<int> idxs = project(tree, st);

        int src_idx = project(tree, src);
        int tgt_idx = project(tree, tgt);
        return dijkstra(src_idx, tgt_idx, tree, adj);
    }
    
    vector<Point> dijkstra_points(Point src, Point tgt){
        // project to get index for src and tgt
        vector<Point> st;
        st.push_back(src);
        st.push_back(tgt);
        vector<int> idxs = project(tree, st);

        int src_idx = project(tree, src);
        int tgt_idx = project(tree, tgt);
        vector<int> path = dijkstra(src_idx, tgt_idx, tree, adj);
        vector<Point> ans;
        for (int i = 0; i < path.size(); i++){
            ans.push_back(points[path[i]]);
        }
        return ans;
    }
};

GraphStat eps_cover_graph(vector<Point> points, Params* params){
    Float radius = params->radius;
    cerr << "computing eps-cover, eps=" << radius << "\t...";
    vector<int> assign = eps_cover(points, radius);

    // clusters: store locations of clusters
    // translate assign into a function that maps
    // [points.size()] to [clusters.size()]
    int* dense_id = new int[points.size()];
    vector<Point> clusters;
    for (int i = 0; i < points.size(); i++){
        dense_id[i] = -1;
    }
    for (int i = 0; i < points.size(); i++){
        int v = assign[i];
        if (dense_id[v] >= 0){
            // already in clusters
            assign[i] = dense_id[v];
            continue;
        }
        dense_id[v] = clusters.size();
        assign[i] = dense_id[v];
        clusters.push_back(points[v]);
    }
    int num_cluster = clusters.size();
    GraphStat gs(clusters, params);
    
    cerr << ", #cluster=" << clusters.size();

    delete[] dense_id;
    cerr << "done" << endl;
    return gs;
}

inline vector<Point> remove_fluc(vector<Point> trace, Params* params){
    vector<Point> ans;
    for (int t = 0; t < trace.size(); t++){
        //if (ans.size() >= 1 && sqrt(distance(ans[ans.size()-1], trace[t])) <= params->interval){
        //    continue;
        //}
        if (ans.size() < 2){
            ans.push_back(trace[t]);
            continue;
        }
        if (distance(ans[ans.size()-2], trace[t]) <= 1e-15){
            ans.pop_back();
            continue;
        }
        Point delta1 = ans[ans.size()-1] - ans[ans.size()-2];
        delta1 = delta1 / norm(delta1);
        Point delta2 = trace[t] - ans[ans.size()-1];
        delta2 = delta2 / norm(delta2);
        Float dp = dot(delta1, delta2);
        if (dp <= 0.0){
            ans.pop_back();
        }
        ans.push_back(trace[t]);
    }
    return ans;
}

inline vector<vector<pair<Float, int>>> trim_adj(vector<vector<pair<Float, int>>> adj, unordered_map<int, Float>& map, Params* params){
    int N = adj.size();
    Float r = params->ratio;
    Float t = params->threshold;
    int* indeg = new int[N];
    int* outdeg = new int[N];
    for (int i = 0; i < N; i++){
        indeg[i] = 0;
        outdeg[i] = 0;
    }
    for (int i = 0; i < N; i++){
        vector<pair<Float, int>> next;
        for (int t = 0; t < adj[i].size(); t++){
            int j = adj[i][t].second;
            if (i == j){
                continue;
            }
            Float wij = 0.0;
            unordered_map<int, Float>::iterator it = map.find(i*N+j);
            if (it != map.end()){
                wij = it->second;
            }
            if (wij < params->threshold - 1e-4){
                continue;
            }

            Float wji = 0.0;
            it = map.find(j*N+i);
            if (it != map.end()){
                wji = it->second;
            }
            if (wij < wji * r){
                continue;
            }
            next.push_back(adj[i][t]);
            outdeg[i]++;
            indeg[j]++;
        }
        adj[i] = next;
    }

    int max_iter = params->max_iter;
    for (int iter = 0; iter < max_iter; iter++){
        for (int i = 0; i < N; i++){
            vector<pair<Float, int>> next;
            for (int t = 0; t < adj[i].size(); t++){
                int j = adj[i][t].second;
                if (indeg[i] == 0 || outdeg[j] == 0){
                    indeg[j]--;
                    outdeg[i]--;
                    continue;
                }
                next.push_back(adj[i][t]);
            }
            adj[i] = next;
        }
    }
    return adj;
}

int main(int argc, char** argv){
    srand(time(NULL));
    // Read Traces
    Params params(argc, argv);
    params.dump();
    int T = 1;
    //for (int iter = 0; iter < T; iter++){
    //    string command = "mkdir -p GPS/"+params.output_folder+"/map"+to_string(iter);
    //    cerr << command.c_str() << endl;
    //    int e = system(command.c_str());
    //    while (e != 0){
    //        e = system(command.c_str());
    //    }
    //}
    vector<string> filelist = read_filelist(params.filelist_path);
    vector<Trace*> traces;
    for (int i = 0; i < filelist.size(); i++){
        traces.push_back(NULL);
    }

    #pragma omp parallel for
    for (int i = 0; i < filelist.size(); i++){
        cerr << "reading " << filelist[i];
        if (params.D == 2){
            GPSTrace* trace = new GPSTrace(filelist[i]);
            cerr << ", #frames=" << trace->traces.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
            traces[i] = trace;
        } else {
            HumanMotion* trace = new HumanMotion(filelist[i]);
            cerr << ", #frames=" << trace->traces.size() << " ( " << i << " / " << filelist.size() << " )" << endl;
            traces[i] = trace;
        }
    }

    // retrieve points
    vector<Point> points;
    for (int i = 0; i < traces.size(); i++){
        for (int j = 0; j < traces[i]->traces.size(); j++){
            points.push_back(traces[i]->traces[j]);
        }
    }

    // build original graph
    GraphStat graph(points, &params);
    points = graph.points;

    // build eps-cover graph
    GraphStat cover_graph = eps_cover_graph(points, &params);
    vector<Point>& clusters = cover_graph.points;
    vector<vector<pair<Float, int>>>& adj = cover_graph.adj;
    int num_cluster = cover_graph.N;
    cerr << "#cluster=" << cover_graph.N << " out of " << points.size() << " points." << endl;

    // retrive shortest path, optionally: smooth it
    vector<vector<Point>> smoothed;
    vector<vector<Point>> trace_points;
    for (int i = 0; i < traces.size(); i++){
        vector<Point> si;
        smoothed.push_back(si);
        trace_points.push_back(si);
    }
    #pragma omp parallel for
    for (int i = 0; i < traces.size(); i++){
        if (i % 100 == 0){
            cerr << ".";
        }
        vector<Point> trace_i = traces[i]->traces;
        vector<Point> path_i;
        for (int j = 0; j < trace_i.size() - 1; j++){
            vector<Point> pj = cover_graph.dijkstra_points(trace_i[j], trace_i[j+1]);
            int st = 0;
            if (j != 0){
                st = 1;
            }
            for (int k = st; k < pj.size(); k++){
                path_i.push_back(pj[k]);
            }
        }
        trace_points[i] = unique(path_i);
        
        vector<Point> smooth_i = remove_fluc(trace_points[i], &params);
        //smooth_i = smooth(params.lambda, smooth_i);
        smoothed[i] = smooth_i;
    }
    cerr << endl;
    /////////////////////////////////DEBUG
    int num_fluc = 0;
    int total_length = 0;
    for (int i = 0; i < traces.size(); i++){
        vector<Point> path_i = smoothed[i];
        for (int j = 2; j < path_i.size(); j++){
            if (distance(path_i[j], path_i[j-2]) <= 1e-15){
                num_fluc++;
            }
        }
        total_length += (path_i.size() - 1);
    }
    cerr << "#fluc=" << num_fluc << ", total length=" << total_length << endl;
    ///////////////////////////////////DEBUGEND

    // compute weights
    unordered_map<int, Float> map;
    map.clear();
    vector<pair<int, int>>* paths = new vector<pair<int, int>>[traces.size()];

    #pragma omp parallel for
    for (int i = 0; i < smoothed.size(); i++){
        vector<pair<int, int>>& path_i = paths[i];
        vector<int> project_i = project(cover_graph.tree, smoothed[i]);
        for (int t = 1; t < smoothed[i].size(); t++){
            path_i.push_back(make_pair(project_i[t-1], project_i[t]));
        }
        if (i % 1000 == 0){
            cerr << i << "/" << smoothed.size() << endl;
        }
    }

    for (int i = 0; i < traces.size(); i++){
        vector<pair<int, int>> path_i = paths[i];
        for (int t = 0; t < path_i.size(); t++){
            pair<int, int> p = path_i[t];
            int a = p.first;
            int b = p.second;
            if (a == b){
                continue;
            }
            //if (a > b){
            //    int temp = a; a = b; b = temp;
            //}
            unordered_map<int, Float>::iterator it = map.find(a*num_cluster + b);

            Float weight = 0.0;
            if (it == map.end()){
                map.insert(make_pair(a*num_cluster+b, 1.0));
            } else {
                weight = it->second;
                it->second = weight + 1.0;
            }

            Float next_weight = map.find(a*num_cluster + b)->second;
            assert(abs(next_weight - (weight+1.0))<1e-3);
        }
    }

    Float threshold = params.threshold;
    vector<pair<int, int>> edges;
    int* indeg = new int[num_cluster];
    int* outdeg = new int[num_cluster];
    for (int i = 0; i < num_cluster; i++){
        indeg[i] = 0;
        outdeg[i] = 0;
    }

    adj = trim_adj(adj, map, &params);

    int num_edge = 0;
    for (int i = 0; i < clusters.size(); i++){
        num_edge += adj[i].size();
        for (int j = 0; j < adj[i].size(); j++){
            int a = i;
            int b = adj[i][j].second;
            edges.push_back(make_pair(a,b)); 
            outdeg[a]++;
            indeg[b]++;
        }
    }

    // get reverse adjacency
    vector<vector<pair<Float, int>>> reverse_adj;
    
    for (int i = 0; i < clusters.size(); i++){
        vector<pair<Float, int>> ri;
        reverse_adj.push_back(ri);
    }
    
    for (int i = 0; i < clusters.size(); i++){
        vector<pair<Float, int>>& adj_i = adj[i];
        for (int t = 0; t < adj_i.size(); t++){
            int j = adj_i[t].second;
            reverse_adj[j].push_back(make_pair(adj_i[t].first, j));
        }
    }

    bool* visited = new bool[clusters.size()];
    // find long chains
    for (int i = 0; i < clusters.size(); i++){
        
        if ((!visited[i]) && indeg[i] == 1 && outdeg[i] == 1){
            //cerr << "Found a potential new chain..." << endl;
            vector<int> chain;
            int v = i;
            while (outdeg[v] == 1 && indeg[v] == 1 && !visited[v]){
                chain.push_back(v);
                visited[v] = true;
                v = reverse_adj[v][0].second;
            }
            reverse(chain.begin(), chain.end());
            chain.pop_back();
            v = i;
            visited[v] = false;
            while (outdeg[v] == 1 && indeg[v] == 1 && !visited[v]){
                chain.push_back(v);
                visited[v] = true;
                v = adj[v][0].second;
            }
            if (chain.size() <= 2){
                continue;
            }
            cerr << "Got a new chain with size: " << chain.size() << endl;
            vector<Point> b;
            for (int t = 0; t < chain.size(); t++){
                b.push_back(clusters[chain[t]]);
            }
            b = smooth_strong_head(params.lambda, b, chain.size());
            for (int t = 0; t < chain.size(); t++){
                clusters[chain[t]] = b[t];
            }
        }
    }

    delete[] visited;

    
    int* new_idx = new int[clusters.size()];
    vector<Point> visible_points;
    for (int i = 0; i < clusters.size(); i++){
        if (indeg[i] > 0 || outdeg[i] > 0){
            new_idx[i] = visible_points.size();
            visible_points.push_back(clusters[i]);
        } else {
            new_idx[i] = -1;
        }
    }

    ofstream fout(params.output_folder, fstream::out);
    fout.precision(15);
    fout << "Points:" << endl;
    for (int i = 0; i < visible_points.size(); i++){
        for (int d = 0; d < visible_points[i].size(); d++){
            if (d != 0){
                fout << " ";
            }
            fout << visible_points[i][d];
        }
        fout << endl;
    }
    fout << "Edges:" << endl;
    for (int i = 0; i < edges.size(); i++){
        int a = new_idx[edges[i].first];
        int b = new_idx[edges[i].second];
        assert(a < visible_points.size());
        assert(b < visible_points.size());
        //assert(outdeg[edges[i].first] > 0);
        //assert(indeg[edges[i].second] > 0);
        //cerr << edges[i].first << " " << num_cluster << " " << edges[i].second << endl;
        Float w = map.find(edges[i].first*num_cluster+edges[i].second)->second;
        fout << a << " " << b << " " << w << endl;
    }
    fout.close();

    delete[] indeg;
    delete[] outdeg;
    delete[] new_idx;
    delete[] paths;
}
