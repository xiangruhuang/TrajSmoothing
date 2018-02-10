#ifndef CORRS_H
#define CORRS_H

#include "util.h"
#include "Graph.h"

using namespace std;

inline void euclidean_correspondences(vector<Motion>& motions, unordered_map<int, vector<pair<pair<int,int>, Float>>>& match){
    match.clear();
    Graph graph(motions);
    int K = 20;
    for (int i = 0; i < n; i++){
        graph.search();
    }
}

void compute_correspondences(vector<vector<Point>>& traces, string method){
   if (method == "euclidean"){
        
   }
}

#endif
