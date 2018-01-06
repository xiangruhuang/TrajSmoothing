#ifndef ALGO
#define ALGO
#include <iostream>
#include <vector>
#include "util.h"

using namespace std;

//
// Longest Common Subsequences
//
vector<pair<int, int>> LCS(vector<Point> a, vector<Point> b, Float threshold, int period){
    vector<pair<int, int>> matchings;
    int n = a.size(), m = b.size();
    int** f = new int*[n];
    int** g = new int*[n];
    Float aug=2.0;
    for (int i = 0; i < n; i++){
        f[i] = new int[m];
        g[i] = new int[m];
        //
        // f[i][j] = max{path0:f[i-1][j-1]+1_{a[i] == b[j]}, path1:f[i][j-1], path2:f[i-1][j]}
        //
        for (int j = 0; j < m; j++){
            Float dist = distance(a[i], b[j]);
            int match = 0;
            Float cur_threshold = threshold;
            if (j % (period + 1) == 0 || (i % (period+1) == 0)){
                cur_threshold *= aug;
            }
            if (dist < cur_threshold){
                match = 1;
            }
            if (i == 0 || j == 0){
                f[i][j] = match;
                g[i][j] = 0;
            } else {
                f[i][j] = f[i-1][j-1] + match;
                g[i][j] = 0;
            }
            if ((i > 0) && (f[i][j] < f[i-1][j])){
                f[i][j] = f[i-1][j];
                g[i][j] = 1;
            }
            if ((j > 0) && (f[i][j] < f[i][j-1])){
                f[i][j] = f[i][j-1];
                g[i][j] = 2;
            }
        }
    }
    int i = n-1, j = m-1;
    while (i >= 0 && j >= 0){
        if (g[i][j] == 0){
            Float cur_threshold = threshold;
            if (j % (period + 1) == 0 || (i % (period+1) == 0)){
                cur_threshold *= aug;
            }
            if (distance(a[i], b[j]) < cur_threshold){
                matchings.push_back(make_pair(i, j));
            }
            i--; j--;
            continue;
        }
        if (g[i][j] == 1){
            i--;
            continue;
        }
        if (g[i][j] == 2){
            j--;
        }
    }
    for (int i = 0; i < matchings.size() / 2; i++){
        pair<int, int> temp = matchings[i];
        matchings[i] = matchings[matchings.size() - i - 1];
        matchings[matchings.size() - i - 1] = temp;
    }
    for (int i = 0; i < n; i++){
        delete[] f[i];
        delete[] g[i];
    }
    delete[] f;
    delete[] g;
    return matchings;
}
#endif
