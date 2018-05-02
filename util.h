#ifndef UTIL_H
#define UTIL_H

#include<fstream>
#include<string>
#include<vector>
#include<cassert>
#include<iostream>
#include <iomanip>
#include<cmath>
#include <ANN/ANN.h>
#include <algorithm>
#include <unordered_set>

using namespace std;

typedef double Float;
typedef vector<Float> Point;
class Matrix{
    public:
    Matrix(int n, int m){
        this->n = n;
        this->m = m;
        entry = new Float[n*m];
    }
    ~Matrix(){
        delete entry;
    }
    int n, m;
    Float* entry;
    
};

const int LINE_LEN = 100000000;

inline double* to_double(Point& point){
    int D = point.size();
    double* ans = new double[D];

    for (int d = 0; d < D; d++){
        ans[d] = point[d];
    }
    return ans;
}

inline void print_point(Point & a, string comment){
    cerr << comment;
    cerr << ":\t";
    for (int i = 0; i < a.size(); i++){
        if (i != 0){
            cerr << " ";
        }
        cerr << a[i];
    }
    cerr << endl;
}

inline Point zero_point(int d){
    Point zero_point;
    for (int i = 0; i < d; i++){
        zero_point.push_back(0.0);
    }
    return zero_point;
}

Float dot(Point& a, Point& b){
    Float ans = 0.0;
    for (int i = 0; i < a.size(); i++){
        ans += a[i] * b[i];
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

inline double diffsq(double* a, double* b, int D){
    Float ans = 0.0;
    for (int i = 0; i < D; i++){
        ans += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return ans;
}

inline Point operator + (Point const &a, Point const &b){
    vector<Float> c;
    for (int i = 0; i < a.size(); i++){
        c.push_back(a[i] + b[i]);
    }
    return c;
}

inline Point operator + (Point const &a, Float const &c){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] + c);
    }
    return ans;
}

inline Point operator - (Point const &a, Point const &b){
    Point c;
    for (int i = 0; i < a.size(); i++){
        c.push_back(a[i] - b[i]);
    }
    return c;
}

inline Point operator - (Point const &a, double* b){
    Point c;
    for (int i = 0; i < a.size(); i++){
        c.push_back(a[i] - b[i]);
    }
    return c;
}

inline Point operator * (Point const &a, Float const &c){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] * c);
    }
    return ans;
}

inline Point operator / (Point const &a, Float const &c){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] / c);
    }
    return ans;
}

inline Point operator / (Point const &a, Point const &b){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] / b[i]);
    }
    return ans;
}

inline Point operator * (Point const &a, Point const &b){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(a[i] * b[i]);
    }
    return ans;
}

inline Point sqrt_point(Point const &a){
    vector<Float> ans;
    for (int i = 0; i < a.size(); i++){
        ans.push_back(sqrt(a[i]));
    }
    return ans;
}

inline Float normsq(Point const &a){
    Float ans = 0.0;
    for (int i = 0; i < a.size(); i++){
        ans += a[i] * a[i];
    }
    return ans;
}

inline Float norm(Point const &a){
    Float ans = 0.0;
    for (int i = 0; i < a.size(); i++){
        ans += a[i] * a[i];
    }
    return sqrt(ans);
}

inline Float distance(vector<Float> a, vector<Float> b){
    Float dist = 0.0;
    for (int i = 0; i < a.size(); i++){
        Float diff = a[i] - b[i];
        dist += diff * diff;
    }
    return dist;
};

inline Float distance(Float* a, Float* b, int D){
    Float dist = 0.0;
    for (int i = 0; i < D; i++){
        Float diff = a[i] - b[i];
        dist += diff * diff;
    }
    return dist;
};

inline Float distance(pair<Float, Float> a, pair<Float, Float> b){
    Float diff1 = a.first - b.first, diff2 = a.second - b.second;
    Float dist = diff1 * diff1 + diff2 * diff2;
    return dist;
};

vector<string> split(string str, string pattern){
    vector<string> str_split;
    size_t i=0;
    size_t index=0;
    while( index != string::npos ){

        index = str.find(pattern,i);
        str_split.push_back(str.substr(i,index-i));

        i = index+1;
    }

    if( str_split.back()=="" )
        str_split.pop_back();

    return str_split;
};

vector<string> readlines(ifstream& fin){
    vector<string> ans;
    char* line = new char[LINE_LEN];
    while (!fin.eof()){
        fin.getline(line, LINE_LEN);
        string line_str(line);
        ans.push_back(line_str);
    }
    return ans;
}

vector<vector<Point>> read_trace(string filename, int D){
    ifstream fin(filename, fstream::in);
    char* line = new char[LINE_LEN];
    vector<vector<Point>> traces;
    int line_count = 0;
    while (!fin.eof()){
        fin.getline(line, LINE_LEN);
        string line_str(line);
        if (line_str.length() < 2 && fin.eof()){
            break;
        }
        vector<string> tokens = split(line_str, " ");
        vector<Point> trace;
        int len = tokens.size() / D;

        assert(tokens.size() % D == 0);
        for (int i = 0; i < len; i++){
            Point new_point;
            for (int d = 0; d < D; d++){
                new_point.push_back(stof(tokens[i*D + d]));
            }
            trace.push_back(new_point);
        }
        traces.push_back(trace);
    }
    delete[] line;

    return traces;
};

void augment_trace(vector<vector<Point>>& traces, int n){
    for (int i = 0; i < traces.size(); i++){
        vector<Point> trace = traces[i];
        vector<Point> new_trace;
        new_trace.push_back(trace[0]);
        for (int j = 1; j < trace.size(); j++){
            for (int t = 0; t < n; t++){
                Float frac = t * 1.0 / n;
                Point middle_t = trace[j-1] * (1.0-frac) + trace[j] * frac;
                new_trace.push_back(middle_t);
            }
            new_trace.push_back(trace[j]);
        }
        traces[i] = new_trace;
    }
}

inline Point project_onto_segment(Point& p, Point& a, Point& b){
    Point ap = p - a;
    Point ab = b - a;
    Point ab_hat = ab / norm(ab);
    Float dp = dot(ap, ab_hat);
    if (dp <= 0.0){
        return a;
    }
    if (dp >= norm(ab)){
        return b;
    }
    return a + ab_hat * dp;
}

inline Point project(Point& p, vector<Point>& traj){
    Float min_dist = 1e10;
    Point ans = traj[0];
    for (int i = 0; i < traj.size()-1; i++){
        Point proj = project_onto_segment(p, traj[i], traj[i+1]);
        Float d = distance(proj, p);
        if (d < min_dist){
            min_dist = d;
            ans = proj;
        }
    }
    return ans;
}

inline string dir(string fullname){
    bool absolute = (fullname.find("/") == 0);
    vector<string> tokens = split(fullname, "/");
    string ans = "";
    if (absolute){
        ans += "/";
    }
    for (int i = 0; i < tokens.size()-1; i++){
        if (i != 0){
            ans += "/";
        }
        ans += tokens[i];
    }
    return ans;
}

vector<string> read_filelist(string filename){
    ifstream fin(filename, fstream::in);
    vector<string> filelist;
    char* line = new char[LINE_LEN];
    while (!fin.eof()){
        fin.getline(line, LINE_LEN);
        string line_str(line);
        if (line_str.size() < 2 && fin.eof()){
            break;
        }
        filelist.push_back(line_str);
    }
    fin.close();
    delete line;
    return filelist;
}

string replace(string a, string before, string after){
    int m = a.find(before);
    if (m < 0){
        return a;
    }
    string ans = a.replace(m, before.size(), after);
    return ans;
}

inline double Frechet(double* _a, double* _b, int D, int D0){
    vector<Point> a;
    int n = D / D0;
    for (int i = 0; i < n; i++){
        Point a_i(D0);
        for (int j = 0; j < D0; j++){
            a_i[j] = _a[i*D0+j];
        }
        a.push_back(a_i);
    }
    vector<Point> b;
    for (int i = 0; i < n; i++){
        Point b_i(D0);
        for (int j = 0; j < D0; j++){
            b_i[j] = _b[i*D0+j];
        }
        b.push_back(b_i);
    }
    Float** f = new Float*[2];
    f[0] = new Float[n];
    f[1] = new Float[n];
    for (int i = 0; i < n; i++){
        Float* last_f = f[0];
        last_f[i] = norm(a[i]-b[0]);
        if (i > 0 && last_f[i] < last_f[i-1]){
            last_f[i] = last_f[i-1];
        }
    }
    for (int j = 1; j < n; j++){
        Float* last_f = f[(j + 1) % 2];
        Float* cur_f = f[j % 2];
        for (int i = 0; i < n; i++){
            Float dij = norm(a[i] - b[j]); 

            cur_f[i] = last_f[i];
            if (i > 0 && cur_f[i-1] < cur_f[i]){
                cur_f[i] = cur_f[i-1];
            }
            if (i > 0 && j > 0 && last_f[i-1] < cur_f[i]){
                cur_f[i] = last_f[i-1];
            }
            if (cur_f[i] < dij){
                cur_f[i] = dij;
            }
            //cerr << "i=" << i << ", j=" << j << ", f=" << cur_f[i] << endl;
        }
    }
    Float ans = f[(n-1) % 2][n-1];
    delete[] f[0];
    delete[] f[1];
    delete[] f;
    return ans;
}

//shift up, maintain reverse index, min heap
inline void siftUp(pair<Float, int>* heap, int index, int* rev_index){
	pair<Float, int> cur = heap[index];
	while (index > 0){
		int parent = (index-1) >> 1;
		if (cur < heap[parent]){
			heap[index] = heap[parent];
			rev_index[heap[parent].second] = index;
			index = parent;
		} else {
			break;
		}
	}
	rev_index[cur.second] = index;
	heap[index] = cur;
}

//shift down, maintain reverse index, min heap
inline void siftDown(pair<Float, int>* heap, int index, int* rev_index, int size_heap){
	pair<Float, int> cur = heap[index];
	int lchild = index * 2 + 1;
	int rchild = lchild+1;
	while (lchild < size_heap){
		int next_index = index;
		if (heap[lchild] < heap[index]){
			next_index = lchild;
		}
		if (rchild < size_heap && heap[rchild] < heap[next_index]){
			next_index = rchild;
		}
		if (index == next_index)
			break;
		heap[index] = heap[next_index];
		rev_index[heap[index].second] = index;
		heap[next_index] = cur;
		index = next_index;
		lchild = index * 2 +1; rchild = lchild+1;
	}
	rev_index[cur.second] = index;
}

vector<int> eps_cover(vector<Point> points, double eps){
    int N = points.size();
    assert(N > 0);
    int D = points[0].size();
    int num_covered = 0;
    bool* covered = new bool[N];
    vector<int> uncovered;
    Float* d = new Float[N];
    vector<int> assign;
    for (int i = 0; i < N; i++){
        covered[i] = false;
        uncovered.push_back(i);
        d[i] = 1e10;
        assign.push_back(-1);
    }
    double** _points = new double*[N];
    for (int i = 0; i < N; i++){
        _points[i] = new double[D];
        double* _pi = _points[i];
        Point& pi = points[i];
        for (int d = 0; d < D; d++){
            _pi[d] = pi[d];
        }
    }
    
    ANNkd_tree* knn = new ANNkd_tree(_points, N, D);
    int K = 1;
    vector<Point> cluster;
    while (uncovered.size() > 0){
        int idx = rand() % uncovered.size();
        int v = uncovered[idx];
        if (covered[v]){
            uncovered[idx] = uncovered[uncovered.size()-1];
            uncovered.pop_back();
            continue;
        }
        int tK = K;
        vector<int> ans;
        
        while(ans.size() == 0){
            int* nn_idx = new int[tK];
            double* dists = new double[tK];
            knn->annkSearch(_points[v], tK, nn_idx, dists, 1e-6);
            for (int i = 0; i < tK; i++){
                if (dists[i] > eps * eps){
                    continue;
                }
                ans.push_back(nn_idx[i]);
            }
            delete[] dists;
            delete[] nn_idx;
            if (ans.size() == tK){
                ans.clear();
                tK *= 2;
            }
        }
        cluster.push_back(points[idx]);
        for (int i = 0; i < ans.size(); i++){
            covered[ans[i]] = true;
            Float dist = distance(points[ans[i]], points[idx]);
            if (dist < d[ans[i]]){
                d[ans[i]] = dist;
                assign[ans[i]] = idx;
            }
        }
    }

    for (int i = 0; i < N; i++){
        delete[] _points[i];
    }
    delete[] _points;
    delete[] covered;
    delete[] d;
    return assign;
}
    
inline vector<int> dijkstra(int src, int tgt, ANNkd_tree* tree, vector<vector<pair<Float, int>>>& adj){
    vector<int> ans;
    //if (src == tgt){
    //    ans.push_back(src);
    //    return ans;
    //}
    int n = tree->nPoints();
    assert(src < n);
    assert(tgt < n);
    double** points = tree->thePoints();
    int D = tree->theDim();
    Float* d = new Float[n];
    int* rec_prev = new int[n];
    bool* visited = new bool[n];
    pair<Float, int>* heap = new pair<Float, int>[n];
    int* heap_index = new int[n];
    int heap_size = n;
    for (int i = 0; i < n; i++){
        d[i] = 1e100;
        rec_prev[i] = -1;
        visited[i] = false;
        heap[i] = make_pair(d[i], i);
        heap_index[i] = i;
    }

    d[src] = 0.0;
    rec_prev[src] = src;
    heap[src].first = 0.0;
    siftUp(heap, src, heap_index);

    for (int t = 0; t < n - 1; t++){
        pair<Float, int> p = heap[0];
        heap[0] = heap[heap_size-1];
        heap_size--;
        if (heap_size > 0){
            siftDown(heap, 0, heap_index, heap_size);
        }
        int min_index = p.second;
        visited[min_index] = true;

        if (min_index == tgt){
            break;
        }

        vector<pair<Float, int>>& v = adj[min_index];

        for (int t = 0; t < v.size(); t++){
            int i = v[t].second;
            Float dist = v[t].first;
            assert(i <= n);
            if (visited[i]){
                continue;
            }
            //Float dist = distance(points[i], points[min_index], D);
            if (d[i] > d[min_index] + dist){
                d[i] = d[min_index] + dist;
                rec_prev[i] = min_index;
                heap[heap_index[i]].first = d[i];
                siftUp(heap, heap_index[i], heap_index);
            }
        }
    }
    if (d[tgt] >= 1e10){
        cerr << points[tgt][0] << " " << points[tgt][1] << " " << points[src][0] << " " << points[src][1] << endl;
    }
    assert(d[tgt] < 1e10);
    int cur = tgt;
    ans.clear();
    ans.push_back(tgt);
    while (rec_prev[cur] != cur){
        assert(rec_prev[cur] != -1);
        cur = rec_prev[cur];
        ans.push_back(cur);
    }
    reverse(ans.begin(), ans.end());
    assert(ans[0] == src);
    assert(ans[ans.size()-1] == tgt);

    delete[] rec_prev;
    delete[] visited;
    delete[] d;
    delete[] heap;
    delete[] heap_index;
    return ans;
}
        
inline vector<int> search_points(int K, double* query, ANNkd_tree* knn){
    int* nn_idx = new int[K];
    double* dists = new double[K];
    knn->annkSearch(query, K, nn_idx, dists, 1e-6);
    vector<int> ans;
    for (int i = 0; i < K; i++){
        ans.push_back(nn_idx[i]);
    }

    delete[] nn_idx;
    delete[] dists;

    return ans;
}

inline bool connected(vector<vector<pair<Float, int>>>& adj, int D, double** points){
    int N = adj.size();
    bool* visited = new bool[N];
    for (int i = 0; i < N; i++){
        visited[i] = false;
    }
    visited[0] = true;
    vector<int> q;
    q.push_back(0);
    int l = 0;
    while (l != q.size()){
        //expand q[l];
        int i = q[l];
        vector<pair<Float, int>>& adj_i = adj[i];
        for (int t = 0; t < adj_i.size(); t++){
            int j = adj_i[t].second;
            if (!visited[j]){
                visited[j] = true;
                q.push_back(j);
            }
        }
        l++;
    }
    for (int i = 0; i < N; i++){
        if (!visited[i]){
            Float min = 1e10;
            for (int j = 0; j < N; j++){
                if (visited[j]){
                    Float dij = sqrt(distance(points[i], points[j], D));
                    if (dij < min){
                        min = dij;
                    }
                }
            }
            cerr << "min=" << min << endl;
            break;
        }
    }
    if (q.size() < N){
        cerr << "visited " << q.size() << " out of " << N << endl;
        return false;
    }
    delete[] visited;
    return true;
}

inline vector<vector<pair<Float, int>>> adjacency_l2(int K, ANNkd_tree* knn, Float r){
    double** points = knn->thePoints();
    int n = knn->nPoints();
    vector<vector<pair<Float, int>>> adj;
    for (int i = 0; i < n; i++){
        vector<pair<Float, int>> adj_i;
        adj.push_back(adj_i);
    }
    int D = knn->theDim();
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        vector<int> candidates = search_points(K, points[i], knn);
        vector<pair<Float, int>>& adj_i = adj[i];
        for (int t = 0; t < candidates.size(); t++){
            int j = candidates[t];
            if (i == j){
                continue;
            }
            Float dij = sqrt(distance(points[i], points[j], D));
            if (dij <= r){
                adj_i.push_back(make_pair(dij, j));
                assert(i != j);
            }
        }
    }
    unordered_set<int> h;
    for (int i = 0; i < n; i++){
        assert(adj[i].size() > 0);
        for (int t = 0; t < adj[i].size(); t++){
            int j = adj[i][t].second;
            h.insert(i*n+j);
        }
    }
    for (int i = 0; i < n; i++){
        for (int t = 0; t < adj[i].size(); t++){
            int j = adj[i][t].second;
            Float dij = adj[i][t].first;
            assert(i != j);
            if (h.find(j*n+i) == h.end()){
                // j --> i is absent
                adj[j].push_back(make_pair(dij, i));
            }
        }
    }

    int num_edge = 0;
    for (int i = 0; i < n; i++){
        num_edge += adj[i].size();
    }
    cerr << "average #adjacent edge=" << (num_edge * 1.0 / n) << endl;
    //assert(connected(adj, D, points));
    return adj;
}
        
inline vector<vector<pair<Float, int>>> adjacency(int K, ANNkd_tree* knn, int r){
    double** points = knn->thePoints();
    int n = knn->nPoints();
    vector<vector<pair<Float, int>>> adj;
    for (int i = 0; i < n; i++){
        vector<pair<Float, int>> adj_i;
        adj.push_back(adj_i);
    }
    int D = knn->theDim();
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        vector<int> candidates = search_points(K, points[i], knn);
        vector<pair<Float, int>>& adj_i = adj[i];
        for (int t = 0; t < candidates.size(); t++){
            int j = candidates[t];
            if (i == j){
                continue;
            }
            bool flag = true;
            Float dij = distance(points[i], points[j], D);
            for (int t2 = 0; t2 < candidates.size(); t2++){
                int k = candidates[t2];
                if (j == k || i == k){
                    continue;
                }

                Float dik = distance(points[i], points[k], D);
                Float dkj = distance(points[k], points[j], D);
                if (dij > dik + dkj){
                    flag = false;
                    break;
                }
            }
            if (flag){
                adj_i.push_back(make_pair(dij, j));
                assert(i != j);
            }
        }
    }
    unordered_set<int> h;
    for (int i = 0; i < n; i++){
        for (int t = 0; t < adj[i].size(); t++){
            int j = adj[i][t].second;
            h.insert(i*n+j);
        }
    }
    for (int i = 0; i < n; i++){
        for (int t = 0; t < adj[i].size(); t++){
            int j = adj[i][t].second;
            Float dij = adj[i][t].first;
            assert(i != j);
            if (h.find(j*n+i) == h.end()){
                // j --> i is absent
                adj[j].push_back(make_pair(dij, i));
            }
        }
    }

    int num_edge = 0;
    for (int i = 0; i < n; i++){
        num_edge += adj[i].size();
    }
    cerr << "average #adjacent edge=" << (num_edge * 1.0 / n) << endl;
    return adj;
}

ANNkd_tree* build_ann_tree(vector<Point> points){
    int N = points.size();
    int D = points[0].size();
    double** _points = new double*[N];
    for (int i = 0; i < N; i++){
        _points[i] = new double[D];
        double* _pi = _points[i];
        for (int d = 0; d < D; d++){
            _pi[d] = points[i][d];
        }
    }
    
    ANNkd_tree* tree = new ANNkd_tree(_points, N, D);
    return tree;
}

// M * ans = b
// M is square, tri-diagonal (l, diag, r) matrix
inline vector<Point> tri_diag_solve(vector<Point> b, vector<Float> diag, vector<Float> l, vector<Float> r){
    int N = diag.size();
    assert(l.size() == N-1);
    assert(r.size() == N-1);
    int D = b[0].size();
    vector<Point> ans;
    for (int i = 0; i < N; i++){
        ans.push_back(Point(D, 0.0));
    }
    for (int i = 1; i < N; i++){
        Float ratio = - l[i-1] / diag[i-1];
        diag[i] += r[i-1] * ratio;
        b[i] = b[i] +  b[i-1] * ratio;
    }
    // Now l is all zero
    for (int i = N-1; i >= 0; i--){
        ans[i] = b[i] / diag[i];
        if (i > 0){
            b[i-1] = b[i-1] - ans[i] * r[i-1];
        }
    }
    
    return ans;
}

inline vector<Point> smooth_strong_head(Float lambda, vector<Point> b, Float stubborn){
    int N = b.size();
    assert(N > 0);
    int D = b[0].size();
    // a = (I + lambda * L)^{-1} b
    vector<Float> diag;
    vector<Float> l;
    vector<Float> r;
    for (int i = 0; i < N; i++){
        Float di = 1.0;
        if (i == 0 || i == N-1){
            di += lambda / stubborn;
        } else {
            di += 2 * lambda;
        }
        diag.push_back(di);
    }
    for (int i = 0; i < N-1; i++){
        l.push_back(-lambda);
        r.push_back(-lambda);
    }
    r[0] /= stubborn;
    l[l.size()-1] /= stubborn;
    return tri_diag_solve(b, diag, l, r);
}

inline vector<Point> smooth(Float lambda, vector<Point> b){
    int N = b.size();
    assert(N > 0);
    int D = b[0].size();
    // a = (I + lambda * L)^{-1} b
    vector<Float> diag;
    vector<Float> l;
    vector<Float> r;
    for (int i = 0; i < N; i++){
        Float di = 1.0;
        if (i == 0 || i == N-1){
            di += lambda;
        } else {
            di += 2 * lambda;
        }
        diag.push_back(di);
    }
    for (int i = 0; i < N-1; i++){
        l.push_back(-lambda);
        r.push_back(-lambda);
    }
    return tri_diag_solve(b, diag, l, r);
}

// Gauss Elimination A x = b, for sparse A, get x, A looks like (I + 2 lambda Laplacian)
//inline vector<Point> sparse_ge(vector<vector<pair<int, Float>>> A, vector<Point> b){
//    vector<Point> ans;
//    for (int i = 1; i < N; i++){
//        vector<pair<int, Float>>& Ai = A[i];
//        while (Ai[0].first != i){
//            int j = Ai[0].first;
//            Float ratio = Ai[0].second / A[j].second;
//        }
//        for (vector<pair<int, Float>>::iterator it = Ai.begin(); it != Ai.end(); it++){
//            
//        }
//    }
//    return ans;
//}

//inline vector<Point> structural_smooth(Float lambda, vector<vector<pair<Float, int>>>& adj, vector<Point> b){
//    vector<vector<pair<int, Float>>> sparse_rows;
//    int N = b.size();
//    for (int i = 0; i < N; i++){
//        vector<pair<int, Float>> sr;
//        int deg = 0;
//        for (int t = 0; t < adj[i].size(); t++){
//            int j = adj[i][t].second;
//            if (i == j){
//                continue;
//            }
//            deg++;
//            sr.push_back(make_pair(j, -1.0 * lambda));
//        }
//        sr.push_back(make_pair(i, 1.0 + deg * 2.0 * lambda));
//    }
//    
//    delete[] deg;
//}

inline vector<int> project(ANNkd_tree* tree, vector<Point> points){
    vector<int> ans;
    for (int i = 0; i < points.size(); i++){
        int* nn_idx = new int[1];
        Float* dists = new Float[1];
        int D = points[i].size();
        assert(D == 2);
        Float* pi = new Float[D];
        for (int d = 0; d < D; d++){
            pi[d] = points[i][d];
        }
        tree->annkSearch(pi, 1, nn_idx, dists, 1e-6);
        ans.push_back(nn_idx[0]);

        delete[] nn_idx;
        delete[] dists;
        delete[] pi;
    }
    return ans;
}

inline int project(ANNkd_tree* tree, Point point){
    int* nn_idx = new int[1];
    Float* dists = new Float[1];
    int D = point.size();
    Float* pi = new Float[D];
    for (int d = 0; d < D; d++){
        pi[d] = point[d];
    }
    tree->annkSearch(pi, 1, nn_idx, dists, 1e-6);
    int ans = nn_idx[0];

    delete[] nn_idx;
    delete[] dists;
    delete[] pi;
    return ans;
}

inline vector<Point> unique(vector<Point> p){
    vector<Point> ans;
    ans.push_back(p[0]);
    for (int t = 1; t < p.size(); t++){
        if (distance(ans[ans.size()-1], p[t]) <= 1e-15){
            continue;
        }
        ans.push_back(p[t]);
    }
    return ans;
}

#endif


