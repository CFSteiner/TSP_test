//
//  tsp.h
//  solver_tsp
//
//  Created by Christoph Steiner on 07.02.21.
//  Copyright © 2021 Christoph Steiner. All rights reserved.
//

#ifndef tsp_h
#define tsp_h
struct Cell {
    int nmax;
    std::vector<int> jz;
};

int read_x_y(std::ifstream &temp_file,int &nodeCount
             ,std::vector<float> &x, std::vector<float>& y);
float length(float xi, float yi, float xj, float yj);
void initialize(std::vector<std::vector<Cell>> &cells);
float calc_distance(int &n, std::vector<int> &sno, std::vector<std::vector<float>> &dist);
void swap3opt(const int nodeCount, std::vector<int>& sno, const std::vector<std::vector<float>> &dist, const int i, const int k, const int l);
void swap3opt_rand(const int nodeCount, std::vector<int>& sno, const std::vector<std::vector<float>> &dist, const int i, const int k, const int l);
void swap3opt_rand_p(int nodeCount, int* sno, std::vector<std::vector<float>> dist, int i, int k, int l);
std::vector<std::vector<float> > dist_vector(int &i, int &j, int nodeCount, std::vector<float> &x, std::vector<float> &y);
void cell_creator(std::vector<std::vector<int> > &cell_id, std::vector<std::vector<Cell> > &cells, int &i, int &ix, int &iy, int nodeCount, int nx, int ny, std::vector<float> &x, std::vector<float> &xb1, std::vector<float> &xb2, std::vector<float> &y, std::vector<float> &yb1, std::vector<float> &yb2);
#endif /* tsp_h */
