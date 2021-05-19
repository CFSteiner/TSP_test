// tsp_cpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <math.h>       /* pow */
#include <stdio.h>
#include "tsp.h"


int read_x_y(std::ifstream &temp_file,int &nodeCount
    ,std::vector<float> &x, std::vector<float>& y)
{
    int ierr;
    float xi, yi;
    std::string CARD;

    if (temp_file.is_open()) {
        std::getline(temp_file, CARD);
        nodeCount = std::stoi(CARD);
        while (std::getline(temp_file, CARD)) {
            std::stringstream str_ss(CARD);
            str_ss >> xi >> yi;
            x.push_back(xi);
            y.push_back(yi);
        }
    }
    ierr = 0;
    return ierr;
}

std::vector<std::vector<float> > dist_vector(int &i, int &j, int nodeCount, std::vector<float> &x, std::vector<float> &y) {
    std::vector<std::vector<float>> dist(nodeCount, std::vector<float>(nodeCount, 0));
    for (i = 0; i < nodeCount; i++) {
        for (j = 0; j < nodeCount; j++) {
            dist[i][j] = length(x[i], y[i], x[j], y[j]);
        }
    }
    return dist;
}


float length(float xi, float yi, float xj, float yj) {
    float len = 0;
    len = sqrt(pow((xi - xj),2.0) + pow((yi - yj),2.0));
    return len;
}

void initialize(std::vector<std::vector<Cell>> &cells) {
    auto nx = cells.size();
    auto ny = cells[0].size();
    for (int ix = 0; ix < nx; ix++) {  for (int iy = 0; iy < ny; iy++) {
            cells[ix][iy].nmax = 0;
        }  }
}



/*-------------------------------------------------------------- - !
!###############################################################*/
float calc_distance(int &n, std::vector<int> &sno, std::vector<std::vector<float>> &dist) {
    float disti = 0;
    disti = dist[sno[n - 1]][sno[0]];
    for (int i = 0; i < n - 1; i++) disti += dist[sno[i]][sno[i + 1]];
    return disti;
}
/*##############################################################!
!-------------------------------------------------------------- */


void swap3opt(const int nodeCount, std::vector<int> &sno, const std::vector<std::vector<float>> &dist, const int i, const int k, const int l) {
    int a=0, b, c, d, e, f;
    float d0, d1, d2, d3, d4;
    if (i > 0)  a = sno[i - 1];
    if (i == 0) a = sno[nodeCount - 1];
    b = sno[i];
    c = sno[k - 1]; d = sno[k];
    e = sno[l - 1]; f = sno[l];

    d0 = dist[a][b] + dist[c][d] + dist[e][f];
    d1 = dist[a][c] + dist[b][d] + dist[e][f];
    d2 = dist[a][b] + dist[c][e] + dist[d][f];
    d3 = dist[a][d] + dist[e][b] + dist[c][f];
    d4 = dist[f][b] + dist[c][d] + dist[e][a];


    if (d0 > d1) {
        auto i0 = sno.begin();
        auto i1 = i0 + i;
        auto i2 = i0 + k;
        std::reverse(i1, i2);
    }
    else if (d0 > d2) {
        auto i0 = sno.begin();
        auto i1 = i0 + k;
        auto i2 = i0 + l;
        std::reverse(i1, i2);
    }
    else if (d0 > d4) {
        auto i0 = sno.begin();
        auto i1 = i0 + i;
        auto i2 = i0 + l;
        std::reverse(i1, i2);
    }
    else if (d0 > d3) {
        auto i0 = sno.begin();
        auto k0 = i0 + k;
        auto l0 = i0 + l;
        i0 = i0 + i;
        std::vector<int> temp(0);
        auto t0 = temp.begin();
        temp.insert(t0,k0,l0);
        temp.insert(temp.end(), i0, k0);
        std::copy(temp.begin(), temp.end(), i0);
    }

}
void swap3opt_rand(const int nodeCount, std::vector<int>& sno, const std::vector<std::vector<float>> &dist, const int i, const int k, const int l) {
    int a=0, b, c, d, e, f;
    float d0, d1, d2, d3, d4;
    std::vector<int>::iterator i0,i1,k0,l0,i2,t0;
    std::vector<int> temp(0);
    int dif1, dif2, dif3,iii;
    if (i > 0)  a = sno[i - 1];
    if (i == 0) a = sno[nodeCount - 1];
    b = sno[i];
    c = sno[k - 1]; d = sno[k];
    e = sno[l - 1]; f = sno[l];

    d0 = dist[a][b] + dist[c][d] + dist[e][f];
    d1 = dist[a][c] + dist[b][d] + dist[e][f];
    d2 = dist[a][b] + dist[c][e] + dist[d][f];
    d3 = dist[a][d] + dist[e][b] + dist[c][f];
    d4 = dist[f][b] + dist[c][d] + dist[e][a];

    int iopt = 0;
    if (d0 > d1) {
        iopt = 1;
        i0 = sno.begin();
        i1 = i0 + i;
        i2 = i0 + k;
        std::reverse(i1, i2);
    }
    else if (d0 > d2) {
        iopt = 1;
        i0 = sno.begin();
        i1 = i0 + k;
        i2 = i0 + l;
        std::reverse(i1, i2);
    }
    else if (d0 > d4) {
        iopt = 1;
        i0 = sno.begin();
        i1 = i0 + i;
        i2 = i0 + l;
        std::reverse(i1, i2);
    }
    else if (d0 > d3) {
        iopt = 1;
        i0 = sno.begin();
        k0 = i0 + k;
        l0 = i0 + l;
        i0 = i0 + i;
        dif1 = l - k;
        dif2 = k - i;
        dif3 = i - l + 1;
        iii = 0;
        t0 = temp.begin();
        temp.insert(t0, k0, l0);
        iii = 0;
        temp.insert(temp.end(), i0, k0);
        std::copy(temp.begin(), temp.end(), i0);
        iii = 0;
    }

    if (iopt == 0) {
        int iselect = rand() % 4 + 1;
        switch (iselect) {
        case 1:
            i0 = sno.begin();
            i1 = i0 + i;
            i2 = i0 + k;
            std::reverse(i1, i2);
            break;
        case 2:
            i0 = sno.begin();
            i1 = i0 + k;
            i2 = i0 + l;
            std::reverse(i1, i2);
            break;
        case 3:
            i0 = sno.begin();
            i1 = i0 + i;
            i2 = i0 + l;
            std::reverse(i1, i2);
            break;
        case 4:
            iopt = 1;
            i0 = sno.begin();
            k0 = i0 + k;
            l0 = i0 + l;
            i0 = i0 + i;
            dif1 = l - k;
            dif2 = k - i;
            dif3 = i - l + 1;
            iii = 0;
            t0 = temp.begin();
            temp.insert(t0, k0, l0);
            iii = 0;
            temp.insert(temp.end(), i0, k0);
            std::copy(temp.begin(), temp.end(), i0);
            iii = 0;
            break;
        default:  
            iii = 0;
        }
    }
}

void swap3opt_rand_p(int nodeCount, int *sno, std::vector<std::vector<float>> dist, int i, int k, int l) {
    float a, b, c, d, e, f;
    float d0, d1, d2, d3, d4;

}

void cell_creator(std::vector<std::vector<int> > &cell_id, std::vector<std::vector<Cell> > &cells, int &i, int &ix, int &iy, int nodeCount, int nx, int ny, std::vector<float> &x, std::vector<float> &xb1, std::vector<float> &xb2, std::vector<float> &y, std::vector<float> &yb1, std::vector<float> &yb2) {
    ix = 0; iy = 0;
    for (i = 0; i < nodeCount; i++) {
        for (ix = 0; ix < nx; ix++) if (x[i] > xb1[ix] && x[i] < xb2[ix]) break;
        for (iy = 0; iy < ny; iy++) if (y[i] > yb1[iy] && y[i] < yb2[iy]) break;
        cells[ix][iy].nmax += 1;
        cells[ix][iy].jz.push_back(i);
        cell_id[i][0] = ix;
        cell_id[i][1] = iy;
    }
}
