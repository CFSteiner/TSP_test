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


static void check_if_continue(int &icount, float &t_init) {
    if (icount == 0) {
        std::cout << "Continue how many iterations? Negative number fo rdoubling T_INIT again\n";
        std::cin >> icount;
        if (icount < 0) {
            icount = abs(icount);
            t_init = 2. * t_init;
        }
    }
}

//!###################################################################
int main()
{
    int iii, i, j, k, l, nodeCount, ierr;
    int nx, ny, ixi, iyi,iz,jz,nj;
    int imin;
    float dx, dy;
    float disti;
    std::ifstream temp_file;
    std::ofstream out_file;
    std::string CARD;
    std::stringstream str_ss(CARD);
    std::vector<float> x;
    std::vector<float> y;

    /*---------------------------------------------------------------- - !
     !################        S T A R T       ##########################
    !---------------------------------------------------------------- - */

    /* Read input data from tmp.data file */
    temp_file.open("tmp.data");
    if (!temp_file) {
        std::cout << "Cannot open file \n";
        return (1);
    }
    else {
        ierr = read_x_y(temp_file, nodeCount, x, y);
        temp_file.close();
    }
    iii = 0;


    /* Create distance array with pointers */
    std::vector<std::vector<float> > dist = dist_vector(i, j, nodeCount, x, y);

    float xmin = *std::min_element(x.begin(), x.end());
    float xmax = *std::max_element(x.begin(), x.end());
    float ymin = *std::min_element(y.begin(), y.end());
    float ymax = *std::max_element(y.begin(), y.end());
    xmin -= 0.01; xmax += 0.01;
    ymin -= 0.01; ymax += 0.01;

    nx = 20; ny = 20;
    dx = (xmax - xmin) / nx;
    dy = (ymax - ymin) / ny;


    std::vector<float> xb1(nx), xb2(nx), yb1(ny), yb2(ny);
    for (int ix = 0; ix < nx; ix++) {
        xb1[ix] = xmin + ix * dx;
        xb2[ix] = xmin + (ix + 1) * dx;
    }
    for (int iy = 0; iy < ny; iy++) {
        yb1[iy] = ymin + iy * dy;
        yb2[iy] = ymin + (iy + 1) * dy;
    }

    /* Create Cells */
    std::vector<std::vector<Cell>> cells(nx, std::vector<Cell>(ny));
    std::vector<std::vector<int>> cell_id(nodeCount, std::vector<int>(2));
    initialize(cells);
    iii = 0;
    int ix;
    int iy;
    cell_creator(cell_id, cells, i, ix, iy, nodeCount, nx, ny, x, xb1, xb2, y, yb1, yb2);
    /* end create cells */

    /* ---------------------------------------------------------------- - !
      !################       INITIAL LIST    ##########################
     !---------------------------------------------------------------- - */

    std::vector<int> sno_0(nodeCount), sno(nodeCount), sno_new(nodeCount);
    iii = 0;
    int i0 = 0;
    jz = 0;
    float distmin = 999999999.99;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = std::max(0,ix-1); ix < std::min(ny,ix+1); iy++) {
            distmin = 999999999.99;
            for (int nj = 0; nj < cells[ix][iy].nmax; nj++) {
                jz = cells[ix][iy].jz[nj];
                disti = length(0.0, 0.0, x[jz], y[jz]);
                if (disti < distmin) {
                    distmin = disti;
                    i0 = jz;
                }
            }
            if (distmin < 999999999.9) goto found_starting_point;
        }
    }
    found_starting_point:
    int istart = i0;

    int db = 2;
    sno_0[0] = i0;
    int ncount = nodeCount - 1;
    i = 0;
    std::vector<bool> ytaken(nodeCount,false);
    ytaken[sno_0[0]] = true;
    bool yfound = false;
    for (;;) {
        if (ncount < 1) break;
        ixi = cell_id[sno_0[i]][0];
        iyi = cell_id[sno_0[i]][1];
        iz = sno_0[i];
        distmin = 999999999.99;
        imin = -1;
        yfound = false;
        for (ix = std::max(0, ixi - 2 * db); ix < std::min(nx, ixi + db); ix++) {
            for (iy = std::max(0, iyi - 2 * db); iy < std::min(ny, iyi + db); iy++) {
                for (nj = 0; nj < cells[ix][iy].nmax; nj++) {
                    jz = cells[ix][iy].jz[nj];
                    if (ytaken[jz]) continue;
                    if (dist[iz][jz] < distmin) {
                        imin = jz;
                        distmin = dist[iz][jz];
                    }
                }
            }
        }

        if (imin > -1) {
            i++;
            sno_0[i] = imin;
            ytaken[imin] = true;
            ncount--;
        }
        else {
            db++;
        }
    }
    iii = 0;
    cells.clear();
    //------------------------------ - END CREATE INITIAL SNO LIST
    float current_dist;
    current_dist = calc_distance(nodeCount, sno_0, dist);
    int iter = 0;
    //CALL WRITE_OUT(NODECOUNT, SNO_0, ITER, DIST) */

    /* ---------------------------------------------------------------- - !
      !################       3 -OPT    ##########################
     !---------------------------------------------------------------- - */
    bool ynotre = true,ynot;
    int itmax = 5;
    if (nodeCount < 500) itmax = 10;
    sno = sno_0;
    sno_0.clear();
    int it = 0;
    int ixk, iyk, ixl, iyl;
    float new_dist;

    for (;;) {
        it++;
        if (!ynotre) break;
        if (it > itmax) break;
        current_dist = calc_distance(nodeCount, sno, dist);
        printf("%d %f \n", it, current_dist);

        for (i = 0; i < nodeCount - 2; i++) {
            for (k = i+1; k < nodeCount - 1; k++) {
                for (l = k + 1; l < nodeCount - 1; l++) {
                    ixi = cell_id[i][0]; iyi = cell_id[i][1];
                    ixk = cell_id[k][0]; iyk = cell_id[k][1];
                    ixl = cell_id[l][0]; iyl = cell_id[l][1];
                    ynot = ixi > ixk + 1 || ixi < ixk - 3;
                    ynot = ynot || ixi > ixl + 1 || ixi < ixl - 3;
                    ynot = ynot || iyi > iyk + 1 || iyi < iyk - 3;

                    sno_new = sno;
                    swap3opt(nodeCount, sno_new, dist, i, k, l);
                    new_dist = calc_distance(nodeCount, sno_new, dist);
                    if (new_dist + 0.00001 < current_dist) {
                        sno = sno_new;
                        current_dist = new_dist;
                        goto end_l1;
                    }
                }
            }
        }
        end_l1:
        iii = 0;
    }


    /*--------------------------------------------------------------- - !
    !---------------------------------------------------------------- - !
    !################      Simulated Annealing  ######################
    !----------------------------------------------------------------  */
    std::vector<int> sno_best(nodeCount);
    int itmax_ini = 70000;
    float t_init = 5 * sqrt(nodeCount);
    float alpha = 12.0;
    float temp0 = 0.01;
    if (nodeCount > 1000) itmax_ini = 50000;
    db = 3;
    sno_best = sno;
    float dist_best = current_dist;
    int icount = 10;
    bool yaccept;
    float temp,delta,prop,frand;
    for (iter =1;;iter++) {
        printf("%d \n", iter);
        printf("Current best solution %f \n", dist_best);
        icount--;
        check_if_continue(icount, t_init);
        if (icount == 0) break;

        t_init = t_init * 0.98;
        itmax = itmax_ini + int(itmax_ini * (iter) / 8.);
        temp = t_init * exp(-alpha*float(it/itmax))+temp0;

        for (it = 0; it < itmax; it++) {
            if (it % 1000 == 0) printf("%d %f %f \n", it, temp, current_dist);
            i = rand() % (nodeCount - 3);
            k = rand() % (nodeCount - 2 - i) + i+1;
            l = rand() % (nodeCount - 1 - k) + k+1;
            sno_new = sno;
            iii = 0;
            swap3opt_rand(nodeCount, sno_new, dist, i, k, l);
            new_dist = calc_distance(nodeCount, sno_new, dist);
            yaccept = false;
            if (new_dist < current_dist) {
                yaccept = true;
            }
            else if (new_dist > current_dist + 0.0001) {
                delta = new_dist - current_dist;
                prop = exp(-delta / temp);
                frand = (float) rand()/(RAND_MAX);
                if (prop > frand) yaccept = true;
            }

            if (yaccept) {
                current_dist = new_dist;
                sno = sno_new;
            }

            if(it % 100 == 0) temp = t_init * exp(-alpha * float(it) / float(itmax)) + temp0;

        }

        if (current_dist < dist_best) {
            sno_best = sno;
            dist_best = current_dist;
        }

    }

    cell_id.clear();
    dist.clear();
    iii = 0;
    return 0;
}


/*-------------------------------------------------------------- - !
!###############################################################*/
