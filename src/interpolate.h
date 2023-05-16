#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

#include <vector>
#include <cstdlib>
#include "input.h"
#include <omp.h>
#include <cmath>
#include "mytimer.h"

class Interpolator{
    public:
    dist_data *dist;
    double* M;
    Interpolator(dist_data *dist):dist(dist){
        mytimer::start("itp calc");
        std::cout << "[itp] calculating coefficents\n";

        double h = dist->dr;
        double *f = dist->f;
        int n = dist->mesh - 1;
        
        M = (double*)malloc(sizeof(double) * dist->mesh);
        M[0] = M[n] = 0.0;

        if(n == 2){
            M[1] = 1.5 * (f[2] - f[1] - f[1] + f[0]) / h / h;
            return;
        }else if(n == 3){
            double c1_div3 = (f[2] - f[1] - f[1] + f[0]) / h / h; // c1 / 3;
            double c2_div3 = (f[3] - f[2] - f[2] + f[1]) / h / h; // c2 / 3;
            M[1] = 1.6 * c1_div3 - 0.4 * c2_div3;
            M[2] = 1.6 * c2_div3 - 0.4 * c1_div3;
            return;
        }

        std::vector<double> u(n - 2);
        std::vector<double> a(n - 2);
        // calc u, a
        a[0] = 0.25, u[0] = 1.875;
        for(int i = 1;i < n - 2;i++){
            a[i] = 0.5 / u[i - 1];
            u[i] = 2.0 - 0.5 * a[i];
        }

        std::vector<double> c(n - 1);
        double coe = 3.0 / h / h;

        omp_set_num_threads(omp_get_max_threads());
        #pragma omp parallel
        {
            for(int i = 0;i <= n - 2;i++){
                c[i] = coe * (f[i + 2] - f[i + 1] - f[i + 1] + f[i]);
            }
        }

        for(int i = 1;i <= n - 2;i++){
            c[i] -= a[i - 1] * c[i - 1];
        }
        
        M[n - 1] = c[n - 2] / u[n - 3];
        for(int i = n - 2;i >= 2;i--) M[i] = (c[i - 1] - 0.5 * M[i + 1]) / u[i - 2];
        M[1] = 0.5 * c[0] - 0.25 * M[2];

        std::cout << "[itp] coeffients calculated\n";
        mytimer::end("itp calc");
    }

    void unload(){
        free(M);
    }

    double F(double x){
        if(x > dist->cutoff) return 0.0;
        double h = dist->dr;
        int n0 = (int)floor(x / h);
        // [n0, n0 + 1]
        double M0 = M[n0], M1 = M[n0 + 1];
        double f0 = dist->f[n0], f1 = dist->f[n0 + 1];
        double x0 = x - n0 * h; double x1 = x0 - h;
        
        double ans = M1 * x1 * x1 * x1 - M0 * x0 * x0 * x0 + x0 * (6 * f1 - M1 * h * h) - x1 * (6 * f0 - M0 * h * h);
        return ans / 6 / h;
    }
};

#endif