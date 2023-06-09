#include "calc.h"
#include "interpolate.h"
#include "input.h"
#include <vector>
#include "myerror.h"
#include <cstring>
#include <algorithm>
#include <iostream>
#include "omp.h"
#include "mytimer.h"
using std::max;
using std::min;

using u64 = unsigned long long;

const u64 MAX_TOTAL_BUFSIZE = 200000000ULL;
inline double S2(const double& x){return x * x;}

void calculator(double* H, int N, Interpolator& itp, v_data &V,
    std::vector<point_data> &points, double dx, double dy, double dz){

    mytimer::start("calc");
    std::cout << "[calculator] start\n";

    memset(H, 0, sizeof(double) * N * N);
    
    u64 n_total = 1ULL * V.nx * V.ny * V.nz;
    u64 fb_size = min(n_total, MAX_TOTAL_BUFSIZE / N);
    double** fbuf = (double**)malloc(sizeof(double*) * N);
    for(int i = 0;i < N;i++){
        fbuf[i] = (double*)malloc(sizeof(double) * fb_size);
    }
    
    u64 mod_yz = 1ULL * V.ny * V.nz; // to get id;

    double dV = dx * dy * dz;

    for(u64 l = 0ULL, r_ = fb_size;l < n_total;l += fb_size, r_ += fb_size){
        u64 r = min(n_total, r_);
        
        std::cout << "[calculator] calculating f in [" << l << "," << r << ")\n";
        
        // calc f_{ai}
        #pragma omp parallel for collapse(2)
        for(int i = 0;i < N;i++){
            for(u64 p = l;p < r;p++){
                point_data& pt = points[i];
                int ix = p / mod_yz, iy = (p / V.nz) % V.ny, iz = p % V.nz;
                double d_2 = S2(pt.x-ix*dx)+S2(pt.y-iy*dy)+S2(pt.z-iz*dz);
                fbuf[i][p - l] = itp.F(sqrt(d_2));
            }
        }

        #pragma omp barrier

        std::cout << "[calculator] reducing in [" << l << "," << r << ")\n";

        for(int i = 0;i < N;i++){
            for(int j = i;j < N;j++){
                double s = 0.0;

                #pragma omp parallel for reduction(+:s)
                for(u64 p = l;p < r;p++){
                    s += fbuf[i][p - l] * V.d[p] * fbuf[j][p - l] * dV;
                }
                
                H[i * N + j] += s;
            }
        }

    }

    for(int i = 0;i < N;i++) for(int j = 0;j < i;j++)
        H[i * N + j] = H[j * N + i];

    std::cout << "[calculator] end\n";
    mytimer::end("calc");
}