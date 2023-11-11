#include "fluid-Euler.h"

#include <iostream>
#include <math.h>
//#include <cuda_runtime.h>

int main() {
    int n1 = 20;
    int n2 = 20;
    int n3 = 20;
    float l = 0.01f;
    std::vector<float> phi;
    std::vector<int> solid;
    phi.resize(n1 * n2 * n3);
    solid.resize(n1 * n2 * n3, 0);
    // set solid
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            solid[(0) * n2 * n3 + (i)*n3 + (j)] = 1;
            solid[(n1 - 1) * n2 * n3 + (i)*n3 + (j)] = 1;
            solid[(i)*n2 * n3 + (0) * n3 + (j)] = 1;
            solid[(i)*n2 * n3 + (n2 - 1) * n3 + (j)] = 1;
            solid[(i)*n2 * n3 + (j)*n3 + (0)] = 1;
            solid[(i)*n2 * n3 + (j)*n3 + (n3 - 1)] = 1;
        }
    }
    // set phi
    float d[6] = { 0 };
    float dis;
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
        d[0] = -(15.0f - i) * l - 0.5 * l;
        d[1] = -(i - 5.0f) * l - 0.5 * l;
        d[2] = -(15.0f - j) * l - 0.5 * l;
        d[3] = -(j - 5.0f) * l - 0.5 * l;
        d[4] = -(15.0f - k) * l - 0.5 * l;
        d[5] = -(k - 5.0f) * l - 0.5 * l;
        dis = d[0];
        for (int idx = 1; idx < 6; idx++) {
            if (d[idx] < 0) {
                dis = d[idx] > dis ? d[idx] : dis;
            }
            else {
                if (dis < 0) {
                    dis = d[idx];
                }
                else {
                    dis = d[idx] < dis ? d[idx] : dis;
                }
            }
        }
        phi[i * n2 * n3 + j * n3 + k] = dis;
    }

    if (0) {
        for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
            if (i < 4 || i>16 || j < 4 || j>16 || k < 4 || k>16) {
                phi[i * n2 * n3 + j * n3 + k] = 100.0f;
            }
            else {
                phi[i * n2 * n3 + j * n3 + k] = -100.0f;
            }
        }
    }

    Fluid_Euler fluid_Euler(n1, n2, n3, l, phi, solid);
    fluid_Euler.show(0.001f, 1);
    return 0;
}
