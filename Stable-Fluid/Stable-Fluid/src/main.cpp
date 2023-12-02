#include "fluid-Euler.h"
#include "geometry.h"

#include <iostream>
#include <math.h>
//#include <cuda_runtime.h>

int main() {
    int n1 = 30;
    int n2 = 30;
    int n3 = 30;
    float l = 0.01f;      

    std::vector<float> phi;
    std::vector<float> solid_phi;
    obj_2_SDF(n1, n2, n3, 16, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test.obj", phi);    
    obj_2_SDF(n1, n2, n3, 28, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test.obj", solid_phi, true);

    Fluid_Euler fluid_Euler(n1, n2, n3, l, phi, solid_phi);
    fluid_Euler.show(0.002f, 1);
    return 0;
}
