#include "fluid-Euler.h"
#include "geometry.h"

#include <iostream>
#include <math.h>
//#include <cuda_runtime.h>

int main() {
    int n1 = 80;
    int n2 = 80;
    int n3 = 80;
    float l = 0.01f;      

    std::vector<float> phi;
    std::vector<float> solid_phi;
    //obj_2_SDF(n1, n2, n3, 16, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test.obj", phi);    
    obj_2_SDF_py(n1, n2, n3, 60, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/bunny.obj", phi);
    obj_2_SDF(n1, n2, n3, 76, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test.obj", solid_phi, true);

    Fluid_Euler fluid_Euler(n1, n2, n3, l, phi, solid_phi);
    //fluid_Euler.show(0.002f, 1);
    fluid_Euler.run(0.002f, 1);    
    //fluid_Euler.outputPLY(30, 20, 0.002f, 1);
    //fluid_Euler.simple_render("C:/Users/11862/Desktop/output_video.mp4");
    return 0;
}
