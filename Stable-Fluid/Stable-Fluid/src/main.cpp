#include "fluid-Euler.h"
#include "geometry.h"

#include <iostream>
#include <math.h>
//#include <cuda_runtime.h>

int main() {
    omp_set_num_threads(12);
    int n1 = 40;
    int n2 = 40;
    int n3 = 80;
    float l = 0.01f;

    std::vector<float> phi;
    std::vector<float> phi2;
    std::vector<float> solid_phi;
    //obj_2_SDF(n1, n2, n3, 16, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test.obj", phi);    
    //obj_2_SDF_py(n1, n2, n3, 60, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/bunny3.obj", phi, false);
    //obj_2_SDF_py(n1, n2, n3, l, "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test.obj", phi, 0.3f, 0.2f, 0.2f, 0.56, false);

    std::string cube_path = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test.obj";
    std::string cuboid_path = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/cuboid.obj";
    std::string boundary_path = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/OBJ/test2.obj";

    obj_2_SDF_py(n1, n2, n3, 36, l, boundary_path, solid_phi, true);
    obj_2_SDF_py(n1, n2, n3, l, cuboid_path,phi, 0.16f, 0.2f, 0.2f, 0.64, false);
    //obj_2_SDF_py(n1, n2, n3, l, cuboid_path, phi2, 0.16f, 0.2f, 0.2f, 0.24, false);
    //phi = union_phi(phi, phi2);
    Fluid_Euler fluid_Euler(n1, n2, n3, l, phi, solid_phi);
    //fluid_Euler.set_velocity(glm::vec3(0.0f, -2.0f, 0.0f));    
    fluid_Euler.add_rigidBody(
        cube_path,
        0.13f,
        glm::vec3(0.2f, 0.085f, 0.15f),
        0.0012f,
        glm::mat3(
            3e-5f, 0.0f, 0.0f,
            0.0f, 3e-5f, 0.0f,
            0.0f, 0.0f, 3e-5f),
        glm::vec3(0.2f, 0.85f, 0.15f),
        glm::quat(1.0f, 0.0f, 0.0f, 0.0f),
        glm::vec3(0.0f),
        glm::vec3(0.0f, 0.0f, 0.0f)
    );
    fluid_Euler.set_fixedOperation();

    int fps = 30;
    int t = 5;
    //fluid_Euler.show(0.001f, 1);
    //fluid_Euler.run(0.001f, 1);    
    fluid_Euler.outputPLY(fps, t, 0.001f, 5);
    //fluid_Euler.outputXYZ(3, 1, 0.001f, 3);    
    fluid_Euler.pbrt_render("C:/Users/11862/Desktop/check_3.mp4", fps, t);
    return 0;
}
