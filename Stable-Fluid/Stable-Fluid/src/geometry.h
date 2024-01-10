#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "tiny_obj_loader.h"
#include "tinyply.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

struct Face {
	int v0;
	int v1;
	int v2;
	Face() {}
	Face(int v0, int v1, int v2) :v0(v0), v1(v1), v2(v2) {}
};

aiVector3D Cross(const aiVector3D& e1, const aiVector3D& e2);
float Dot(const aiVector3D& e1, const aiVector3D& e2);

float computeDis(const aiVector3D& v0, const aiVector3D& v1, const aiVector3D& v2, const aiVector3D& x);
float computeSDF(const aiScene* scene, const aiVector3D& x);

void read_from_OBJ(std::string obj_path, std::vector<glm::vec3>& vertices, std::vector<Face>& faces, float scale = 1.0f, glm::vec3 translate = glm::vec3(0.0f));
void write_to_PLY(std::string ply_path, const std::vector<glm::vec3>& vertices, const std::vector<Face>& faces, float scale = 1.0f);

bool obj_2_SDF(int N1, int N2, int N3, float size, float l_, std::string obj_path, std::vector<float>& phi, bool inverse = false);
bool obj_2_SDF_py(int N1, int N2, int N3, int size, float l, std::string obj_path, std::vector<float>& phi, bool inverse = false);
bool obj_2_SDF_py(int N1, int N2, int N3, float l, std::string obj_path, std::vector<float>& phi, float scale = 1.0f, float translate_x = 0.0f, float translate_y = 0.0f, float translate_z = 0.0f, bool inverse = false);

#endif