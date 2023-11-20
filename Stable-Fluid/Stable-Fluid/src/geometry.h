#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <cmath>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

aiVector3D Cross(const aiVector3D& e1, const aiVector3D& e2);
float Dot(const aiVector3D& e1, const aiVector3D& e2);

float computeDis(const aiVector3D& v0, const aiVector3D& v1, const aiVector3D& v2, const aiVector3D& x);
float computeSDF(const aiScene* scene, const aiVector3D& x);

bool obj_2_SDF(int N1, int N2, int N3, float size, float l_, std::string obj_path, std::vector<float>& phi);

#endif