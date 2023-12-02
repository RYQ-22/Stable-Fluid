#include "geometry.h"

aiVector3D Cross(const aiVector3D& e1, const aiVector3D& e2) {
	return aiVector3D(e1.y * e2.z - e1.z * e2.y, e1.z * e2.x - e1.x * e2.z, e1.x * e2.y - e1.y * e2.x);
}

float Dot(const aiVector3D& e1, const aiVector3D& e2) {
	return e1.x * e2.x + e1.y * e2.y + e1.z * e2.z;
}

float computeDis(const aiVector3D& v0, const aiVector3D& v1, const aiVector3D& v2, const aiVector3D& x) {
	aiVector3D rc = (v0 + v1 + v2) / 3.0f;
	aiVector3D r = x - v0;
	aiVector3D e1 = v1 - v0;
	aiVector3D e2 = v2 - v1;
	aiVector3D e3 = v0 - v2;
	aiVector3D n = Cross(e1, e2).Normalize();
	aiVector3D xp = x - Dot(r, n) * n;// in triangle plane
	aiVector3D d01 = Cross(e1, xp - v0);
	aiVector3D d12 = Cross(e2, xp - v1);
	aiVector3D d20 = Cross(e3, xp - v2);
	bool isInTriangle = (Dot(d01, n) >= 0 && Dot(d12, n) >= 0 && Dot(d20, n) >= 0) || (Dot(d01, n) <= 0 && Dot(d12, n) <= 0 && Dot(d20, n) <= 0);
	if (isInTriangle)
		return Dot(r, n);
	float dis = (r - rc).Length();
	return dis;
}

float computeSDF(const aiScene* scene, const aiVector3D& x) {
	float phi;
	bool isFirst = true;
	for (int i = 0; i < scene->mNumMeshes; i++) {
		const aiMesh* mesh = scene->mMeshes[i];
		for (int j = 0; j < mesh->mNumFaces; ++j) {
			const aiFace& face = mesh->mFaces[j];
			if (face.mNumIndices == 3) {
				aiVector3D v0 = mesh->mVertices[face.mIndices[0]];
				aiVector3D v1 = mesh->mVertices[face.mIndices[1]];
				aiVector3D v2 = mesh->mVertices[face.mIndices[2]];
				if (isFirst) {
					phi = computeDis(v0, v1, v2, x);
					isFirst = false;
				}
				else {
					if (abs(computeDis(v0, v1, v2, x)) < abs(phi)) phi = computeDis(v0, v1, v2, x);
				}
			}
		}
	}
	//std::cout << phi << std::endl;
	return phi;
}

bool obj_2_SDF(int N1, int N2, int N3, float size, float l, std::string obj_path, std::vector<float>& phi, bool inverse) {
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(obj_path, aiProcess_Triangulate | aiProcess_FlipUVs);
	if (!scene) {
		std::cerr << "Failed to load the model." << std::endl;
		return false;
	}
	else
	{
		phi.resize(N1 * N2 * N3);
		float x_max, x_min;
		bool isFirst = true;
		for (int i = 0; i < scene->mNumMeshes; i++) {
			const aiMesh* mesh = scene->mMeshes[i];
			for (int j = 0; j < mesh->mNumFaces; ++j) {
				const aiFace& face = mesh->mFaces[j];
				if (face.mNumIndices == 3) {
					aiVector3D v0 = mesh->mVertices[face.mIndices[0]];
					aiVector3D v1 = mesh->mVertices[face.mIndices[1]];
					aiVector3D v2 = mesh->mVertices[face.mIndices[2]];
					if (isFirst) {
						x_max = v0.x;
						x_min = v0.x;
						x_max = std::max(x_max, std::max(v1.x, v2.x));
						x_min = std::min(x_min, std::min(v1.x, v2.x));
						isFirst = false;
					}
					else {
						x_max = std::max(x_max, std::max(v0.x, std::max(v1.x, v2.x)));
						x_min = std::min(x_min, std::min(v0.x, std::min(v1.x, v2.x)));
					}
				}
			}
		}
		for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
			aiVector3D x(i + 0.5f, j + 0.5f, k + 0.5f);
			x = x - aiVector3D(N1 / 2, N2 / 2, N3 / 2);// move to the center
			x = x * (x_max - x_min) / size;
			phi[i * N2 * N3 + j * N3 + k] = (size * l) / (x_max - x_min) * computeSDF(scene, x);
			if (inverse) phi[i * N2 * N3 + j * N3 + k] *= -1;
		}
	}
}