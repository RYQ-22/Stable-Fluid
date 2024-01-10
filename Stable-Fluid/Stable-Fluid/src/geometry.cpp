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

void read_from_OBJ(std::string obj_path, std::vector<glm::vec3>& vertices, std::vector<Face>& faces, float scale, glm::vec3 translate) {
	std::cout << "Reading " + obj_path << "\n";
	/*
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(obj_path, aiProcess_Triangulate | aiProcess_FlipUVs);
	if (!scene) {
		std::cerr << "Failed to load the model." << std::endl;
		abort();
	}
	else {
		vertices.clear();
		faces.clear();
		for (unsigned int i = 0; i < scene->mNumMeshes; i++) {
			aiMesh* mesh = scene->mMeshes[i];
			// vertices
			for (unsigned int j = 0; j < mesh->mNumVertices; j++) {
				aiVector3D vertex = mesh->mVertices[j];
				vertices.push_back(glm::vec3(vertex.x * scale + translate.x, vertex.y * scale + translate.y, vertex.z * scale + translate.z));
			}
			// faces
			for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
				aiFace face = mesh->mFaces[j];
				if (face.mNumIndices == 3) {
					faces.push_back(Face(face.mIndices[0], face.mIndices[1], face.mIndices[2]));
				}				
			}
		}
	}
	*/
	vertices.clear();
	faces.clear();

	/*
	std::ifstream file(obj_path);
	std::string line;
	if (!file.is_open()) {
		std::cerr << "Unable to open OBJ file: " << obj_path << std::endl;
		return;
	}
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string prefix;
		iss >> prefix;

		if (prefix == "v") {
			glm::vec3 vertex;
			iss >> vertex.x >> vertex.y >> vertex.z;
			vertex = vertex * scale + translate;
			vertices.push_back(vertex);
		}
		else if (prefix == "f") {
			Face face;
			iss >> face.v0 >> face.v1 >> face.v2;		
			face.v0 -= 1;
			face.v1 -= 1;
			face.v2 -= 1;
			faces.push_back(face);
		}
	}
	*/
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string warn;
	std::string err;
	if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, obj_path.c_str())) {
		throw std::runtime_error(err);
	}

	vertices.reserve(attrib.vertices.size() / 3);
	faces.reserve(shapes.size());
	// load vertices
	for (size_t v = 0; v < attrib.vertices.size(); v += 3) {
		vertices.push_back(glm::vec3(attrib.vertices[v], attrib.vertices[v + 1], attrib.vertices[v + 2]) * scale + translate);
	}
	// load faces
	for (size_t s = 0; s < shapes.size(); s++) {
		size_t index_offset = 0;		
		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
			int fv = shapes[s].mesh.num_face_vertices[f];		
			if (fv == 3) { // for triangles
				Face face;
				for (size_t v = 0; v < fv; v++) {
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
					if (v == 0) face.v0 = idx.vertex_index;
					if (v == 1) face.v1 = idx.vertex_index;
					if (v == 2) face.v2 = idx.vertex_index;
				}
				faces.push_back(face);
			}
			index_offset += fv;
		}
	}

	return;
}

void write_to_PLY(std::string ply_path, const std::vector<glm::vec3>& vertices, const std::vector<Face>& faces, float scale) {	
	/*
	std::ofstream file(ply_path);
	// write header
	file << "ply\n";
	file << "format binary_little_endian 1.0\n";
	file << "element vertex " << vertices.size() << "\n";
	file << "property float x\n";
	file << "property float y\n";
	file << "property float z\n";
	file << "element face " << faces.size() << "\n";
	file << "property list uchar int vertex_indices\n";
	file << "end_header\n";

	glm::vec3 vertex_;
	for (const auto& vertex : vertices) {
		vertex_ = vertex * scale;
		file.write((char*)&vertex_.x, sizeof(vertex_.x));
		file.write((char*)&vertex_.y, sizeof(vertex_.y));
		file.write((char*)&vertex_.z, sizeof(vertex_.z));
	}

	unsigned char n = 3;
	for (const auto& face : faces) {
		file.write((char*)&n, sizeof(n));
		file.write((char*)&face.v0, sizeof(face.v0));		
		file.write((char*)&face.v1, sizeof(face.v1));
		file.write((char*)&face.v2, sizeof(face.v2));
	}
	file.close();
	return;
	*/
	tinyply::PlyFile plyFile;
	std::vector<float> vertexData;
	for (int i = 0; i < vertices.size(); i++) {
		vertexData.push_back(vertices[i].x * scale);
		vertexData.push_back(vertices[i].y * scale);
		vertexData.push_back(vertices[i].z * scale);
	}
	plyFile.add_properties_to_element("vertex", { "x", "y", "z" }, tinyply::Type::FLOAT32, vertices.size(), reinterpret_cast<uint8_t*>(vertexData.data()), tinyply::Type::INVALID, 0);
	std::vector<int32_t> faceData;
	faceData.reserve(faces.size() * 3);
	for (int i = 0; i < faces.size(); i++) {		
		faceData.push_back(faces[i].v0);
		faceData.push_back(faces[i].v1);
		faceData.push_back(faces[i].v2);
	}
	plyFile.add_properties_to_element("face", { "vertex_indices" }, tinyply::Type::INT32, faces.size(), reinterpret_cast<uint8_t*>(faceData.data()), tinyply::Type::UINT8, 3);
	std::filebuf fb;
	fb.open(ply_path, std::ios::out | std::ios::binary);
	std::ostream outputStream(&fb);
	if (outputStream.fail())
		throw std::runtime_error("failed to open " + ply_path);
	plyFile.write(outputStream, true);

	return;
}


bool obj_2_SDF(int N1, int N2, int N3, float size, float l, std::string obj_path, std::vector<float>& phi, bool inverse) {
	std::cout << "Loading " + obj_path << "\n";
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
	std::cout << "Model is loaded.\n";
	return true;
}

bool obj_2_SDF_py(int N1, int N2, int N3, int size, float l, std::string obj_path, std::vector<float>& phi, bool inverse) {
	std::cout << "Loading " + obj_path << "\n";
	std::string file_path = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/sdf/temp.txt";
	// call python script	
	std::string command = "D:/Anaconda/envs/myenv/python.exe C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/obj2sdf.py"
		+ (" -N1 " + std::to_string(N1))
		+ " -N2 " + std::to_string(N2)
		+ " -N3 " + std::to_string(N3)
		+ " -size " + std::to_string(size)
		+ " -voxel_size " + std::to_string(l)
		+ " -input_file " + obj_path
		+ " -output_file " + file_path;
	system(command.c_str());
	// read sdf data
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_path << std::endl;
		return false;
	}
	file >> N1 >> N2 >> N3 >> l;
	phi.resize(N1 * N2 * N3);
	float value;
	for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
		if (!(file >> value)) {
			abort();
		}
		phi[i * N2 * N3 + j * N3 + k] = value;		
		if (inverse) phi[i * N2 * N3 + j * N3 + k] *= -1;
	}
	file.close();
	std::cout << "Model is loaded.\n";
	// delet temp.txt	
	try {
		if (std::filesystem::exists(file_path)) { // check whether file exists
			std::filesystem::remove(file_path); // delete file
		}
		else {
			std::cout << "File does not exist." << std::endl;
		}
	}
	catch (const std::filesystem::filesystem_error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}

	return true;
}

bool obj_2_SDF_py(int N1, int N2, int N3, float l, std::string obj_path, std::vector<float>& phi, float scale, float translate_x, float translate_y, float translate_z, bool inverse) {
	std::cout << "Loading " + obj_path << "\n";
	std::string file_path = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/sdf/temp.txt";
	// call python script
	std::string command = "D:/Anaconda/envs/myenv/python.exe C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/obj2sdf.py"
		+ (" -N1 " + std::to_string(N1))
		+ " -N2 " + std::to_string(N2)
		+ " -N3 " + std::to_string(N3)
		+ " -scale " + std::to_string(scale)
		+ " -translate_x " + std::to_string(translate_x)
		+ " -translate_y " + std::to_string(translate_y)
		+ " -translate_z " + std::to_string(translate_z)
		+ " -voxel_size " + std::to_string(l)
		+ " -input_file " + obj_path
		+ " -output_file " + file_path;
	system(command.c_str());
	// read sdf data
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_path << std::endl;
		return false;
	}
	file >> N1 >> N2 >> N3 >> l;
	phi.resize(N1 * N2 * N3);
	float value;
	for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
		if (!(file >> value)) {
			abort();
		}
		phi[i * N2 * N3 + j * N3 + k] = value;
		if (inverse) phi[i * N2 * N3 + j * N3 + k] *= -1;
	}
	file.close();
	std::cout << "Model is loaded.\n";
	// delet temp.txt	
	try {
		if (std::filesystem::exists(file_path)) { // check whether file exists
			std::filesystem::remove(file_path); // delete file
		}
		else {
			std::cout << "File does not exist." << std::endl;
		}
	}
	catch (const std::filesystem::filesystem_error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}

	return true;
}
