#ifndef SHADER_H
#define SHADER_H

#include <glad/glad.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class Shader {
public:
	unsigned int id_;

	Shader(const char* vertexPath, const char* fragmentPath);
	~Shader();
	void use();
	//set uniform
	void setValue(const std::string& name, bool value) const;
	void setValue(const std::string& name, int value) const;
	void setValue(const std::string& name, float value) const;
	void setValue(const std::string& name, glm::mat4 value) const;
};

#endif