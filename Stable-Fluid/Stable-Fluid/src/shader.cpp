#include "shader.h"

Shader::Shader(const char* vertexPath, const char* fragmentPath) {
	//1. Read from files
	std::string vertexCode;
	std::string fragmentCode;
	std::ifstream vShaderFile;
	std::ifstream fShaderFile;
	const char* vShaderCode;//Address of ShaderCode
	const char* fShaderCode;
	//avoid error
	vShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try {
		//Open files
		vShaderFile.open(vertexPath);
		fShaderFile.open(fragmentPath);
		std::stringstream vShaderStream, fShaderStream;
		//Read buffer to stream
		vShaderStream << vShaderFile.rdbuf();
		fShaderStream << fShaderFile.rdbuf();
		vShaderFile.close();
		fShaderFile.close();
		//Translate stream to string
		vertexCode = vShaderStream.str();
		fragmentCode = fShaderStream.str();
	}
	catch (std::ifstream::failure error) {
		std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
		std::cout << vertexCode.length() << std::endl;
		std::cout << fragmentCode.length() << std::endl;
		std::cout << "Vertex Shader Path: " << vertexPath << std::endl;
		std::cout << "Fragment Shader Path: " << fragmentPath << std::endl;

	}
	vShaderCode = vertexCode.c_str();
	fShaderCode = fragmentCode.c_str();

	//2. Complile shader
	unsigned int vertex, fragment;
	int success;
	char infoLog[512];

	//vertexShader
	vertex = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex, 1, &vShaderCode, NULL);
	glCompileShader(vertex);
	//if error
	glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(vertex, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
	}
	//fragmentShader
	fragment = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragment, 1, &fShaderCode, NULL);
	glCompileShader(fragment);
	//if error
	glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(fragment, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
	}

	//ShaderProgram
	id_ = glCreateProgram();
	glAttachShader(id_, vertex);
	glAttachShader(id_, fragment);
	glLinkProgram(id_);
	//if error
	glGetProgramiv(id_, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(id_, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
	}

	//Delete shader
	glDeleteShader(vertex);
	glDeleteShader(fragment);
}

void Shader::use() {
	glUseProgram(id_);
}

void Shader::setValue(const std::string& name, bool value) const {
	glUniform1i(glGetUniformLocation(id_, name.c_str()), (int)value);
}

void Shader::setValue(const std::string& name, int value) const {
	glUniform1i(glGetUniformLocation(id_, name.c_str()), value);
}

void Shader::setValue(const std::string& name, float value) const {
	glUniform1f(glGetUniformLocation(id_, name.c_str()), value);
}

void Shader::setValue(const std::string& name, glm::mat4 value) const {
	glUniformMatrix4fv(glGetUniformLocation(id_, name.c_str()), 1, GL_FALSE, &value[0][0]);
}

Shader::~Shader() {
	glDeleteProgram(id_);
}