#ifndef SHADER_GEOMETRY_H
#define SHADER_GEOMETRY_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>


class ShaderGeometry{
public:
    // build and compiles our shader
    ShaderGeometry(const GLchar* vertexPath, const GLchar* geometryPath, const GLchar* fragmentPath);
    ~ShaderGeometry();
    // use shader
    void Use();
    GLuint ThisProgram();

private:
    GLuint Program;
};

#endif
