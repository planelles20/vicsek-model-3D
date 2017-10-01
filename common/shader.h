#ifndef SHADER_H
#define SHADER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>


class Shader{
public:
    // build and compiles our shader
    Shader(const GLchar* vertexPath, const GLchar* fragmentPath);
    ~Shader();
    // use shader
    void Use();
    GLuint ThisProgram();

private:
    GLuint Program;
};

#endif
