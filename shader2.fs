#version 330 core

in vec3 pos;

out vec4 color;

void main()
{
    float r = 0.0f;
    float g = 1.0f;
    float b = 0.0f;
    float a = 1.0f;

    color = vec4(r, g, b, a);
}
