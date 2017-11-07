#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 velocity;

out vec3 pass_vel;

void main()
{
    gl_Position = vec4(2.0*position-1.0, 1.0f);
	pass_vel = velocity;
}
