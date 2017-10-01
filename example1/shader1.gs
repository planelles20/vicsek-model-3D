#version 330 core

layout (points) in;
//layout (triangle_strip, max_vertices = 24) out; //cube
layout (triangle_strip, max_vertices = 60) out; //arrow

in vec3 pass_vel[];

out vec3 finalVel;
out vec3 finalColour;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

// change for cube size
const float size = 0.001;
// light
const vec3 lightDirection = normalize(vec3(0.4, -1.0, 0.8));
const float ambientLighting = 0.3;

mat3 Rx(const float radians){
  return mat3(1.0,          0.0,         0.0,
              0.0, cos(radians),sin(radians),
              0.0,-sin(radians),cos(radians));
}

mat3 Ry(const float radians){
  return mat3(cos(radians), 0.0, -sin(radians),
                       0.0, 1.0,           0.0,
              sin(radians), 0.0,  cos(radians));
}

mat3 Rz(const float radians){
  return mat3(cos(radians), sin(radians), 0.0,
             -sin(radians), cos(radians), 0.0,
                       0.0,          0.0, 1.0);
}


void createVertex(vec3 offset,  vec3 colour) {

    float yaw = atan(pass_vel[0].y, pass_vel[0].x);
    float pitch = atan(pass_vel[0].z/sqrt(pass_vel[0].x*pass_vel[0].x+pass_vel[0].y*pass_vel[0].y));
    float roll = 0.0;

    offset = Rz(yaw)*Ry(-pitch)*Rx(roll)*offset;
    vec4 actualOffset = vec4(offset * size, 0.0);
    vec4 worldPosition = gl_in[0].gl_Position + actualOffset;
    gl_Position = projection * view * model * worldPosition;
    finalVel = pass_vel[0];
    finalColour = colour;
    EmitVertex();
}

vec3 calculateLighting(vec3 faceNormal){
	float brightness = max(dot(-lightDirection, faceNormal), ambientLighting);
	return abs(pass_vel[0]) * brightness;
}

// cube geometry
void cube();
void arrow();

void main() {
    //cube();
    arrow();
}


void cube(){

    	vec3 faceNormal = vec3(0.0, 0.0, 1.0);
    	vec3 colour = calculateLighting(faceNormal);
    	createVertex(vec3(-1.0, 1.0, 1.0), colour);
    	createVertex(vec3(-1.0, -1.0, 1.0), colour);
    	createVertex(vec3(1.0, 1.0, 1.0), colour);
    	createVertex(vec3(1.0, -1.0, 1.0), colour);
    	EndPrimitive();

    	faceNormal = vec3(1.0, 0.0, 0.0);
    	colour = calculateLighting(faceNormal);
    	createVertex(vec3(1.0, 1.0, 1.0), colour);
    	createVertex(vec3(1.0, -1.0, 1.0), colour);
    	createVertex(vec3(1.0, 1.0, -1.0), colour);
    	createVertex(vec3(1.0, -1.0, -1.0), colour);
    	EndPrimitive();

    	faceNormal = vec3(0.0, 0.0, -1.0);
    	colour = calculateLighting(faceNormal);
    	createVertex(vec3(1.0, 1.0, -1.0), colour);
    	createVertex(vec3(1.0, -1.0, -1.0), colour);
    	createVertex(vec3(-1.0, 1.0, -1.0), colour);
    	createVertex(vec3(-1.0, -1.0, -1.0), colour);
    	EndPrimitive();

    	faceNormal = vec3(-1.0, 0.0, 0.0);
    	colour = calculateLighting(faceNormal);
    	createVertex(vec3(-1.0, 1.0, -1.0), colour);
    	createVertex(vec3(-1.0, -1.0, -1.0), colour);
    	createVertex(vec3(-1.0, 1.0, 1.0), colour);
    	createVertex(vec3(-1.0, -1.0, 1.0), colour);
    	EndPrimitive();

    	faceNormal = vec3(0.0, 1.0, 0.0);
    	colour = calculateLighting(faceNormal);
    	createVertex(vec3(1.0, 1.0, 1.0), colour);
    	createVertex(vec3(1.0, 1.0, -1.0), colour);
    	createVertex(vec3(-1.0, 1.0, 1.0), colour);
    	createVertex(vec3(-1.0, 1.0, -1.0), colour);
    	EndPrimitive();

    	faceNormal = vec3(0.0, -1.0, 0.0);
    	colour = calculateLighting(faceNormal);
    	createVertex(vec3(-1.0, -1.0, 1.0), colour);
    	createVertex(vec3(-1.0, -1.0, -1.0), colour);
    	createVertex(vec3(1.0, -1.0, 1.0), colour);
    	createVertex(vec3(1.0, -1.0, -1.0), colour);
    	EndPrimitive();
}

void arrow(){
    //NOTE: review faceNormal vectors
    //Heat
    vec3 faceNormal = vec3(-1.0, 0.0, 0.0);
    vec3 colour = calculateLighting(faceNormal);
    createVertex(vec3(0.0,  -1.0,  1.0), colour);
    createVertex(vec3(0.0, -1.0,  -1.0), colour);
    createVertex(vec3(0.0,  1.0, 1.0), colour);
    createVertex(vec3(0.0, 1.0, -1.0), colour);
    EndPrimitive();

    faceNormal = normalize(cross(vec3(-2.0, 1.0, 1.0), vec3(-2.0, -1.0, 1.0)));
    colour = calculateLighting(faceNormal);
    createVertex(vec3(0.0,  1.0,  1.0), colour);
    createVertex(vec3(0.0,  -1.0,  1.0), colour);
    createVertex(vec3(2.0,  0.0,  0.0), colour);
    EndPrimitive();

    faceNormal = normalize(cross(vec3(-2.0, 1.0, -1.0), vec3(-2.0, -1.0, -1.0)));
    colour = calculateLighting(faceNormal);
    createVertex(vec3(0.0,  1.0, -1.0), colour);
    createVertex(vec3(0.0, -1.0, -1.0), colour);
    createVertex(vec3(2.0,  0.0,  0.0), colour);
    EndPrimitive();

    faceNormal = normalize(cross(vec3(-2.0,  1.0, 1.0), vec3(-2.0, 1.0, -1.0)));
    colour = calculateLighting(faceNormal);
    createVertex(vec3(0.0,  1.0,  1.0), colour);
    createVertex(vec3(0.0,  1.0, -1.0), colour);
    createVertex(vec3(2.0,  0.0,  0.0), colour);
    EndPrimitive();

    faceNormal = normalize(cross(vec3(-2.0, -1.0, 1.0), vec3(-2.0, -1.0, -1.0)));
    colour = calculateLighting(faceNormal);
    createVertex(vec3(0.0, -1.0,  1.0), colour);
    createVertex(vec3(0.0, -1.0, -1.0), colour);
    createVertex(vec3(2.0,  0.0,  0.0), colour);
    EndPrimitive();

    //body
    faceNormal = vec3(-1.0, 0.0, 0.0);
    colour = calculateLighting(faceNormal);
    createVertex(vec3(-4.0, -0.5,  0.5), colour);
    createVertex(vec3(-4.0, -0.5, -0.5), colour);
    createVertex(vec3(-4.0,  0.5,  0.5), colour);
    createVertex(vec3(-4.0,  0.5, -0.5), colour);
    EndPrimitive();

    faceNormal = vec3(0.0, 1.0, 0.0);
    colour = calculateLighting(faceNormal);
    createVertex(vec3( 0.0,  0.5, -0.5), colour);
    createVertex(vec3( 0.0,  0.5,  0.5), colour);
    createVertex(vec3(-4.0,  0.5, -0.5), colour);
    createVertex(vec3(-4.0,  0.5,  0.5), colour);
    EndPrimitive();

    faceNormal = vec3(0.0, 0.0, 1.0);
    colour = calculateLighting(faceNormal);
    createVertex(vec3( 0.0,  0.5,  0.5), colour);
    createVertex(vec3( 0.0, -0.5,  0.5), colour);
    createVertex(vec3(-4.0,  0.5,  0.5), colour);
    createVertex(vec3(-4.0, -0.5,  0.5), colour);
    EndPrimitive();

    faceNormal = vec3(0.0, -1.0, 0.0);
    colour = calculateLighting(faceNormal);
    createVertex(vec3( 0.0,  -0.5,  0.5), colour);
    createVertex(vec3( 0.0,  -0.5, -0.5), colour);
    createVertex(vec3(-4.0,  -0.5,  0.5), colour);
    createVertex(vec3(-4.0,  -0.5, -0.5), colour);
    EndPrimitive();

    faceNormal = vec3(0.0, 0.0, -1.0);
    colour = calculateLighting(faceNormal);
    createVertex(vec3( 0.0,  0.5, -0.5), colour);
    createVertex(vec3( 0.0, -0.5, -0.5), colour);
    createVertex(vec3(-4.0,  0.5, -0.5), colour);
    createVertex(vec3(-4.0, -0.5, -0.5), colour);
    EndPrimitive();

}
