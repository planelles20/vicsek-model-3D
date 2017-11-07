
//-----------------------------------------------------------------------
// Main file of Vicsek model simulation using Smoothed-particle hydrodynamics (SPH)
// method.
//
// Licensing: This code is distributed under the Apache License 2.0
// Author: Carlos Planelles Alemany, planelles20(at)gmail(dot)com
//-----------------------------------------------------------------------

// GLM Mathematics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// our class includes
#include "./common/display.h"
#include "./common/shaderGeometry.h"
#include "./common/shader.h"
#include "systemSPH.h"

//particles = THREADS*BLOCKS
#define THREADS 512 // 2^9
#define BLOCKS 64 // 2^15   16384
//mesh
// dR =  min(1/X_MESH, 1/Y_MESH, 1/Z_MESH)
#define X_MESH 100
#define Y_MESH 100
#define Z_MESH 100


// Window dimensions
const GLuint WIDTH = 1800, HEIGHT = 1000;

// Camera
glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f,  3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);
GLfloat yaw    = -90.0f;	// Yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right (due to how Eular angles work) so we initially rotate a bit to the left.
GLfloat pitch  =  0.0f;
GLfloat lastX  =  WIDTH  / 2.0;
GLfloat lastY  =  HEIGHT / 2.0;
GLfloat fov =  -45.0f;
bool keys[1024];

// Deltatime
GLfloat deltaTime = 0.0f;	// Time between current frame and last frame
GLfloat lastFrame = 0.0f;  	// Time of last frame

//openGL functions
void printInstructions(void);
// Function prototypes
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void button_mouse_callback(GLFWwindow* window, int button, int action, int mode);
void do_movement();

int main(void){
    // Init GLFW
    glfwInit();

    // create our window
    Display display = Display(WIDTH, HEIGHT, "Vicsek 3D simulation, Carlos Planelles",
                              key_callback, mouse_callback, scroll_callback, button_mouse_callback);

    // Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
    glewExperimental = GL_TRUE;

    // Initialize GLEW to setup the OpenGL Function pointers
    glewInit();

    // Define the viewport dimensions
    glViewport(0, 0, WIDTH, HEIGHT);

    // enable z-buffer
    glEnable(GL_DEPTH_TEST);

    // Build and compile our shader programs
    ShaderGeometry ourShader1 = ShaderGeometry("shader1.vs", "shader1.gs", "shader1.fs");
    Shader ourShader2 = Shader("shader2.vs", "shader2.fs");

    glEnable(GL_PROGRAM_POINT_SIZE);

    // create SPH sistem
    SystemSPH systemSPH = SystemSPH(BLOCKS, THREADS, X_MESH, Y_MESH, Z_MESH);

    int i = 0;

    while(i<1000000 && display.State()){
        // Calculate deltatime of current frame
        GLfloat currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        //start simulation
        //sistem.SeedUpdate(i); //update seed each iteration
        systemSPH.Calculate();
        //sistem.Particle_print();
        //sistem.Indices_print();
        //sistem.Calc_print();

        /*
        // save current data
        if(i%10000==0){
            sistem.Save("save0.csv");
        }
        */


        display.PollEvents();
        do_movement();

        systemSPH.BackGround(0.0, 0.0, 0.0, 1.0);

        // draw particles
        // use our shader 1
        ourShader1.Use();
        // Camera/View transformation
        glm::mat4 view;
        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
        glm::mat4 model;
        glm::mat4 projection;
        model = glm::rotate(model, 90.0f, glm::vec3(1.0f, 0.0f, 0.0f));
        view = glm::translate(view, glm::vec3(0.0f, 0.0f, -3.0f));
        projection = glm::perspective(45.0f, (GLfloat)WIDTH / (GLfloat)HEIGHT, 0.1f, 100.0f);
        // Get their uniform location shader 1
        GLint modelLoc1 = glGetUniformLocation(ourShader1.ThisProgram(), "model");
        GLint viewLoc1 = glGetUniformLocation(ourShader1.ThisProgram(), "view");
        GLint projLoc1 = glGetUniformLocation(ourShader1.ThisProgram(), "projection");
        GLint camerajLoc1 = glGetUniformLocation(ourShader1.ThisProgram(), "camera");
        // Pass them to the shader 1
        glUniformMatrix4fv(modelLoc1, 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(viewLoc1, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(projLoc1, 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(camerajLoc1, 1, glm::value_ptr(cameraPos));
        //draw particles
        systemSPH.DrawParticles();

        ////// draw boundaries
        // use our shader 2
        ourShader2.Use();
        // Get their uniform location shader 2
        GLint modelLoc2 = glGetUniformLocation(ourShader2.ThisProgram(), "model");
        GLint viewLoc2 = glGetUniformLocation(ourShader2.ThisProgram(), "view");
        GLint projLoc2 = glGetUniformLocation(ourShader2.ThisProgram(), "projection");
        // Pass them to the shader 2
        glUniformMatrix4fv(modelLoc2, 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(viewLoc2, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(projLoc2, 1, GL_FALSE, glm::value_ptr(projection));
        //draw boundary
        systemSPH.DrawBoundary();

        //Swap buffers
        display.SwapWindows();

        i++;
    }

    glfwTerminate();
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
///////////////////////// Function prototypes source /////////////////////////
//////////////////////////////////////////////////////////////////////////////
void printInstructions(void){
    std::cout << "w: move forward" << '\n';
    std::cout << "s: move backwards" << '\n';
    std::cout << "d: move right" << '\n';
    std::cout << "a: move left" << '\n';
    std::cout << "i: pitch up" << '\n';
    std::cout << "k: pitch down" << '\n';
    std::cout << "l: yaw right" << '\n';
    std::cout << "j: yaw left" << '\n';
    std::cout << "use mouse too" << '\n';
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
    if (key >= 0 && key < 1024)
    {
        if (action == GLFW_PRESS)
            keys[key] = true;
        else if (action == GLFW_RELEASE)
            keys[key] = false;
    }
}

void button_mouse_callback(GLFWwindow* window, int button, int action, int mode) {
    if (button >= 0 && button < 1024)
    {
        if (action == GLFW_PRESS)
            keys[button] = true;
        else if (action == GLFW_RELEASE)
            keys[button] = false;
    }
}

void do_movement() {
    // Camera controls
    GLfloat cameraSpeed = 3.0f * deltaTime;
    GLfloat sensitivity = 0.5;
    if (keys[GLFW_KEY_W])
        cameraPos += cameraSpeed * cameraFront;
    if (keys[GLFW_KEY_S])
        cameraPos -= cameraSpeed * cameraFront;
    if (keys[GLFW_KEY_A])
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (keys[GLFW_KEY_D])
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (keys[GLFW_KEY_I]) {
        pitch += sensitivity;

        glm::vec3 front;
        front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraFront = glm::normalize(front);
    }
    if (keys[GLFW_KEY_K]) {
        pitch -= sensitivity;

        glm::vec3 front;
        front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraFront = glm::normalize(front);
    }
    if (keys[GLFW_KEY_L]) {
        yaw += sensitivity;

        glm::vec3 front;
        front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraFront = glm::normalize(front);
    }
    if (keys[GLFW_KEY_J]) {
        yaw -= sensitivity;

        glm::vec3 front;
        front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraFront = glm::normalize(front);
    }
    // print position
    //std::cout << "cameraPos: " << cameraPos.x << "," << cameraPos.y << "," << cameraPos.z << std::endl;
    //std::cout << "cameraFront: " << cameraFront.x << "," << cameraFront.y << "," << cameraFront.z << std::endl;
    //std::cout << "cameraUp: " << cameraUp.x << "," << cameraUp.y << "," << cameraUp.z << std::endl;
    //glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f,  3.0f);
    //glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
    //glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

}

bool firstMouse = true;
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (firstMouse || !keys[GLFW_MOUSE_BUTTON_LEFT]) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    if(keys[GLFW_MOUSE_BUTTON_LEFT]){

        GLfloat xoffset = xpos - lastX;
        GLfloat yoffset = lastY - ypos; // Reversed since y-coordinates go from bottom to left
        lastX = xpos;
        lastY = ypos;

        GLfloat sensitivity = 0.05;	// Change this value to your liking
        xoffset *= sensitivity;
        yoffset *= sensitivity;

        yaw   += xoffset;
        pitch += yoffset;

        // Make sure that when pitch is out of bounds, screen doesn't get flipped
        if (pitch > 89.0f)
            pitch = 89.0f;
        if (pitch < -89.0f)
            pitch = -89.0f;

        glm::vec3 front;
        front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraFront = glm::normalize(front);
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    if (fov >= 1.0f && fov <= 45.0f)
        fov -= yoffset;
    if (fov <= 1.0f)
        fov = 1.0f;
    if (fov >= 45.0f)
        fov = 45.0f;
}
