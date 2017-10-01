#include "display.h"

Display::Display(const unsigned int WIDTH, const unsigned int HEIGHT, const std::string& title,
                 void (*key_callback)(GLFWwindow*,int,int,int,int),
                 void (*mouse_callback)(GLFWwindow*,double,double),
                 void (*scroll_callback)(GLFWwindow*,double,double),
                 void (*button_mouse_callback)(GLFWwindow*,int,int,int)) {
    // Init GLFW
    glfwInit();

    // Set all the required options for GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    // Create a GLFWwindow object that we can use for GLFW's functions
    this->window = glfwCreateWindow(WIDTH, HEIGHT, title.c_str(), nullptr, nullptr);
    glfwMakeContextCurrent(this->window);

    // Set the required callback functions
    glfwSetKeyCallback(this->window, key_callback);
    glfwSetCursorPosCallback(this->window, mouse_callback);
    glfwSetScrollCallback(this->window, scroll_callback);
    glfwSetMouseButtonCallback(this->window, button_mouse_callback);
}

Display::~Display(){
    glfwTerminate();
}

void Display::PollEvents(){
    // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
    glfwPollEvents();
}

void Display::SwapWindows(){
    // Swap the screen buffers
    glfwSwapBuffers(this->window);
}

bool Display::State() {
    return !glfwWindowShouldClose(this->window);
}

void Display::Instructions(){
    std::cout << "Press 'esc' to exit" << std::endl;
}
