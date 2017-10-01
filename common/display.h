#ifndef DISPLAY_H
#define DISPLAY_H

#include <iostream>
#include <string>

#include <GLFW/glfw3.h>

// GLM Mathematics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class Display {
    public:
        //build our window
        Display(const unsigned int WIDTH, const unsigned int HEIGHT, const std::string& title,
                void (*key_callback)(GLFWwindow*,int,int,int,int) = NULL,
                void (*mouse_callback)(GLFWwindow*,double,double) = NULL,
                void (*scroll_callback)(GLFWwindow*,double, double) = NULL,
                void (*button_mouse_callback)(GLFWwindow*,int, int, int) = NULL);

        virtual ~Display();
        void PollEvents();
        void SwapWindows();
        bool State();
        void Instructions();

    private:
        GLFWwindow* window;
};

#endif // DISPLAY_H
