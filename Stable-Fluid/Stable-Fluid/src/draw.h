#ifndef DRAW_H
#define DRAW_H
#include "shader.h"
#include "camera.h"

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <GLFW/glfw3.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>

void mouse_callback(GLFWwindow* window, double xpos, double ypos);

void processInput(GLFWwindow* window);

template <class T>
int draw(T* obj, float dt, int n = 1);

#endif