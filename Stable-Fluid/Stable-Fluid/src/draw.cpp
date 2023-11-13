#include "draw.h"
#include "fluid-Euler.h"

#define SCR_WIDTH 1600
#define SCR_HEIGHT 1600

//camera
Camera camera(glm::vec3(0.0f, -5.0f, -3.0f), glm::vec3(0.0f, 1.0f, 0.0f), 90.0f, -30.0f);

//mouse
bool firstMouse = true;
float lastX = (float)SCR_WIDTH / 2.0f;
float lastY = (float)SCR_HEIGHT / 2.0f;

//timing
float deltaTime = 0.0f;
float lastTime = 0.0f;
float deltaTime_tot = 0.0f;

//calculate frame rate
float frameRate = 0;
int count = 0;

//pause
bool isPause = true;
bool keyWasPressed1 = false;
bool keyWasPressed2 = false;
bool keyWasPressed3 = false;

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = ypos - lastY;
    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

void processInput(GLFWwindow* window)
{
    //esc
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    //pause
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && !keyWasPressed1) {
        keyWasPressed1 = true;
        isPause = !isPause;
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_RELEASE && keyWasPressed1)
        keyWasPressed1 = false;

    //change speed
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS && !keyWasPressed2) {
        keyWasPressed2 = true;
        camera.movementSpeed_ += 2;
    }
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_RELEASE && keyWasPressed2)
        keyWasPressed2 = false;
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS && !keyWasPressed3) {
        keyWasPressed3 = true;
        camera.movementSpeed_ = camera.movementSpeed_ > 2 ? camera.movementSpeed_ - 1 : 1;
    }
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_RELEASE && keyWasPressed3)
        keyWasPressed3 = false;

    //control
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
}

//___________________________________________________________________draw()______________________________________________________________________________
template <class T>
int draw(T* obj, float dt, int n) {
    //init GLFW
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        //return -1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, u8"Learn OpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);//set mouse invisible

    glfwSwapInterval(0); // 0表示禁用垂直同步，1表示启用.

    //init GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);

    Shader ourShader("C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/src/shader/vertexShader.txt", "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/src/shader/fragmentShader.txt");

    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glBindVertexArray(VAO);
    //VBO. Buffer vertex array
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * obj->vertices().size(), obj->vertices().data(), GL_STREAM_DRAW);
    //EBO. Buffer index array
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * obj->indices().size(), obj->indices().data(), GL_STREAM_DRAW);
    //Configure vertex attributes 
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    //Unbind
    //glBindBuffer(GL_ARRAY_BUFFER, 0);
    //glBindVertexArray(0);

    //set depth test
    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(window))
    {
        //get deltaTime
        float currentTime = glfwGetTime();
        deltaTime = currentTime - lastTime;
        lastTime = currentTime;

        //print frameRate
        deltaTime_tot += deltaTime;
        frameRate += 1.0 / deltaTime;
        count++;
        if (deltaTime_tot >= 0.5f) {
            //std::cout << "fps: " << frameRate / 60 << std::endl;
            std::stringstream ss;
            ss << "Stable Fluid (" << std::setprecision(4) << frameRate / count << " FPS)";
            std::string title = ss.str();
            glfwSetWindowTitle(window, title.c_str());
            frameRate = 0;
            count = 0;
            deltaTime_tot = 0.0f;
        }

        //update vertices
        if (!isPause) {
            for (int i = 0; i < n; i++) {
                obj->update(dt);
            }
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, obj->vertices().size() * sizeof(float), obj->vertices().data(), GL_STREAM_DRAW);
        }

        //检查并调用事件.
        glfwPollEvents();
        //处理输入事件.
        processInput(window);
        //处理渲染指令.
        glClearColor(0.16f, 0.3f, 0.3f, 0.9f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        //使用着色器程序对象.
        ourShader.use();

        glm::mat4 model = glm::mat4(1.0f);
        glm::mat4 view = glm::mat4(1.0f);
        glm::mat4 projection = glm::mat4(1.0f);
        model = glm::rotate(model, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
        view = glm::lookAt(camera.position_, camera.front_ + camera.position_, camera.up_);
        projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 200.0f);
        //pass them to shaders
        ourShader.setValue("model", model);
        ourShader.setValue("view", view);
        ourShader.setValue("projection", projection);

        glBindVertexArray(VAO);

        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, obj->indices().size(), GL_UNSIGNED_INT, 0);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        //glBindVertexArray(0); just unbind once 
        //交换缓冲区.
        glfwSwapBuffers(window);

    }
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);

    glfwTerminate();
    return 0;
}
