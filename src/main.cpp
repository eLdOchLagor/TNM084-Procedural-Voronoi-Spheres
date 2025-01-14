#pragma once

#define _USE_MATH_DEFINES

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <quickhull/QuickHull.hpp>
#include <quickhull/Structs/Vector3.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <math.h>
#include <vector>
#include <chrono>



void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
std::string readShaderFile(const char* filePath);
unsigned int compileShader(const char* shaderSource, GLenum shaderType);
double generateRandomValue(double min = 0, double max = 1);
std::vector<float> generateRandomPointsOnSphere(int n, float r);
std::vector<quickhull::Vector3<float>> generateGridPointsOnSphere(int n, float r);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
std::vector<unsigned int> generateAndUploadBuffers(unsigned int& VAO, unsigned int& VBO, unsigned int& EBO);

glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

float deltaTime = 0.0f;	// Time between current frame and last frame
float lastFrame = 0.0f; // Time of last frame

bool firstMouse = true;
float lastX = 400.f, lastY = 300.f, yaw = -90.f, pitch = 0.f;

int numberOfPoints = 100;
int previousNumberOfPoints = numberOfPoints;
float radius = 1;

int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
   
    
    // 
    GLFWwindow* window = glfwCreateWindow(800, 600, "Voronoi Spheres", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // Read the vertex shader source from the file
    //std::string vertexShaderSource = readShaderFile("C:\\TNM084 project\\VoronoiSphere\\src\\VertexShader.vert");
    std::string vertexShaderSource = readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\VertexShader.vert");
    if (vertexShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    // Read the fragment shader source from the file
    //std::string fragmentShaderSource = readShaderFile("C:\\TNM084 project\\VoronoiSphere\\src\\FragmentShader.frag");
    std::string fragmentShaderSource = readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\FragmentShader.frag");
    if (fragmentShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    // Read the geometry shader source from the file
    //std::string geometryShaderSource = readShaderFile("C:\\TNM084 project\\VoronoiSphere\\src\\GeometryShader.geom");
    std::string geometryShaderSource = readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\GeometryShader.geom");
    if (geometryShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    unsigned int vertexShader = compileShader(vertexShaderSource.c_str(), GL_VERTEX_SHADER);
    unsigned int fragmentShader = compileShader(fragmentShaderSource.c_str(), GL_FRAGMENT_SHADER);
    unsigned int geometryShader = compileShader(geometryShaderSource.c_str(), GL_GEOMETRY_SHADER);

    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    //glAttachShader(shaderProgram, geometryShader);
    glLinkProgram(shaderProgram);

    // Check for linking errors
    int success;
    char infoLog[512];
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cerr << "ERROR::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }

    // Delete the shaders as it is no longer needed
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glDeleteShader(geometryShader);

    glm::mat4 model = glm::mat4(1.0f);
    model = glm::rotate(model, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));

    glm::mat4 view = glm::mat4(1.0f);

    glm::mat4 projection;
    projection = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.1f, 100.0f);

    unsigned int modelLoc = glGetUniformLocation(shaderProgram, "model");
    unsigned int viewLoc = glGetUniformLocation(shaderProgram, "view");
    unsigned int projectionLoc = glGetUniformLocation(shaderProgram, "projection");
    unsigned int numberOfPointsLoc = glGetUniformLocation(shaderProgram, "numberOfPoints");

    glViewport(0, 0, 800, 600);

    // Updates glViewport when window is resized
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    glfwSetCursorPosCallback(window, mouse_callback);

    // OpenGL buffers
    unsigned int VAO, VBO, EBO;

    // Initialize OpenGL buffers
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    std::vector<unsigned int> castedIndexBuffer = generateAndUploadBuffers(VAO, VBO, EBO);

    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);



    // Rendering loop
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        // Process user input
        processInput(window);

        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

        // Clear the screen
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Pass uniforms
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
        glUniform1i(numberOfPointsLoc, numberOfPoints);

        //numberOfPoints = abs(sin(currentFrame)) * 100;
        
        // Detect change in numberOfPoints
        if (numberOfPoints != previousNumberOfPoints) {
            previousNumberOfPoints = numberOfPoints;
            generateAndUploadBuffers(VAO, VBO, EBO); // Regenerate points and upload buffers
        }

        

        // Use the shader program
        glUseProgram(shaderProgram);
        
        glLineWidth(2.0f);
        glPointSize(5.0f);

        // Bind the VAO and draw
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, castedIndexBuffer.size(), GL_UNSIGNED_INT, 0);

        // Unbind VAO
        glBindVertexArray(0);

        // Swap the buffers
        glfwSwapBuffers(window);

        // Poll for events
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.05f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;

    glm::vec3 direction;
    direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    direction.y = sin(glm::radians(pitch));
    direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(direction);
}

void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    const float cameraSpeed = 2.5 * deltaTime;
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos += cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos -= cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
}

std::string readShaderFile(const char* filePath)
{
    std::ifstream shaderFile;
    std::stringstream shaderStream;

    // Open the file
    shaderFile.open(filePath);
    if (!shaderFile.is_open())
    {
        std::cerr << "Failed to open shader file: " << filePath << std::endl;
        return "";
    }

    // Read the file into a string stream
    shaderStream << shaderFile.rdbuf();

    // Close the file
    shaderFile.close();

    // Convert the stream into a string and return
    return shaderStream.str();
}

unsigned int compileShader(const char* shaderSource, GLenum shaderType)
{
    // Create the shader
    unsigned int shader = glCreateShader(shaderType);

    // Attach the source code and compile
    glShaderSource(shader, 1, &shaderSource, NULL);
    glCompileShader(shader);

    // Check for compilation errors
    int success;
    char infoLog[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(shader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    return shader;
}

double generateRandomValue(double min, double max) {
    static std::random_device rd;   // Obtain a random seed
    static std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(min, max);

    return dis(rd);
}

// Genererar random seed points på sfär, kan användas om voronoi beräknas på cpu
std::vector<float> generateRandomPointsOnSphere(int n, float r) {
    std::vector<float> randomPoints;
    randomPoints.reserve(3*n);

    for (size_t i = 0; i < n; i++)
    {
        float randomNumberAz = generateRandomValue();
        float randomNumberInc = generateRandomValue();

        float inclinationAngle = acos(2 * randomNumberInc - 1);
        float azimuthAngle = 2 * M_PI * randomNumberAz;

        float x = r * sin(inclinationAngle) * cos(azimuthAngle);
        float y = r * sin(inclinationAngle) * sin(azimuthAngle);
        float z = r * cos(inclinationAngle);

        randomPoints.push_back(x);
        randomPoints.push_back(y);
        randomPoints.push_back(z);
    }

    return randomPoints;
}

std::vector<quickhull::Vector3<float>> generateGridPointsOnSphere(int n, float r) {
    std::vector<quickhull::Vector3<float>> points;

    float goldenRatio = (1.0f + sqrt(5.0f)) / 2.0f; // Golden ratio

    for (int i = 0; i < n; i++) {
        float z = 1.0f - (2.0f * i) / (n - 1); // z goes from 1 to -1
        float radiusAtZ = sqrt(1.0f - z * z); // Radius at this z

        float theta = 2.0f * M_PI * i / goldenRatio; // Golden angle increment

        float x = r * radiusAtZ * cos(theta);
        float y = r * radiusAtZ * sin(theta);
        
        glm::vec3 pos{ x, y, z };
        
        // Generate a random direction vector for offset
        glm::vec3 randomOffset{
            generateRandomValue(-1.0, 1.0),
            generateRandomValue(-1.0, 1.0),
            generateRandomValue(-1.0, 1.0)
        };
        randomOffset = glm::normalize(randomOffset);

        double offsetMagnitude = generateRandomValue(0.0, 0.1); // Adjust the range as needed
        randomOffset *= offsetMagnitude;

        pos += randomOffset;
        pos = glm::normalize(pos);

        points.push_back(quickhull::Vector3<float>{pos.x, pos.y, pos.z});
    }

    return points;
}

std::vector<unsigned int> generateAndUploadBuffers(unsigned int& VAO, unsigned int& VBO, unsigned int& EBO) {
    auto start = std::chrono::high_resolution_clock::now(); // TIMER START

    // Generate grid points
    std::vector<quickhull::Vector3<float>> points = generateGridPointsOnSphere(numberOfPoints, radius);

    quickhull::QuickHull<float> qh;
    auto hull = qh.getConvexHull(points, true, false);
    const auto& indexBuffer = hull.getIndexBuffer();
    const auto& vertexBuffer = hull.getVertexBuffer();

    // Flatten the vertex buffer
    std::vector<float> flattenedVertexBuffer;
    flattenedVertexBuffer.reserve(vertexBuffer.size() * 3);
    for (const auto& vertex : vertexBuffer) {
        flattenedVertexBuffer.push_back(vertex.x);
        flattenedVertexBuffer.push_back(vertex.y);
        flattenedVertexBuffer.push_back(vertex.z);
    }

    std::vector<unsigned int> castedIndexBuffer(indexBuffer.begin(), indexBuffer.end());

    auto end = std::chrono::high_resolution_clock::now(); // TIMER STOP
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "The generation took " << duration << " milliseconds.\n";

    // Update OpenGL buffers
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, flattenedVertexBuffer.size() * sizeof(float), flattenedVertexBuffer.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, castedIndexBuffer.size() * sizeof(unsigned int), castedIndexBuffer.data(), GL_STATIC_DRAW);

    // Reconfigure vertex attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

    return castedIndexBuffer;
}


