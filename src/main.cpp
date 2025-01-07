#pragma once

#define _USE_MATH_DEFINES

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <math.h>
#include <vector>



void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
std::string readShaderFile(const char* filePath);
unsigned int compileShader(const char* shaderSource, GLenum shaderType);
double generateRandomValue(double min = 0, double max = 1);
std::vector<float> generateRandomPointsOnSphere(int n, float r);
std::vector<float> generateGridPointsOnSphere(int n, float r);


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

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // Read the vertex shader source from the file
    std::string vertexShaderSource = readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\VertexShader.vert");
    if (vertexShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    // Read the fragment shader source from the file
    std::string fragmentShaderSource = readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\FragmentShader.frag");
    if (fragmentShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    // Read the geometry shader source from the file
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
    glAttachShader(shaderProgram, geometryShader);
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
    view = glm::translate(view, glm::vec3(0.0f, 0.0f, -3.0f));

    glm::mat4 projection;
    projection = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.1f, 100.0f);

    unsigned int modelLoc = glGetUniformLocation(shaderProgram, "model");
    unsigned int viewLoc = glGetUniformLocation(shaderProgram, "view");
    unsigned int projectionLoc = glGetUniformLocation(shaderProgram, "projection");

    glViewport(0, 0, 800, 600);

    // Updates glViewport when window is resized
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    float vertices[] = {
        // Positions          // Texture Coords
        -1.f, -1.f, 0.0f,     0.0f, 0.0f, // Bottom-left
         1.f, -1.f, 0.0f,     1.0f, 0.0f, // Bottom-right
        -1.f,  1.f, 0.0f,     0.0f, 1.0f, // Top-left
         1.f, -1.f, 0.0f,     1.0f, 0.0f, // Bottom-right
         1.f,  1.f, 0.0f,     1.0f, 1.0f, // Top-right
        -1.f,  1.f, 0.0f,     0.0f, 1.0f  // Top-left
    };

    int numberOfPoints = 2000;
    float radius = 1;
    std::vector<float> points = generateGridPointsOnSphere(numberOfPoints, radius);
    
    unsigned int VAO;
    glGenVertexArrays(1, &VAO);

    unsigned int VBO;
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    // Copy vertices into VBO
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(float), points.data(), GL_STATIC_DRAW);

    // Position attribute (location = 0)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Texture coordinate attribute (location = 1)
    //glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));

    //glEnableVertexAttribArray(0);
    //glEnableVertexAttribArray(1);

    glEnable(GL_DEPTH_TEST);

    // Rendering loop
    while (!glfwWindowShouldClose(window))
    {
        // Process user input
        processInput(window);

        // Clear the screen
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Pass uniforms
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));

        // Use the shader program
        glUseProgram(shaderProgram);

        glPointSize(5.0f);

        // Bind the VAO and draw the triangle
        glBindVertexArray(VAO);
        glDrawArrays(GL_POINTS, 0, numberOfPoints);

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

void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
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

// Genererar random seed points på sfär, kan användas om voronoi beräknas på cpu
std::vector<float> generateGridPointsOnSphere(int n, float r) {
    std::vector<float> randomPoints;
    randomPoints.reserve(3 * n);

    int pointsPerAxis = sqrt(n); // This will cause the nested loop to generate n points since sqrt(n)*sqrt(n)=n

    for (size_t i = 0; i < pointsPerAxis; i++) // Iterates around sphere
    {
        float stepAz = (float)(i + 1) / (float)(pointsPerAxis + 1);
        float azimuthAngle = 2 * M_PI * stepAz;

        for (size_t j = 0; j < pointsPerAxis; j++) // Iterates up along sphere
        {
            float stepInc = (float)(j + 1) / (float)(pointsPerAxis + 1);

            float inclinationAngle = M_PI * stepInc;
            
            float x = r * sin(inclinationAngle) * cos(azimuthAngle);
            float y = r * sin(inclinationAngle) * sin(azimuthAngle);
            float z = r * cos(inclinationAngle);

            randomPoints.push_back(x);
            randomPoints.push_back(y);
            randomPoints.push_back(z);
        }
    }

    return randomPoints;
}


