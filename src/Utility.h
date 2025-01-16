#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <glad/glad.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

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

bool renderTriangles = false;
bool renderPoints = false;

std::vector<float> triangleStripVertices{};

namespace Utility {

    void exportToOBJ(const std::vector<float>& meshData, const std::string& filename) {
        std::ofstream outFile(filename);

        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
            return;
        }

        // Track indices
        int vertexIndex = 1; // OBJ indices start from 1

        // Step 1: Write vertices and normals
        for (size_t i = 0; i < meshData.size(); i += 6) {
            // Vertex positions
            outFile << "v " << meshData[i] << " " << meshData[i + 1] << " " << meshData[i + 2] << "\n";

            // Normals
            outFile << "vn " << meshData[i + 3] << " " << meshData[i + 4] << " " << meshData[i + 5] << "\n";
        }

        // Step 2: Write faces
        for (size_t i = 0; i < meshData.size() / 6; i += 3) {
            // Faces (referencing vertex and normal indices)
            outFile << "f "
                << vertexIndex << "//" << vertexIndex << " "
                << vertexIndex + 1 << "//" << vertexIndex + 1 << " "
                << vertexIndex + 2 << "//" << vertexIndex + 2 << "\n";

            // Increment vertex index by 3 (one triangle = 3 vertices)
            vertexIndex += 3;
        }

        outFile.close();
        std::cout << "Mesh exported to " << filename << std::endl;
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
        static float lastTime = 0.0f; // Last time the E key was processed
        float currentTime = glfwGetTime();

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

        const float keyDelay = 0.01f; // Delay in seconds
        const float keyDelayExport = 3.0f; // Delay in seconds

        if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS && (currentTime - lastTime > keyDelay))
        {
            numberOfPoints += 2;
            lastTime = currentTime; // Update the last processed time
        }

        if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS && (currentTime - lastTime > keyDelay))
        {
            numberOfPoints -= 2;
            lastTime = currentTime; // Update the last processed time
        }

        if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS && (currentTime - lastTime > keyDelay))
        {
            renderTriangles = !renderTriangles;
            lastTime = currentTime; // Update the last processed time
        }

        if (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS && (currentTime - lastTime > keyDelayExport))
        {
            exportToOBJ(triangleStripVertices, "mesh.obj");
            lastTime = currentTime; // Update the last processed time
        }
    }

    double generateRandomValue(double min = 0.0, double max = 1.0) {
        static std::random_device rd;   // Obtain a random seed
        static std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(min, max);

        return dis(rd);
    }
};
