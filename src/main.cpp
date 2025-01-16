#pragma once

#define _USE_MATH_DEFINES

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Quickhull for creating convex hull, CURRENTLY NOT USED
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
#include <functional>

#include "Utility.h"


std::vector<float> generateRandomPointsOnSphere(int n, float r);
std::vector<quickhull::Vector3<float>> generateGridPointsOnSphere(int n, float r);

std::vector<unsigned int> generateAndUploadBuffers(unsigned int& VAO, unsigned int& VBO, unsigned int& EBO);
std::vector<std::pair<glm::vec3, glm::vec3>> computeVoronoiEdges(const std::vector<float>& vertices, const std::vector<glm::vec3>& circumcenters, const std::vector<unsigned int>& indices);
std::vector<glm::vec3> computeCircumcenters(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

float viewportWidth = 800.0;
float viewportHeight = 600.0;

struct hash_pair {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        auto hash1 = std::hash<T1>()(pair.first);  // Hash the first element
        auto hash2 = std::hash<T2>()(pair.second); // Hash the second element
        // Combine the two hashes using bitwise XOR and a shift
        return hash1 ^ (hash2 << 1);
    }
};

int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
   
    
    // 
    GLFWwindow* window = glfwCreateWindow(viewportWidth, viewportHeight, "Voronoi Spheres", NULL, NULL);
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
    std::string vertexShaderSource = Utility::readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\VertexShader.vert");
    if (vertexShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    // Read the fragment shader source from the file
    //std::string fragmentShaderSource = readShaderFile("C:\\TNM084 project\\VoronoiSphere\\src\\FragmentShader.frag");
    std::string fragmentShaderSource = Utility::readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\FragmentShader.frag");
    if (fragmentShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    // Read the geometry shader source from the file
    //std::string geometryShaderSource = readShaderFile("C:\\TNM084 project\\VoronoiSphere\\src\\GeometryShader.geom");
    std::string geometryShaderSource = Utility::readShaderFile("C:\\TNM084\\VoronoiSpheres\\src\\GeometryShader.geom");
    if (geometryShaderSource.empty())
    {
        return -1; // Exit if the file couldn't be read
    }

    unsigned int vertexShader = Utility::compileShader(vertexShaderSource.c_str(), GL_VERTEX_SHADER);
    unsigned int fragmentShader = Utility::compileShader(fragmentShaderSource.c_str(), GL_FRAGMENT_SHADER);
    unsigned int geometryShader = Utility::compileShader(geometryShaderSource.c_str(), GL_GEOMETRY_SHADER);

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

    glm::mat4 projection;
    projection = glm::perspective(glm::radians(45.0f), viewportWidth / viewportHeight, 0.1f, 100.0f);

    unsigned int modelLoc = glGetUniformLocation(shaderProgram, "model");
    unsigned int viewLoc = glGetUniformLocation(shaderProgram, "view");
    unsigned int projectionLoc = glGetUniformLocation(shaderProgram, "projection");
    unsigned int numberOfPointsLoc = glGetUniformLocation(shaderProgram, "numberOfPoints");

    glViewport(0, 0, viewportWidth, viewportHeight);

    // Updates glViewport when window is resized
    glfwSetFramebufferSizeCallback(window, Utility::framebuffer_size_callback);

    glfwSetCursorPosCallback(window, Utility::mouse_callback);

    // OpenGL buffers
    unsigned int VAO, VBO, EBO;

    // Initialize OpenGL buffers
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    std::vector<unsigned int> castedIndexBuffer = generateAndUploadBuffers(VAO, VBO, EBO);

    glEnable(GL_DEPTH_TEST);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);



    // Rendering loop
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        // Process user input
        Utility::processInput(window);

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
            castedIndexBuffer = generateAndUploadBuffers(VAO, VBO, EBO); // Regenerate points and upload buffers
        }

        

        // Use the shader program
        glUseProgram(shaderProgram);
        
        glLineWidth(2.0f);
        glPointSize(5.0f);

        // Bind the VAO and draw
        glBindVertexArray(VAO);

        // For drawing triangulation
        if (renderTriangles)
        {
            glDrawElements(GL_TRIANGLES, castedIndexBuffer.size(), GL_UNSIGNED_INT, 0);
        }
        else
        {
            glDrawArrays(GL_LINES, 0, 10000);
        }
        

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



// Genererar random seed points på sfär, kan användas om voronoi beräknas på cpu
std::vector<float> generateRandomPointsOnSphere(int n, float r) {
    std::vector<float> randomPoints;
    randomPoints.reserve(3*n);

    for (size_t i = 0; i < n; i++)
    {
        float randomNumberAz = Utility::generateRandomValue();
        float randomNumberInc = Utility::generateRandomValue();

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
            Utility::generateRandomValue(-1.0, 1.0),
            Utility::generateRandomValue(-1.0, 1.0),
            Utility::generateRandomValue(-1.0, 1.0)
        };
        randomOffset = glm::normalize(randomOffset);

        double offsetMagnitude = Utility::generateRandomValue(0.0, 0.2); // Adjust the range as needed
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

    // Compute voronoi edges
    std::vector<glm::vec3> circumcenters = computeCircumcenters(flattenedVertexBuffer, castedIndexBuffer);
    std::vector<std::pair<glm::vec3, glm::vec3>> voronoiEdges = computeVoronoiEdges(flattenedVertexBuffer, circumcenters, castedIndexBuffer);

    // Flatten Voronoi edges for OpenGL
    std::vector<float> voronoiEdgeVertices;

    for (const auto& edge : voronoiEdges) {
        voronoiEdgeVertices.push_back(edge.first.x);
        voronoiEdgeVertices.push_back(edge.first.y);
        voronoiEdgeVertices.push_back(edge.first.z);

        voronoiEdgeVertices.push_back(edge.second.x);
        voronoiEdgeVertices.push_back(edge.second.y);
        voronoiEdgeVertices.push_back(edge.second.z);
    }

    auto end = std::chrono::high_resolution_clock::now(); // TIMER STOP
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "The generation took " << duration << " milliseconds.\n";

    if (renderTriangles)
    {
        // Update OpenGL buffers
        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, flattenedVertexBuffer.size() * sizeof(float), flattenedVertexBuffer.data(), GL_DYNAMIC_DRAW);

        // For drawing the triangulation
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, castedIndexBuffer.size() * sizeof(unsigned int), castedIndexBuffer.data(), GL_DYNAMIC_DRAW);
        
        // configure vertex attributes
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
    }
    else
    {
        // Update OpenGL buffers
        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, voronoiEdgeVertices.size() * sizeof(float), voronoiEdgeVertices.data(), GL_DYNAMIC_DRAW);

        // configure vertex attributes
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
    }


    return castedIndexBuffer;
}

std::vector<glm::vec3> computeCircumcenters(const std::vector<float>& vertices,
    const std::vector<unsigned int>& indices) {
    std::vector<glm::vec3> circumcenters;

    for (size_t i = 0; i < indices.size(); i += 3) {
        // Get the three vertices of the triangle
        glm::vec3 A{ vertices[indices[i] * 3], vertices[indices[i] * 3 + 1], vertices[indices[i] * 3 + 2] };
        glm::vec3 B{ vertices[indices[i + 1] * 3], vertices[indices[i + 1] * 3 + 1], vertices[indices[i + 1] * 3 + 2] };
        glm::vec3 C{ vertices[indices[i + 2] * 3], vertices[indices[i + 2] * 3 + 1], vertices[indices[i + 2] * 3 + 2] };
        
        glm::vec3 a = A - C;
        glm::vec3 b = B - C;

        glm::vec3 circumcenter = glm::cross((glm::length(a) * glm::length(a) * b - glm::length(b) * glm::length(b) * a), glm::cross(a, b)) / (2.0f * glm::length(glm::cross(a, b)) * glm::length(glm::cross(a, b))) + C;

        // Project the circumcenter onto the sphere surface
        // circumcenter = glm::normalize(circumcenter);

        circumcenters.push_back(circumcenter);
    }

    return circumcenters;
}

std::vector<std::pair<glm::vec3, glm::vec3>> computeVoronoiEdges(const std::vector<float>& vertices,
    const std::vector<glm::vec3>& circumcenters,
    const std::vector<unsigned int>& indices) {

    std::unordered_map<std::pair<int, int>, glm::vec3, hash_pair> edgeToCircumcenter;
    std::vector<std::pair<glm::vec3, glm::vec3>> voronoiEdges;

    for (size_t i = 0; i < indices.size(); i += 3) { // loop through every triangle
        int a = indices[i];
        int b = indices[i + 1];
        int c = indices[i + 2];

        glm::vec3 A{ vertices[a], vertices[a + 1], vertices[a + 2] };
        glm::vec3 B{ vertices[b], vertices[b + 1], vertices[b + 2] };
        glm::vec3 C{ vertices[c], vertices[c + 1], vertices[c + 2] };

        glm::vec3 AB = B - A;
        glm::vec3 AC = C - A;
        glm::vec3 normal = glm::cross(AB, AC);
        if (glm::length(normal) < 0.001) {
            // Skip degenerate triangle
            continue;
        }

        glm::vec3 circumcenter = circumcenters[i / 3];

        if (circumcenter == glm::vec3(0.0))
        {
            std::cout << "Do not create triangle \n";
            continue;
        }

        auto insertEdge = [&](int u, int v) {
            if (u > v) std::swap(u, v);
            auto edge = std::make_pair(u, v);

            if (edgeToCircumcenter.find(edge) != edgeToCircumcenter.end()) { // If already added
                voronoiEdges.emplace_back(edgeToCircumcenter[edge], circumcenter);
                edgeToCircumcenter.erase(edge);
            }
            else {
                edgeToCircumcenter[edge] = circumcenter;
            }
            };

        insertEdge(a, b);
        insertEdge(b, c);
        insertEdge(c, a);
    }

    return voronoiEdges;
}

