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

// CGAL for computing delaunay triangulation
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Point_3                                Point_3;
typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;


void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
std::string readShaderFile(const char* filePath);
unsigned int compileShader(const char* shaderSource, GLenum shaderType);
double generateRandomValue(double min = 0, double max = 1);
std::vector<float> generateRandomPointsOnSphere(int n, float r);
std::vector<Point_3> generateGridPointsOnSphere(int n, float r);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
std::vector<unsigned int> generateAndUploadBuffers(unsigned int& VAO, unsigned int& VBO, unsigned int& EBO);
std::vector<std::pair<glm::vec3, glm::vec3>> computeVoronoiEdges(const std::vector<float>& vertices, const std::vector<glm::vec3>& circumcenters, const std::vector<unsigned int>& indices);
std::vector<glm::vec3> computeCircumcenters(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

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

        // For drawing triangulation
        //glDrawElements(GL_TRIANGLES, castedIndexBuffer.size(), GL_UNSIGNED_INT, 0);
        glDrawArrays(GL_LINES, 0, 1764);

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

    const float keyDelay = 0.2f; // Delay in seconds
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS && (currentTime - lastTime > keyDelay))
    {
        numberOfPoints += 2;
        lastTime = currentTime; // Update the last processed time
    }
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

std::vector<Point_3> generateGridPointsOnSphere(int n, float r) {
    std::vector<Point_3> points;

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

        double offsetMagnitude = generateRandomValue(0.0, 0.0); // Adjust the range as needed
        randomOffset *= offsetMagnitude;

        pos += randomOffset;
        pos = glm::normalize(pos);

        points.push_back(Point_3{pos.x, pos.y, pos.z});
    }

    return points;
}

std::vector<unsigned int> generateAndUploadBuffers(unsigned int& VAO, unsigned int& VBO, unsigned int& EBO) {
    auto start = std::chrono::high_resolution_clock::now(); // TIMER START


    // Generate grid points
    std::vector<Point_3> points = generateGridPointsOnSphere(numberOfPoints, radius);

    // define polyhedron to hold convex hull
    Polyhedron_3 poly;

    // compute convex hull of non-collinear points
    CGAL::convex_hull_3(points.begin(), points.end(), poly);

    std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;

    auto end = std::chrono::high_resolution_clock::now(); // TIMER STOP
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "The generation took " << duration << " milliseconds.\n";

    return std::vector<unsigned int>{};
}

std::vector<glm::vec3> computeCircumcenters(const std::vector<float>& vertices,
    const std::vector<unsigned int>& indices) {
    std::vector<glm::vec3> circumcenters;

    for (size_t i = 0; i < indices.size(); i += 3) {
        // Get the three vertices of the triangle
        glm::vec3 A{ vertices[indices[i] * 3], vertices[indices[i] * 3 + 1], vertices[indices[i] * 3 + 2] };
        glm::vec3 B{ vertices[indices[i + 1] * 3], vertices[indices[i + 1] * 3 + 1], vertices[indices[i + 1] * 3 + 2] };
        glm::vec3 C{ vertices[indices[i + 2] * 3], vertices[indices[i + 2] * 3 + 1], vertices[indices[i + 2] * 3 + 2] };

        // Calculate the vectors for two triangle edges
        glm::vec3 AB = B - A;
        glm::vec3 AC = C - A;

        // Calculate the normal of the plane formed by the triangle
        glm::vec3 normal = glm::normalize(glm::cross(AB, AC));

        // Calculate midpoints of AB and AC
        glm::vec3 midpointAB = (A + B) * 0.5f;
        glm::vec3 midpointAC = (A + C) * 0.5f;

        // Vectors perpendicular to AB and AC, lying on the plane
        glm::vec3 perpAB = glm::cross(normal, AB);
        glm::vec3 perpAC = glm::cross(normal, AC);

        // Solve for intersection of lines (midpointAB + t1 * perpAB) and (midpointAC + t2 * perpAC)
        glm::vec3 t = glm::cross(perpAC, midpointAB - midpointAC) / glm::dot(perpAC, perpAB);

        // Circumcenter in 3D
        glm::vec3 circumcenter = midpointAB + t * perpAB;

        // Project the circumcenter onto the sphere surface
        circumcenter = glm::normalize(circumcenter);

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

