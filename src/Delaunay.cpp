#include "Delaunay.h"

float Delaunay::squaredDistance(const glm::vec3& a, const glm::vec3& b) {
    // Return the squared distance between two points
    return glm::dot(a - b, a - b);
}

float Delaunay::determinant(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& d, const glm::vec3& p) {
    // Compute squared distances from point p to points a, b, c, and d
    float A_w = Delaunay::squaredDistance(a, p);
    float B_w = Delaunay::squaredDistance(b, p);
    float C_w = Delaunay::squaredDistance(c, p);
    float D_w = Delaunay::squaredDistance(d, p);

    // Create the 4x4 matrix
    glm::mat4 mat(
        glm::vec4(a.x, a.y, a.z, A_w),
        glm::vec4(b.x, b.y, b.z, B_w),
        glm::vec4(c.x, c.y, c.z, C_w),
        glm::vec4(d.x, d.y, d.z, D_w)
    );

    // Return the determinant of the 4x4 matrix
    return glm::determinant(mat);
}

bool Delaunay::isPointInCircumsphere(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& d, const glm::vec3& p) {
    double det = Delaunay::determinant(a, b, c, d, p);
    return det > 0; // Positive determinant means the point is inside the circumsphere
}

Delaunay::Tetrahedron Delaunay::createSuperTetrahedron(const std::vector<glm::vec3>& points) {
    float minX = points[0].x, minY = points[0].y, minZ = points[0].z;
    float maxX = points[0].x, maxY = points[0].y, maxZ = points[0].z;

    // Find bounding box of points
    for (const auto& point : points) {
        minX = std::min(minX, point.x);
        minY = std::min(minY, point.y);
        minZ = std::min(minZ, point.z);
        maxX = std::max(maxX, point.x);
        maxY = std::max(maxY, point.y);
        maxZ = std::max(maxZ, point.z);
    }
    
    double dx = maxX - minX;
    double dy = maxY - minY;
    double dz = maxZ - minZ;

    double deltaMax = std::max({ dx, dy, dz });

    // Create the super-tetrahedron points
    glm::vec3 p1 = { minX - 10 * deltaMax, minY - 10 * deltaMax, minZ - 10 * deltaMax };
    glm::vec3 p2 = { minX + 10 * deltaMax, minY - 10 * deltaMax, minZ - 10 * deltaMax };
    glm::vec3 p3 = { minX + 10 * deltaMax, minY + 10 * deltaMax, minZ + 10 * deltaMax };
    glm::vec3 p4 = { minX - 10 * deltaMax, minY + 10 * deltaMax, minZ + 10 * deltaMax };

    // Add these points to the points vector
    return Tetrahedron{ static_cast<int>(points.size()), static_cast<int>(points.size() + 1),
                       static_cast<int>(points.size() + 2), static_cast<int>(points.size() + 3) };
}

void Delaunay::removeDuplicateEdges(std::vector<Delaunay::Edge>& edges) {
    std::vector<Delaunay::Edge> uniqueEdges;
    for (const auto& edge : edges) {
        auto it = std::find(uniqueEdges.begin(), uniqueEdges.end(), edge);
        if (it == uniqueEdges.end()) {
            uniqueEdges.push_back(edge);
        }
        else {
            uniqueEdges.erase(it); // Remove duplicate
        }
    }
    edges = uniqueEdges;
}

std::vector<Delaunay::Tetrahedron> Delaunay::bowyerWatson3D(std::vector<glm::vec3>& points) {
    std::vector<Tetrahedron> tetrahedra;

    // Add the super-tetrahedron
    Tetrahedron superTetrahedron = createSuperTetrahedron(points);
    tetrahedra.push_back(superTetrahedron);

    // Insert each point into the triangulation
    for (int i = 0; i < points.size(); ++i) {
        const glm::vec3& point = points[i];

        // Step 1: Find all tetrahedra whose circumsphere contains the point
        std::vector<Tetrahedron> badTetrahedra;
        for (const auto& tet : tetrahedra) {
            if (isPointInCircumsphere(points[tet.p1], points[tet.p2], points[tet.p3], points[tet.p4], point)) {
                badTetrahedra.push_back(tet);
            }
        }

        // Step 2: Find the boundary of the cavity (edges not shared by two tetrahedra)
        std::vector<Edge> cavityEdges;
        for (const auto& tet : badTetrahedra) {
            // Add edges of the bad tetrahedra to the cavity list
            cavityEdges.push_back(Edge{ tet.p1, tet.p2 });
            cavityEdges.push_back(Edge{ tet.p2, tet.p3 });
            cavityEdges.push_back(Edge{ tet.p3, tet.p1 });
            cavityEdges.push_back(Edge{ tet.p1, tet.p4 });
            cavityEdges.push_back(Edge{ tet.p2, tet.p4 });
            cavityEdges.push_back(Edge{ tet.p3, tet.p4 });
        }
        removeDuplicateEdges(cavityEdges);

        // Step 3: Remove bad tetrahedra
        tetrahedra.erase(std::remove_if(tetrahedra.begin(), tetrahedra.end(),
            [&badTetrahedra](const Tetrahedron& tet) {
                return std::find(badTetrahedra.begin(), badTetrahedra.end(), tet) != badTetrahedra.end();
            }),
            tetrahedra.end());

        // Step 4: Re-triangulate the cavity
        for (const auto& edge : cavityEdges) {
            tetrahedra.push_back(Tetrahedron{ edge.p1, edge.p2, i });
        }
    }

    // Step 5: Remove tetrahedra that share vertices with the super-tetrahedron
    tetrahedra.erase(std::remove_if(tetrahedra.begin(), tetrahedra.end(),
        [&superTetrahedron](const Tetrahedron& tet) {
            return (tet.p1 == superTetrahedron.p1 || tet.p1 == superTetrahedron.p2 || tet.p1 == superTetrahedron.p3 || tet.p1 == superTetrahedron.p4 ||
                tet.p2 == superTetrahedron.p1 || tet.p2 == superTetrahedron.p2 || tet.p2 == superTetrahedron.p3 || tet.p2 == superTetrahedron.p4 ||
                tet.p3 == superTetrahedron.p1 || tet.p3 == superTetrahedron.p2 || tet.p3 == superTetrahedron.p3 || tet.p3 == superTetrahedron.p4 ||
                tet.p4 == superTetrahedron.p1 || tet.p4 == superTetrahedron.p2 || tet.p4 == superTetrahedron.p3 || tet.p4 == superTetrahedron.p4);
        }),
        tetrahedra.end());

    return tetrahedra;
}
