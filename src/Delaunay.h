#pragma once

#include <vector>
#include <cmath>
#include <tuple>
#include <glm/glm.hpp>

class Delaunay
{
public:

	struct Tetrahedron {
		int p1, p2, p3, p4; // Indices of the vertices in the points vector
	};

	struct Edge {
		int p1, p2; // Indices of the edge's endpoints

		bool operator==(const Edge& other) const {
			return (p1 == other.p1 && p2 == other.p2) || (p1 == other.p2 && p2 == other.p1);
		}
	};

	float squaredDistance(const glm::vec3& a, const glm::vec3& b);
	float determinant(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& d, const glm::vec3& p);

	void removeDuplicateEdges(std::vector<Edge>& edges);

	bool isPointInCircumsphere(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& d, const glm::vec3& p);

	Tetrahedron createSuperTetrahedron(const std::vector<glm::vec3>& points);

	std::vector<Tetrahedron> bowyerWatson3D(std::vector<glm::vec3>& points);


private:

};

