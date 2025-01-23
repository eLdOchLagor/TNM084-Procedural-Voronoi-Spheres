# Procedurally Generate Voronoi Spheres
This is a project made in the course TNM084 at Linköping University. 
The project is a program that allows the user to generate different sphere meshes using spherical Voronoi noise.
It is accomplished by generating a bunch of random seed points on the surface of a sphere using Fibonacci lattice, calculating the convex hull of the points to get the Delaunay triangulation, and then extracting the Voronoi edges since it is the dual of the Delaunay. A very simple mesh is then generated for exporting.
The convex hull is calculated using the Quickhull algorithm as implemented here: https://github.com/akuukka/quickhull

![final](https://github.com/user-attachments/assets/0aa877af-5774-4941-85b4-4c80701d3aaf)
![fibonacciRed](https://github.com/user-attachments/assets/38f7dff7-17b8-4dad-a95c-33051351cef4)
![user interface](https://github.com/user-attachments/assets/2e399b46-ced1-404d-bb76-467a532af99e)
