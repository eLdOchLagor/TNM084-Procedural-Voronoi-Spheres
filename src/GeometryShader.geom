#version 330 core

layout(points) in;                 // Input: single point
layout(points, max_vertices = 1) out; // Output: single point

void main()
{
    gl_Position = gl_in[0].gl_Position; // Pass-through position
    EmitVertex();                       // Emit the single vertex
    EndPrimitive();                     // End the primitive
}