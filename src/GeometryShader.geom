#version 330 core

layout(points) in;                 // Input: single point
layout(points, max_vertices = 1) out; // Output: single point

in ivec2 gridIndex[];

flat out ivec2 outGridIndex;

// Random3 function taken from Lab4 geometry shader
vec3 random3(vec3 st)
{
    st = vec3( dot(st,vec3(127.1,311.7, 543.21)),
              dot(st,vec3(269.5,183.3, 355.23)),
              dot(st,vec3(846.34,364.45, 123.65)) ); // Haphazard additional numbers by IR
    return -1.0 + 2.0*fract(sin(st)*43758.5453123);
}

void main()
{
    outGridIndex = gridIndex[0];
    gl_Position = gl_in[0].gl_Position; // Pass-through position
    EmitVertex();                       // Emit the single vertex
    EndPrimitive();                     // End the primitive
}