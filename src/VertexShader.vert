#version 330 core
layout (location = 0) in vec3 aPos;
layout(location = 1) in ivec2 inGridIndex; // Grid indices (i, j)

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out ivec2 gridIndex; // Pass indices to the next stage (optional)

void main()
{
    gridIndex = inGridIndex;
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}