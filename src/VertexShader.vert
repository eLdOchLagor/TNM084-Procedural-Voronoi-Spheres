#version 330 core
layout (location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 normal;
out vec3 fragPos; // Pass fragment position to fragment shader

void main()
{
    // Transform vertex position to world space
    vec4 worldPos = model * vec4(aPos, 1.0);
    fragPos = vec3(worldPos);

    // Transform normal to world space
    normal = mat3(transpose(inverse(model))) * aNormal;

    // Final position
    gl_Position = projection * view * worldPos;
}