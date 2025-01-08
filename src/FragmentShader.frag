#version 330 core

out vec4 FragColor;   // Output fragment color

flat in ivec2 outGridIndex;

void main()
{
    FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}