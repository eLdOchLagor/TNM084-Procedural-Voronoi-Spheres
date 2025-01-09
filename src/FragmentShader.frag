#version 330 core

uniform int numberOfPoints;

out vec4 FragColor;   // Output fragment color

flat in float outIndex;

void main()
{
    FragColor = vec4(outIndex/float(numberOfPoints), 0.0, 0.0, 1.0);
}