#version 330 core

#define M_PI 3.1415926535897932384626433832795

layout(points) in;
layout(line_strip, max_vertices = 2) out;

uniform int numberOfPoints;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

flat in vec2 index[];
flat out vec2 outIndex;

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
    /*
    for( int i=-1; i<=1; i++ ) {

        float inc = 2.0f * M_PI * (index[0].x + i) / sqrt(numberOfPoints);

        for( int j=-1; j<=1; j++ ) {
            float az = M_PI * (index[0].y + j) / sqrt(numberOfPoints);

            vec3 neighborPoint = vec3(1.0 * sin(az) * cos(inc), 1.0 * sin(az) * sin(inc), 1.0 * cos(az));
        }
    }
    */

    outIndex = index[0];

    vec4 spherePoint = gl_in[0].gl_Position;

    gl_Position = projection * view * model * spherePoint;
    EmitVertex();

    vec4 endPoint = spherePoint * vec4(1.1, 1.1, 1.1, 1.0);
    gl_Position = projection * view * model * endPoint;
    EmitVertex();

    EndPrimitive();
}