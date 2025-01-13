#version 330 core
layout (location = 0) in vec3 aPos;
//layout(location = 1) in vec2 inIndex; // Grid index i

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

//flat out vec2 index;

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

    float jitterLevel = 0.05;

    vec3 randomOffset = random3(aPos.xyz);

    //index = inIndex;
    vec4 newPosition = vec4(normalize(aPos + randomOffset * jitterLevel), 1.0);
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}