#version 330 core

#define M_PI 3.1415926535897932384626433832795

layout(lines) in;
layout(triangle_strip, max_vertices = 4) out;

uniform int numberOfPoints;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;



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
    vec4 quadCL = gl_in[0].gl_Position;
    vec4 quadCR = gl_in[1].gl_Position;

    vec3 offsetCL = vec3(normalize(quadCL)) * 0.02;
    vec3 offsetCR = vec3(normalize(quadCR)) * 0.02;
    
    vec4 quadTL = quadCL + vec4(-offsetCL,0.0);
    vec4 quadTR = quadCR + vec4(-offsetCR,0.0);
    vec4 quadBL = quadCL + vec4(offsetCL,0.0);    
    vec4 quadBR = quadCR + vec4(offsetCR,0.0);

    //Push the quad

    //Quad
    //0: TL
    gl_Position = projection * view * model * quadTL;
    EmitVertex();
    //1: BL
    gl_Position = projection * view * model * quadBL;
    EmitVertex();
    //2: TR
    gl_Position = projection * view * model * quadTR;
    EmitVertex();
    //3: BR
    gl_Position = projection * view * model * quadBR;
    EmitVertex();

    EndPrimitive(); 

}