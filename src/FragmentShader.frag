#version 330 core
out vec4 FragColor;

in vec3 normal; // from vertex shader
in vec3 fragPos; // add this to the vertex shader output

uniform vec3 cameraPos;

vec3 lightPos = vec3(0.0, 2.0, 2.0);
vec3 lightColor = vec3(1.0, 1.0, 1.0);
vec3 objectColor = vec3(0.2, 0.2, 0.2);

void main()
{
    // Normalize the interpolated normal
    vec3 norm = normalize(normal);

    // Calculate the light direction
    vec3 lightDir = normalize(lightPos - fragPos);

    // Ambient component
    float ambientStrength = 0.2;
    vec3 ambient = ambientStrength * lightColor;

    // Diffuse component
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // Specular component (Blinn-Phong)
    float specularStrength = 0.5;
    vec3 viewDir = normalize(cameraPos - fragPos); // Direction to the viewer
    vec3 halfwayDir = normalize(lightDir + viewDir); // Blinn-Phong halfway vector
    float spec = pow(max(dot(norm, halfwayDir), 0.0), 32.0); // 32 is the shininess factor
    vec3 specular = specularStrength * spec * lightColor;

    // Combine components
    vec3 result = (ambient + diffuse + specular) * objectColor;
    FragColor = vec4(result, 1.0);
}