#version 330 core

out vec4 FragColor;

in vec2 TexCoords;

int numTiles = 5; // Number of tiles along one axis

float random(vec2 st) {
    st = vec2(dot(st, vec2(127.1, 311.7)), dot(st, vec2(269.5, 183.3)));
    return fract(sin(st.x * st.y) * 43758.5453123);
}

vec3 randomColor(vec2 seed) {
    return vec3(random(seed), random(seed + vec2(1.0)), random(seed + vec2(2.0)));
}

void main() {
    float tileSize = 1.0 / float(numTiles);
    vec2 gridIndex = floor(TexCoords / tileSize);
    vec2 localUV = fract(TexCoords / tileSize);

    // Jittered seed position
    vec2 seed = vec2(
        random(gridIndex + vec2(0.1, 0.2)),
        random(gridIndex + vec2(0.3, 0.4))
    ) * tileSize;
    float jitterAmount = 3; // Reduce jitter to avoid overlaps
    seed += jitterAmount * (random(gridIndex + vec2(1.0)) - 0.5) * tileSize;

    vec2 seedPosition = (gridIndex + seed) * tileSize;
    float minDist = length(TexCoords - seedPosition);
    vec2 closestSeedPosition = seedPosition;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            vec2 neighborIndex = gridIndex + vec2(dx, dy);

            // Consistent neighbor seed jitter
            vec2 neighborSeed = vec2(
                random(neighborIndex + vec2(0.1, 0.2)),
                random(neighborIndex + vec2(0.3, 0.4))
            ) * tileSize;
            neighborSeed += jitterAmount * (random(neighborIndex + vec2(1.0)) - 0.5) * tileSize;

            vec2 neighborSeedPosition = (neighborIndex + neighborSeed) * tileSize;
            float dist = length(TexCoords - neighborSeedPosition);
            if (dist < minDist) {
                minDist = dist;
                closestSeedPosition = neighborSeedPosition;
            }
        }
    }

    vec3 cellColor = randomColor(closestSeedPosition);
    FragColor = vec4(cellColor, 1.0);
}
