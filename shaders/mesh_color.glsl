#version 330

#ifdef VERTEX_SHADER

layout(location = 0) in vec3 position;
layout(location = 3) in vec4 color;

uniform mat4 mvpMatrix;

out vec4 vColor;

void main()
{
    gl_Position = mvpMatrix * vec4(position, 1.0);
    vColor = color;
}

#endif


#ifdef FRAGMENT_SHADER

in vec4 vColor;

uniform vec3 baseColor;
uniform int useVertexColor;

out vec4 fragment_color;

void main()
{
    if(useVertexColor == 1)
        fragment_color = vColor;
    else
        fragment_color = vec4(baseColor, 1.0);
}

#endif