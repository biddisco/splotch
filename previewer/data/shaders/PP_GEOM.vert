#version 120

uniform mat4 MVP; //ModelViewProjectionMatrix
varying vec3 normal; 

void main()
{	

    //pass through colour and position after multiplying pos by matrices
    gl_FrontColor = vec4(gl_Color.x, gl_Color.y, gl_Color.z, 1.0);

	normal = gl_Normal;

    gl_Position = MVP * gl_Vertex;
}


