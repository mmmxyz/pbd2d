#version 150 core

in vec4 position;
in vec4 color;
out vec4 vcolor;

void main()
{
		gl_Position = position;
		vcolor = color;
}
