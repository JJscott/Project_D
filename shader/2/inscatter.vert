#version 120

varying vec4 fragpos_p;

void main() {
	gl_Position = gl_Vertex;
	fragpos_p = gl_Position;
}
