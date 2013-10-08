#version 120

varying vec4 p;
varying vec3 n;

varying float logz;

uniform float far;

void main()
{
	p = gl_ModelViewMatrix * gl_Vertex;
	n = normalize(gl_NormalMatrix * gl_Normal);
	gl_Position = gl_ProjectionMatrix * p;

	// http://outerra.blogspot.com/2012/11/maximizing-depth-buffer-range-and.html
	logz = log(gl_Position.w * 0.01 + 1.0) / log(far * 0.01 + 1.0);
	// the depth of all frags is written in the frag shader, so i think this isnt needed? seems to break things anyway
	// gl_Position.z = (2.0 * logz - 1.0) * gl_Position.w;
}
