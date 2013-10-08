#version 120

// distance to the far plane
uniform float far;

// inverse view transform
uniform mat4 view_inv;

// magic to resolve depth buffer precision issues
varying float logz;

varying vec3 fragpos_v;
varying vec3 fragpos_w;

void main() {
	vec4 p = gl_ModelViewMatrix * gl_Vertex;
	vec4 n = vec4(gl_NormalMatrix * gl_Normal, 0.0);
	fragpos_w = (view_inv * p).xyz;
	fragpos_v = p.xyz;
	gl_Position = gl_ProjectionMatrix * p;
	// http://outerra.blogspot.com/2012/11/maximizing-depth-buffer-range-and.html
	logz = log(gl_Position.w * 0.01 + 1.0) / log(far * 0.01 + 1.0);
}