#version 120

const vec3 planetpos_w = vec3(0.0);

// sun direction
uniform vec3 sunnorm_w;

// distance to the far plane
uniform float far;

// inverse view transform
uniform mat4 view;
uniform mat4 view_inv;

// magic to resolve depth buffer precision issues
varying float logz;

varying vec3 fragpos_w;
varying vec3 fragpos_v;
varying vec3 fragnorm_w;
varying vec3 fragnorm_v;
varying vec3 planetpos_v;

void main() {
	vec4 p = gl_ModelViewMatrix * gl_Vertex;
	vec4 n = vec4(gl_NormalMatrix * gl_Normal, 0.0);
	fragpos_v = p.xyz;
	fragnorm_v = normalize(n.xyz);
	fragpos_w = (view_inv * p).xyz;
	fragnorm_w = normalize(view_inv * n).xyz;
	// planet in view space
	planetpos_v = (view * vec4(planetpos_w, 1.0)).xyz;
	gl_Position = gl_ProjectionMatrix * p;
	// http://outerra.blogspot.com/2012/11/maximizing-depth-buffer-range-and.html
	logz = log(gl_Position.w * 0.01 + 1.0) / log(far * 0.01 + 1.0);
}