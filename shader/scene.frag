#version 120
// #include this

// magic to resolve depth buffer precision issues
varying float logz;

varying vec3 fragpos_v;
varying vec3 fragnorm_v;

void fragment_data(vec3 pos_v, inout float tag, inout vec3 norm_v, inout vec4 diffuse);

void main() {
	gl_FragDepth = logz;
	// the 4th component of the normal texture is used for the tag
	float tag = 0.0;
	vec3 n = normalize(fragnorm_v);
	vec4 kd = vec4(0.0);
	fragment_data(fragpos_v, tag, n, kd);
	gl_FragData[0].xyz = fragpos_v;
	gl_FragData[1] = vec4(n, tag);
	gl_FragData[2] = kd;
}
