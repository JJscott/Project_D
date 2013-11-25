#version 120

#include "scene.frag"

void fragment_data(vec3 pos_v, inout float tag, inout vec3 norm_v, inout vec4 diffuse) {
	diffuse.rgb = vec3(0.0);
	tag = -1.0;
}
