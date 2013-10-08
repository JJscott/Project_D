#version 120

#include "atmos.frag"

// position texture
uniform sampler2D sampler_pos;

varying vec4 fragpos_p;

void main() {
	// grab original position
	vec2 coords = (fragpos_p.xy * 0.5) + vec2(0.5);
	vec3 pos_v = texture2D(sampler_pos, coords).xyz;
	vec3 L, att;
	float scatterdist;
	float type; // 0 => frag, < 0 => planet, > 0 => sky
	inscatter(pos_v, L, att, scatterdist, type);

	// need to do incident sunlight in a post-process shader, should be per pixel
	
	gl_FragData[0] = vec4(L, scatterdist);
	gl_FragData[1] = vec4(att, type);
}
