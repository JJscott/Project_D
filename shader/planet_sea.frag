#version 120

#include "planet2.frag"

vec3 fragment_L1(vec3 Lsun) {
	return Lsun * vec3(0.0, 0.0, 0.05) * max(dot(fragnorm_v, sunnorm_v), 0.0);
}
