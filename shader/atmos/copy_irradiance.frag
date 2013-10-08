#version 120

#include "common.frag"

uniform sampler2D deltaESampler;

void main() {
	vec2 uv = gl_FragCoord.xy / vec2(SKY_W, SKY_H);
	gl_FragColor = texture2D(deltaESampler, uv);
}