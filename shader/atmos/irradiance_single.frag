#version 120

#include "common.frag"

void main() {
	float r, muS;
	getIrradianceRMuS(r, muS);
	gl_FragColor = vec4(transmittance(r, muS) * max(muS, 0.0), 0.0);
}
