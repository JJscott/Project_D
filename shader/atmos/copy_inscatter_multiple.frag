#version 120

#include "common.frag"

// radius
uniform float r;

// magic
uniform vec4 dhdH;

// current layer of inscatter texture
uniform int layer;

// deltaS texture - for multiple scattering both forms are combined
uniform sampler3D deltaSSampler;

void main() {
	float mu, muS, nu;
	getMuMuSNu(r, dhdH, mu, muS, nu);
	vec3 uvw = vec3(gl_FragCoord.xy, float(layer) + 0.5) / vec3(ivec3(RES_MU_S * RES_NU, RES_MU, RES_R));
	gl_FragColor = vec4(texture3D(deltaSSampler, uvw).rgb / phaseFunctionR(nu), 0.0);
}
