#version 120

#include "common.frag"

// radius
uniform float r;

// magic
uniform vec4 dhdH;

// current layer of inscatter texture
uniform int layer;

void integrand(float r, float mu, float muS, float nu, float t, out vec3 ray, out vec3 mie) {
	ray = vec3(0.0);
	mie = vec3(0.0);
	float ri = sqrt(r * r + t * t + 2.0 * r * mu * t);
	float muSi = (nu * t + muS * r) / ri;
	//ri = max(Rg, ri);
	if (ri <= Rg) {
		ray = vec3(0.0);
		mie = vec3(0.0);
	}
	if (muSi >= -sqrt(1.0 - Rg * Rg / (ri * ri))) {
		vec3 ti = transmittance(r, mu, t) * transmittance(ri, muSi);
		ray = exp(-(ri - Rg) / HR) * ti;
		mie = exp(-(ri - Rg) / HM) * ti;
	}
}

void inscatter(float r, float mu, float muS, float nu, out vec3 ray, out vec3 mie) {
	ray = vec3(0.0);
	mie = vec3(0.0);
	float dx = limit(r, mu) / float(INSCATTER_INTEGRAL_SAMPLES);
	//float xi = 0.0;
	//vec3 rayi;
	//vec3 miei;
	//integrand(r, mu, muS, nu, 0.0, rayi, miei);
	for (int i = 0; i <= INSCATTER_INTEGRAL_SAMPLES; ++i) {
		float xj = (float(i) + 0.5) * dx;
		vec3 rayj;
		vec3 miej;
		integrand(r, mu, muS, nu, xj, rayj, miej);
		//ray += (rayi + rayj) / 2.0 * dx;
		//mie += (miei + miej) / 2.0 * dx;
		ray += rayj * dx;
		mie += miej * dx;
		//xi = xj;
		//rayi = rayj;
		//miei = miej;
	}
	ray *= betaR;
	mie *= betaM;
}

void main() {
	vec3 ray;
	vec3 mie;
	float mu, muS, nu;
	getMuMuSNu(r, dhdH, mu, muS, nu);
	inscatter(r, mu, muS, nu, ray, mie);
	// store separately Rayleigh and Mie contributions, WITHOUT the phase function factor
	// (cf "Angular precision")
	gl_FragData[0].rgb = ray;
	gl_FragData[1].rgb = mie;
}
