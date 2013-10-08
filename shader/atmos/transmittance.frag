#version 120

#include "common.frag"

// intersection of a ray with an origin centered sphere
bool intersect(in float radius, in vec3 p0, in vec3 n, out vec3 p) {
	p = vec3(0.0, 0.0, 0.0);
	n = normalize(n);
	// a = dot(n, n) == 1
	float b = 2.0 * dot(n, p0);
	float c = dot(p0, p0) - pow(radius, 2);
	float disc = pow(b, 2) - 4.0 * c;
	if (disc < 0) return false;
	disc = sqrt(disc);
	float t0 = (-b - disc) / (2.0);
	float t1 = (-b + disc) / (2.0);
	if (t1 < 0.0) return false;
	if (t0 < 0.0) {
		p = p0 + t1 * n;
		return true;
	}
	p = p0 + t0 * n;
	return true;
}

// density ratio integral for 2 scale heights at the same time
// assumes r is less than Rt
vec2 densityRatioIntegral(in vec2 h0, in float r, in float mu) {
	vec3 p0 = vec3(0.0, r, 0.0);
	// ray direction: rotate j acos(mu) towards i
	vec3 n = vec3(-sin(acos(mu)), mu, 0);
	// test if we hit the planet
	vec3 p1;
	if (intersect(Rg, p0, n, p1)) {
		// its a hit! return arbitrary big value
		return vec2(1e22, 1e22);
	}
	// find exit from atmosphere
	intersect(Rt, p0, n, p1);
	// ray delta
	vec3 dray = n * (distance(p0, p1) * 1.2 / float(TRANSMITTANCE_INTEGRAL_SAMPLES));
	// length of a segment
	float dx = length(dray);
	// result of the integral
	vec2 dri = vec2(0.0, 0.0);
	vec3 pa = p0;
	vec3 pb;
	for (int i = 0; i < TRANSMITTANCE_INTEGRAL_SAMPLES; i++) {
		// exp(-h / h0)
		pb = pa + dray;
		float ri = length(0.5 * (pa + pb));
		if (ri < 0.999 * Rg) return vec2(1.0e11);
		vec2 ddri = exp((Rg - ri) / h0) * dx;
		if (ddri == ddri) dri += ddri; // nan guard
		pa = pb;
	}
	return dri;
}

void main() {
	float r, mu;
	getTransmittanceRMu(r, mu);
	vec2 dri = densityRatioIntegral(vec2(HR, HM), r, mu);
	vec3 t = betaR * dri.x + betaM * dri.y;
	gl_FragColor = vec4(exp(-t), 0.0);
}

