#version 120

// this is sensitive in single precision (1000 looks worse than 500, 10000 is really bad)
// TODO whyyyyyyyyyyyyyy?
#define SAMPLES 500

#define H0_MAX 0.005
#define RT 1.1
// dont use RC > 1.0 when working in float (1.02 breaks earth mie scatter)
#define RC 1.0005

// for 1024 mu samples, using 1.0 >= mu >= -0.05 doesnt seem to be an improvement over 1.0 >= mu >= -1.0

uniform int RES_H0;
uniform int RES_MU;

// ray-sphere intersection, n must be normalized
bool intersect(float radius, vec3 p0, vec3 n, out vec3 p) {
	float b = 2.0 * dot(n, p0);
	float c = dot(p0, p0) - radius * radius;
	float disc = b * b - 4.0 * c;
	if (disc < 0) return false;
	float sqrt_disc = sqrt(disc);
	float t0 = (-b - sqrt_disc) * 0.5;
	float t1 = (-b + sqrt_disc) * 0.5;
	if (t1 < 0.0) return false;
	if (t0 < 0.0) {
		p = p0 + t1 * n;
		return true;
	}
	p = p0 + t0 * n;
	return true;
}

float densityRatioIntegralLog1p(float rg, float rt, float h0, float r, float theta) {
	vec3 p0 = vec3(0, r, 0);
	// rotate j theta rads towards i
	vec3 n = vec3(-sin(theta), cos(theta), 0);
	vec3 p1;
	intersect(rt, p0, n, p1);

	// ray delta
	vec3 dp = (p1 - p0) / float(SAMPLES);
	float lx = log(length(dp));

	// _dont_ treat planet intersection specially

	// we're gonna compute this thing in log space for precision benefits

	// log identities:
	// 1. log(x * y) = log(x) + log(y)
	// 2. log(x / y) = log(x) - log(y)
	// 3. log(x + y) = log(x) + log(1 + exp(log(y) - log(x))) [if 0<y<x]
	// 4. log(x - y) = log(x) + log(1 - exp(log(y) - log(x))) [if 0<y<x]

	// this starts dri = 1.0, not 0.0
	// but that's ok, cause we're gonna store log(dri + 1) in the table
	float ldri = 0.0;
	vec3 pa = p0;
	vec3 pb = vec3(0.0);

	for (int i = 0; i < SAMPLES; i++) {
		pb = pa + dp;
		float ri = (length(pa) + length(pb)) * 0.5;
		// use log identity 1
		float lddri = ((rg - ri) / h0) + lx;
		if (lddri == lddri) {
			// we get nans when r == rg and h0 == 0
			
			// use log identity 3
			if (ldri > lddri) {
				ldri = ldri + log(1.0 + exp(lddri - ldri));
			} else {
				ldri = lddri + log(1.0 + exp(ldri - lddri));
			}

		}
		pa = pb;
	}

	return ldri;
}

void main() {
	float h0 = pow(gl_FragCoord.x / float(RES_H0), 3.0) * H0_MAX;
	float mu = gl_FragCoord.y / float(RES_MU) * 2.0 - 1.0;
	float log_dri = densityRatioIntegralLog1p(1.0, RT, h0, RC, acos(mu));
	// TODO this is to catch nans - is it actually needed ?
	if (log_dri != log_dri) log_dri = 0;
	gl_FragColor.r = log_dri;
}
