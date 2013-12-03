#version 120

#define SAMPLES 1000

#define H0_MAX 0.005
#define RT 1.1
#define RC 1.02

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

float densityRatioIntegral(float rg, float rt, float h0, float r, float theta) {
	vec3 p0 = vec3(0, r, 0);
	// rotate j theta rads towards i
	vec3 n = vec3(-sin(theta), cos(theta), 0);
	vec3 p1;
	intersect(rt, p0, n, p1);

	// ray delta
	vec3 dp = (p1 - p0) / float(SAMPLES);

	// _dont_ treat planet intersection specially

	float dri = 0.0;
	vec3 pa = p0;
	vec3 pb = vec3(0.0);

	for (int i = 0; i < SAMPLES; i++) {
		pb = pa + dp;
		float ri = (length(pa) + length(pb)) * 0.5;
		float ddri = exp((rg - ri) / h0) * length(dp);
		if (ddri == ddri) {
			// we get nans when r == rg and h0 == 0
			dri += ddri;
		}
		pa = pb;
	}

	return dri;
}

void main() {
	float h0 = pow(gl_FragCoord.x / float(RES_H0), 3.0) * H0_MAX;
	float mu = gl_FragCoord.y / float(RES_MU) * 2.0 - 1.0;
	float log_dri = log(densityRatioIntegral(1.0, RT, h0, RC, acos(mu)));
	gl_FragColor.r = -1.0e10;
	if (log_dri == log_dri) gl_FragColor.r = log_dri;
}