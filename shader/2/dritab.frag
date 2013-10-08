#version 120

// scale heights
uniform float h0_r;
uniform float h0_m;

varying vec4 fragpos_p;

const float PI = 3.1415926536;

const int SAMPLES = 500;

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

float densityRatioIntegral(float rg, float rt, float h0, float r, float theta) {
	
	vec3 p0 = vec3(0.0, r, 0.0);
	vec3 n = vec3(-sin(theta), cos(theta), 0);
	vec3 p1;
	if (!intersect(rt, p0, n, p1)) return 0.0;
	
	vec3 dp = (p1 - p0) / float(SAMPLES);
	float dp_mag = length(dp);
	
	float dri = 0.0;
	vec3 pa = p0;
	vec3 pb;
	
	for (int i = 0; i < SAMPLES; i++) {
		pb = pa + dp;
		float ri = length(0.5 * (pa + pb));
		if (ri < 0.999 * rg) return 1.0e11;
		float ddri = exp((rg - ri) / h0) * dp_mag;
		if (ddri == ddri) dri += ddri; // nan guard
		pa = pb;
	}
	
	return dri;
}

void main() {
	// equiv texture coords
	vec2 coords = vec2(fragpos_p.x * 0.5 + 0.5, fragpos_p.y * 0.5 + 0.5);

	float rt = 1.0 + 20.0 * h0_r;
	float r = pow(coords.t, 5.0) * (rt - 1.0) + 1.0;
	float theta = PI * coords.s; //(pow(uv.x * 2.0 - 1.0, 3.0) * 0.5 + 0.5);
	
	// this mapping from Bruneton's shaders
	//float mu = -0.15 + tan(1.5 * coords.s) / tan(1.5) * (1.0 + 0.15);
	//float theta = acos(mu);

	// adjust r slightly for some error margin
	float r_r = r + 0.01 * h0_r;
	float r_m = r + 0.01 * h0_m;

	float dri_r = densityRatioIntegral(1.0, rt, h0_r, r_r, theta);
	float dri_m = densityRatioIntegral(1.0, rt, h0_m, r_m, theta);
	float dr_r = exp((1.0 - r_r) / h0_r);
	float dr_m = exp((1.0 - r_m) / h0_m);
	
	gl_FragColor.rgba = vec4(dri_r, dri_m, dr_r, dr_m);
}
















