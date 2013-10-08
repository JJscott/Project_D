// #version 120
// #include this

// TODO -> uniforms
const vec3 beta0_sc_r = vec3(5.47e-6, 1.28e-5, 3.12e-5);
const vec3 beta0_sc_m = vec3(0.00001);
const float h0_r = 7994.0;
const float h0_m = 1200.0;

// sun 'intensity', this includes the solid angle
uniform vec3 Isun;

// view space sun normal
uniform vec3 sunnorm_v;

// view space planet position
uniform vec3 planetpos_v;

// planet radius
uniform float rg;

// table used for optical depth
uniform sampler2D sampler_dritab;

#define PI 3.1415926536
#define SAMPLES 30

// get density ratio integral from non-normalized values
// returns (dri_r, dri_m, dr_r, dr_m)
vec4 densityRatioIntegral(float r, float theta) {
	// all these values must match the 'encoding' in the dritab shader and the app
	// except the radius adjustment in the dritab shader
	theta = abs(theta); // TODO big angles?
	vec2 coords;
	coords.s = theta / PI;
	// Bruneton's mapping
	// coords.s = atan((cos(theta) + 0.15) / (1.0 + 0.15) * tan(1.5)) / 1.5;
	coords.t = pow((r / rg - 1.0) / (20.0 * h0_r / rg), 1.0 / 5.0);
	vec4 t = texture2D(sampler_dritab, coords);
	return vec4(rg * t.xy, t.zw);
}

vec4 densityRatioIntegral(float r, float theta, bool invert) {
	return densityRatioIntegral(r, invert ? (PI - theta) : theta);
}

// density ratio integral over finite path
vec2 densityRatioIntegral(vec3 pp, vec3 pa, vec3 pb) {
	pa -= pp;
	pb -= pp;
	if (length(pb) < length(pa)) {
		// we want the lookup rays to point up if possible
		vec3 temp = pa;
		pa = pb;
		pb = temp;
	}
	float pa_mag = length(pa);
	float pb_mag = length(pb);
	vec3 n = normalize(pb - pa);
	float theta_a = acos(dot(pa, n) / pa_mag);
	float theta_b = acos(dot(pb, n) / pb_mag);
	return densityRatioIntegral(pa_mag, theta_a).xy - densityRatioIntegral(pb_mag, theta_b).xy;
}

vec3 attenuation(float r, float theta) {
	vec2 dri = densityRatioIntegral(r, theta).xy;
	vec3 att_r = exp(-beta0_sc_r * dri.x);
	vec3 att_m = exp(-beta0_sc_m * 1.1 * dri.y);
	return min(att_r * att_m, vec3(1.0));
}

vec3 attenuation(vec3 pp, vec3 pa, vec3 pb) {
	vec2 dri = densityRatioIntegral(pp, pa, pb);
	vec3 att_r = exp(-beta0_sc_r * dri.x);
	vec3 att_m = exp(-beta0_sc_m * 1.1 * dri.y);
	return min(att_r * att_m, vec3(1.0)); // TODO investigate?
}

// intersection of ray with arbitrary sphere
bool intersect(vec3 ps, float radius, vec3 p0, vec3 n, out vec3 p) {
	p0 -= ps;
	n = normalize(n);
	// a = dot(n, n) == 1
	float b = 2.0 * dot(n, p0);
	float c = dot(p0, p0) - radius * radius;
	float disc = b * b - 4.0 * c;
	if (disc < 0) return false;
	disc = sqrt(disc);
	float t0 = (-b - disc) * 0.5;
	float t1 = (-b + disc) * 0.5;
	if (t1 < 0.0) return false;
	if (t0 < 0.0) {
		p = ps + p0 + t1 * n;
		return true;
	}
	p = ps + p0 + t0 * n;
	return true;
}

float phase_r(float theta) {
	float mu = cos(theta);
	return 3.0 / (16.0 * PI) * (1 + mu * mu);
}

float phase_m(float theta) {
	float g = 0.76; // TODO -> uniform
	float mu = cos(theta);
	return 3 * (1 - g * g) * (1 + mu * mu) / (8 * PI * (2 + g * g) * pow(1 + g * g - 2 * g * mu, 1.5));
}

// type: == 0 => frag, < 0 => planet, > 0 => sky
void inscatter(vec3 pos_v, out vec3 L, out vec3 att, out float scatterdist, out float type) {
	// 'top' of atmosphere
	float rt = rg + 20.0 * h0_r; // reduced from 20 * h0 for optimisation

	// distance from camera to fragment
	float fragdist = length(pos_v);

	// ray direction
	vec3 dp = normalize(pos_v);
	// view vector (to camera)
	vec3 v = -dp;

	// find interval
	vec3 p0 = vec3(0.0); // camera
	vec3 p1;

	// init; we might not enter the atmosphere
	type = 0.0;
	L = vec3(0.0);
	att = vec3(1.0);
	scatterdist = fragdist;

	int samples = 30;

	if (length(planetpos_v) > rt) {
		// camera not in atmosphere, find entrance
		if (intersect(planetpos_v, rt, vec3(0.0), dp, p0)) {
			// is fragment closer?
			if (fragdist <= length(p0)) {
				samples = 0;
			}
		} else {
			// didnt enter
			samples = 0;
		}
	}

	if (samples > 0) {
		type = -1.0;
		// we entered the atmosphere
		if (!intersect(planetpos_v, rg, p0, dp, p1)) {
			// missed planet, find exit from atmosphere
			// this should always succeed...
			intersect(planetpos_v, rt, p0 + 100.0 * dp, dp, p1); // TODO adjust epsilon
			type = 1.0;
		}
		// is fragment closer?
		if (fragdist <= length(p1)) {
			p1 = pos_v;
			type = 0.0;
		}
	}

	// distance to where inscatter calculations end
	scatterdist = samples == 0 ? fragdist : length(p1);

	// ray delta
	float dp_mag = length(p1 - p0) / float(samples);
	dp *= dp_mag;

	// TODO sample 'rate' proportional to density <-- this is proving really hard
	
	// now trace the path!
	float P_r = phase_r(acos(dot(v, -sunnorm_v))); // yay for parallel sunlight
	float P_m = phase_m(acos(dot(v, -sunnorm_v)));
	vec3 pa = p0;
	float ra = max(length(pa - planetpos_v), rg);

	for (int i = 0; i < samples; i++) {
		vec3 pb = pa + dp;
		vec3 p = 0.5 * pa + 0.5 * pb;

		float r = max(length(p - planetpos_v), rg);
		float rb = max(length(pb - planetpos_v), rg);

		att *= attenuation(planetpos_v, pa, pb); // TODO less texture lookups?
		
		// single inscatter
		vec3 att_sun = attenuation(r, acos(dot(normalize(p - planetpos_v), sunnorm_v)));
		vec3 S_r = beta0_sc_r * exp((rg - r) / h0_r) * P_r;
		vec3 S_m = beta0_sc_m * exp((rg - r) / h0_m) * P_m;
		L += Isun * att_sun * (S_r + S_m) * att * dp_mag;
		
		pa = pb;
		ra = rb;
	}
}