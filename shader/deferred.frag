#version 120

#define H0_MAX 0.005
// RC == 1.01 for dritab4 v1006
// RC == 1.0 for dritab4 v1007 (which uses a log1p table)
// shader-computed table now stores log(dri + 1)
#define RC 1.0

#define PI 3.14159265

const vec3 ISun = vec3(1.5);

const vec3 beta0_r = vec3(7.196e-6, 1.461e-5, 2.983e-5);
const vec3 beta0_m = vec3(0.00002);

const float h0_r = 7994.0;
const float h0_m = 1200.0;

uniform float time;

// planet radius
uniform float Rg;
uniform float Rt;

// view space sun normal
uniform vec3 sunnorm_v;

// view space planet position
uniform vec3 planetpos_v;

// exposure for HDR
uniform float exposure;

// original fragment position texture
uniform sampler2D sampler_position;

// surface normal texture
uniform sampler2D sampler_normal;

// diffuse material texture
uniform sampler2D sampler_diffuse;

// dritab texture
uniform sampler2D sampler_dritab;

varying vec4 fragpos_p;

vec3 hdr(vec3 L) {
	return vec3(1.0) - exp(-exposure * L);
}

// intersection of ray with arbitrary sphere, n must be normalized
bool intersect(vec3 ps, float radius, vec3 p0, vec3 n, out vec3 p) {
	p0 -= ps;
	float b = 2.0 * dot(n, p0);
	float c = dot(p0, p0) - radius * radius;
	float disc = b * b - 4.0 * c;
	if (disc < 0) return false;
	float sqrt_disc = sqrt(disc);
	float t0 = (-b - sqrt_disc) * 0.5;
	float t1 = (-b + sqrt_disc) * 0.5;
	if (t1 < 0.0) return false;
	if (t0 < 0.0) {
		p = ps + p0 + t1 * n;
		return true;
	}
	p = ps + p0 + t0 * n;
	return true;
}

// intersection of ray with arbitrary sphere (no output arg), n must be normalized
bool intersect(vec3 ps, float radius, vec3 p0, vec3 n) {
	p0 -= ps;
	float b = 2.0 * dot(n, p0);
	float c = dot(p0, p0) - radius * radius;
	float disc = b * b - 4.0 * c;
	if (disc < 0) return false;
	float t1 = (-b + sqrt(disc)) * 0.5;
	return t1 >= 0.0;
}

float dritab_eval(float h0, float r, float mu) {
	float u = pow(h0 / (Rg * H0_MAX), 1.0 / 3.0);
	float v = (mu + 1.0) * 0.5;
	float t = texture2D(sampler_dritab, vec2(u, v)).r;
	// t = min(t, 85.0); // infinities can cause problems, 87 ~= log(1e38), close to max ieee 32-bit float value
	// texture stores log(dri + 1)
	float k = (RC - r / Rg) / (h0 / Rg);
	return Rg * (exp(t + k) - exp(k));
}

// density ratio integral over finite path
float dri(float h0, vec3 pa_v, vec3 pb_v) {
	float ra = distance(pa_v, planetpos_v);
	float rb = distance(pb_v, planetpos_v);
	if (rb < ra) {
		// need pa to be lower
		vec3 temp0 = pa_v;
		pa_v = pb_v;
		pb_v = temp0;
		float temp1 = ra;
		ra = rb;
		rb = temp1;
	}
	vec3 n_v = normalize(pb_v - pa_v);
	return abs(dritab_eval(h0, ra, dot(normalize(pa_v - planetpos_v), n_v)) - dritab_eval(h0, rb, dot(normalize(pb_v - planetpos_v), n_v)));
}

// density ratio integral over infinite path
float dri_n(float h0, vec3 pa_v, vec3 n_v) {
	float ra = distance(pa_v, planetpos_v);
	return dritab_eval(h0, ra, dot(normalize(pa_v - planetpos_v), n_v));
}

// optical thickness over finite path
vec3 thickness(vec3 pa_v, vec3 pb_v) {
	return beta0_r * dri(h0_r, pa_v, pb_v) + 1.1 * beta0_m * dri(h0_m, pa_v, pb_v);
}

// optical thickness over infinite path
vec3 thickness_n(vec3 pa_v, vec3 n_v) {
	return beta0_r * dri_n(h0_r, pa_v, n_v) + 1.1 * beta0_m * dri_n(h0_m, pa_v, n_v);
}

float phase_r(float mu) {
	return (3.0 / (16.0 * PI)) * (1.0 + mu * mu);
}

float phase_m(float mu) {
	const float g = 0.75;
	return 3.0 * (1.0 - g * g) * (1.0 + mu * mu) / (8.0 * PI * (2.0 + g * g) * pow(1.0 + g * g - 2.0 * g * mu, 1.5));
}

void main() {
	vec4 temp;

	// get scene data
	vec2 coords = (fragpos_p.xy * 0.5) + vec2(0.5);
	vec3 pos_v = texture2D(sampler_position, coords).xyz;
	temp = texture2D(sampler_normal, coords);
	vec3 norm_v = temp.xyz;
	float tag = temp.w;
	vec4 diffuse = texture2D(sampler_diffuse, coords);
	
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
	vec3 L = vec3(0.0);
	vec3 att = vec3(1.0);
	
	bool hit_atmos = true;
	if (length(planetpos_v) > Rt + 5.0) {
		// camera not in atmosphere, find entrance
		if (intersect(planetpos_v, Rt, vec3(0.0), dp, p0)) {
			// is fragment closer?
			if (fragdist <= length(p0)) {
				hit_atmos = false;
			}
		} else {
			// didnt enter
			hit_atmos = false;
		}
	}
	
	//hit_atmos = false;
	
	bool hit_planet = false;
	if (hit_atmos) {
		// entered the atmosphere
		if (intersect(planetpos_v, Rg, p0, dp, p1)) {
			// hit planet
			hit_planet = true;
			if (fragdist > length(p1)) { // TODO adjust epsilon
				// fragment behind planet surface, replace with planet surface
				// TODO properly?
				tag = 0.0;
				diffuse = vec4(0.0, 0.05, 0.0, 0.0);
				pos_v = p1;
				norm_v = normalize(pos_v - planetpos_v);
			}
		} else {
			// missed planet, find exit from atmosphere
			p1 = p0; // just in case...
			intersect(planetpos_v, Rt, p0 + 50.0 * dp, dp, p1); // TODO adjust epsilon
		}
		if (fragdist <= length(p1)) {
			// fragment is closer
			p1 = pos_v;
		}

		// integrate inscatter over path...

		vec3 pa = p0;
		vec3 pb = p0;

		float ra = distance(pa, planetpos_v);
		float rb = ra;

		float min_x = 50.0 + distance(p0, p1) / 65.0;

		float xx = 0.0;
		float xx_max = distance(p0, p1) - 10.0;
		vec3 th_cam = vec3(0.0);

		float mu = dot(sunnorm_v, dp);
		float phr = phase_r(mu);
		float phm = phase_m(mu);

		// potential optimisation for thickness lookups:
		// do reverse lookups until horizon transition, forward lookups after (2 loops?)
		// still have to do 2 lookups because varying coefs, but can get rid of branches

		for(; xx < xx_max;) {
			// end point of this segment
			float x = min_x; // + 2.5 * h0_r * (1.0 - exp(0.5 * (Rg - ra) / h0_r));
			x = min(xx_max - xx + 5.0, x);
			xx += x;
			pb = pa + dp * x;
			float rb = distance(pb, planetpos_v);

			vec3 Lsun = ISun * exp(-thickness_n(pa, sunnorm_v));

			L += x * att * Lsun * ((beta0_r * exp((Rg - ra) / h0_r) * phr) + (beta0_m * exp((Rg - ra) / h0_m) * phm));
			
			att *= exp(-thickness(pa, pb));
			pa = pb;
			ra = rb;
		}

		// L = att;
	}

	// reflected light
	vec3 L1 = vec3(0.0);
	if (tag == 1.0) {
		// eg nice water
	} else {
		// normal operation
		vec3 Lsun = ISun * exp(-thickness_n(pos_v, sunnorm_v));
		vec3 L1a = vec3(0.0) * diffuse.rgb;
		float dot_d = dot(sunnorm_v, norm_v);
		vec3 L1d = Lsun * max(dot_d, 0.0) * diffuse.rgb;
		vec3 L1s = vec3(0.0);
		if (dot_d > 0) {
			// TODO specular
		}
		L1 = L1a + L1d + L1s;
	}

	// add inscattered and reflected light
	L += L1 * att;
	
	// direct light from the sun
	if (tag < 0.0 && !hit_planet && dot(sunnorm_v, dp) > 0.9999905) {
		L += ISun * att;
	}

	// holy shit this function is brilliant
	// http://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
	float junk = fract(sin(dot(gl_FragCoord.xy + time, vec2(12.9898, 78.233))) * 43758.5453);
	// gl_FragColor.rgb = vec3(junk);

	gl_FragColor.rgb = hdr(L) + vec3((junk - 0.5) / 127.0);
}
