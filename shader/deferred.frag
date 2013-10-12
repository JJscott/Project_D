#version 120

// this has the sampler for the transmittance texture
#include "atmos/common.frag"

const vec3 ISun = vec3(100.0);

// view space sun normal
uniform vec3 sunnorm_v;

// view space planet position
uniform vec3 planetpos_v;

// exposure for HDR
uniform float exposure;

// inscatter texture (4D)
uniform sampler3D sampler_inscatter;

// irradiance texture
uniform sampler2D sampler_irradiance;

// original fragment position texture
uniform sampler2D sampler_position;

// surface normal texture
uniform sampler2D sampler_normal;

// diffuse material texture
uniform sampler2D sampler_diffuse;

varying vec4 fragpos_p;

vec3 hdr(vec3 L) {
	return vec3(1.0) - exp(-exposure * L);
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

float magic_term(float a, float b, float c, float x) {
	return a * exp(-pow((x - b) / c, 2));
}

// analytic density ratio integral (for Rg = 6360000)
float dri(float H, float r, float theta) {
	// THIS ACTUALLY WORKS! praise the maou.
	theta = abs(theta);
	float magic = 0.0;
	// piecewise function fitted in matlab (in halves) using the 'fit' function with model 'gauss3'
	if (theta < M_PI * 0.5) {
		//magic += magic_term(1.303e13, 5.627, 0.7397, theta);
		//magic += magic_term(5.03, 2.238, 0.7728, theta);
		//magic += magic_term(0.2276, 0.9542, 0.4843, theta);
		magic += magic_term(+2.4906166e+09, +4.4987157e+00, +6.3150025e-01, theta);
		magic += magic_term(+4.9114098e+00, +2.2356545e+00, +7.8121914e-01, theta);
		magic += magic_term(+2.1232654e-01, +9.3667162e-01, +4.7475181e-01, theta);
	} else {
		theta -= M_PI * 0.5;
		//magic += magic_term(-56.31, 0.4758, 0.1257, theta);
		//magic += magic_term(13.99, 0.3213, 0.0009591, theta);
		//magic += magic_term(6.829e16, 8.512, 1.39, theta);
		magic += magic_term(-5.6851890e+01, +4.7739437e-01, +1.2632559e-01, theta);
		magic += magic_term(+0.0000000e+00, +6.2822074e-01, +9.5788390e-04, theta);
		magic += magic_term(+9.7339527e+15, +8.0646437e+00, +1.3527238e+00, theta);
	}
	return H * exp((Rg - r) / H) * exp(magic);
}

vec3 transmittance_Mk3(float r, float mu) {
	float theta = acos(mu);
	return exp(-betaR * dri(HR, r, theta) - 1.1 * betaM * dri(HM, r, theta));
}

vec3 transmittance_Mk3(vec3 pp, vec3 pa, vec3 pb) {
	float d = distance(pa, pb);
	if (d < 50.0) return vec3(1.0);
	vec3 n = normalize(pb - pa);
	pa -= pp;
	pb -= pp;
	
	//vec3 foo = vec3(1.0, 0.0, 0.0);
	
	if (length(pb) < length(pa)) {
		// we want the lookup rays to point up if possible
		vec3 temp = pa;
		pa = pb;
		pb = temp;
		n = -n;
		
		//foo = vec3(0.0, 0.0, 1.0);
	}
	float mu_a = dot(normalize(pa), n);
	float mu_b = dot(normalize(pb), n);
	//if (d < 100.0) return analyticTransmittance(length(pa) + 5.0, mu_a, d);
	vec3 trans_a = clamp(transmittance_Mk3(length(pa) + 5.0, mu_a), vec3(0.0), vec3(1.0));
	vec3 trans_b = clamp(transmittance_Mk3(length(pb) + 5.0, mu_b), vec3(0.0), vec3(1.0));
	return clamp(trans_a / trans_b, vec3(0.0), vec3(1.0));
}

vec3 transmittance(vec3 pp, vec3 pa, vec3 pb) {
	float d = distance(pa, pb);
	if (d < 10.0) return vec3(1.0);
	vec3 n = normalize(pb - pa);
	pa -= pp;
	pb -= pp;
	
	//vec3 foo = vec3(1.0, 0.0, 0.0);
	
	if (length(pb) < length(pa)) {
		// we want the lookup rays to point up if possible
		vec3 temp = pa;
		pa = pb;
		pb = temp;
		n = -n;
		
		//foo = vec3(0.0, 0.0, 1.0);
	}
	float mu_a = dot(normalize(pa), n);
	float mu_b = dot(normalize(pb), n);
	//if (d < 100.0) return analyticTransmittance(length(pa) + 5.0, mu_a, d);
	vec3 trans_a = clamp(transmittance(length(pa) + 5.0, mu_a), vec3(0.0), vec3(1.0));
	vec3 trans_b = clamp(transmittance(length(pb) + 5.0, mu_b), vec3(0.0), vec3(1.0));
	return clamp(trans_a / trans_b, vec3(0.0), vec3(1.0));
}

// this works in pseudo-world space
vec3 transmittance_naive(vec3 pa, vec3 pb, vec3 n) {
	float mu_a = dot(normalize(pa), n);
	float mu_b = dot(normalize(pb), n);
	vec3 trans_a = clamp(transmittance(length(pa) + 5.0, mu_a), vec3(0.0), vec3(1.0));
	vec3 trans_b = clamp(transmittance(length(pb) + 5.0, mu_b), vec3(0.0), vec3(1.0));
	return clamp(trans_a / trans_b, vec3(0.0), vec3(1.0));
}

vec3 transmittance_Mk2(vec3 pp, vec3 pa, vec3 pb) {
	float d = distance(pa, pb);
	if (d < 10.0) return vec3(1.0);
	vec3 n = normalize(pb - pa);
	pa -= pp;
	pb -= pp;
	
	// lookup both ways, then take the biggest
	vec3 trans_ab = transmittance_naive(pa, pb, n);
	vec3 trans_ba = transmittance_naive(pb, pa, -n);
	
	return max(trans_ab, trans_ba);
	
	//if (distance(trans_ab, vec3(0.8)) < distance(trans_ba, vec3(0.8))) return trans_ab;
	//else return trans_ba;
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
	vec3 L1sun = ISun;
	vec3 L1irr = vec3(0.0);

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
	
	if (hit_atmos) {
		// entered the atmosphere
		if (intersect(planetpos_v, Rg, p0, dp, p1)) {
			// hit planet
			if (fragdist > length(p1)) { // TODO adjust epsilon
				// fragment behind planet surface
				// TODO properly?
				tag = 0;
				diffuse = vec4(0.0, 0.05, 0.0, 0.0);
				pos_v = p1;
				norm_v = normalize(pos_v - planetpos_v);
			}
		} else {
			// missed planet, find exit from atmosphere
			// this should always succeed...
			intersect(planetpos_v, Rt, p0 + 100.0 * dp, dp, p1); // TODO adjust epsilon
		}
		if (fragdist <= length(p1)) {
			// fragment is closer
			p1 = pos_v;
		}

		float mu_vs = dot(dp, sunnorm_v);
		float Pr = phaseFunctionR(mu_vs);
		float Pm = phaseFunctionM(mu_vs);

		float r0 = max(length(p0 - planetpos_v), Rg + 5.0);
		float r1 = max(length(p1 - planetpos_v), Rg + 5.0);
		
		float mu0_vx = dot(dp, normalize(p0 - planetpos_v));
		float mu0_sx = dot(sunnorm_v, normalize(p0 - planetpos_v));
		
		float mu1_vx = dot(dp, normalize(p1 - planetpos_v));
		float mu1_sx = dot(sunnorm_v, normalize(p1 - planetpos_v));
		
		att = transmittance_Mk3(planetpos_v, p0, p1);

		vec4 insc0 = max(texture4D(sampler_inscatter, r0, mu0_vx, mu0_sx, mu_vs), vec4(0.0));
		vec4 insc1 = max(texture4D(sampler_inscatter, r1, mu1_vx, mu1_sx, mu_vs), vec4(0.0));
		
		vec4 insc = max(insc0 - att.rgbr * insc1, vec4(0.0));
		
		// TODO shadow test in transmittance here
		L1sun = ISun * clamp(transmittance_Mk3(r1, mu1_sx), vec3(0.0), vec3(1.0));
		L1irr = ISun * irradiance(sampler_irradiance, r1, mu1_sx);
		
		L = max(insc.rgb * Pr + getMie(insc) * Pm, vec3(0.0)) * ISun;
	}

	// reflected light
	vec3 L1 = vec3(0.0);
	if (tag == 1.0) {
		// eg nice water
	} else {
		// normal operation
		vec3 L1a = L1irr * diffuse.rgb;
		float dot_d = dot(sunnorm_v, norm_v);
		vec3 L1d = L1sun * max(dot_d, 0.0) * diffuse.rgb;
		vec3 L1s = vec3(0.0);
		if (dot_d > 0) {
			// TODO specular
		}
		L1 = L1a + L1d + L1s;
	}

	// add inscattered and reflected light
	L += L1 * att;

	gl_FragColor = vec4(hdr(L), 1.0);
}
