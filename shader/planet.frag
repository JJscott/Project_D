#version 120

#include "atmos/common.frag"

// view space sun normal
uniform vec3 sunnorm_v;

// view space planet position
uniform vec3 planetpos_v;

// camera position and sun normal, world space
uniform vec3 c;
uniform vec3 s;

// exposure for HDR
uniform float exposure;

// inscatter and irradiance textures (transmittanceSampler is in common)
uniform sampler3D inscatterSampler;
uniform sampler2D irradianceSampler;

// depth buffer precision magic
varying float logz;

// frag position
varying vec3 fragpos_v;
varying vec3 fragpos_w;

const vec3 ISun = vec3(100.0, 100.0, 100.0);

//#define FIX

// inscattered light along ray x+tv, when sun in direction s (=S[L]-T(x,x0)S[L]|x0)
vec3 inscatter(inout vec3 x, inout float t, vec3 v, vec3 s, out float r, out float mu, out vec3 attenuation) {
	vec3 result;
	r = length(x);
	mu = dot(x, v) / r;
	float d = -r * mu - sqrt(r * r * (mu * mu - 1.0) + Rt * Rt);
	if (d > 0.0) {
		// x in space and ray intersects atmosphere
		// move x to nearest intersection of ray with top atmosphere boundary
		x += d * v;
		t -= d;
		mu = (r * mu + d) / Rt;
		r = Rt;
	}
	if (r <= Rt) {
		// ray intersects atmosphere
		float nu = dot(v, s);
		float muS = dot(x, s) / r;
		float phaseR = phaseFunctionR(nu);
		float phaseM = phaseFunctionM(nu);
		vec4 inscatter = max(texture4D(inscatterSampler, r, mu, muS, nu), 0.0);
		if (t > 0.0) {
			vec3 x0 = x + t * v;
			float r0 = length(x0);
			float rMu0 = dot(x0, v);
			float mu0 = rMu0 / r0;
			float muS0 = dot(x0, s) / r0;
#ifdef FIX
			// avoids imprecision problems in transmittance computations based on textures
			attenuation = analyticTransmittance(r, mu, t);
#else
			attenuation = transmittance(r, mu, v, x0);
#endif
			//r0 = max(r0, Rg + 10.0);
			if (r0 >= Rg + 10.0) {
				// computes S[L]-T(x,x0)S[L]|x0
				inscatter = max(inscatter - attenuation.rgbr * texture4D(inscatterSampler, r0, mu0, muS0, nu), 0.0);
#ifdef FIX
				// avoids imprecision problems near horizon by interpolating between two points above and below horizon
				const float EPS = 0.004;
				float muHoriz = -sqrt(1.0 - (Rg / r) * (Rg / r));
				if (abs(mu - muHoriz) < EPS) {
					float a = ((mu - muHoriz) + EPS) / (2.0 * EPS);

					mu = muHoriz - EPS;
					r0 = sqrt(r * r + t * t + 2.0 * r * t * mu);
					mu0 = (r * mu + t) / r0;
					vec4 inScatter0 = texture4D(inscatterSampler, r, mu, muS, nu);
					vec4 inScatter1 = texture4D(inscatterSampler, r0, mu0, muS0, nu);
					vec4 inScatterA = max(inScatter0 - attenuation.rgbr * inScatter1, 0.0);

					mu = muHoriz + EPS;
					r0 = sqrt(r * r + t * t + 2.0 * r * t * mu);
					mu0 = (r * mu + t) / r0;
					inScatter0 = texture4D(inscatterSampler, r, mu, muS, nu);
					inScatter1 = texture4D(inscatterSampler, r0, mu0, muS0, nu);
					vec4 inScatterB = max(inScatter0 - attenuation.rgbr * inScatter1, 0.0);

					inscatter = mix(inScatterA, inScatterB, a);
				}
#endif
			} else {
				inscatter = vec4(1.0, 0.0, 0.0, 0.0);
			}
		}
#ifdef FIX
		// avoids imprecision problems in Mie scattering when sun is below horizon
		inscatter.w *= smoothstep(0.00, 0.02, muS);
#endif
		result = max(inscatter.rgb * phaseR + getMie(inscatter) * phaseM, 0.0);
	} else { // x in space and ray looking in space
		result = vec3(0.0);
	}
	return result * ISun;
	//return vec3(0.0);
}

// ground radiance at end of ray x+tv, when sun in direction s
// attenuated bewteen ground and viewer (=R[L0]+R[L*])
vec3 groundColor(vec3 x, float t, vec3 v, vec3 s, float r, float mu, vec3 attenuation) {
	vec3 result;
	if (t > 0.0) {
		// ray hits ground surface
		// ground reflectance at end of ray, x0
		vec3 x0 = x + t * v;
		float r0 = length(x0);
		vec3 n = x0 / r0;
		
		// direct sun light (radiance) reaching x0
		float muS = dot(n, s);
		vec3 sunLight = transmittanceWithShadow(max(r0, Rg + 10.0), muS);

		// precomputed sky light (irradiance) (=E[L*]) at x0
		vec3 groundSkyLight = irradiance(irradianceSampler, max(r0, Rg + 10.0), muS);

		// light reflected at x0 (=(R[L0]+R[L*])/T(x,x0))
		vec3 groundColor = vec3(0.01, 0.1, 0.0) * (max(muS, 0.0) * sunLight + groundSkyLight) * ISun / M_PI;

		//vec3 groundColor = vec3(0f);

		// water specular color due to sunLight
		vec3 h = normalize(s - v);
		float fresnel = 0.02 + 0.98 * pow(1.0 - dot(-v, h), 5.0);
		float waterBrdf = fresnel * pow(max(dot(h, n), 0.0), 150.0);
		//groundColor += max(waterBrdf, 0.0) * sunLight * ISun;

		//groundColor = vec3(0f, 0.5, 0f);

		result = attenuation * groundColor; //=R[L0]+R[L*]
	} else {
		// ray looking at the sky
		result = vec3(0.0);
	}
	//return result;
	return vec3(0.1, 0.2, 0.0);
}

// direct sun light for ray x+tv, when sun in direction s (=L0)
vec3 sunColor(vec3 x, float t, vec3 v, vec3 s, float r, float mu) {
	if (t > 0.0) {
		return vec3(0.0);
	} else {
		vec3 transmittance = r <= Rt ? transmittanceWithShadow(r, mu) : vec3(1.0); // T(x,xo)
		vec3 isun = step(cos(M_PI / 180.0), dot(v, s)) * ISun; // Lsun
		return transmittance * isun; // Eq (9)
	}
}

vec3 HDR(vec3 L) {
	// L = L * exposure;
	// L.r = L.r < 1.413 ? pow(L.r * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.r);
	// L.g = L.g < 1.413 ? pow(L.g * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.g);
	// L.b = L.b < 1.413 ? pow(L.b * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.b);
	// return L;
	// why not just this?
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

void main_not() {
	gl_FragDepth = logz;

	vec3 x = c;
	vec3 v = normalize(fragpos_w - c);

	float r = length(x);
	float mu = dot(x, v) / r;
	
	//float t = -r * mu - sqrt(r * r * (mu * mu - 1.0) + Rg * Rg);
	//
	//vec3 g = x - vec3(0.0, 0.0, Rg + 10.0);
	//float a = v.x * v.x + v.y * v.y - v.z * v.z;
	//float b = 2.0 * (g.x * v.x + g.y * v.y - g.z * v.z);
	//float c = g.x * g.x + g.y * g.y - g.z * g.z;
	//float d = -(b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
	//bool cone = d > 0.0 && abs(x.z + d * v.z - Rg) <= 10.0;
	//
	//if (t > 0.0) {
	//	if (cone && d < t) {
	//		t = d;
	//	}
	//} else if (cone) {
	//	t = d;
	//}

	vec3 tv;
	bool hit = intersect(vec3(0.0), Rg + 10.0, x, v, tv);
	float t = hit ? length(tv - c) : -1.0e9;

	vec3 attenuation;
	vec3 inscatterColor = inscatter(x, t, v, s, r, mu, attenuation); //S[L]-T(x,xs)S[l]|xs
	vec3 groundColor = groundColor(x, t, v, s, r, mu, attenuation); //R[L0]+R[L*]
	vec3 sunColor = sunColor(x, t, v, s, r, mu); //L0

	gl_FragColor = vec4(HDR(sunColor + groundColor + inscatterColor), 1.0); // Eq (16)
}

vec3 transmittance(vec3 pp, vec3 pa, vec3 pb) {
	vec3 n = normalize(pb - pa);
	pa -= pp;
	pb -= pp;
	if (length(pb) < length(pa)) {
		// we want the lookup rays to point up if possible
		vec3 temp = pa;
		pa = pb;
		pb = temp;
		n = -n;
	}
	float mu_a = dot(normalize(pa), n);
	float mu_b = dot(normalize(pb), n);
	vec3 trans_a = clamp(transmittance(length(pa) + 5.0, mu_a), vec3(0.0), vec3(1.0));
	vec3 trans_b = clamp(transmittance(length(pb) + 5.0, mu_b), vec3(0.0), vec3(1.0));
	return clamp(trans_a / trans_b, vec3(0.0), vec3(1.0));
}

void main() {
	gl_FragDepth = logz;
	
	// distance from camera to fragment
	float fragdist = length(fragpos_v);

	// ray direction
	vec3 dp = normalize(fragpos_v);
	// view vector (to camera)
	vec3 v = -dp;

	// find interval
	vec3 p0 = vec3(0.0); // camera
	vec3 p1;

	// init; we might not enter the atmosphere
	vec3 L1 = vec3(0.0);
	vec3 Ls = vec3(0.0);
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

	vec3 L1sun = ISun;
	vec3 L1irr = vec3(0.0);
	float mu1_sx = 1.0; // TODO from fragnorm
	
	if (hit_atmos) {
		// we entered the atmosphere
		L1 = vec3(0.2, 0.5, 0.0);
		if (!intersect(planetpos_v, Rg, p0, dp, p1)) {
			// missed planet, find exit from atmosphere
			// this should always succeed...
			intersect(planetpos_v, Rt, p0 + 100.0 * dp, dp, p1); // TODO adjust epsilon
			L1 = vec3(0.0);
		}
		// is fragment closer?
		if (fragdist <= length(p1)) {
			p1 = fragpos_v;
			L1 = vec3(0.0);
		}

		float mu_vs = dot(dp, sunnorm_v);
		float Pr = phaseFunctionR(mu_vs);
		float Pm = phaseFunctionM(mu_vs);

		float r0 = max(length(p0 - planetpos_v), Rg + 5.0);
		float r1 = max(length(p1 - planetpos_v), Rg + 5.0);
		
		float mu0_vx = dot(dp, normalize(p0 - planetpos_v));
		float mu0_sx = dot(sunnorm_v, normalize(p0 - planetpos_v));
		
		float mu1_vx = dot(dp, normalize(p1 - planetpos_v));
		mu1_sx = dot(sunnorm_v, normalize(p1 - planetpos_v));
		
		att = transmittance(planetpos_v, p0, p1);

		vec4 insc0 = max(texture4D(inscatterSampler, r0, mu0_vx, mu0_sx, mu_vs), vec4(0.0));
		vec4 insc1 = max(texture4D(inscatterSampler, r1, mu1_vx, mu1_sx, mu_vs), vec4(0.0));
		
		vec4 insc = max(insc0 - att.rgbr * insc1, vec4(0.0));
		
		L1sun = ISun * clamp(transmittance(r1, mu1_sx), vec3(0.0), vec3(1.0));
		L1irr = ISun * irradiance(irradianceSampler, r1, mu1_sx);
		
		Ls = max(insc.rgb * Pr + getMie(insc) * Pm, vec3(0.0)) * ISun;
	}
	
	L1 *= (L1sun * max(mu1_sx, 0.0) + L1irr) * L1 * att;

	gl_FragColor = vec4(HDR(L1 + Ls), 1.0);
}































