#version 120

// dimensions of the transmittance table
uniform int TRANSMITTANCE_W;
uniform int TRANSMITTANCE_H;

// dimensions of the irradiance table
uniform int SKY_W;
uniform int SKY_H;

// dimensions of the inscatter tables
uniform int RES_R;
uniform int RES_MU;
uniform int RES_MU_S;
uniform int RES_NU;

// rayleigh and mie scattering coefs
uniform vec3 betaR;
uniform vec3 betaM;
uniform float mieG;

// rayleigh and mie scale heights
uniform float HR;
uniform float HM;

// radius of planet and atmosphere
uniform float Rg;
uniform float Rt;

// some other radius for black magic
uniform float RL;

// transmittance texture
uniform sampler2D transmittanceSampler;

// ground reflectance for multiple scattering
const float AVERAGE_GROUND_REFLECTANCE = 0.1;

// constants
const int TRANSMITTANCE_INTEGRAL_SAMPLES       = 500;
// inscatter samples increased to reduce sunset banding
// do not increase spherical integral samples further; my video driver exceeded Windows timeout with 64
const int INSCATTER_INTEGRAL_SAMPLES           = 300; // 50
const int INSCATTER_SPHERICAL_INTEGRAL_SAMPLES = 32; // 16
const int IRRADIANCE_INTEGRAL_SAMPLES          = 32;
const float M_PI = 3.1415926535898; // ha, Bruneton was off by 4 in his last (9th) dp. not that it matters, only about 7sf in a float.

vec2 getTransmittanceUV(float r, float mu) {
	float uR, uMu;
	uR = sqrt((r - Rg) / (Rt - Rg));
	uMu = atan((mu + 0.15) / (1.0 + 0.15) * tan(1.5)) / 1.5;
	return vec2(uMu, uR);
}

void getTransmittanceRMu(out float r, out float mu) {
	r = gl_FragCoord.y / float(TRANSMITTANCE_H);
	mu = gl_FragCoord.x / float(TRANSMITTANCE_W);
	r = Rg + (r * r) * (Rt - Rg);
	mu = -0.15 + tan(1.5 * mu) / tan(1.5) * (1.0 + 0.15);
}

vec2 getIrradianceUV(float r, float muS) {
	float uR = (r - Rg) / (Rt - Rg);
	float uMuS = (muS + 0.2) / (1.0 + 0.2);
	return vec2(uMuS, uR);
}

void getIrradianceRMuS(out float r, out float muS) {
	r = Rg + (gl_FragCoord.y - 0.5) / (float(SKY_H) - 1.0) * (Rt - Rg);
	muS = -0.2 + (gl_FragCoord.x - 0.5) / (float(SKY_W) - 1.0) * (1.0 + 0.2);
}

vec4 texture4D(sampler3D table, float r, float mu, float muS, float nu) {
	float H = sqrt(Rt * Rt - Rg * Rg);
	float rho = sqrt(r * r - Rg * Rg);
	float rmu = r * mu;
	float delta = rmu * rmu - r * r + Rg * Rg;
	vec4 cst = rmu < 0.0 && delta > 0.0 ? vec4(1.0, 0.0, 0.0, 0.5 - 0.5 / float(RES_MU)) : vec4(-1.0, H * H, H, 0.5 + 0.5 / float(RES_MU));
	float uR = 0.5 / float(RES_R) + rho / H * (1.0 - 1.0 / float(RES_R));
	float uMu = cst.w + (rmu * cst.x + sqrt(delta + cst.y)) / (rho + cst.z) * (0.5 - 1.0 / float(RES_MU));
	// paper formula
	// float uMuS = 0.5 / float(RES_MU_S) + max((1.0 - exp(-3.0 * muS - 0.6)) / (1.0 - exp(-3.6)), 0.0) * (1.0 - 1.0 / float(RES_MU_S));
	// better formula
	float uMuS = 0.5 / float(RES_MU_S) + (atan(max(muS, -0.1975) * tan(1.26 * 1.1)) / 1.1 + (1.0 - 0.26)) * 0.5 * (1.0 - 1.0 / float(RES_MU_S));
	float lerp = (nu + 1.0) / 2.0 * (float(RES_NU) - 1.0);
	float uNu = floor(lerp);
	lerp = lerp - uNu;
	
	vec4 t3da = texture3D(table, vec3((uNu + uMuS) / float(RES_NU), uMu, uR));
	vec4 t3db = texture3D(table, vec3((uNu + uMuS + 1.0) / float(RES_NU), uMu, uR));

	// __HACK__ to stop nans coming out and breaking things
	// TODO fix this in inscatter_multiple_a, somehow...
	// this may have been caused by a nan bug in my transmittance table?
	float nan_test = dot(t3da, t3db);
	if (!(nan_test == nan_test)) return vec4(0f);

	return t3da * (1.0 - lerp) + t3db * lerp;
}

void getMuMuSNu(float r, vec4 dhdH, out float mu, out float muS, out float nu) {
	float x = gl_FragCoord.x - 0.5;
	float y = gl_FragCoord.y - 0.5;
	if (y < float(RES_MU) / 2.0) {
		float d = 1.0 - y / (float(RES_MU) / 2.0 - 1.0);
		d = min(max(dhdH.z, d * dhdH.w), dhdH.w * 0.999);
		mu = (Rg * Rg - r * r - d * d) / (2.0 * r * d);
		mu = min(mu, -sqrt(1.0 - (Rg / r) * (Rg / r)) - 0.001);
	} else {
		float d = (y - float(RES_MU) / 2.0) / (float(RES_MU) / 2.0 - 1.0);
		d = min(max(dhdH.x, d * dhdH.y), dhdH.y * 0.999);
		mu = (Rt * Rt - r * r - d * d) / (2.0 * r * d);
	}
	muS = mod(x, float(RES_MU_S)) / (float(RES_MU_S) - 1.0);
	// paper formula
	// muS = -(0.6 + log(1.0 - muS * (1.0 -  exp(-3.6)))) / 3.0;
	// better formula
	muS = tan((2.0 * muS - 1.0 + 0.26) * 1.1) / tan(1.26 * 1.1);
	nu = -1.0 + floor(x / float(RES_MU_S)) / (float(RES_NU) - 1.0) * 2.0;
}

// nearest intersection of ray r,mu with ground or top atmosphere boundary
// mu=cos(ray zenith angle at ray origin)
float limit(float r, float mu) {
	float dout = -r * mu + sqrt(r * r * (mu * mu - 1.0) + RL * RL);
	float delta2 = r * r * (mu * mu - 1.0) + Rg * Rg;
	if (delta2 >= 0.0) {
		float din = -r * mu - sqrt(delta2);
		if (din >= 0.0) {
			dout = min(dout, din);
		}
	}
	return dout;
}

// transmittance(=transparency) of atmosphere for infinite ray (r,mu)
// (mu=cos(view zenith angle)), intersections with ground ignored
vec3 transmittance(float r, float mu) {
	vec2 uv = getTransmittanceUV(r, mu);
	return texture2D(transmittanceSampler, uv).rgb;
}

// transmittance(=transparency) of atmosphere for infinite ray (r,mu)
// (mu=cos(view zenith angle)), or zero if ray intersects ground
vec3 transmittanceWithShadow(float r, float mu) {
	return mu < -sqrt(1.0 - (Rg / r) * (Rg / r)) ? vec3(0.0) : transmittance(r, mu);
}

// transmittance(=transparency) of atmosphere between x and x0
// assume segment x,x0 not intersecting ground
// r=||x||, mu=cos(zenith angle of [x,x0) ray at x), v=unit direction vector of [x,x0) ray
vec3 transmittance(float r, float mu, vec3 v, vec3 x0) {
	vec3 result;
	float r1 = length(x0);
	float mu1 = dot(x0, v) / r;
	if (mu > 0.0) {
		result = min(transmittance(r, mu) / transmittance(r1, mu1), 1.0);
	} else {
		result = min(transmittance(r1, -mu1) / transmittance(r, -mu), 1.0);
	}
	return result;
}

// optical depth for ray (r,mu) of length d, using analytic formula
// (mu=cos(view zenith angle)), intersections with ground ignored
// H=height scale of exponential density function
float opticalDepth(float H, float r, float mu, float d) {
	float a = sqrt((0.5/H)*r);
	vec2 a01 = a*vec2(mu, mu + d / r);
	vec2 a01s = sign(a01);
	vec2 a01sq = a01*a01;
	float x = a01s.y > a01s.x ? exp(a01sq.x) : 0.0;
	vec2 y = a01s / (2.3193*abs(a01) + sqrt(1.52*a01sq + 4.0)) * vec2(1.0, exp(-d/H*(d/(2.0*r)+mu)));
	return sqrt((6.2831*H)*r) * exp((Rg-r)/H) * (x + dot(y, vec2(1.0, -1.0)));
}

// transmittance(=transparency) of atmosphere for ray (r,mu) of length d
// (mu=cos(view zenith angle)), intersections with ground ignored
// uses analytic formula instead of transmittance texture
vec3 analyticTransmittance(float r, float mu, float d) {
	return exp(- betaR * opticalDepth(HR, r, mu, d) - (betaM * 1.1) * opticalDepth(HM, r, mu, d));
}

// transmittance(=transparency) of atmosphere between x and x0
// assume segment x,x0 not intersecting ground
// d = distance between x and x0, mu=cos(zenith angle of [x,x0) ray at x)
vec3 transmittance(float r, float mu, float d) {
	vec3 result;
	float r1 = sqrt(r * r + d * d + 2.0 * r * mu * d);
	float mu1 = (r * mu + d) / r1;
	if (mu > 0.0) {
		result = min(transmittance(r, mu) / transmittance(r1, mu1), 1.0);
	} else {
		result = min(transmittance(r1, -mu1) / transmittance(r, -mu), 1.0);
	}
	return result;
}

vec3 irradiance(sampler2D sampler, float r, float muS) {
	vec2 uv = getIrradianceUV(r, muS);
	return texture2D(sampler, uv).rgb;
}

// Rayleigh phase function
float phaseFunctionR(float mu) {
	return (3.0 / (16.0 * M_PI)) * (1.0 + mu * mu);
}

// Mie phase function
float phaseFunctionM(float mu) {
	return 1.5 * 1.0 / (4.0 * M_PI) * (1.0 - mieG*mieG) * pow(1.0 + (mieG*mieG) - 2.0*mieG*mu, -3.0/2.0) * (1.0 + mu * mu) / (2.0 + mieG*mieG);
}

// approximated single Mie scattering (cf. approximate Cm in paragraph "Angular precision")
vec3 getMie(vec4 rayMie) { // rayMie.rgb=C*, rayMie.w=Cm,r
	return rayMie.rgb * rayMie.w / max(rayMie.r, 1e-4) * (betaR.r / betaR);
}
