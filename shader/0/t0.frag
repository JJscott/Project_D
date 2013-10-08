#version 120

varying vec4 p;
varying vec3 n;

varying float logz;

void main()
{
	gl_FragDepth = logz;
	vec4 l = gl_LightSource[0].position;
	l -= l.w * p;
	// ambient
	vec4 c = gl_FrontMaterial.ambient * gl_LightSource[0].ambient;
	// diffuse
	float dot_d = dot(normalize(l.xyz), normalize(n));
	c += gl_FrontMaterial.diffuse * clamp(dot_d, 0.0, 1.0) * gl_LightSource[0].diffuse;
	// specular
	vec3 v = normalize(-p.xyz);
	vec3 lu = normalize(l.xyz);
	vec3 r = normalize((2.0 * dot(lu, n) * n) - lu);
	if (dot_d > 0.0) {
		c += gl_FrontMaterial.specular * pow(clamp(dot(r, v), 0.0, 1.0), gl_FrontMaterial.shininess) * gl_LightSource[0].specular;
	}
	
	gl_FragColor = c;
}
