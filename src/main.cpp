
#include "GLee.h"
#include <GL/glu.h>
#include <GLFW/glfw3.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

#include "glfw_helper.h"
#include "initial3d.h"
#include "hrtime.h"
#include "shader.h"
#include "atmos.h"
#include "camera.h"

using namespace std;
using namespace initial3d;

using atmos::Rg;

// forward declarations
void display();
void reshape(GLFWwindow *, int, int);
void keyboard(GLFWwindow *, unsigned int);

// shaders
ShaderManager *shaderman = NULL;
GLuint prog_scene_space;
GLuint prog_scene_sea;
GLuint prog_scene_test;
GLuint prog_deferred;

// fps counter variables
int fps;
hrtime_t fps_start;

// window / projection
GLFWwindow *window;
int win_width = 512, win_height = 512;
double fov_y = 30;
float zfar = 20000000;

float exposure = 1;

// deferred shading fbo + textures
GLuint fbo_scene = 0;
GLuint rbuf_scene_depth = 0; // depth
GLuint tex_scene_pos = 0;
GLuint tex_scene_norm = 0;
GLuint tex_scene_diffuse = 0;

// sun
vec3d sun = vec3d::j();

// camera
Camera *camera;

void initShaders() {
	if (shaderman == NULL) shaderman = new ShaderManager("./shader");
	shaderman->unloadAll();
	
	prog_scene_space = shaderman->getProgram("scene.vert;scene_space.frag");
	prog_scene_sea = shaderman->getProgram("scene.vert;scene_sea.frag");
	prog_scene_test = shaderman->getProgram("scene.vert;scene_test.frag");
	
	prog_deferred = shaderman->getProgram("deferred.vert;deferred.frag");
}

class main_event_handler : public glfwpp::KeyListener {
public:
	virtual void handleKey(GLFWwindow *window, int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS) {
			if (key == GLFW_KEY_F12) initShaders();
		}
	}
	
	virtual void handleChar(GLFWwindow *window, unsigned int c) {
		//if (c == 'w') cam_pos.z() -= Rg / 100;
		//if (c == 'W') cam_pos.z() -= Rg / 10000;
		//if (c == 's') cam_pos.z() += Rg / 100;
		//if (c == 'S') cam_pos.z() += Rg / 10000;
		if (c == 'r') sun = quatd::axisangle(vec3d::i(), math::pi() / 180) * sun;
		if (c == 'f') sun = quatd::axisangle(vec3d::i(), -math::pi() / 180) * sun;
		if (c == '=') exposure *= 1.2;
		if (c == '-') exposure /= 1.2;
		//cout << +cam_pos << endl;
	}
};

void wait_exit(int code) {
	string line;
	getline(cin, line);
	exit(code);
}

void do_fps() {
	hrtime_t now;
	getHighResTime(&now);
	if (highResTimeToSec(&fps_start, &now) < 1) {
		fps++;
	} else {
		char strbuf[42];
		sprintf(strbuf, "COMP308 x1, now with GLFW (%d FPS)", fps);
		//glfwSetWindowTitle(window, strbuf);
		fps = 0;
		getHighResTime(&fps_start);
	}
}

void check_gl() {
	GLenum err = glGetError();
	if (err != GL_NO_ERROR) {
		printf("%s\n", gluErrorString(err));
		throw runtime_error("BOOM!");
	}
}

void check_fbo() {
	if (glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		cerr << "YOU BROKE THE FRAMEBUFFER!" << endl;
		throw runtime_error("OH NOES! THE FRAMEBUFFER IS BROKED");
	}
}

void init_fbo_scene(int w, int h) {

	// make new framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	if (fbo_scene) glDeleteFramebuffers(1, &fbo_scene);
	glGenFramebuffers(1, &fbo_scene);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo_scene);
	
	// delete existing textures
	if (rbuf_scene_depth) glDeleteRenderbuffers(1, &rbuf_scene_depth);
	if (tex_scene_pos) glDeleteTextures(1, &tex_scene_pos);
	if (tex_scene_norm) glDeleteTextures(1, &tex_scene_norm);
	if (tex_scene_diffuse) glDeleteTextures(1, &tex_scene_diffuse);
	
	// create new texture handles
	glGenRenderbuffers(1, &rbuf_scene_depth);
	glGenTextures(1, &tex_scene_pos);
	glGenTextures(1, &tex_scene_norm);
	glGenTextures(1, &tex_scene_diffuse);	
	
	// depth
	glBindRenderbuffer(GL_RENDERBUFFER, rbuf_scene_depth);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32, w, h);
	glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbuf_scene_depth);
	
	// position (could reconstruct from depth, potential precision problems...)
	glBindTexture(GL_TEXTURE_2D, tex_scene_pos);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, w, h, 0, GL_RGB, GL_FLOAT, NULL);
	glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_scene_pos, 0);
	
	// normal
	glBindTexture(GL_TEXTURE_2D, tex_scene_norm);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, NULL);
	glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, tex_scene_norm, 0);
	
	// diffuse
	glBindTexture(GL_TEXTURE_2D, tex_scene_diffuse);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_FLOAT, NULL);
	glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, tex_scene_diffuse, 0);
	
	check_fbo();
}

int main(int argc, char *argv[]) {

	if (!glfwInit()) exit(1);

	glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 1);
	//glfwWindowHint(GLFW_SAMPLES, 16);
	window = glfwCreateWindow(win_width, win_height, "DAVE", NULL, NULL);
	glfwMakeContextCurrent(window);

	glfwpp::setCallbacks(window);
	glfwpp::addKeyListener(window, new main_event_handler());

	glfwSetWindowSizeCallback(window, reshape);

	if (!GLeeInit()) {
		cerr << "GLee init failed: " << GLeeGetErrorString() << endl;
		wait_exit(1);
	}

	if (!GLEE_VERSION_2_1) {
		cerr << "GL 2.1 not available." << endl;
		wait_exit(1);
	}

	if (!GLEE_ARB_framebuffer_object) {
		cerr << "ARB_framebuffer_object extension not available." << endl;
		wait_exit(1);
	}

	if (!GLEE_ARB_texture_float) {
		cerr << "ARB_texture_float extension not available." << endl;
		wait_exit(1);
	}

	if (!GLEE_ARB_draw_buffers) {
		cerr << "ARB_draw_buffers extension not available." << endl;
		wait_exit(1);
	}

	cout << "GLee init (GL 2.1) success." << endl;

	initShaders();
	
	// init camera
	camera = new FPSCamera(window, vec3d(0, Rg + 100, 0));

	fps = 0;
	getHighResTime(&fps_start);

	atmos::makeTables();

	reshape(window, win_width, win_height);

	while (!glfwWindowShouldClose(window)) {
		display();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();

	return 0;
}


void draw_fullscreen_quad() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	double p = 0.0;

	glBegin(GL_POLYGON);
	glTexCoord3d(0, 1, p);
	glVertex3d(-1, 1, 0);
	glTexCoord3d(0, 0, p);
	glVertex3d(-1, -1, 0);
	glTexCoord3d(1, 0, p);
	glVertex3d(1, -1, 0);
	glTexCoord3d(1, 1, p);
	glVertex3d(1, 1, 0);
	glEnd();
}

void test_tex2d(GLuint tex) {
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glViewport(0, 0, win_width, win_height);
	glUseProgram(0);
	glDisable(GL_LIGHTING);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex);
	glEnable(GL_TEXTURE_2D);
	glColor3d(1, 1, 1);
	draw_fullscreen_quad();
	glFinish();
	glfwSwapBuffers(window);
	wait_exit(0);
}

void display() {
	
	camera->update();

	//
	// FIRST PASS : assemble scene buffer
	//
	
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo_scene);
	// pos, norm, diffuse
	GLenum bufs_scene[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
	glDrawBuffers(3, bufs_scene);
	glViewport(0, 0, win_width, win_height);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	check_fbo();

	// setup projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov_y, double(win_width) / double(win_height), 1, zfar);

	// setup view matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	mat4d view = camera->getViewMatrix();

	// begin drawing
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	// setup terrain modelview matrix
	mat4d modelview_terrain = view * mat4d::translate(0, Rg, 0);
	glPushMatrix();
	glMultMatrixd(~modelview_terrain);

	// draw terrain here
	
	glUseProgram(prog_scene_test);
	glUniform1f(glGetUniformLocation(prog_scene_test, "far"), zfar);
	
	GLUquadric *quad = gluNewQuadric();
	gluSphere(quad, 10000, 100, 100);
	gluDeleteQuadric(quad);
	
	glPopMatrix();

	// now any world space stuff
	glPushMatrix();
	glMultMatrixd(~view);

	// draw here

	glPopMatrix();

	// pseudo-skybox
	glUseProgram(prog_scene_space);
	glUniform1f(glGetUniformLocation(prog_scene_space, "far"), zfar);
	glBegin(GL_POLYGON);
	glNormal3d(0, 0, 1);
	glVertex3d(-2 * zfar, -2 * zfar, -0.95 * zfar);
	glVertex3d(2 * zfar, -2 * zfar, -0.95 * zfar);
	glVertex3d(2 * zfar, 2 * zfar, -0.95 * zfar);
	glVertex3d(-2 * zfar, 2 * zfar, -0.95 * zfar);
	glEnd();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glFinish();
	check_gl();

	//
	// SECOND PASS : shadows (not yet)
	//
	
	glFinish();
	check_gl();

	//
	// THIRD PASS : deferred shading
	//
	
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glDrawBuffer(GL_BACK);
	glViewport(0, 0, win_width, win_height);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	// bind textures
	enum { TU_TRANSMITTANCE = 0, TU_INSCATTER, TU_IRRADIANCE, TU_POSITION, TU_NORMAL, TU_DIFFUSE };
	glActiveTexture(GL_TEXTURE0 + TU_TRANSMITTANCE);
	glBindTexture(GL_TEXTURE_2D, atmos::tex_transmittance);
	glActiveTexture(GL_TEXTURE0 + TU_INSCATTER);
	glBindTexture(GL_TEXTURE_3D, atmos::tex_inscatter);
	glActiveTexture(GL_TEXTURE0 + TU_IRRADIANCE);
	glBindTexture(GL_TEXTURE_2D, atmos::tex_irradiance);
	glActiveTexture(GL_TEXTURE0 + TU_POSITION);
	glBindTexture(GL_TEXTURE_2D, tex_scene_pos);
	glActiveTexture(GL_TEXTURE0 + TU_NORMAL);
	glBindTexture(GL_TEXTURE_2D, tex_scene_norm);
	glActiveTexture(GL_TEXTURE0 + TU_DIFFUSE);
	glBindTexture(GL_TEXTURE_2D, tex_scene_diffuse);

	// shader config
	glUseProgram(prog_deferred);
	atmos::setCommon(prog_deferred);
	//glUniform1f(glGetUniformLocation(prog_deferred, "far"), zfar);
	//glUniformMatrix4fv(glGetUniformLocation(prog_deferred, "view_inv"), 1, true, mat4f(!view));
	glUniform1i(glGetUniformLocation(prog_deferred, "transmittanceSampler"), TU_TRANSMITTANCE);
	glUniform1i(glGetUniformLocation(prog_deferred, "sampler_inscatter"), TU_INSCATTER);
	glUniform1i(glGetUniformLocation(prog_deferred, "sampler_irradiance"), TU_IRRADIANCE);
	glUniform1i(glGetUniformLocation(prog_deferred, "sampler_position"), TU_POSITION);
	glUniform1i(glGetUniformLocation(prog_deferred, "sampler_normal"), TU_NORMAL);
	glUniform1i(glGetUniformLocation(prog_deferred, "sampler_diffuse"), TU_DIFFUSE);
	glUniform1f(glGetUniformLocation(prog_deferred, "exposure"), exposure);
	glUniform3fv(glGetUniformLocation(prog_deferred, "sunnorm_v"), 1, (view * vec4d(sun, 0)).xyz<float>());
	glUniform3fv(glGetUniformLocation(prog_deferred, "planetpos_v"), 1, (view * vec3d::zero()).xyz<float>());

	draw_fullscreen_quad();

	// unbind textures
	glActiveTexture(GL_TEXTURE0 + TU_TRANSMITTANCE);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0 + TU_INSCATTER);
	glBindTexture(GL_TEXTURE_3D, 0);
	glActiveTexture(GL_TEXTURE0 + TU_IRRADIANCE);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0 + TU_POSITION);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0 + TU_NORMAL);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0 + TU_DIFFUSE);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0);

	glFinish();
	check_gl();

	// clean up
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glUseProgram(0);
	glFinish();
	do_fps();
}

void reshape(GLFWwindow *window, int w, int h) {
	cout << "RESHAPE" << endl;
	if (h == 0) h = 1;
	win_width = w;
	win_height = h;
	init_fbo_scene(win_width, win_height);
}

void test_tex3d(GLuint tex) {
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glViewport(0, 0, win_width, win_height);
	glUseProgram(0);
	glDisable(GL_LIGHTING);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_3D, tex);
	glEnable(GL_TEXTURE_3D);
	glColor3d(1, 1, 1);
	draw_fullscreen_quad();
	glFinish();
	//glutSwapBuffers();
	wait_exit(0);
}

