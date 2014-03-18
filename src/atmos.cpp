/*
 * Atmospheric Scattering precomputation
 *
 */

#include "GLee.h"
#include <GL/glu.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <cassert>

#include "atmos.h"
#include "shader.h"
#include "hrtime.h"
#include "initial3d.h"

using namespace std;
using namespace initial3d;

namespace atmos {
	
	ShaderManager *shaderman = NULL;
	
	// planet radius
	float Rg = 6360000;
	float Rt = 6420000;
	
	// dritab texture
	GLuint tex_dritab = 0;

	void check_gl() {
		GLenum err = glGetError();
		if (err != GL_NO_ERROR) {
			printf("%s\n", gluErrorString(err));
			throw runtime_error("GL!!!");
		}
	}

	void check_fbo() {
		if (glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
			cerr << "YOU BROKE THE FRAMEBUFFER!" << endl;
			throw runtime_error("FBO!!!");
		}
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

	void setCommon(GLuint prog) {
		glUniform1f(glGetUniformLocation(prog, "Rg"), Rg);
		glUniform1f(glGetUniformLocation(prog, "Rt"), Rt);
	}
	
	void load_dritab() {
		ifstream ifs("./dritab4");
		if (ifs.is_open()) {
			float *data = new float[1024 * 1024];
			
			string version;
			getline(ifs, version);
			if (version != "version_1004") {
				cout << "wrong dritab version" << endl;
				return;
			}
			
			cout << "loading dritab from file..." << endl;
			
			for (size_t c = 0; c < 1024; c++) {
				for (size_t r = 0; r < 1024; r++) {
					float *p = data + r * 1024 + c;
					string token;
					ifs >> token;
					if (token == "Infinity") {
						*p = math::inf<float>();
					} else if (token == "NaN") {
						*p = math::nan<float>();
					} else {
						istringstream iss(token);
						iss >> *p;
					}
				}
			}
			
			glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1024, 1024, 0, GL_RED, GL_FLOAT, data);
			
			delete[] data;
			
			cout << "dritab loaded" << endl;
			
		} else {
			cout << "failed to open dritab file" << endl;
		}
		
	}
	
	void make_dritab() {
		if (tex_dritab != 0) glDeleteTextures(1, &tex_dritab);
		glGenTextures(1, &tex_dritab);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, tex_dritab);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1024, 1024, 0, GL_RGB, GL_FLOAT, nullptr);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_dritab, 0);
		check_fbo();
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glViewport(0, 0, 1024, 1024);
		GLuint prog = shaderman->getProgram("precomp.vert;dritab.frag");
		glUseProgram(prog);
		glUniform1i(glGetUniformLocation(prog, "RES_H0"), 1024);
		glUniform1i(glGetUniformLocation(prog, "RES_MU"), 1024);
		draw_fullscreen_quad();
		glUseProgram(0);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
		glFinish();

		// load precomputed table
		// load_dritab();
		
		glFinish();
		check_gl();
	}

	void makeTables() {
		// init shader manager
		if (shaderman == NULL) shaderman = new ShaderManager("./shader/atmos");
		// make a framebuffer
		GLuint fbo;
		glGenFramebuffers(1, &fbo);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
		// make the tables
		glDisable(GL_DEPTH_TEST);
		make_dritab();
		// clean up
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glDeleteFramebuffers(1, &fbo);
	}
	
}
