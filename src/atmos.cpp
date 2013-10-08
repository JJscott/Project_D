/*
 * Atmospheric Scattering precomputation
 *
 */

#include "GLee.h"
#include <GL/glu.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <cstdlib>

#include "atmos.h"
#include "shader.h"
#include "hrtime.h"
#include "initial3d.h"

using namespace std;
using namespace initial3d;

namespace atmos {
	
	ShaderManager *shaderman = NULL;
	
	// scattering constants
	color betaR(5.47e-6, 1.28e-5, 3.12e-5); // normal
	//color betaR(3.12e-5, 1.28e-5, 5.47e-6); // red sky, blue sunset
	color betaM(0.000021, 0.000021, 0.000021);
	float h0R  = 7994;
	float h0M  = 1200;
	// somewhere between radii of 600000 and 6000000 something breaks;
	// blending in copy_inscatter_multiple breaks (for high altitudes)
	// a bug somewhere in inscatter_multiple_a outputs nans, hackfixed in texture4D
	float Rg   = 6360000;
	float Rt   = 6420000;
	float RL   = 6421000;
	float mieG = 0.8;

	// texture units to bind during precomputation
	enum { TU_TRANSMITTANCE = 0, TU_DELTA_SM, TU_DELTA_SR, TU_DELTA_J, TU_DELTA_E };

	// lookup table textures
	GLuint tex_transmittance = 0;
	GLuint tex_inscatter = 0;
	GLuint tex_irradiance = 0;

	void check_gl() {
		GLenum err = glGetError();
		if (err != GL_NO_ERROR) {
			printf("%s\n", gluErrorString(err));
			exit(1);
		}
	}

	void check_fbo() {
		if (glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
			cerr << "YOU BROKE THE FRAMEBUFFER!" << endl;
			exit(1);
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
		// transmittance table size
		glUniform1i(glGetUniformLocation(prog, "TRANSMITTANCE_W"), RES_TRANSMITTANCE_MU);
		glUniform1i(glGetUniformLocation(prog, "TRANSMITTANCE_H"), RES_TRANSMITTANCE_H);
		// irradiance table size
		glUniform1i(glGetUniformLocation(prog, "SKY_W"), RES_IRRADIANCE_SX);
		glUniform1i(glGetUniformLocation(prog, "SKY_H"), RES_IRRADIANCE_H);
		// inscatter table size
		glUniform1i(glGetUniformLocation(prog, "RES_R"), RES_INSCATTER_H);
		glUniform1i(glGetUniformLocation(prog, "RES_MU"), RES_INSCATTER_VX);
		glUniform1i(glGetUniformLocation(prog, "RES_MU_S"), RES_INSCATTER_SX);
		glUniform1i(glGetUniformLocation(prog, "RES_NU"), RES_INSCATTER_VS);
		// constants
		glUniform3fv(glGetUniformLocation(prog, "betaR"), 1, betaR);
		glUniform3fv(glGetUniformLocation(prog, "betaM"), 1, betaM);
		glUniform1f(glGetUniformLocation(prog, "HR"), h0R);
		glUniform1f(glGetUniformLocation(prog, "HM"), h0M);
		glUniform1f(glGetUniformLocation(prog, "Rg"), Rg);
		glUniform1f(glGetUniformLocation(prog, "Rt"), Rt);
		glUniform1f(glGetUniformLocation(prog, "RL"), RL);
		glUniform1f(glGetUniformLocation(prog, "mieG"), mieG);
		// transmittance texture
		glUniform1i(glGetUniformLocation(prog, "transmittanceSampler"), TU_TRANSMITTANCE);
	}

	void set_layer(GLuint prog, int layer) {
		// im getting sick of all this black magic
		double r = layer / (RES_INSCATTER_H - 1.0);
		r = r * r;
		r = sqrt(Rg * Rg + r * (Rt * Rt - Rg * Rg)) + (layer == 0 ? 10 : (layer == RES_INSCATTER_H - 1 ? -1 : 0.0));
		double dmin = Rt - r;
		double dmax = sqrt(r * r - Rg * Rg) + sqrt(Rt * Rt - Rg * Rg);
		double dminp = r - Rg;
		double dmaxp = sqrt(r * r - Rg * Rg);
		glUniform1f(glGetUniformLocation(prog, "r"), float(r));
		glUniform4f(glGetUniformLocation(prog, "dhdH"), float(dmin), float(dmax), float(dminp), float(dmaxp));
		glUniform1i(glGetUniformLocation(prog, "layer"), layer);
	}

	void test_tex2d(GLuint tex) {
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glViewport(0, 0, 512, 512);
		glUseProgram(0);
		glDisable(GL_LIGHTING);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, tex);
		glEnable(GL_TEXTURE_2D);
		glColor3d(1, 1, 1);
		draw_fullscreen_quad();
		glFinish();
		//glutSwapBuffers();
		exit(0);
	}

	void test_tex3d(GLuint tex) {
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glViewport(0, 0, 512, 512);
		glUseProgram(0);
		glDisable(GL_LIGHTING);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, tex);
		glEnable(GL_TEXTURE_3D);
		glColor3d(1, 1, 1);
		draw_fullscreen_quad();
		glFinish();
		//glutSwapBuffers();
		exit(0);
	}

	void make_transmittance() {
		cout << "Making transmittance lookup table..." << endl;

		// shader
		GLuint prog = shaderman->getProgram("precomp.vert;transmittance.frag");
		glUseProgram(prog);

		// delete existing texture
		if (tex_transmittance != 0) glDeleteTextures(1, &tex_transmittance);

		// make texture
		glGenTextures(1, &tex_transmittance);
		glBindTexture(GL_TEXTURE_2D, tex_transmittance);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, RES_TRANSMITTANCE_MU, RES_TRANSMITTANCE_H, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		// attach texture
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_transmittance, 0);
		check_fbo();

		// draw to the texture
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glViewport(0, 0, RES_TRANSMITTANCE_MU, RES_TRANSMITTANCE_H);
		setCommon(prog);
		draw_fullscreen_quad();
		glFinish();
		check_gl();
		cout << "Making transmittance lookup table done." << endl;

		//test_tex2d(tex_transmittance);
	}

	GLuint make_inscatter_texture() {
		// sliced by h, subdivided by vs
		GLuint tex;
		glGenTextures(1, &tex);
		glBindTexture(GL_TEXTURE_3D, tex);
		glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA16F, RES_INSCATTER_SX * RES_INSCATTER_VS, RES_INSCATTER_VX, RES_INSCATTER_H, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
		return tex;
	}

	GLuint make_irradiance_texture() {
		GLuint tex;
		glGenTextures(1, &tex);
		glBindTexture(GL_TEXTURE_2D, tex);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, RES_IRRADIANCE_SX, RES_IRRADIANCE_H, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		return tex;
	}

	void do_irradiance_single(GLuint tex_deltaE) {
		cout << "Irradiance single..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;irradiance_single.frag");
		glUseProgram(prog);
		setCommon(prog);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_deltaE, 0);
		check_fbo();
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glViewport(0, 0, RES_IRRADIANCE_SX, RES_IRRADIANCE_H);
		draw_fullscreen_quad();
		check_gl();
	}

	void do_inscatter_single(GLuint tex_deltaSR, GLuint tex_deltaSM) {
		cout << "Inscatter single..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;inscatter_single.frag");
		glUseProgram(prog);
		setCommon(prog);
		GLuint bufs[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
		glDrawBuffers(2, bufs);
		for (int sample_h = 0; sample_h < RES_INSCATTER_H; sample_h++) {
			glFramebufferTexture3D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, tex_deltaSR, 0, sample_h);
			glFramebufferTexture3D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_3D, tex_deltaSM, 0, sample_h);
			glViewport(0, 0, RES_INSCATTER_VS * RES_INSCATTER_SX, RES_INSCATTER_VX);
			check_fbo();
			set_layer(prog, sample_h);
			draw_fullscreen_quad();
			check_gl();
		}
		// detach color1 so it doesnt interfere (by limiting viewport size)
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, 0, 0);

		// test_tex3d(tex_deltaSR);
	}

	void do_inscatter_multiple_a(bool first, GLuint tex_deltaJ) {
		cout << "Inscatter multiple A..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;inscatter_multiple_a.frag");
		glUseProgram(prog);
		setCommon(prog);
		glUniform1i(glGetUniformLocation(prog, "deltaESampler"), TU_DELTA_E);
		glUniform1i(glGetUniformLocation(prog, "deltaSRSampler"), TU_DELTA_SR);
		glUniform1i(glGetUniformLocation(prog, "deltaSMSampler"), TU_DELTA_SM);
		glUniform1i(glGetUniformLocation(prog, "first"), first);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		for (int sample_h = 0; sample_h < RES_INSCATTER_H; sample_h++) {
			glFramebufferTexture3D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, tex_deltaJ, 0, sample_h);
			glViewport(0, 0, RES_INSCATTER_VS * RES_INSCATTER_SX, RES_INSCATTER_VX);
			check_fbo();
			set_layer(prog, sample_h);
			draw_fullscreen_quad();
			check_gl();
		}

		// test_tex3d(tex_deltaJ);
	}

	void do_irradiance_multiple(bool first, GLuint tex_deltaE) {
		cout << "Irradiance multiple..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;irradiance_multiple.frag");
		glUseProgram(prog);
		setCommon(prog);
		glUniform1i(glGetUniformLocation(prog, "deltaSRSampler"), TU_DELTA_SR);
		glUniform1i(glGetUniformLocation(prog, "deltaSMSampler"), TU_DELTA_SM);
		glUniform1i(glGetUniformLocation(prog, "first"), first);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_deltaE, 0);
		check_fbo();
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glViewport(0, 0, RES_IRRADIANCE_SX, RES_IRRADIANCE_H);
		draw_fullscreen_quad();
		check_gl();
	}

	void do_inscatter_multiple_b(GLuint tex_deltaSR) {
		cout << "Inscatter multiple B..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;inscatter_multiple_b.frag");
		glUseProgram(prog);
		setCommon(prog);
		glUniform1i(glGetUniformLocation(prog, "deltaJSampler"), TU_DELTA_J);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		for (int sample_h = 0; sample_h < RES_INSCATTER_H; sample_h++) {
			glFramebufferTexture3D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, tex_deltaSR, 0, sample_h);
			glViewport(0, 0, RES_INSCATTER_VS * RES_INSCATTER_SX, RES_INSCATTER_VX);
			check_fbo();
			set_layer(prog, sample_h);
			draw_fullscreen_quad();
			check_gl();
		}

		// test_tex3d(tex_deltaSR);
	}

	void copy_inscatter_single() {
		cout << "Copy inscatter single..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;copy_inscatter_single.frag");
		glUseProgram(prog);
		setCommon(prog);
		glUniform1i(glGetUniformLocation(prog, "deltaSRSampler"), TU_DELTA_SR);
		glUniform1i(glGetUniformLocation(prog, "deltaSMSampler"), TU_DELTA_SM);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		for (int sample_h = 0; sample_h < RES_INSCATTER_H; sample_h++) {
			glFramebufferTexture3D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, tex_inscatter, 0, sample_h);
			glViewport(0, 0, RES_INSCATTER_VS * RES_INSCATTER_SX, RES_INSCATTER_VX);
			check_fbo();
			set_layer(prog, sample_h);
			draw_fullscreen_quad();
		}
		check_gl();
	}

	void copy_inscatter_multiple() {
		cout << "Copy inscatter multiple..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;copy_inscatter_multiple.frag");
		glUseProgram(prog);
		setCommon(prog);
		glUniform1i(glGetUniformLocation(prog, "deltaSSampler"), TU_DELTA_SR);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glEnable(GL_BLEND);
		glBlendEquation(GL_FUNC_ADD);
		glBlendFunc(GL_ONE, GL_ONE);
		for (int sample_h = 0; sample_h < RES_INSCATTER_H; sample_h++) {
			glFramebufferTexture3D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, tex_inscatter, 0, sample_h);
			glViewport(0, 0, RES_INSCATTER_VS * RES_INSCATTER_SX, RES_INSCATTER_VX);
			check_fbo();
			set_layer(prog, sample_h);
			draw_fullscreen_quad();
		}
		glDisable(GL_BLEND);
		check_gl();
	}

	void copy_irradiance() {
		cout << "Copy irradiance..." << endl;
		GLuint prog = shaderman->getProgram("precomp.vert;copy_irradiance.frag");
		glUseProgram(prog);
		setCommon(prog);
		glUniform1i(glGetUniformLocation(prog, "deltaESampler"), TU_DELTA_E);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_irradiance, 0);
		check_fbo();
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glViewport(0, 0, RES_IRRADIANCE_SX, RES_IRRADIANCE_H);
		glEnable(GL_BLEND);
		glBlendEquation(GL_FUNC_ADD);
		glBlendFunc(GL_ONE, GL_ONE);
		draw_fullscreen_quad();
		glDisable(GL_BLEND);
		check_gl();
	}

	void make_inscatter_irradiance() {
		cout << "Making inscatter and irradiance lookup tables..." << endl;

		// delete existing textures
		if (tex_inscatter != 0) glDeleteTextures(1, &tex_inscatter);
		if (tex_irradiance != 0) glDeleteTextures(1, &tex_irradiance);

		// make and bind textures

		// make these first because they dont need to be bound atm
		tex_inscatter = make_inscatter_texture();
		tex_irradiance = make_irradiance_texture();

		// clear the irradiance texture
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_irradiance, 0);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glClearColor(0, 0, 0, 0);
		glClear(GL_COLOR_BUFFER_BIT);

		// bind the transmittance texture
		glActiveTexture(GL_TEXTURE0 + TU_TRANSMITTANCE);
		glBindTexture(GL_TEXTURE_2D, tex_transmittance);
		// make the temp textures, bound to the right unit
		glActiveTexture(GL_TEXTURE0 + TU_DELTA_SM);
		GLuint tex_deltaSM = make_inscatter_texture();
		glActiveTexture(GL_TEXTURE0 + TU_DELTA_SR);
		GLuint tex_deltaSR = make_inscatter_texture();
		glActiveTexture(GL_TEXTURE0 + TU_DELTA_J);
		GLuint tex_deltaJ = make_inscatter_texture();
		glActiveTexture(GL_TEXTURE0 + TU_DELTA_E);
		GLuint tex_deltaE = make_irradiance_texture();

		// do single scattering / irradiance first
		do_irradiance_single(tex_deltaE);
		do_inscatter_single(tex_deltaSR, tex_deltaSM);
		// we dont copy irradiance here (its only direct sunlight)
		// this step in Bruneton's source just clears the irradiance texture
		// (the cheeky ******* - that confused me for a while)
		copy_inscatter_single();
		glFinish();
	
		//test_tex3d(tex_inscatter);
	
		// now multiple scattering
		for (int order = 2; order <= ATMOS_SCATTER_ORDERS; order++) {
			cout << "Multiple scattering order " << order << "..." << endl;
			// this _is_ the right order
			do_inscatter_multiple_a(order == 2, tex_deltaJ);
			do_irradiance_multiple(order == 2, tex_deltaE);
			do_inscatter_multiple_b(tex_deltaSR);
			copy_irradiance();
			copy_inscatter_multiple();
			glFinish();
		}

		//test_tex3d(tex_inscatter);

		// delete temp textures
		glDeleteTextures(1, &tex_deltaSM);
		glDeleteTextures(1, &tex_deltaSR);
		glDeleteTextures(1, &tex_deltaJ);
		glDeleteTextures(1, &tex_deltaE);
		glActiveTexture(GL_TEXTURE0);
	
		glFinish();
	
		check_gl();

		cout << "Making inscatter and irradiance lookup tables done." << endl;
	}
	
	void makeTables() {
		// init shader manager
		if (shaderman == NULL) shaderman = new ShaderManager("./shader/atmos");
		hrtime_t time_start;
		getHighResTime(&time_start);
		cout << "Making lookup tables..." << endl;
		// make a framebuffer
		GLuint fbo;
		glGenFramebuffers(1, &fbo);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
		// make the tables
		glDisable(GL_DEPTH_TEST);
		make_transmittance();
		make_inscatter_irradiance();
		// clean up
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glDeleteFramebuffers(1, &fbo);
		hrtime_t time_end;
		getHighResTime(&time_end);
		cout << "Making lookup tables done in " << highResTimeToSec(&time_start, &time_end) << " seconds." << endl;
	}
	
}
