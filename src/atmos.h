/*
 * Atmospheric Scattering header
 *
 */

#ifndef ATMOS_H
#define ATMOS_H

#include "GLee.h"
#include "initial3d.h"

namespace atmos {

	extern float Rg;
	extern float Rt;

	// lookup table textures
	extern GLuint tex_dritab;

	// function to (re)make the tables
	void makeTables();
	
	// function to set 'common' uniforms in atmos shaders (that include 'atmos/common.frag')
	// this sets the transmittance sampler as used by precomputation, re-set as necessary
	void setCommon(GLuint prog);
	
}

#endif // ATMOS_H
