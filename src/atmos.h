/*
 * Atmospheric Scattering precomputation header
 *
 */

#ifndef ATMOS_H
#define ATMOS_H

#include "GLee.h"
#include "initial3d.h"

// lookup table sizes
#define RES_TRANSMITTANCE_MU	1024
#define RES_TRANSMITTANCE_H		512
#define RES_INSCATTER_VS		8
#define RES_INSCATTER_SX		32
#define RES_INSCATTER_VX		128
#define RES_INSCATTER_H			32
#define RES_IRRADIANCE_SX		1024
#define RES_IRRADIANCE_H		512

// number of orders of multiple scattering
#define ATMOS_SCATTER_ORDERS	4

namespace atmos {

	// scattering constants
	extern initial3d::color betaR;
	extern initial3d::color betaM;
	extern float h0R;
	extern float h0M;
	extern float Rg;
	extern float Rt;
	extern float RL;
	extern float mieG;

	// lookup table textures
	extern GLuint tex_transmittance;
	extern GLuint tex_inscatter;
	extern GLuint tex_irradiance;

	// function to (re)make the tables
	void makeTables();
	
	// function to set 'common' uniforms in atmos shaders (that include 'atmos/common.frag')
	// this sets the transmittance sampler as used by precomputation, re-set as necessary
	void setCommon(GLuint prog);
	
}

#endif // ATMOS_H
