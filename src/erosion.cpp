
//pragma once
#include <cmath> 
#include <algorithm>
#include <cstdlib>
#include <iostream> 
#include "perlin.h"
#include "erosion.h"
#include "nnntt.h"
#include "initial3d.h"

// #include <fenv.h>

#include <stdexcept>
//throw runtime_error("MOO!");


#define TIME_INCR 0.02				// ∆t Time increment [0; 0.05] 0.02
#define RAIN_RATE 0.012				// Kr Rain rate [0; 0.05] 0.012 
#define WATER_EVAP_RATE 0.015		// Ke Water evaporation rate [0; 0.05] 0.015 
#define VIRTUAL_PIPE_AREA 20.0		// A Virtual pipe cross section area [0.1; 60] 20 
#define VIRTUAL_PIPE_LENGTH 1.0		// l Virtual pipe arraySize between sections (any value)
#define GRAVITY 9.81				// g Gravity [0.1; 20] 9.81 
#define SEDIMENT_CAPACITY 1.0		// Kc Sediment capacity [0.1; 3] 1 
#define THERMAL_EROSION_RATE 0.15	// Kt Thermal erosion rate [0; 3] 0.15

#define SOIL_SUSPENSION_RATE 0.5// Ks Soil suspension rate [0.1; 2] 0.5
#define SEDEMENT_DEPO_RATE 1.0	// Kd Sediment deposition rate [0.1; 3] 1
#define SEDEMENT_SOFT_RATE 5.0	// Kh Sediment softening rate [0; 10] 5
#define MAX_EROSION_DEPTH 10.0	// Kdmax Maximal erosion depth [0; 40] 10
#define TALUS_TAN_COEFF 0.8		// Ka Talus angle tangent coeff. [0; 1] 0.8
#define TALUS_TAN_BIAS 0.1		// Ki Talus angle tangent bias [0; 1] 0.1

// +ve y-axis
//
//
//-----------------//
//  6  |  5  |  4  //
//-----------------//
//  7  |  n  |  3  //
//-----------------//
//  0  |  1  |  2  //
//-----------------// -----> +ve x-axis

//imagining x+ is right and x- is left, y+ is up, and y- is down


using namespace std;
using namespace initial3d;


ErosionSimulation::ErosionSimulation(int w, int h) {
	width = w;
	height = h;
	state = new float[w*h*4];
	outFlux = new float[w*h*4];
	velocity = new float[w*h*4];

	//initialise array
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			for (int i = 0; i < 4; i++) {
				state[y*width*4 + x*4 + i] = 0.0f;
				outFlux[y*width*4 + x*4 + i] = 0.0f;
				velocity[y*width*4 + x*4 + i] = 0.0f;
			}
		}
	}
}

ErosionSimulation::~ErosionSimulation() {
	delete[] state;
	delete[] outFlux;
	delete[] velocity;
}



void ErosionSimulation::loadHeightMap(float *hmap) {
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			ground(x,y) = hmap[y*width + x];
			hardness(x,y) = 1.0f;
		}
	}
}


float * ErosionSimulation::getHeightMap() {
	float * heightMap = new float[height*width];
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			heightMap[y*width + x] = ground(x,y);
		}
	}
	return heightMap;
}



void ErosionSimulation::update() {
	waterIncrement();
	waterFlowCalculation();
	// soilFlowCalculation();
	erosionDeposition();
	sedimentTransportation();
	// soilFlowCalculation();
	waterEvaporation();
}



// 1. Water incrementation due to rain or water sources.
void ErosionSimulation::waterIncrement() {
	cout << "1 - Water incrementation" << endl;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			water(x, y) += TIME_INCR * RAIN_RATE;
		}
	}
}


// 2. Flow simulation using shallow-water model. Computation of velocity ﬁeld and water height changes.
void ErosionSimulation::waterFlowCalculation() {
	cout << "2 - Flow simulation" << endl;

	float fluxFactor = TIME_INCR * VIRTUAL_PIPE_AREA * GRAVITY / VIRTUAL_PIPE_LENGTH;

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

			if (lFlux(x,y) < 0.0f || rFlux(x,y) < 0.0f || tFlux(x,y) < 0.0f || bFlux(x,y) < 0.0f)
				throw runtime_error("Erosion error: Negative outflux value in previous itertation.");

			float height = ground(x,y) + water(x,y);

			//left flux
			if (x != 0)
			{
				float dh = height - ground(x-1,y) + water(x-1,y);
				float flux = lFlux(x,y) + fluxFactor * dh;
				lFlux(x,y) = max(0.0f, flux);

			} else {
				lFlux(x,y) = 0.0f;
			}

			//right flux
			if (x != width-1)
			{
				float dh = height - ground(x+1,y) + water(x+1,y);
				float flux = rFlux(x,y) + fluxFactor * dh;
				rFlux(x,y) = max(0.0f, flux);
			} else {
				rFlux(x,y) = 0.0f;
			}

			//top flux
			if (y != height-1)
			{
				float dh = height - ground(x,y+1) + water(x,y+1);
				float flux = tFlux(x,y) + fluxFactor * dh;
				tFlux(x,y) = max(0.0f, flux);
			} else {
				tFlux(x,y) = 0.0f;
			}

			//bottom flux
			if (y != 0)
			{
				float dh = height - ground(x,y-1) + water(x,y-1);
				float flux = bFlux(x,y) + fluxFactor * dh;
				bFlux(x,y) = max(0.0f, flux);
			} else {
				bFlux(x,y) = 0.0f;
			}

			float sumFlux = lFlux(x,y) + rFlux(x,y) + tFlux(x,y) + bFlux(x,y);

			if (sumFlux < 0.0f)
				throw runtime_error("Erosion error: Negative outflow sum.");

			assertValid(sumFlux);

			if (sumFlux > water(x,y)) {
				// float scaleF = min(1.0f, (d1_water[index] / sumFlux));
				float scaleF = (water(x,y) * VIRTUAL_PIPE_LENGTH * VIRTUAL_PIPE_LENGTH) / (sumFlux * TIME_INCR);
				scaleF = math::clamp(scaleF, 0.0f, 1.0f);

				lFlux(x,y) *= scaleF;
				rFlux(x,y) *= scaleF;
				tFlux(x,y) *= scaleF;
				bFlux(x,y) *= scaleF;
			}

			assertValid(lFlux(x,y));
			assertValid(rFlux(x,y));
			assertValid(tFlux(x,y));
			assertValid(bFlux(x,y));


		}
	}

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

			float outFlux = lFlux(x, y) + rFlux(x, y) + tFlux(x, y) + bFlux(x, y);

			float inFlux = 0.0f;
			if (x != 0) inFlux += rFlux(x-1, y);
			if (x != width-1) inFlux += lFlux(x+1, y);
			if (y != height-1) inFlux += bFlux(x, y+1);
			if (y != 0) inFlux += tFlux(x, y-1);

			float deltaFlux = TIME_INCR * (inFlux - outFlux);
			float oldWater = water(x,y);
			water(x,y) += deltaFlux / (VIRTUAL_PIPE_LENGTH * VIRTUAL_PIPE_LENGTH);
			water(x,y) = max(water(x,y), 0.0f);

			if (water(x,y) < 0.0f)
				throw runtime_error("Erosion error: Negative water height after update.");

			float meanWater = 0.5 * (oldWater + water(x,y));

			if (meanWater == 0.0f) {
				xVel(x,y) = yVel(x,y) = 0.0f;
			} else {
				xVel(x,y) += rFlux(x,y) - lFlux(x,y);
				yVel(x,y) += tFlux(x,y) - bFlux(x,y);
			
				if (x != 0) 		xVel(x,y) += rFlux(x-1, y);
				if (x != width-1) 	xVel(x,y) -= lFlux(x+1, y);
				if (y != height-1) 	yVel(x,y) -= bFlux(x, y+1);
				if (y != 0) 		yVel(x,y) += tFlux(x, y-1);

				//scale non-edge cases
				if(x != 0 && x != width-1) xVel(x,y) *= 0.5f;
				if(y != 0 && y != height-1) yVel(x,y) *= 0.5f;
			}

			// if(xVel(x,y) < -1.0 || xVel(x,y) > 1.0){
			// 	cout << "-------------------------------" << endl;
			// 	cout << "postition "  << x << "," << y << endl;
			// 	cout << "values  "  << xVel(x,y) << "," << yVel(x,y) << endl;
			// 	cout << "sqrt" << sqrt(xVel(x,y)*xVel(x,y) + yVel(x,y)*yVel(x,y)) << endl;
			// 	if (x != 0) 		cout << "flux r " << rFlux(x-1, y) <<  endl;
			// 	if (x != width-1) 	cout << "flux l " << lFlux(x+1, y) <<  endl;
			// 	if (y != height-1) 	cout << "flux b " << bFlux(x, y+1) <<  endl;
			// 	if (y != 0) 		cout << "flux t " << tFlux(x, y-1) <<  endl;
			// 	cout << "-------------------------------" << endl;
			// }

			assertValid(xVel(x,y));
			assertValid(yVel(x,y));
		}
	}
}


// 3. Soil ﬂow calculation with outﬂow in virtual pipesof thermal erosion model.
void ErosionSimulation::soilFlowCalculation() {
	// 	/*
// 		3. Soil ﬂow calculation with outﬂow in virtual pipes of thermal erosion model.
// 	*/
// 	// ∆St+∆t = a·∆t·Kt·Rt(x,y)·H/2
// 	// α = tan((b − bi)/d)
// 	// A = {bi,b − bi < 0 ∧tan(α) > (R(x,y) * Ka +Ki),i = 1,...8}
// 	// ∆Si = ∆S * (bi / ∑∀ bk ∈ Abk)

// 	// cout << "3 - Soil ﬂow calculation" << endl;

// 	// //thermal sedminent
// 	// float d_thermal[arraySize];
// 	// bool a[8];

// 	// for (int x = 0; x < width; x++) {
// 	// 	for (int y = 0; y < height; y++) {
// 	// 		int index = x * height + y;
// 	// 		float maxHieghtDif = 0;
// 	// 		for (int i = 0; i < 8; i++) {
// 	// 			maxHieghtDif = max(maxHieghtDif, groundDif[index] i)); //get the maximum height difference
// 	// 		}
			
// 	// 		if (maxHieghtDif > 0) {
// 	// 			float thermalVolume = VIRTUAL_PIPE_AREA * TIME_INCR * THERMAL_EROSION_RATE * oldTex->hardness(x, y) * maxHieghtDif / 2;
// 	// 			float totalHeightDff = 0.0;
// 	// 			for (int i = 0; i < 8; i++) {
// 	// 				float alpha = tan(groundDif[index] i) / VIRTUAL_PIPE_LENGTH);
					
// 	// 				if (alpha < 0 && tan(alpha) > oldTex->hardness(x, y) * TALUS_TAN_COEFF + TALUS_TAN_BIAS) {
// 	// 					a[i] = true;
// 	// 					totalHeightDff += groundDif[index] i);
// 	// 				} else {
// 	// 					a[i] = false;
// 	// 				}
// 	// 			}

// 	// 			for (int i = 0; i < 8; i++) {
// 	// 				//TODO change this to the actual square...
// 	// 				if (a[i])
// 	// 					d_thermal[i] = thermalVolume * (groundDif[index] i) / totalHeightDff);
// 	// 			}

// 	// 		} else {
// 	// 			d_thermal[index] = -1;
// 	// 		}
// 	// 	}
// 	// }
}


// 4. Simulation of erosion-deposition process.
void ErosionSimulation::erosionDeposition() {
	cout << "4 - Simulation of erosion-deposition process" << endl;

	float *tempState = new float[height * width * 4];

	// for (int x = 0; x < width; x++) {
	// 	for (int y = 0; y < height; y++) {
	// 		for (int i = 0; i < 4; i++) {
	// 			tempState[y*width*4 + x*4 + i] = 0.0f;
	// 		}
	// 	}
	// }

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			// int index = x * height + y;

			float normX = 0;
			float normY = 0;
			if (x != 0) normX -= ground(x-1, y);
			else normX -= ground(x, y);
			if (x != width-1) normX += ground(x+1, y);
			else normX += ground(x, y);
			if (y != height-1) normY += ground(x, y+1);
			else normY += ground(x, y);
			if (y != 0) normY -= ground(x, y-1);
			else normY -= ground(x, y);

			vec3f normal(normX, normY, 2.0f);
			float cosa = (~normal) * vec3f::j();
			float sinAlpha = max(sin(acos(cosa)), 0.2f); //or 0.1
			float velLength = sqrt(xVel(x,y)*xVel(x,y) + yVel(x,y)*yVel(x,y));

			float capacity = SEDIMENT_CAPACITY * velLength * sinAlpha * lmax(water(x,y));
			float delta = capacity - sediment(x,y);

			if (capacity < 0.0f || math::isnan(capacity))
				throw runtime_error("Erosion error: Invalid sediment capacity.");
			if (sediment(x,y) < 0.0f)
				throw runtime_error("Erosion error: Negative sediment bef1ore update.");

			// //if the capacity larger than the sediment, pick up
			// if (delta > 0.0f) {
			// 	float val = TIME_INCR * hardness(x,y) * SOIL_SUSPENSION_RATE * delta;
			// 	ground(x,y) -= val;
			// 	sediment(x,y) += val;
			// 	water(x,y) += val;

			// //if the capacity smaller than the sediment, drop off
			// } else if (delta < 0.0f) {
			// 	float val = TIME_INCR * SEDEMENT_DEPO_RATE * delta;
			// 	ground(x,y) -= val;
			// 	sediment(x,y) += val;
			// 	water(x,y) += val;
			// }



			int index = y * width * 4 + x * 4;

			assertValid(ground(x,y));
			assertValid(delta);

			//if the capacity larger than the sediment, pick up
			if (delta > 0.0f) {
				float val = TIME_INCR * hardness(x,y) * SOIL_SUSPENSION_RATE * delta;
				tempState[index] = ground(x,y) - val;
				tempState[index + 2] = sediment(x,y) + val;
				tempState[index + 1] = water(x,y) + val;

			//if the capacity smaller than the sediment, drop off
			} else if (delta < 0.0f) {
				float val = TIME_INCR * SEDEMENT_DEPO_RATE * delta;
				tempState[index] = ground(x,y) - val;
				tempState[index + 2] = sediment(x,y) + val;
				tempState[index + 1] = water(x,y) + val;
			}

			tempState[index + 3] = hardness(x,y);

			assertValid(tempState[index]);
			assertValid(tempState[index + 1]);
			assertValid(tempState[index + 2]);
			assertValid(tempState[index + 3]);

			//modifiy hardness HERERERERERERERE
		}
	}
	delete[] state;
	state = tempState;
}


// 5. Transportation of suspended sediment by the velocityﬁeld.
void ErosionSimulation::sedimentTransportation() {
	cout << "5 - Transportation of suspended sediment" << endl;

	float* sedimentTemp = new float[width * height];
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int index = y * width + x;

			float fromPosX = x - xVel(x,y) * TIME_INCR;
			float fromPosY = y - yVel(x,y) * TIME_INCR;

			// integer coordinates
			int x0 = floor(fromPosX);
			int y0 = floor(fromPosY);
			int x1 = x0+1;
			int y1 = y0+1;

			// interpolation factors
			float fX = fromPosX - x0;
			float fY = fromPosY - y0;

			// clamp to grid borders
			x0 = math::clamp(x0, 0, width-1);
			x1 = math::clamp(x1, 0, width-1);
			y0 = math::clamp(y0, 0, height-1);
			y1 = math::clamp(y1, 0, height-1);

			sedimentTemp[index] = math::lerp<float, float>( 
									math::lerp<float, float>(sediment(x0,y0), sediment(x1,y0),fX), 
									math::lerp<float, float>(sediment(x0,y1), sediment(x1,y1),fX), fY);

			assertValid(sedimentTemp[index]);
		}
	}

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			sediment(x,y) = sedimentTemp[y * width + x];
		}
	}
}


// 6. Thermal erosion material amount calculation.
void ErosionSimulation::thermalErosion() {
	
// 	// cout << "6 - Thermal erosion material amount calculation" << endl;

// 	// for (int x = 0; x < width; x++)
// 	// 	for (int y = 0; y < height; y++)
// 	// 		newTex->ground(x, y) += d_thermal[x * height + y];
}


// 7. Water evaporation.
void ErosionSimulation::waterEvaporation() {
	cout << "7 - Water evaporation" << endl;

	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			water(x,y) = water(x,y) * (1 - WATER_EVAP_RATE * TIME_INCR);

}

//MAX_EROSION_DEPTH
float lmax (float x){
	if (x <= 0) return 0;
	if (x >= MAX_EROSION_DEPTH) return 1;
	return 1-(MAX_EROSION_DEPTH - x) / MAX_EROSION_DEPTH;
}

void assertValid(float v) {
	if (math::isnan(v))
		throw runtime_error("ERROR :: NaN float detected");
	if (math::isinf(v))
		throw runtime_error("ERROR :: Infinite float detected");
}


			// if (x == X_CH && y == Y_CH) {
			// 	cout << "X=10, Y=10" << endl;
			// 	cout << "capacity " << capacity << endl;
			// 	cout << "delta " << delta << endl;
			// 	cout << "ground " << ground(x,y) << endl;
			// 	cout << "water " << water(x,y) << endl;
			// 	cout << "sediment " << sediment(x,y) << endl << endl;

			// 	cout << "normal " << normal << endl;
			// 	cout << "sinAlpha " << sinAlpha << endl << endl;


			// 	cout << "olFlux " << lFlux(x,y) << endl;
			// 	cout << "orFlux " << rFlux(x,y) << endl;
			// 	cout << "otFlux " << tFlux(x,y) << endl;
			// 	cout << "obFlux " << bFlux(x,y) << endl << endl;

			// 	cout << "ilFlux " << rFlux(x-1,y) << endl;
			// 	cout << "irFlux " << lFlux(x+1,y) << endl;
			// 	cout << "itFlux " << bFlux(x+1,y) << endl;
			// 	cout << "ibFlux " << tFlux(x-1,y) << endl << endl;

			// 	cout << "xVel " << xVel(x,y) << endl;
			// 	cout << "yVel " << yVel(x,y) << endl << endl;

			// 	cout << "value" << TIME_INCR * hardness(x,y) * SOIL_SUSPENSION_RATE * delta << endl;
			// 	cout << "-------------------------" << endl;
			// }