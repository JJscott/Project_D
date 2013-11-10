
#include <cmath> 
#include <algorithm>
#include <cstdlib>
#include "nnntt.h"

float lmax (float);
void assertValid(float);

class ErosionSimulation {

private:

	int width;
	int height;
	float* state;
	float* outFlux;
	float* velocity;

	inline float & ground (int x, int y)         { return state[y*width*4 + x*4 + 0]; }
	inline float   ground (int x, int y)   const { return state[y*width*4 + x*4 + 0]; }
	inline float & water (int x, int y)          { return state[y*width*4 + x*4 + 1]; }
	inline float   water (int x, int y)    const { return state[y*width*4 + x*4 + 1]; }
	inline float & sediment (int x, int y)       { return state[y*width*4 + x*4 + 2]; }
	inline float   sediment (int x, int y) const { return state[y*width*4 + x*4 + 2]; }
	inline float & hardness (int x, int y)       { return state[y*width*4 + x*4 + 3]; }
	inline float   hardness (int x, int y) const { return state[y*width*4 + x*4 + 3]; }

	inline float & lFlux (int x, int y)       { return outFlux[y*width*4 + x*4 + 0]; }
	inline float   lFlux (int x, int y) const { return outFlux[y*width*4 + x*4 + 0]; }
	inline float & rFlux (int x, int y)       { return outFlux[y*width*4 + x*4 + 1]; }
	inline float   rFlux (int x, int y) const { return outFlux[y*width*4 + x*4 + 1]; }
	inline float & tFlux (int x, int y)       { return outFlux[y*width*4 + x*4 + 2]; }
	inline float   tFlux (int x, int y) const { return outFlux[y*width*4 + x*4 + 2]; }
	inline float & bFlux (int x, int y)       { return outFlux[y*width*4 + x*4 + 3]; }
	inline float   bFlux (int x, int y) const { return outFlux[y*width*4 + x*4 + 3]; }

	inline float & xVel (int x, int y)       { return velocity[y*width*4 + x*4 + 0]; }
	inline float   xVel (int x, int y) const { return velocity[y*width*4 + x*4 + 0]; }
	inline float & yVel (int x, int y)       { return velocity[y*width*4 + x*4 + 1]; }
	inline float   yVel (int x, int y) const { return velocity[y*width*4 + x*4 + 1]; }

public:

	ErosionSimulation(int, int);
	~ErosionSimulation();
	
	void loadHeightMap(float*);
	float * getHeightMap();
	void update();

	void waterIncrement();			// 1. Water incrementation due to rain or water sources.
	void waterFlowCalculation();	// 2. Flow simulation using shallow-water model. Computation of velocity ﬁeld and water height changes.
	void soilFlowCalculation();		// 3. Soil ﬂow calculation with outﬂow in virtual pipesof thermal erosion model.
	void erosionDeposition();		// 4. Simulation of erosion-deposition process.
	void sedimentTransportation();	// 5. Transportation of suspended sediment by the velocityﬁeld.
	void thermalErosion();		// 6. Thermal erosion material amount calculation.
	void waterEvaporation();		// 7. Water evaporation.
};