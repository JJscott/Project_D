//perlin.h
#pragma once

#define GRADIENT_TABLE_SIZE 256

#include "initial3d.h"


inline double pi() {
	return 3.1415926535897932384626433832795028841971693993751058209749445923078;
}

class perlin {

private:
	initial3d::random m_random;
	double m_gradients[GRADIENT_TABLE_SIZE*3];

	double smooth(double);
	double lattice(int, int, int, double, double, double);
	int getIndex(int, int, int);
	int permutate(int);
	double lerp(double, double, double);

public:
	perlin();
	perlin(long);
	~perlin();

	double getNoise(double, double, double, int);
	double getNoise(double, double, double);
};