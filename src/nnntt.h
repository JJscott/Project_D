#pragma once

//mapType
#define MAP_FLAT 1
#define MAP_PERLIN 2
#define MAP_M_PERLIN 3
#define MAP_EROSION 4

#include "initial3d.h"
#include "GLee.h"


struct heightArray{
	int w, h, scale;
	double *heightMap;

	heightArray(int ws, int hs) : w(ws), h(hs), heightMap(new double[ws*hs]) {}

	inline double & y(int x, int z) { 
		return heightMap[x*h + z];
	}
	
	inline double y(int x, int z) const {
		return heightMap[x*h + z];
	}
	
	initial3d::vec3d pos(int x, int z) const;
	initial3d::vec3d norm(int x, int z) const;
	
	~heightArray() {
		delete[] heightMap;
	}
};

class NewNewNewTerrainTest{

private:

	static heightArray * flatMap(int);
	static heightArray * perlinMap(int, long);
	static heightArray * modifiedPerlinMap(int, long);
	static heightArray * hydroThermErosion(int, int, double);
	static void scale(heightArray *ha, double s);
	static void compile(GLuint, const heightArray *);

public:

	static void map(GLuint, int, double);

};


