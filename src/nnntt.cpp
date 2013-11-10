#include "perlin.h"
#include "nnntt.h"
#include "erosion.h"
#include "terrainsim.h"
#include "GLee.h"
#include "initial3d.h"

#undef min
#undef max

using namespace initial3d;

double m_size;

vec3d heightArray::pos(int x, int z) const { 
	return vec3d(x-(w/2.0), y(x,z), z-(h/2.0));
}


vec3d heightArray::norm(int x, int z) const { 
	//norm calc here
	try {
		vec3d norm;
		if (x > 0) norm   += ~(vec3d(y(x,z)-y(x-1,z), 1, 0));
		if (z > 0) norm   += ~(vec3d(0, 1, y(x,z)-y(x,z-1)));
		if (x < w-1) norm += ~(vec3d(y(x,z)-y(x+1,z), 1, 0));
		if (z < h-1) norm += ~(vec3d(0, 1, y(x,z)-y(x,z+1)));
		return ~norm;
	} catch (nan_error &e) {
		return vec3d::j();
	}
}


void NewNewNewTerrainTest::map(GLuint id, int size, double height){
	//std::cout << "Started map creation " << size << "x" << size << std::endl;
	//int t = glutGet(GLUT_ELAPSED_TIME);
	heightArray *ha;
	// size = size * 2; //for tapering the edges
	// switch(mapType){
	// case MAP_FLAT:
	// 	ha = flatMap(size); break;
	// case MAP_PERLIN:
	// 	ha = perlinMap(size, seed); break;
	// case MAP_M_PERLIN:
	// 	ha = modifiedPerlinMap(size, seed); break;
	// case MAP_EROSION:
		ha = hydroThermErosion(size, size, height);//erosion e;
		//ha = e.hydroThermErosion(size, size, height); 
	// 	break;
	// }
	ha->scale = height;
	//scale(ha, height);
	compile(id, ha);
	//std::cout << "Completed map creation int " <<glutGet(GLUT_ELAPSED_TIME)-t << "ms" << std::endl;
}

heightArray * NewNewNewTerrainTest::flatMap(int size){
	heightArray *ha = new heightArray(size, size);
	for (int x = 0; x < ha->w; x++)
		for (int z = 0; z < ha->h; z++)
			ha->y(x,z) = 0;
	return ha;
}

heightArray * NewNewNewTerrainTest::perlinMap(int size, long seed){
	heightArray *ha = new heightArray(size, size);
	perlin p(seed);
	for (int x = 0; x < ha->w; x++)
		for (int z = 0; z < ha->h; z++)
			ha->y(x,z) = p.getNoise(x/(double)ha->w, z/(double)ha->h, 0, 2);
	return ha;
}

heightArray * NewNewNewTerrainTest::modifiedPerlinMap(int size, long seed){
	heightArray *ha = new heightArray(size, size);
	perlin p1(seed);
	perlin p2(seed+1);
	for (int x = 0; x < ha->w; x++)
		for (int z = 0; z < ha->h; z++)
			ha->y(x,z) = p1.getNoise(x/(double)ha->w, z/(double)ha->h, 0, 3) * 
						p2.getNoise(x/(double)ha->w, z/(double)ha->h, 0, 3) * 2;
	return ha;
}


heightArray * NewNewNewTerrainTest::hydroThermErosion(int width, int height, double scale) {
	// ErosionSimulation *es = new ErosionSimulation(width, height);
	// float *heightMap = new float[width * height];
	// perlin p(32);

	// float midX = width/2;
	// float midY = height/2;
	// float maxDis = sqrt(midX * midX + midY * midY);
	// for (int y = 0; y < height; y++) {
	// 	for (int x = 0; x < width; x++) {
	// 		// float island = (maxDis - sqrt((midX-x)*(midX-x)+(midY-y)*(midY-y))) / maxDis;
	// 		heightMap[y*height + x] = 0.5 + p.getNoise(x/(double)width, y/(double)height, 0, 3);
	// 	}
	// }

	// es->loadHeightMap(heightMap);

	// for (int i = 0; i < 500; i++){
	// 	std::cout << std::endl << "Iteration  " << i << std::endl;
	// 	es->update();
	// }


	// float *heightMap2 = es->getHeightMap();
	// heightArray* ha = new heightArray(width, height);

	// float diff = 0;
	// for (int x = 0; x < ha->w; x++)
	// 	for (int z = 0; z < ha->h; z++) {
	// 		diff += abs(heightMap[x*ha->h+z] - heightMap2[x*ha->h+z]);
	// 		ha->y(x, z) = scale * heightMap2[x*ha->h+z];
	// 	}

	// std::cout << "DIFERENCE ::  " << diff/(width*height) << std::endl;


	heightArray *ha = new heightArray(width, height);
	perlin p1(1);
	perlin p2(2);
	initial3d::random r(4);

	double midX = width/2;
	double midZ = height/2;
	float maxDis = sqrt(midX * midX + midZ * midZ);
	for (int x = 0; x < ha->w; x++)
		for (int z = 0; z < ha->h; z++){
			double difX = abs(x-midX);
			double difZ = abs(z-midZ);
			double island = (1 - sqrt(difX*difX+difZ*difZ)/maxDis) * 0.4;

			ha->y(x,z) = island + 0.4 * (ha->y(x,z) = p1.getNoise(x/(double)ha->w, z/(double)ha->h, 0, 8) * 
						p2.getNoise(x/(double)ha->w, z/(double)ha->h, 0, 16) * 2);
			ha->y(x,z) += 0.05 * r.nextDouble(); 
		}

	generateErosionMap(ha->heightMap, width, height, scale, 500);

	for (int x = 0; x < ha->w; x++){
		for (int z = 0; z < ha->h; z++) {
			ha->y(x, z) = math::max(-5.0, ha->y(x, z)-scale*0.25);
		}
	}


	return ha;
}




void NewNewNewTerrainTest::scale(heightArray *ha, double s){
	for (int x = 0; x < ha->w; x++)
		for (int z = 0; z < ha->h; z++)
			ha->y(x,z) *= s;
}

void NewNewNewTerrainTest::compile(GLuint id, const heightArray *ha){
	glNewList(id, GL_COMPILE);
		glBegin(GL_TRIANGLES);
		for (int x = 0; x < ha->w - 1; x++) {
			for (int z = 0; z < ha->h - 1; z++) {
				//triangle 1
				glColor3f(1.0, ha->y(x,z)/(ha->scale/2) + 0.25, ha->y(x,z)/ha->scale + 0.5);
				glNormal3dv(ha->norm(x,z));
				glVertex3dv(ha->pos(x,z));
				
				glColor3f(1.0, ha->y(x,z)/(ha->scale/2) + 0.25, ha->y(x,z)/ha->scale + 0.5);
				glNormal3dv(ha->norm(x,z+1));
				glVertex3dv(ha->pos(x,z+1));

				glColor3f(1.0, ha->y(x,z)/(ha->scale/2) + 0.25, ha->y(x,z)/ha->scale + 0.5);
				glNormal3dv(ha->norm(x+1,z));
				glVertex3dv(ha->pos(x+1,z));

				//triangle 2
				glColor3f(1.0, ha->y(x,z)/(ha->scale/2) + 0.25, ha->y(x,z)/ha->scale + 0.5);
				glNormal3dv(ha->norm(x+1,z+1));
				glVertex3dv(ha->pos(x+1,z+1));

				glColor3f(1.0, ha->y(x,z)/(ha->scale/2) + 0.25, ha->y(x,z)/ha->scale + 0.5);
				glNormal3dv(ha->norm(x+1,z));
				glVertex3dv(ha->pos(x+1,z));

				glColor3f(1.0, ha->y(x,z)/(ha->scale/2) + 0.25, ha->y(x,z)/ha->scale + 0.5);
				glNormal3dv(ha->norm(x,z+1));
				glVertex3dv(ha->pos(x,z+1));
			}
		}
		glEnd();
	glEndList();
}
