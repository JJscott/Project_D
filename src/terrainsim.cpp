
#include <iostream>
#include <math.h>
#include <cstring>
#include <sstream>

#include "initial3d.h"
#include "terrainsim.h"

using initial3d::math::max;
using initial3d::math::min;

using namespace std;
using namespace initial3d;


struct GGen_OutflowValues{
    double top;
    double left;
    double right;
    double bottom;

	GGen_OutflowValues() : top(0), left(0), right(0), bottom(0) { }
};

struct GGen_Vec3{
    double x, y, z;

    GGen_Vec3() : x(0), y(0), z(0) { }
};

class ErosionSimulator {

private:
    //variable data
    int width;
    int height;
    int length;

    //constant data
    double deltaT;
    double pipeLength;
    double pipeCrossectionArea;
    double graviationalAcceleration;
    double sedimentCapacityConstant;
    double dissolvingConstant;
    double depositionConstant;
    double minimumComputedSurfaceTilt;
    double talusAngle;

public:
    ErosionSimulator(int width, int height);

    void ApplyWaterSources(double* waterMap, double waterAmount);
    void ApplyEvaporation(double* waterMap);
    void ApplyFlowSimulation(double* heightMap, double* waterMap, GGen_OutflowValues* outflowFluxMap, GGen_Vec3* velocityVectorMap);
    void ApplyErosion(double* heightMap, double* waterMap, GGen_Vec3* velocityVectorMap, double* sedimentMap);
    void ApplyThermalWeathering(double* heightMap, double powerMultiplier);

    //helper methods
    int getLength();

};




void generateErosionMap(double *heightMap, int width, int height, int scale, int iterations)
{
    // if(width < 2 || height < 2 || scale < 1 || iterations < 1)
    //     throw runtime_error("ERROR : invalid values");

    // ErosionSimulator * simulator = new ErosionSimulator(width, height);

    // //declare other values
    // // heightMap
    // double *waterMap = new double[simulator->getLength()];
    // double *sedimentMap  = new double[simulator->getLength()];
    // GGen_OutflowValues *outflowFluxMap = new GGen_OutflowValues[simulator->getLength()];
    // GGen_Vec3 *velocityVectorMap = new GGen_Vec3[simulator->getLength()];


    // for (int i = 0; i < iterations; i++){
    //     simulator->ApplyWaterSources(waterMap, 0.012);
    //     simulator->ApplyFlowSimulation(heightMap, waterMap, outflowFluxMap, velocityVectorMap);
    //     simulator->ApplyErosion(heightMap, waterMap, velocityVectorMap, sedimentMap);
    //     simulator->ApplyThermalWeathering(heightMap, 1.0);
    //     simulator->ApplyEvaporation(waterMap);
    // }


    if(width < 2 || height < 2 || scale < 1 || iterations < 1)
        throw runtime_error("ERROR : invalid values");

    // cout << "Begin the END!!!!!!" << endl;

    ErosionSimulator * simulator = new ErosionSimulator(width, height);

    //declare other values
    // heightMap
    double *waterMap = new double[simulator->getLength()];
    double *sedimentMap  = new double[simulator->getLength()];
    GGen_OutflowValues *outflowFluxMap = new GGen_OutflowValues[simulator->getLength()];
    GGen_Vec3 *velocityVectorMap = new GGen_Vec3[simulator->getLength()];

	memset(waterMap, 0, simulator->getLength() * sizeof(double));
	memset(sedimentMap, 0, simulator->getLength() * sizeof(double));

    for (int i = 0; i < simulator->getLength(); i++){
        heightMap[i] *= scale;
    }


    for (int i = 0; i < iterations; i++){
        // cout << "iteration : " << i << endl;
        simulator->ApplyWaterSources(waterMap, 0.012);
        simulator->ApplyFlowSimulation(heightMap, waterMap, outflowFluxMap, velocityVectorMap);
        simulator->ApplyErosion(heightMap, waterMap, velocityVectorMap, sedimentMap);
        simulator->ApplyThermalWeathering(heightMap, 4.0);
        simulator->ApplyEvaporation(waterMap);

        for(int y = 0; y < height; y++){
            for(int x = 0; x < width; x++){
                if ((waterMap[y*width + x] + heightMap[y*width + x]) < 0.0) {
                    waterMap[y*width + x] = 0.0 - heightMap[y*width + x];
                }
            }
        }
    }
}









ErosionSimulator::ErosionSimulator(int width, int height):width(width), height(height), length(width * height)
{
    this->deltaT = 0.005;
    this->pipeLength = 1;
    this->pipeCrossectionArea = 20;
    this->graviationalAcceleration = 9.7;
    this->sedimentCapacityConstant = 1;
    this->dissolvingConstant = 0.5;
    this->depositionConstant = 1;
    this->minimumComputedSurfaceTilt = 0.1;
    this->talusAngle = 0.5;
}


int ErosionSimulator::getLength()
{
    return length;
}

// struct GGen_ErosionPointData{
//         GGen_Height height;
//         double waterAmount;
//         double sedimentAmount;
//         double vectorX;
//         double vectorY;
// };


//Notes
//is this thing even needed?
//NEED?
//to clean up
//TODO



// double* ErosionSimulator::ImportHeightMap( double * heightMap)
// {
//         double* returnData = new double[this->length];

//         for(GGen_Index i = 0; i < this->length; i++){
//         if(heightMap.data[i] < 0) {
//             //cout << "E";
//         }

//                 returnData[i] = (double) heightMap.data[i] * 100. / (double) GGEN_MAX_HEIGHT;
//         }

//         return returnData;
// }

// double ErosionSimulator::ExportHeightMap( double* heightMap, GGen_Data_2D& ggenHeightMap )
// {
//         double max = 0;
//         double min = 1000000000;

//         for(GGen_Index i = 0; i < this->length; i++){
//                 if(heightMap[i] > max) max  = heightMap[i];
//                 if(heightMap[i] < min) min  = heightMap[i];
//         }

//         //cout << "Min << " << min << "  ---- MAX: " << max << endl;

//         for(GGen_Index i = 0; i < this->length; i++){
//                 ggenHeightMap.data[i] = (GGen_Height) (heightMap[i] * GGEN_MAX_HEIGHT / max);

//         if(heightMap[i] <= -100.){
//             ggenHeightMap.data[i] = GGEN_MIN_HEIGHT;
//         }
//         }

//     return GGEN_MAX_HEIGHT / max;
// }



    /*
        1. Water incrementation due to rain or water sources.
    */

void ErosionSimulator::ApplyWaterSources(double* waterMap, double waterAmount)
{
    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < this->width; x++){
            waterMap[x + y * this->width] += waterAmount * this->deltaT;
        }
    }
}

    /*
        7. Water evaporation.
    */

void ErosionSimulator::ApplyEvaporation( double* waterMap)
{
    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < this->width; x++){                    
            waterMap[x + y * this->width] *= 0.985;
        }
    }
}

    /*
        2. Flow simulation using shallow-water model. Computation of velocity ﬁeld and water height changes.
    */

void ErosionSimulator::ApplyFlowSimulation(double* heightMap, double* waterMap, GGen_OutflowValues* outflowFluxMap, GGen_Vec3* velocityVectorMap)
{
    double* waterMapCopy  = new double[this->length];
    memcpy(waterMapCopy, waterMap, this->length * sizeof(double));

    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < this->width; x++){
            int currentIndex = x + this->width * y;
            double currentHeight = heightMap[currentIndex];
                    
            if(outflowFluxMap[currentIndex].left < 0 || outflowFluxMap[currentIndex].right < 0 || outflowFluxMap[currentIndex].top < 0 || outflowFluxMap[currentIndex].bottom < 0){
                throw runtime_error("Erosion error: Negative outflux value.");
            }

            // Calculate outflow values for individual directions.
            if(x > 0)
            {
                double heightDifference = currentHeight + waterMapCopy[currentIndex] - heightMap[currentIndex - 1] - waterMapCopy[currentIndex - 1];
                double flowValue = outflowFluxMap[currentIndex].left + this->deltaT * this->pipeCrossectionArea * (this->graviationalAcceleration * heightDifference) / this->pipeLength;
                outflowFluxMap[currentIndex].left = max(0.0, flowValue);
            }

            if(x + 1 < this->width)
            {
                double heightDifference = currentHeight + waterMapCopy[currentIndex] - heightMap[currentIndex + 1] - waterMapCopy[currentIndex + 1];
                double flowValue = outflowFluxMap[currentIndex].right + this->deltaT * this->pipeCrossectionArea * (this->graviationalAcceleration * heightDifference) / this->pipeLength;
                outflowFluxMap[currentIndex].right = max(0.0, flowValue);
            }

            if(y > 0)
            {
                double heightDifference = currentHeight + waterMapCopy[currentIndex] - heightMap[currentIndex - this->width] - waterMapCopy[currentIndex - this->width];
                double flowValue = outflowFluxMap[currentIndex].top + this->deltaT * this->pipeCrossectionArea * (this->graviationalAcceleration * heightDifference) / this->pipeLength;
                outflowFluxMap[currentIndex].top = max(0.0, flowValue);
            }

            if(y + 1 < this->height)
            {
                double heightDifference = currentHeight + waterMapCopy[currentIndex] - heightMap[currentIndex + this->width] - waterMapCopy[currentIndex + this->width];
                double flowValue = outflowFluxMap[currentIndex].bottom + this->deltaT * this->pipeCrossectionArea * (this->graviationalAcceleration * heightDifference) / this->pipeLength;
                outflowFluxMap[currentIndex].bottom = max(0.0, flowValue);
            }

            // Scale the outflow values so sum(outflow) < amount of water in this tile.
            double sumOutflow = 
                    outflowFluxMap[currentIndex].left + 
                    outflowFluxMap[currentIndex].right + 
                    outflowFluxMap[currentIndex].top + 
                    outflowFluxMap[currentIndex].bottom;

            if((int)sumOutflow < 0)
                    throw runtime_error("Erosion error: Negative outflow sum.");


            if(sumOutflow > waterMap[currentIndex]){
                double factor = min(1.0, waterMap[currentIndex] * this->pipeLength * this->pipeLength / (sumOutflow * this->deltaT));

                outflowFluxMap[currentIndex].left *= factor;
                outflowFluxMap[currentIndex].right *= factor;
                outflowFluxMap[currentIndex].top *= factor;
                outflowFluxMap[currentIndex].bottom *= factor;

                sumOutflow = //NEED?
                        outflowFluxMap[currentIndex].left + 
                        outflowFluxMap[currentIndex].right + 
                        outflowFluxMap[currentIndex].top + 
                        outflowFluxMap[currentIndex].bottom;
            }
        }               
    }

    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < this->width; x++){
            int currentIndex = x + this->width * y;

            // Update water level using outflow and inflow information from surrounding cells.
            double sumOutflow = 
                    outflowFluxMap[currentIndex].left + 
                    outflowFluxMap[currentIndex].right + 
                    outflowFluxMap[currentIndex].top + 
                    outflowFluxMap[currentIndex].bottom;

            double sumInflow = 0;
            if(x > 0)           sumInflow += outflowFluxMap[currentIndex - 1].right;
            if(x + 1 < this->width) sumInflow += outflowFluxMap[currentIndex + 1].left;
            if(y > 0)           sumInflow += outflowFluxMap[currentIndex - this->width].bottom;
            if(y + 1 < this->height)sumInflow += outflowFluxMap[currentIndex + this->width].top;
            
            double oldWaterLevel =  waterMapCopy[currentIndex];
            double newWaterLevel = waterMapCopy[currentIndex] + (this->deltaT * (sumInflow - sumOutflow)) / (this->pipeLength * this->pipeLength);

            waterMap[currentIndex] = max(0.0, newWaterLevel);

            if(velocityVectorMap == NULL) {
                continue;
            }

            double waterAverage = (oldWaterLevel + waterMap[currentIndex]) / 2;

            if(waterMap[currentIndex] < 0){
                    throw runtime_error("Erosion error: Negative water height after update.");
            }

            // The velocity field must be updated

            // Horizontal (x-axis) velocity field vector component
            if(x == 0) // Left border.
            {
                velocityVectorMap[currentIndex].x = 
                        outflowFluxMap[x + this->width * y].right -
                        outflowFluxMap[x + 1 + this->width * y].left;
            }
            else if(x + 1 == this->width) // Right border.
            {
                velocityVectorMap[currentIndex].x = 
                        outflowFluxMap[x - 1 + this->width * y].right -
                        outflowFluxMap[x + this->width * y].left;
            }
            else // The rest.
            {
                velocityVectorMap[currentIndex].x = 
                    (
                        outflowFluxMap[x - 1 + this->width * y].right -
                        outflowFluxMap[x + this->width * y].left +
                        outflowFluxMap[x + this->width * y].right -
                        outflowFluxMap[x + 1 + this->width * y].left
                    ) / 2;
            }

            // Vertical (y-axis) velocity field vector component
            if(y == 0) // Top border.
            {
                velocityVectorMap[currentIndex].y = 
                        outflowFluxMap[x + this->width * y].bottom -
                        outflowFluxMap[x + 1 + this->width * (y + 1)].top;
            }
            else if(y + 1 == this->height) // Bottom border.
            {
                velocityVectorMap[currentIndex].y = 
                        outflowFluxMap[x + this->width * (y - 1)].bottom -
                        outflowFluxMap[x + this->width * y].top;
            }
            else // The rest.
            {
                velocityVectorMap[currentIndex].y = 
                    (
                        outflowFluxMap[x + this->width * (y - 1)].bottom -
                        outflowFluxMap[x + this->width * y].top +
                        outflowFluxMap[x + this->width * y].bottom -
                        outflowFluxMap[x + this->width * (y + 1)].top
                    ) / 2;
            }

            //TODO NEED?
        ///velocityVectorMap[x + width * y].x /= waterAverage * this->pipeLength;
        //velocityVectorMap[x + width * y].y /= waterAverage * this->pipeLength;

            // if(velocityVectorMap[x + width * y].x * this->deltaT > 1 || velocityVectorMap[x + width * y].y * this->deltaT > 1){
            //     throw runtime_error("Erosion error: Too long velocity vector.");
            // }
        }
    }

    delete [] waterMapCopy;
}


void ErosionSimulator::ApplyErosion( double* heightMap, double* waterMap, GGen_Vec3* velocityVectorMap, double* sedimentMap)
{

    /*
        4. Simulation of erosion-deposition process.
    */

    double* sedimentToMoveMap = new double[this->length];

    // Update sediment amount carried by water in the current tile  
    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < this->width; x++){
            int currentIndex = x + this->width * y;

            sedimentToMoveMap[currentIndex] = 0;

            GGen_Vec3 currentVelocityVector = velocityVectorMap[currentIndex];

            //currentVelocityVector.x *= -1;
                        //currentVelocityVector.y *= -1;

            int coordinateOffsetX = currentVelocityVector.x < 0 ? -1 : 0;
            int coordinateOffsetY = currentVelocityVector.y < 0 ? -1 : 0;
            int baseCoordinateX = x + coordinateOffsetX;
            int baseCoordinateY = y + coordinateOffsetY;
            double coordinatePartX = currentVelocityVector.x * this->deltaT - coordinateOffsetX;
            double coordinatePartY = currentVelocityVector.y * this->deltaT - coordinateOffsetY;        


            if(coordinatePartX < 0){
                    coordinatePartX += 1;
            }

            if(coordinatePartY < 0){
                    coordinatePartY += 1;
            }

            // Calculate current height of the target point (it will likely have to be interpolated from 4 surrounding points).
            
            // Top left point.
            double heightSum = 0;
            double heightWeightSum = 0;
            if(
                baseCoordinateX < this->width &&
                x >= -coordinateOffsetX &&
                baseCoordinateY < this->height &&
                y >= -coordinateOffsetY
            ){
                double currentWeight = (1 - coordinatePartX) * (1 - coordinatePartY);
                heightSum += heightMap[baseCoordinateX + this->width * baseCoordinateY] * currentWeight;
                heightWeightSum += currentWeight;
            }

            // Top right point.
            if(
                baseCoordinateX + 1 < this->width &&
                x + 1 >= -coordinateOffsetX &&
                baseCoordinateY < this->height &&
                y >= -coordinateOffsetY
            ){
                double currentWeight = coordinatePartX * (1 - coordinatePartY);
                heightSum += heightMap[baseCoordinateX + 1 + this->width * baseCoordinateY] * currentWeight;
                heightWeightSum += currentWeight;
            }

            // Bottom left point.
            if(
                baseCoordinateX < this->width &&
                x >= -coordinateOffsetX &&
                baseCoordinateY + 1 < this->height &&
                y + 1 >= -coordinateOffsetY
            ){
                double currentWeight = (1 - coordinatePartX) * coordinatePartY;
                heightSum += heightMap[baseCoordinateX + this->width * (baseCoordinateY + 1)] * currentWeight;
                heightWeightSum += currentWeight;
            }

            // Bottom right point.
            if(
                baseCoordinateX + 1 < this->width &&
                x + 1 >= -coordinateOffsetX &&
                baseCoordinateY + 1 < this->height &&
                y + 1 >= -coordinateOffsetY
            ){
                double currentWeight = coordinatePartX * coordinatePartY;
                heightSum += heightMap[baseCoordinateX + 1 + this->width * (baseCoordinateY + 1)] * currentWeight;
                heightWeightSum += currentWeight;
            }

            if(heightWeightSum == 0){
                    throw runtime_error("Erosion error: No height for interpolation.");
            }

            //TODO NEED?
            double targetPointHeightDiff = heightSum / heightWeightSum - heightMap[currentIndex];
            //double targetPointDistance = sqrt(coordinatePartX * coordinatePartX + coordinatePartY * coordinatePartY);

            double velocityVectorLength = sqrt(velocityVectorMap[currentIndex].x * velocityVectorMap[currentIndex].x + velocityVectorMap[currentIndex].y * velocityVectorMap[currentIndex].y);// * this->deltaT;            

            double surfaceTilt;
        
            //NEED?
            //TODO: co s timhle if(targetPointHeightDiff < 0) continue;

            if(velocityVectorLength > 0){
                surfaceTilt = atan2(targetPointHeightDiff, velocityVectorLength);//1 - sin(this->GetSurfaceTilt(heightMap, x, y));
            }
            else {
                surfaceTilt = 0;
            }

            surfaceTilt = max(0.2, surfaceTilt);

            double sedimentCapacity = this->sedimentCapacityConstant * surfaceTilt * velocityVectorLength;
            double sedimentToMove = 0;

            if(sedimentCapacity < 0 || sedimentCapacity != sedimentCapacity)
                throw runtime_error("Erosion error: Invalid sediment capacity.");


            if(sedimentMap[currentIndex] == 0 &&  sedimentCapacity <= 0) {
                    sedimentToMoveMap[currentIndex] = 0;
                    continue;
            };

            double originalHeight = heightMap[currentIndex];

            if(sedimentMap[x + this->width * y] < 0)
                    throw runtime_error("Erosion error: Negative sediment bef1ore update.");

            double sedimentChange;
            if(sedimentCapacity > sedimentMap[currentIndex])
            {
                // The water can carry more sediment - some sediment will be picked up
                sedimentChange = this->dissolvingConstant * this->deltaT * (sedimentCapacity - sedimentMap[currentIndex]);

                //TODO NEED?
                // if(sedimentChange > 1 || sedimentChange != sedimentChange){
                        //throw runtime_error("Erosion error: Excessive negative sediment change.");
                // }
                heightMap[currentIndex] -= sedimentChange;
                waterMap[currentIndex] += sedimentChange;
                sedimentToMove = sedimentMap[currentIndex] + sedimentChange;
            }
            else if(sedimentCapacity < sedimentMap[currentIndex])
            {
                // The water is over saturated with sediment - some sediment will be deposited
                sedimentChange = this->depositionConstant * this->deltaT * (sedimentMap[currentIndex] - sedimentCapacity);

                //TODO NEED?
                // if(sedimentChange > 1 || sedimentChange != sedimentChange){
                        //throw runtime_error("Erosion error: Excessive positive sediment change.");
                // }

                heightMap[currentIndex] += sedimentChange;
                waterMap[currentIndex] -= sedimentChange;                
                if(waterMap[currentIndex] < 0){
                    waterMap[currentIndex] = 0;
                }

                sedimentToMove = sedimentMap[currentIndex] - sedimentChange;
            }

            if(sedimentToMove < 0){
                    throw runtime_error("Erosion error: Negative sediment to move.");
            }

            sedimentToMoveMap[currentIndex] = sedimentToMove;
        }
    }


    /*
        5. Transportation of suspended sediment by the velocity ﬁeld.
    */


    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < this->width; x++){
            int currentIndex = x + this->width * y;
            GGen_Vec3 currentVelocityVector = velocityVectorMap[currentIndex];

            // Move the sediment according to the velocity field map
            if(sedimentToMoveMap[currentIndex] == 0) continue;                      

            // we are doing a step backwards in time, so inverse vector has to be used
            currentVelocityVector.x *= -1; // TODO: Vyzkouset to bez toho NEED?
            currentVelocityVector.y *= -1;

            int coordinateOffsetX = currentVelocityVector.x < 0 ? -1 : 0;
            int coordinateOffsetY = currentVelocityVector.y < 0 ? -1 : 0;
            int baseCoordinateX = x + coordinateOffsetX;
            int baseCoordinateY = y + coordinateOffsetY;
            double coordinatePartX = currentVelocityVector.x * this->deltaT - coordinateOffsetX;
            double coordinatePartY = currentVelocityVector.y * this->deltaT - coordinateOffsetY;        

            if(coordinatePartX > 1 || coordinatePartY > 1 || coordinatePartX < 0 || coordinatePartY < 0){
                    cout << "E";
            }

            // Calculate current height of the target point (it will likely to be interpolated from 4 surrounding points).

            // Top left point.
            double sedimentSum = 0;
            double sedimentWeightSum = 0;
            if(
                    baseCoordinateX < this->width &&
                    x >= -coordinateOffsetX &&
                    baseCoordinateY < this->height &&
                    y >= -coordinateOffsetY
                    ){
                            double currentWeight = (1 - coordinatePartX) * (1 - coordinatePartY);
                            sedimentSum += sedimentToMoveMap[baseCoordinateX + this->width * baseCoordinateY] * currentWeight;
                            sedimentWeightSum += currentWeight;
            }

            // Top right point.
            if(
                    baseCoordinateX + 1 < this->width &&
                    x + 1 >= -coordinateOffsetX &&
                    baseCoordinateY < this->height &&
                    y >= -coordinateOffsetY
                    ){
                            double currentWeight = coordinatePartX * (1 - coordinatePartY);
                            sedimentSum += sedimentToMoveMap[baseCoordinateX + 1 + this->width * baseCoordinateY] * currentWeight;
                            sedimentWeightSum += currentWeight;
            }

            // Bottom left point.
            if(
                    baseCoordinateX < this->width &&
                    x >= -coordinateOffsetX &&
                    baseCoordinateY + 1 < this->height &&
                    y + 1 >= -coordinateOffsetY
                    ){
                            double currentWeight = (1 - coordinatePartX) * coordinatePartY;
                            sedimentSum += sedimentToMoveMap[baseCoordinateX + this->width * (baseCoordinateY + 1)] * currentWeight;
                            sedimentWeightSum += currentWeight;
            }

            // Bottom right point.
            if(
                    baseCoordinateX + 1 < this->width &&
                    x + 1 >= -coordinateOffsetX &&
                    baseCoordinateY + 1 < this->height &&
                    y + 1 >= -coordinateOffsetY
                    ){
                            double currentWeight = coordinatePartX * coordinatePartY;
                            sedimentSum += sedimentToMoveMap[baseCoordinateX + 1 + this->width * (baseCoordinateY + 1)] * currentWeight;
                            sedimentWeightSum += currentWeight;
            }

            if(sedimentSum < 0)
                throw runtime_error("Erosion error: Negative interpolated sediment value.");

            sedimentMap[currentIndex] = sedimentWeightSum > 0 ? sedimentSum / sedimentWeightSum : 0;
        }
    }
        
    delete [] sedimentToMoveMap;
}




void ErosionSimulator::ApplyThermalWeathering(double* heightMap, double powerMultiplier)
{
    //GGen_ThermalWeatheringValues* weatheringData = new GGen_ThermalWeatheringValues[this->length];
    //memset(weatheringData, 0, this->length * sizeof(GGen_ThermalWeatheringValues));

    /*
        3. Soil ﬂow calculation with outﬂow in virtual pipes of thermal erosion model.
    */

    double* heightMapCopy = new double[this->length];
    memcpy(heightMapCopy, heightMap, this->length * sizeof(double));

    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < this->width; x++){
            int currentIndex = x + this->width * y;
            double currentHeight = heightMap[currentIndex];

            double heightDiffTopLeft = 0;
            double heightDiffTop = 0;
            double heightDiffTopRight = 0;
            double heightDiffRight = 0;
            double heightDiffBottomLeft = 0;
            double heightDiffBottom = 0;
            double heightDiffBottomRight = 0;
            double heightDiffLeft = 0;

            if(y > 0)
            {
                if(x > 0){
                    heightDiffTopLeft = currentHeight - heightMapCopy[currentIndex - this->width - 1];
                }

                heightDiffTop = currentHeight - heightMapCopy[currentIndex - this->width];

                if(x < this->width - 1){
                    heightDiffTopRight = currentHeight - heightMapCopy[currentIndex - this->width + 1];
                }
            }

            if(x < this->width - 1)
            {
                heightDiffRight = currentHeight - heightMapCopy[currentIndex + 1];
            }

            if(y < this->height - 1)
            {
                if(x > 0) {
                    heightDiffBottomLeft = currentHeight - heightMapCopy[currentIndex + this->width - 1];
                }

                heightDiffBottom = currentHeight - heightMapCopy[currentIndex + this->width];

                if(x < this->width - 1) {
                    heightDiffBottomRight = currentHeight - heightMapCopy[currentIndex + this->width + 1];
                }
            }

            if(x > 0)
            {
                heightDiffLeft = currentHeight - heightMapCopy[currentIndex - 1];
            }

            double maxHeightDiff = 0;
            maxHeightDiff = max(maxHeightDiff, heightDiffTopLeft);
            maxHeightDiff = max(maxHeightDiff, heightDiffTop);
            maxHeightDiff = max(maxHeightDiff, heightDiffTopRight);
            maxHeightDiff = max(maxHeightDiff, heightDiffRight);
            maxHeightDiff = max(maxHeightDiff, heightDiffBottomLeft);
            maxHeightDiff = max(maxHeightDiff, heightDiffBottom);
            maxHeightDiff = max(maxHeightDiff, heightDiffBottomRight);
            maxHeightDiff = max(maxHeightDiff, heightDiffLeft);

            if(maxHeightDiff < this->talusAngle)
            {
                continue;
            }

    /*
        6. Thermal erosion material amount calculation.
    */

            double amountToTransport = maxHeightDiff / 2 * this->deltaT * powerMultiplier;

            double totalTransportableAmount = 0;

            if(heightDiffTopLeft >= this->talusAngle)
                totalTransportableAmount += heightDiffTopLeft;

            if(heightDiffTop >= this->talusAngle)
                totalTransportableAmount += heightDiffTop;

            if(heightDiffTopRight >= this->talusAngle)
                totalTransportableAmount += heightDiffTopRight;

            if(heightDiffRight >= this->talusAngle)
                totalTransportableAmount += heightDiffRight;

            if(heightDiffBottomLeft >= this->talusAngle)
                totalTransportableAmount += heightDiffBottomLeft;

            if(heightDiffBottom >= this->talusAngle)
                totalTransportableAmount += heightDiffBottom;

            if(heightDiffBottomRight >= this->talusAngle)
                totalTransportableAmount += heightDiffBottomRight;

            if(heightDiffLeft >= this->talusAngle)
                totalTransportableAmount += heightDiffLeft;

            heightMap[currentIndex] -= amountToTransport;

            if(y > 0)
            {
                if(x > 0)
                {
                    if(heightDiffTopLeft >= this->talusAngle) heightMap[currentIndex - this->width - 1] += amountToTransport * heightDiffTopLeft / totalTransportableAmount;
                }

                if(heightDiffTop >= this->talusAngle) heightMap[currentIndex - this->width] += amountToTransport * heightDiffTop / totalTransportableAmount;

                if(x < this->width - 1)
                {
                    if(heightDiffTopRight >= this->talusAngle) heightMap[currentIndex - this->width + 1] += amountToTransport * heightDiffTopRight / totalTransportableAmount;
                }
            }

            if(x < this->width - 1)
            {
                if(heightDiffRight >= this->talusAngle) heightMap[currentIndex + 1] += amountToTransport * heightDiffRight / totalTransportableAmount;
            }

            if(y < this->height - 1)
            {
                if(x > 0)
                {
                    if(heightDiffBottomLeft >= this->talusAngle) heightMap[currentIndex + this->width - 1] += amountToTransport * heightDiffBottomLeft / totalTransportableAmount;
                }

                if(heightDiffBottom >= this->talusAngle) heightMap[currentIndex + this->width] += amountToTransport * heightDiffBottom / totalTransportableAmount;

                if(x < this->width - 1)
                {
                    if(heightDiffBottomRight >= this->talusAngle) heightMap[currentIndex + this->width + 1] += amountToTransport * heightDiffBottomRight / totalTransportableAmount;
                }
            }

            if(x > 0)
            {
                if(heightDiffLeft >= this->talusAngle) heightMap[currentIndex - 1] += amountToTransport * heightDiffLeft / totalTransportableAmount;
            }
        }
    }

    //delete [] weatheringData;
    delete [] heightMapCopy;
}
