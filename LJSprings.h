#ifndef LJSPRINGS_H_
#define LJSPRINGS_H_ 

#include "SystemStructures.h"

void ComputeLJSprings(
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs,
    GeneralParams& generalParams);

void AdvanceLJParticle(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs);

struct LJSpringFunctor {
    double Rcutoff;
    double Rmin;
    double epsilon_rep1;
    double epsilon_rep2;

    double LJLocX;
    double LJLocY;
    double LJLocZ;

    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ LJSpringFunctor(
        double& _Rcutoff,
        double& _Rmin,
        double& _epsilon_rep1,
        double& _epsilon_rep2,

        double& _LJLocX,
        double& _LJLocY,
        double& _LJLocZ,

        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,

        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr) :

        Rcutoff(_Rcutoff),
        Rmin(_Rmin),
        epsilon_rep1(_epsilon_rep1),
        epsilon_rep2(_epsilon_rep2),

        LJLocX(_LJLocX),
        LJLocY(_LJLocY),
        LJLocZ(_LJLocZ),

        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),

        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ CVec4 operator()(const int & id) {
        		
        // compute forces
        double xLoc = locXAddr[id];
        double yLoc = locYAddr[id];
        double zLoc = locZAddr[id];
   

        double R = sqrt( 
            (LJLocX - xLoc) * (LJLocX - xLoc) + 
            (LJLocY - yLoc) * (LJLocY - yLoc) + 
            (LJLocZ - zLoc) * (LJLocZ - zLoc) );
        
        double forceX=0.0;
        double forceY=0.0;
        double forceZ=0.0;
        double energy = 0.0;
        if (R < Rmin) {
            
            double magnitude = 2*epsilon_rep1*
                                        (1-exp(-epsilon_rep2*(R-Rmin)))*
                                        (-exp(-epsilon_rep2*(R-Rmin)))*
                                        (epsilon_rep2/R);

                    forceX = -magnitude*(LJLocX - xLoc);//xLoc_LR;
                    forceY = -magnitude*(LJLocY - yLoc);//yLoc_LR;
                    forceZ = -magnitude*(LJLocZ - zLoc);//zLoc_LR;

            //forceX = epsilon * (-6.0 * pow(Rmin, 6.0) * (LJLocX - xLoc)/pow(R, 8.0) + 
            //                3.0 * pow(Rmin, 3.0)  * (LJLocX - xLoc)/pow(R, 5.0) );

            //forceY = epsilon * (-6.0 * pow(Rmin, 6.0) * (LJLocY - yLoc)/pow(R, 8.0) + 
            //                3.0 * pow(Rmin, 3.0) * (LJLocY - yLoc)/pow(R, 5.0) );

            //forceZ = epsilon * (-6.0 * pow(Rmin, 6.0) * (LJLocZ - zLoc)/pow(R, 8.0) + 
            //                3.0 * pow(Rmin, 3.0) * (LJLocZ - zLoc)/pow(R, 5.0) );

            //WARNING: since this function modifies nodeForceX etc, you cannot rewrite entries. 
            forceXAddr[id] += forceX;
            forceYAddr[id] += forceY;
            forceZAddr[id] += forceZ;
            //energy = epsilon + epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0));//force positivity
            energy = epsilon_rep1 * (1-exp(-epsilon_rep2*(R - Rmin))) * (1-exp(-epsilon_rep2*(R - Rmin)));
        }

        

        //double energy = epsilon + epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0));//force positivity

        return thrust::make_tuple(-forceX, -forceY, -forceZ, energy);


    }
};

#endif