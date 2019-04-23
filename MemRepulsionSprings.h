
#ifndef MEMREPULSIONSPRINGS_H_
#define MEMREPULSIONSPRINGS_H_ 

#include "SystemStructures.h"

void ComputeMemRepulsionSprings(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);

struct MemRepulsionSpringFunctor : public thrust::unary_function<U2CVec3,CVec3> {
    double Rmin;
    double abs_Rmin;
    double epsilon_rep1;
    double epsilon_rep2;
    unsigned membraneMaxNode;

    int* nndata1;
    int* nndata2;
    int* nndata3;
    int* nndata4;
    int* nndata5;
    int* nndata6;
    int* nndata7;
    int* nndata8;
    int* nndata9;
    int* nndata10;
    int* nndata11;
    int* nndata12;

    double* membraneNodeXAddr;
    double* membraneNodeYAddr;
    double* membraneNodeZAddr;

    unsigned* id_value_expanded;
    unsigned* keyBegin;
    unsigned* keyEnd;
    
	__host__ __device__ 
    MemRepulsionSpringFunctor(
        double& _Rmin,
        double& _abs_Rmin,
        double& _epsilon_rep1,
        double& _epsilon_rep2,
        unsigned& _membraneMaxNode,
     
        int* _nndata1,
        int* _nndata2,
        int* _nndata3,
        int* _nndata4,
        int* _nndata5,
        int* _nndata6,
        int* _nndata7,
        int* _nndata8,
        int* _nndata9,
        int* _nndata10,
        int* _nndata11,
        int* _nndata12,

        double* _membraneNodeXAddr,
        double* _membraneNodeYAddr,
        double* _membraneNodeZAddr,

        unsigned* _id_value_expanded,
        unsigned* _keyBegin,
        unsigned* _keyEnd):

        Rmin(_Rmin),
        abs_Rmin(_abs_Rmin),
        epsilon_rep1(_epsilon_rep1),
        epsilon_rep2(_epsilon_rep2),
        membraneMaxNode(_membraneMaxNode),

        nndata1(_nndata1),
        nndata2(_nndata2),
        nndata3(_nndata3),
        nndata4(_nndata4),
        nndata5(_nndata5),
        nndata6(_nndata6),
        nndata7(_nndata7),
        nndata8(_nndata8),
        nndata9(_nndata9),
        nndata10(_nndata10),
        nndata11(_nndata11),
        nndata12(_nndata12),
    
        membraneNodeXAddr(_membraneNodeXAddr),
        membraneNodeYAddr(_membraneNodeYAddr),
        membraneNodeZAddr(_membraneNodeZAddr),

        id_value_expanded(_id_value_expanded),
        keyBegin(_keyBegin),
        keyEnd(_keyEnd) {}

	//hand in counting iterator id for membrane nodes
	__device__ 
    CVec3 operator()(const U2CVec3& u2d3) {

        unsigned node_id = thrust::get<0>(u2d3);
        
        unsigned bucketId = thrust::get<1>(u2d3);//bucket containing nodeId
	
		//beginning and end of attempted attachment id's in id_value_expanded
        //these indices are membrane
	
        double xLoc_LR;
        double yLoc_LR;
        double zLoc_LR;

        //extract current force
        double current_forceX = thrust::get<2>(u2d3);
        double current_forceY = thrust::get<3>(u2d3);
        double current_forceZ = thrust::get<4>(u2d3);

		double posX = membraneNodeXAddr[node_id];
        double posY = membraneNodeYAddr[node_id];
        double posZ = membraneNodeZAddr[node_id];

        //for now iterate through all membrane id's
       // unsigned begin = keyBegin[bucketId];
       // unsigned end = keyEnd[bucketId];

        for (unsigned memId_count = 0; memId_count < membraneMaxNode; memId_count++ ){
            
            unsigned memId = memId_count;//id_value_expanded[ memId_count ];
            xLoc_LR = -( posX - membraneNodeXAddr[memId] ); //This is like Rk - Rj, where Rj is the point we are trying to update
            yLoc_LR = -( posY - membraneNodeYAddr[memId] );
            zLoc_LR = -( posZ - membraneNodeZAddr[memId] );

            double R = sqrt( 
                        (xLoc_LR) * (xLoc_LR) + 
                        (yLoc_LR) * (yLoc_LR) + 
                        (zLoc_LR) * (zLoc_LR) );

            if ((memId < membraneMaxNode) && (memId != node_id) && (R < abs_Rmin)) {
                bool neighbor;
                for (unsigned i = 0; i < membraneMaxNode; i++){
                    if (memId == nndata1[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata2[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata3[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata4[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata5[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata6[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata7[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata8[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata9[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata10[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata11[node_id]){
                        neighbor = true;
                    }
                    else if (memId == nndata12[node_id]){
                        neighbor = true;
                    }
                    else {}
                }
                if (neighbor == false){
                    
                    //if (R < Rmin) {
                        double magnitude = 2*epsilon_rep1*
                                            (1-exp(-epsilon_rep2*(R-Rmin)))*
                                            (-exp(-epsilon_rep2*(R-Rmin)))*
                                            (epsilon_rep2/R);

                        current_forceX += -magnitude*xLoc_LR;
                        current_forceY += -magnitude*yLoc_LR;
                        current_forceZ += -magnitude*zLoc_LR;
                    //}
                    
                    //    current_forceX += (epsilon_rep1)*(R - abs_Rmin)*(xLoc_LR)/R;
                     //   current_forceY += (epsilon_rep1)*(R - abs_Rmin)*(yLoc_LR)/R;
                   //     current_forceZ += (epsilon_rep1)*(R - abs_Rmin)*(zLoc_LR)/R;
                    
                }
            }
        }
        return thrust::make_tuple(current_forceX, current_forceY, current_forceZ);
    }
};

#endif