#ifndef BENDINGTRIANGLESENERGY_H_
#define BENDINGTRIANGLESENERGY_H_

#include "SystemStructures.h"
#include <math.h>

double ComputeCosTrianglesEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs);

struct CosBendingEnergyFunctor {
	double spring_constant;
	double angle_0;
	double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    unsigned* idKey;

	unsigned* triangle2Nodes_1Addr;
	unsigned* triangle2Nodes_2Addr;
	unsigned* triangle2Nodes_3Addr;

	__host__ __device__ CosBendingEnergyFunctor(
		double& _spring_constant,
		double& _angle_0,
		double* _locXAddr,
		double* _locYAddr,
		double* _locZAddr,

		unsigned* _idKey,

		unsigned* _triangle2Nodes_1Addr,
		unsigned* _triangle2Nodes_2Addr,
		unsigned* _triangle2Nodes_3Addr) :
		


		spring_constant(_spring_constant),
		angle_0(_angle_0),
		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		
		idKey(_idKey),

		triangle2Nodes_1Addr(_triangle2Nodes_1Addr),
		triangle2Nodes_2Addr(_triangle2Nodes_2Addr),
		triangle2Nodes_3Addr(_triangle2Nodes_3Addr) {}

	//hand in counting iterator and id's of two triangles, and two nodes involved. 
	__device__ double operator()(const Tuuuuu &u5) {
		double energy = 0.0;
		unsigned counter = thrust::get<0>(u5);
		unsigned place = 4 * counter;//represents location in write to vector.

		unsigned id_l, id_j;
		//Id's of elements on sides of a given edge
		//these are the same if we are at the edge of the membrane, so do not compute those.
		unsigned T1 = thrust::get<1>(u5);
		unsigned T2 = thrust::get<2>(u5);

		//these id's are accurate
		unsigned id_k = thrust::get<3>(u5);
		unsigned id_i = thrust::get<4>(u5);
		if (T1 != T2) {
			//we need to compute rl and rj from the two involved triangles. 

			//j is from the element on edge2elem_1
			//l is from the element node in edge2elem_1, not i or k
			unsigned n1T1 = triangle2Nodes_1Addr[T1];
			unsigned n2T1 = triangle2Nodes_2Addr[T1];
			unsigned n3T1 = triangle2Nodes_3Addr[T1];
			if ((n1T1 != id_i) && (n1T1 != id_k)) {
				id_j = n1T1;
			}
			else if ((n2T1 != id_i) && (n2T1 != id_k)) {
				id_j = n2T1;
			}
			else if ((n3T1 != id_i) && (n3T1 != id_k)) {
				id_j = n3T1;
			}



			//l is from the element node in edge2elem_2, not i or k
			//one of these is l, find it
			unsigned n1T2 = triangle2Nodes_1Addr[T2];
			unsigned n2T2 = triangle2Nodes_2Addr[T2];
			unsigned n3T2 = triangle2Nodes_3Addr[T2];
			if ((n1T2 != id_i) && (n1T2 != id_k)) {
				id_l = n1T2;
			}
			else if ((n2T2 != id_i) && (n2T2 != id_k)) {
				id_l = n2T2;
			}
			else if ((n3T2 != id_i) && (n3T2 != id_k)) {
				id_l = n3T2;
			}

			CVec3 ri = thrust::make_tuple<double>(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
			CVec3 rj = thrust::make_tuple<double>(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
			CVec3 rk = thrust::make_tuple<double>(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);
			CVec3 rl = thrust::make_tuple<double>(locXAddr[id_l], locYAddr[id_l], locZAddr[id_l]);

			CVec3 rjk = CVec3_minus(rk, rj);
			CVec3 rji = CVec3_minus(ri, rj);
			CVec3 rli = CVec3_minus(ri, rl);		
			CVec3 rlk = CVec3_minus(rk, rl);		
			CVec3 rki = CVec3_minus(ri, rk);		 
			
			CVec3 N1 = CVec3_cross(rjk,rji);
			CVec3 N2 = CVec3_cross(rli, rlk);
			double nN1 = sqrt(CVec3_dot(N1,N1));
			double nN2 = sqrt(CVec3_dot(N2,N2));

			//N1 = cross(rjk, rji);
			//N2 = cross(rli, rlk);
			//nN1 = sqrt(sum(N1.^2)); %norm of N1
			//nN2 = sqrt(sum(N2.^2)); %norm of N2

			//due to symmetry, we the angle does not matter if it is the convex or concave one. 
			double cosAngle = CVec3_dot(N1, N2) / (nN1*nN2);
			
			if (cosAngle > 1.0) {
				cosAngle = 1.0;
			}
			else if (cosAngle < -1.0){
				cosAngle = -1.0;
			}

			double theta_current = acos( cosAngle );
			
			energy = spring_constant * (1 - cos(theta_current - angle_0) );
			
		}
		return energy;
	}
};





#endif /* BENDINGTRIANGLESENERGY_H_*/