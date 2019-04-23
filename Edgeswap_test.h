#ifndef EDGESWAP_TEST_H_
#define EDGESWAP_TEST_H_

#include "SystemStructures.h"

class Edgeswap {
    std::vector<bool> boundary_node;
    std::vector<unsigned> nndata;

    public:
    Edgeswap(CoordInfoVecs& coordInfoVecs);

	std::vector<bool> DomainBd (CoordInfoVecs& coordInfoVecs);
	std::vector<unsigned> Number_of_Neighbor(CoordInfoVecs& coordInfoVecs);
    int edge_swap_device_vecs (unsigned iedge, 
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        BendingTriangleInfoVecs& bendingTriangleInfoVecs,
        AreaTriangleInfoVecs& areaTriangleInfoVecs);

    void edge_swap_host_vecs (unsigned iedge, 
        GeneralParams& generalParams,
        HostSetInfoVecs& hostSetInfoVecs,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        BendingTriangleInfoVecs& bendingTriangleInfoVecs,
        AreaTriangleInfoVecs& areaTriangleInfoVecs);
	
    void transferHtoD(CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);
    
    void transferDtoH(CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);

};

#endif