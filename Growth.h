#ifndef GROWTH_H_
#define GROWTH_H_

#include "SystemStructures.h"

class Growth {
    public:
    void growth (int ielem, 
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs);
	
    void transferHtoD(CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);
    
    void transferDtoH(CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);

};

#endif