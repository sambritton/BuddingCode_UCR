#include "System.h"
#include "SystemStructures.h"
#include "AreaTriangles.h"

void ComputeAreaTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

        thrust::counting_iterator<int> elemId(0);
        
        thrust::fill(areaTriangleInfoVecs.tempNodeForceXReduced.begin(),areaTriangleInfoVecs.tempNodeForceXReduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceYReduced.begin(),areaTriangleInfoVecs.tempNodeForceYReduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceZReduced.begin(),areaTriangleInfoVecs.tempNodeForceZReduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),areaTriangleInfoVecs.tempNodeForceXUnreduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),areaTriangleInfoVecs.tempNodeForceYUnreduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceZUnreduced.begin(),areaTriangleInfoVecs.tempNodeForceZUnreduced.end(),0.0);
    
        
        areaTriangleInfoVecs.area_triangle_energy = thrust::transform_reduce( 
			thrust::make_zip_iterator(
				thrust::make_tuple(
                    elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
                    elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())) + coordInfoVecs.num_triangles,
            AreaSpringFunctor( 
                areaTriangleInfoVecs.initial_area,
                areaTriangleInfoVecs.spring_constant,
                areaTriangleInfoVecs.spring_constant_weak,
                thrust::raw_pointer_cast(generalParams.triangles_in_upperhem.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeIdUnreduced.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceXUnreduced.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceYUnreduced.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceZUnreduced.data())),
            0.0, thrust::plus<double>() );
        
        
        //now we have un reduced forces. Sort by id and reduce. 
        //key, then value. Each vector returns sorted		
		thrust::sort_by_key(areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), areaTriangleInfoVecs.tempNodeIdUnreduced.begin() + (areaTriangleInfoVecs.factor*coordInfoVecs.num_triangles),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
					areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
					areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<int>());

        int endKey = thrust::get<0>(
            thrust::reduce_by_key(
                areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
                areaTriangleInfoVecs.tempNodeIdUnreduced.begin() + (areaTriangleInfoVecs.factor*coordInfoVecs.num_triangles),
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())),
                areaTriangleInfoVecs.tempNodeIdReduced.begin(),
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
            thrust::equal_to<int>(), CVec3Add())) - areaTriangleInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
       // std::cout<<"endKey_area = "<<endKey<<std::endl;

        //apply reduced force to all nodes. 
        thrust::for_each(
            thrust::make_zip_iterator(//1st begin
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
            thrust::make_zip_iterator(//1st end
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())) + endKey,
            AddForceFunctor (
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));
            
      //          std::cout<<"Force from area on node 36 = "<<coordInfoVecs.nodeForceX[35]<<" "<<coordInfoVecs.nodeForceY[35]<<" "<<coordInfoVecs.nodeForceZ[35]<<std::endl;
};