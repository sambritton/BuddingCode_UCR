
#include "System.h"
#include "SystemStructures.h"
#include "LinearSprings.h"

void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs) {
        std::cout<<"ERROR 30"<<std::endl;
        thrust::fill(linearSpringInfoVecs.tempNodeForceXReduced.begin(),linearSpringInfoVecs.tempNodeForceXReduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceYReduced.begin(),linearSpringInfoVecs.tempNodeForceYReduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceZReduced.begin(),linearSpringInfoVecs.tempNodeForceZReduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceXUnreduced.begin(),linearSpringInfoVecs.tempNodeForceXUnreduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceYUnreduced.begin(),linearSpringInfoVecs.tempNodeForceYUnreduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceZUnreduced.begin(),linearSpringInfoVecs.tempNodeForceZUnreduced.end(),0.0);
    
    
        thrust::counting_iterator<int> edgeIdBegin(0);
        thrust::counting_iterator<int> edgeIdEnd(generalParams.num_of_edges);
        std::cout<<"ERROR 31"<<std::endl;

    //std::cout<<"pre linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
   /* int id =ljInfoVecs.node_id_close[0];
	std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
      */  

    linearSpringInfoVecs.linear_spring_energy=
    thrust::transform_reduce(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                edgeIdBegin,
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin(),
                linearSpringInfoVecs.edge_initial_length.begin())),
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                edgeIdBegin,
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin(), 
                linearSpringInfoVecs.edge_initial_length.begin())) + generalParams.num_of_edges,
        LinearSpringFunctor(
            linearSpringInfoVecs.spring_constant, 
            linearSpringInfoVecs.spring_constant_weak,
            thrust::raw_pointer_cast(generalParams.edges_in_upperhem.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),

            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeIdUnreduced.data()),
            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceXUnreduced.data()),
            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceYUnreduced.data()),
            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceZUnreduced.data()) ),
        0.0, thrust::plus<double>() ); 
        std::cout<<"ERROR 32"<<std::endl;
        //for (int i = 0; i < linearSpringInfoVecs.tempNodeIdUnreduced.size(); i++){
          //  std::cout<<"tempNodeIdUnreduced"<<linearSpringInfoVecs.tempNodeIdUnreduced[i]<<std::endl;
        //}
      //std::cout<<"linear energy from spring.cu: "<< linearSpringInfoVecs.linear_spring_energy<<std::endl;
    //now we have un reduced forces. Sort by id and reduce. 
    //key, then value. Each vector returns sorted		
    thrust::sort_by_key(linearSpringInfoVecs.tempNodeIdUnreduced.begin(), linearSpringInfoVecs.tempNodeIdUnreduced.begin() + (generalParams.num_of_edges*linearSpringInfoVecs.factor),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeForceXUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceYUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<int>());
                std::cout<<"ERROR 33"<<std::endl;
   /* std::cout<<"mid1 linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
    std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
	*/
    int endKey = thrust::get<0>(
        thrust::reduce_by_key(
            linearSpringInfoVecs.tempNodeIdUnreduced.begin(), 
            linearSpringInfoVecs.tempNodeIdUnreduced.end(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeForceXUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceYUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceZUnreduced.begin())),
            linearSpringInfoVecs.tempNodeIdReduced.begin(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeForceXReduced.begin(),
                linearSpringInfoVecs.tempNodeForceYReduced.begin(),
                linearSpringInfoVecs.tempNodeForceZReduced.begin())),
        thrust::equal_to<int>(), CVec3Add())) - linearSpringInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
        std::cout<<"ERROR 34"<<std::endl;
/*
    std::cout<<"mid2 linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
    std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
	*/
    //apply reduced force to all nodes. 
    thrust::for_each(
        thrust::make_zip_iterator(//1st begin
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeIdReduced.begin(),
                linearSpringInfoVecs.tempNodeForceXReduced.begin(),
                linearSpringInfoVecs.tempNodeForceYReduced.begin(),
                linearSpringInfoVecs.tempNodeForceZReduced.begin())),
        thrust::make_zip_iterator(//1st end
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeIdReduced.begin(),
                linearSpringInfoVecs.tempNodeForceXReduced.begin(),
                linearSpringInfoVecs.tempNodeForceYReduced.begin(),
                linearSpringInfoVecs.tempNodeForceZReduced.begin())) + endKey,
        AddForceFunctor (
            generalParams.num_of_nodes,
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));
    /*
    std::cout<<"post linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
    std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
	*/
    std::cout<<"ERROR 35"<<std::endl;
};