#ifndef GROWTH_NEW_H_
#define GROWTH_NEW_H_ 

#include "SystemStructures.h"

void MembraneGrowth(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);
    
struct MembraneGrowthFunctor {
    
    double spring_constant;
    double spring_constant_weak;
    int* edges_in_upperhem;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ LinearSpringFunctor(
        double& _spring_constant,
        double& _spring_constant_weak,
        int* _edges_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr) :
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        edges_in_upperhem(_edges_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const Tuuud &u3d) {
int k = edges_in_upperhem[p];
//int k = p;
k -= triggered_counter;
int iedge = k;
//std::cout<<"node1 of iedge = "<<coordInfoVecs.edges2Nodes_1[iedge]<<std::endl;
//std::cout<<"node2 of iedge = "<<coordInfoVecs.edges2Nodes_2[iedge]<<std::endl;
int elem1 = coordInfoVecs.edges2Triangles_1[iedge];
int elem2 = coordInfoVecs.edges2Triangles_2[iedge];
//std::cout<<"elem1 of iedge = "<<elem1<<std::endl;
//std::cout<<"elem2 of iedge = "<<elem2<<std::endl;

int first_v = coordInfoVecs.triangles2Nodes_1[elem1];
int second_v = coordInfoVecs.triangles2Nodes_2[elem1];
int third_v = coordInfoVecs.triangles2Nodes_3[elem1];
double v1x = coordInfoVecs.nodeLocX[second_v] - coordInfoVecs.nodeLocX[first_v];
double v1y = coordInfoVecs.nodeLocY[second_v] - coordInfoVecs.nodeLocY[first_v];
double v1z = coordInfoVecs.nodeLocZ[second_v] - coordInfoVecs.nodeLocZ[first_v];
double v2x = coordInfoVecs.nodeLocX[third_v] - coordInfoVecs.nodeLocX[first_v];
double v2y = coordInfoVecs.nodeLocY[third_v] - coordInfoVecs.nodeLocY[first_v];
double v2z = coordInfoVecs.nodeLocZ[third_v] - coordInfoVecs.nodeLocZ[first_v];
double This_area_v = sqrt((v1y*v2z - v2y*v1z)*(v1y*v2z - v2y*v1z) + 
                        ((-v1x*v2z) + (v2x*v1z))*((-v1x*v2z) + (v2x*v1z)) +
                        (v1x*v2y - v2x*v1y)*(v1x*v2y - v2x*v1y))/2.0;
int first_w = coordInfoVecs.triangles2Nodes_1[elem2];
int second_w = coordInfoVecs.triangles2Nodes_2[elem2];
int third_w = coordInfoVecs.triangles2Nodes_3[elem2];
double w1x = coordInfoVecs.nodeLocX[second_w] - coordInfoVecs.nodeLocX[first_w];
double w1y = coordInfoVecs.nodeLocY[second_w] - coordInfoVecs.nodeLocY[first_w];
double w1z = coordInfoVecs.nodeLocZ[second_w] - coordInfoVecs.nodeLocZ[first_w];
double w2x = coordInfoVecs.nodeLocX[third_w] - coordInfoVecs.nodeLocX[first_w];
double w2y = coordInfoVecs.nodeLocY[third_w] - coordInfoVecs.nodeLocY[first_w];
double w2z = coordInfoVecs.nodeLocZ[third_w] - coordInfoVecs.nodeLocZ[first_w];
double This_area_w = sqrt((w1y*v2z - w2y*v1z)*(w1y*w2z - w2y*w1z) + 
                        ((-w1x*w2z) + (w2x*w1z))*((-w1x*w2z) + (w2x*w1z)) +
                        (w1x*w2y - w2x*w1y)*(w1x*w2y - w2x*w1y))/2.0;
if ((This_area_v/ areaTriangleInfoVecs.initial_area >= 2.0) && (This_area_w/ areaTriangleInfoVecs.initial_area >= 2.0)){
    //triggered = true;
    triggered_counter += 1;
}
else{continue;}

int t1e1, t1e2, t1e3, t2e1, t2e2, t2e3;

if (coordInfoVecs.triangles2Edges_1[elem1] == iedge){
    t1e1 = coordInfoVecs.triangles2Edges_2[elem1];
    t1e2 = coordInfoVecs.triangles2Edges_3[elem1];
    t1e3 = coordInfoVecs.triangles2Edges_1[elem1];
}
else if (coordInfoVecs.triangles2Edges_2[elem1] == iedge){
    t1e1 = coordInfoVecs.triangles2Edges_3[elem1];
    t1e2 = coordInfoVecs.triangles2Edges_1[elem1];
    t1e3 = coordInfoVecs.triangles2Edges_2[elem1];
} 
else if (coordInfoVecs.triangles2Edges_3[elem1] == iedge){
    t1e1 = coordInfoVecs.triangles2Edges_1[elem1];
    t1e2 = coordInfoVecs.triangles2Edges_2[elem1];
    t1e3 = coordInfoVecs.triangles2Edges_3[elem1];
}

if (coordInfoVecs.triangles2Edges_1[elem2] == iedge){
    t2e1 = coordInfoVecs.triangles2Edges_2[elem2];
    t2e2 = coordInfoVecs.triangles2Edges_3[elem2];
    t2e3 = coordInfoVecs.triangles2Edges_1[elem2];
}
else if (coordInfoVecs.triangles2Edges_2[elem2] == iedge){
    t2e1 = coordInfoVecs.triangles2Edges_3[elem2];
    t2e2 = coordInfoVecs.triangles2Edges_1[elem2];
    t2e3 = coordInfoVecs.triangles2Edges_2[elem2];
} 
else if (coordInfoVecs.triangles2Edges_3[elem2] == iedge){
    t2e1 = coordInfoVecs.triangles2Edges_1[elem2];
    t2e2 = coordInfoVecs.triangles2Edges_2[elem2];
    t2e3 = coordInfoVecs.triangles2Edges_3[elem2];
}


int n1, n2, n3, n4;
if ((coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]) || (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]) ){
    n1 = coordInfoVecs.edges2Nodes_1[t1e1];
    n2 = coordInfoVecs.edges2Nodes_2[t1e1];
    if (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]){
        n3 = coordInfoVecs.edges2Nodes_2[iedge];
    }
    else if (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]){
        n3 = coordInfoVecs.edges2Nodes_1[iedge];
    }
}
else if ((coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]) || (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]) ){
    n1 = coordInfoVecs.edges2Nodes_2[t1e1];
    n2 = coordInfoVecs.edges2Nodes_1[t1e1];
    if (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]){
        n3 = coordInfoVecs.edges2Nodes_2[iedge];
    }
    else if (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]){
        n3 = coordInfoVecs.edges2Nodes_1[iedge];
    }
}

if (coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
    n4 = coordInfoVecs.edges2Nodes_2[t2e1];
}
else if (coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
    n4 = coordInfoVecs.edges2Nodes_1[t2e1];
}

//std::cout<<"n1 = "<<n1<<std::endl;
//std::cout<<"n2 = "<<n2<<std::endl;
//std::cout<<"n3 = "<<n3<<std::endl;
//std::cout<<"n4 = "<<n4<<std::endl;
//These extract the indices of vertices of the selected triangles "elem1" and "elem2". Now we have n1, n2, n3, n4 in the correct orientation (supposedly).



double newx = (coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
coordInfoVecs.nodeLocX.push_back(newx);
double newy = (coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
coordInfoVecs.nodeLocY.push_back(newy);
double newz = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
coordInfoVecs.nodeLocZ.push_back(newz);
//These are the coordinate of the new vertex. Its index is "coordInfoVecs.nodeLocX.size()-1"
generalParams.maxNodeCount += 1;





int NODESIZE= generalParams.maxNodeCount;//coordInfoVecs.nodeLocX.size();
coordInfoVecs.triangles2Nodes_1.push_back(n1);
coordInfoVecs.triangles2Nodes_2.push_back(n2);
coordInfoVecs.triangles2Nodes_3.push_back(NODESIZE-1);
//This is a new triangle associated with (tn1, tn2, newnode). Its index is "coordinfoVecs.triangles2Nodes_1.size()-4".
coordInfoVecs.triangles2Nodes_1.push_back(n2);
coordInfoVecs.triangles2Nodes_2.push_back(n3);
coordInfoVecs.triangles2Nodes_3.push_back(NODESIZE-1);
//This is a new triangle associated with (tn2, tn3, newnode). Its index is "coordinfoVecs.triangles2Nodes_1.size()-3".
coordInfoVecs.triangles2Nodes_1.push_back(n3);
coordInfoVecs.triangles2Nodes_2.push_back(n4);
coordInfoVecs.triangles2Nodes_3.push_back(NODESIZE-1);
//This is a new triangle associated with (tn3, tn1, newnode). Its index is "coordinfoVecs.triangles2Nodes_1.size()-2".
coordInfoVecs.triangles2Nodes_1.push_back(n4);
coordInfoVecs.triangles2Nodes_2.push_back(n1);
coordInfoVecs.triangles2Nodes_3.push_back(NODESIZE-1);
//This is a new triangle associated with (tn3, tn1, newnode). Its index is "coordinfoVecs.triangles2Nodes_1.size()-1".
generalParams.num_of_triangles += 4;



//Now we add new edges formed by the addition of the new node.
coordInfoVecs.edges2Nodes_1.push_back(generalParams.maxNodeCount-1);
coordInfoVecs.edges2Nodes_2.push_back(n1);
//This is a new edge associated with (newnode, tn1). Its index is "edges2Nodes_1.size()-4".
coordInfoVecs.edges2Nodes_1.push_back(generalParams.maxNodeCount-1);
coordInfoVecs.edges2Nodes_2.push_back(n2);
//This is a new edge associated with (newnode, tn2). Its index is "edges2Nodes_1.size()-3".
coordInfoVecs.edges2Nodes_1.push_back(generalParams.maxNodeCount-1);
coordInfoVecs.edges2Nodes_2.push_back(n3);
//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-2".
coordInfoVecs.edges2Nodes_1.push_back(generalParams.maxNodeCount-1);
coordInfoVecs.edges2Nodes_2.push_back(n4);
//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-1".
generalParams.num_of_edges += 4;



//Now we update the edges2Triangles data structure with new edges.
int edgeindex, a, a1, a2, a3, temp1, temp2;
int TRIANGLESIZE = generalParams.num_of_triangles;//coordInfoVecs.triangles2Nodes_1.size();
if (coordInfoVecs.edges2Triangles_1[t1e1] == elem1){
    coordInfoVecs.edges2Triangles_1[t1e1] = TRIANGLESIZE-4;
}
else if (coordInfoVecs.edges2Triangles_2[t1e1] == elem1){
    coordInfoVecs.edges2Triangles_2[t1e1] = TRIANGLESIZE-4;
}
else{}
if (coordInfoVecs.edges2Triangles_1[t1e2] == elem1){
    coordInfoVecs.edges2Triangles_1[t1e2] = TRIANGLESIZE-3;
}
else if (coordInfoVecs.edges2Triangles_2[t1e2] == elem1){
    coordInfoVecs.edges2Triangles_2[t1e2] = TRIANGLESIZE-3;
}
else{}
if (coordInfoVecs.edges2Triangles_1[t2e1] == elem2){
    coordInfoVecs.edges2Triangles_1[t2e1] = TRIANGLESIZE-2;
}
else if (coordInfoVecs.edges2Triangles_2[t2e1] == elem2){
    coordInfoVecs.edges2Triangles_2[t2e1] = TRIANGLESIZE-2;
}
if (coordInfoVecs.edges2Triangles_1[t2e2] == elem2){
    coordInfoVecs.edges2Triangles_1[t2e2] = TRIANGLESIZE-1;
}
else if (coordInfoVecs.edges2Triangles_2[t2e2] == elem2){
    coordInfoVecs.edges2Triangles_2[t2e2] = TRIANGLESIZE-1;
}
else{}
//The above change the existing edges2Triangles data structure to accomodate new triangles added.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

for (int j = 0; j < 4; j++){
    //Now we will add the new data associated with both new edges and triangles.
    coordInfoVecs.edges2Triangles_1.push_back(TRIANGLESIZE-(4-j));
    if (j == 0){
        coordInfoVecs.edges2Triangles_2.push_back(TRIANGLESIZE-1);
    }
    else if (j == 1){
        coordInfoVecs.edges2Triangles_2.push_back(TRIANGLESIZE-4);
    }
    else if (j == 2){
        coordInfoVecs.edges2Triangles_2.push_back(TRIANGLESIZE-3);
    }
    else if (j == 3){
        coordInfoVecs.edges2Triangles_2.push_back(TRIANGLESIZE-2);
    }
    //Now we check to see if the order of "push_back" done is correct, i.e. are edges2Triangles data in correct orientation.
    //This is crucial in the bendingspring computation.
    edgeindex = (generalParams.num_of_edges - (4-j));
    a = coordInfoVecs.edges2Triangles_1[edgeindex];
    if ((coordInfoVecs.triangles2Nodes_1[a] == coordInfoVecs.edges2Nodes_1[edgeindex]) && (coordInfoVecs.triangles2Nodes_2[a] == coordInfoVecs.edges2Nodes_2[edgeindex])){
        a1 = 1;
    }
    else{
        a1 = 0;
    }
    if ((coordInfoVecs.triangles2Nodes_2[a] == coordInfoVecs.edges2Nodes_1[edgeindex]) && (coordInfoVecs.triangles2Nodes_3[a] == coordInfoVecs.edges2Nodes_2[edgeindex])){
        a2 = 1;
    }
    else{
        a2 = 0;
    }
    if ((coordInfoVecs.triangles2Nodes_3[a] == coordInfoVecs.edges2Nodes_1[edgeindex]) && (coordInfoVecs.triangles2Nodes_1[a] == coordInfoVecs.edges2Nodes_2[edgeindex])){
        a3 = 1;
    }
    else{
        a3 = 0;
    }

    if ((a1+a2+a3) == 0){
        temp1 = coordInfoVecs.edges2Triangles_1[edgeindex];
        temp2 = coordInfoVecs.edges2Triangles_2[edgeindex];
        coordInfoVecs.edges2Triangles_1[edgeindex] = temp2;
        coordInfoVecs.edges2Triangles_2[edgeindex] = temp1;
    }
    else{}
    //This checks if the orientation is correct or not, if not, flip the ordering.
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Now we will take care of the last unedited data structure "triangles2Edges".
//int aa, bb;
int EDGESIZE = generalParams.num_of_edges;//coordInfoVecs.edges2Nodes_1.size();
for (int j = 0; j < 4; j++){
    if (j == 0){
        coordInfoVecs.triangles2Edges_1.push_back(EDGESIZE-4);
        coordInfoVecs.triangles2Edges_2.push_back(t1e1);
        coordInfoVecs.triangles2Edges_3.push_back(EDGESIZE-3);   
    }
    else if (j == 1){
        coordInfoVecs.triangles2Edges_1.push_back(EDGESIZE-3);
        coordInfoVecs.triangles2Edges_2.push_back(t1e2);
        coordInfoVecs.triangles2Edges_3.push_back(EDGESIZE-2);   
    }
    else if (j ==2){
        coordInfoVecs.triangles2Edges_1.push_back(EDGESIZE-2);
        coordInfoVecs.triangles2Edges_2.push_back(t2e1);
        coordInfoVecs.triangles2Edges_3.push_back(EDGESIZE-1);   
    }
    else if (j ==3){
        coordInfoVecs.triangles2Edges_1.push_back(EDGESIZE-1);
        coordInfoVecs.triangles2Edges_2.push_back(t2e2);
        coordInfoVecs.triangles2Edges_3.push_back(EDGESIZE-4);   
    }
    
}


//Finally, we will delete the edge chosen for growth (expansion)

coordInfoVecs.edges2Nodes_1.erase(coordInfoVecs.edges2Nodes_1.begin() + iedge);
coordInfoVecs.edges2Nodes_2.erase(coordInfoVecs.edges2Nodes_2.begin() + iedge);
for (int i = 0; i < generalParams.num_of_triangles; i++){
    if (coordInfoVecs.triangles2Edges_1[i] > iedge){
        coordInfoVecs.triangles2Edges_1[i] = coordInfoVecs.triangles2Edges_1[i] - 1;
    }
    if (coordInfoVecs.triangles2Edges_2[i] > iedge){
        coordInfoVecs.triangles2Edges_2[i] = coordInfoVecs.triangles2Edges_2[i] - 1;
    }
    if (coordInfoVecs.triangles2Edges_3[i] > iedge){
        coordInfoVecs.triangles2Edges_3[i] = coordInfoVecs.triangles2Edges_3[i] - 1;
    }
}
coordInfoVecs.edges2Triangles_1.erase(coordInfoVecs.edges2Triangles_1.begin() + iedge);
coordInfoVecs.edges2Triangles_2.erase(coordInfoVecs.edges2Triangles_2.begin() + iedge);
generalParams.num_of_edges -= 1;

int first_deletion;
int second_deletion;
int ielem;
if (elem1 > elem2){
    first_deletion = elem1;
    second_deletion = elem2;
}
else if (elem2 > elem1){
    first_deletion = elem2;
    second_deletion = elem1;
}
for (int i = 0; i < 2; i++){
    if (i == 0){
        ielem = first_deletion;
    }
    else if (i == 1){
        ielem = second_deletion;
    }
    coordInfoVecs.triangles2Nodes_1.erase(coordInfoVecs.triangles2Nodes_1.begin() + ielem);
    coordInfoVecs.triangles2Nodes_2.erase(coordInfoVecs.triangles2Nodes_2.begin() + ielem);
    coordInfoVecs.triangles2Nodes_3.erase(coordInfoVecs.triangles2Nodes_3.begin() + ielem);
    //Delete the associated vertices information of the selected triangle.
    //Since we delete the chosen triangles, any triangle indexed lower than the deleted one will have its index reduced (or moved up) by 1.
    //Hence, we need to sweep through all data structures using the triangle index to change the index accordingly.
    for (int i = 0; i < generalParams.num_of_edges; i++){
        if (coordInfoVecs.edges2Triangles_1[i] > ielem){
            coordInfoVecs.edges2Triangles_1[i] = coordInfoVecs.edges2Triangles_1[i] - 1;
        }
        if (coordInfoVecs.edges2Triangles_2[i] > ielem){
            coordInfoVecs.edges2Triangles_2[i] = coordInfoVecs.edges2Triangles_2[i] - 1;
        }
    }
    //This completes the sweep. After this, the indices of triangle used in edges2Triangles data structure should be the correct one.
    generalParams.num_of_triangles += -1;
}

for (int i = 0; i < 2; i++){
    if (i == 0){
        ielem = first_deletion;
    }
    else if (i == 1){
        ielem = second_deletion;
    }
    coordInfoVecs.triangles2Edges_1.erase(coordInfoVecs.triangles2Edges_1.begin() + ielem);
    coordInfoVecs.triangles2Edges_2.erase(coordInfoVecs.triangles2Edges_2.begin() + ielem);
    coordInfoVecs.triangles2Edges_3.erase(coordInfoVecs.triangles2Edges_3.begin() + ielem);
}
//Erase the edge infomation related to the deleted triangle. Note the deletion should always start with the largest index.

//Before we delete the edge, determine whether the newly added node is part of nodes_in_upperhem or not.
if (generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[iedge]] == 1 && generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[iedge]] == 1){
    generalParams.nodes_in_upperhem.push_back(1);
}
else{
    generalParams.nodes_in_upperhem.push_back(-1);
}


                //Erase the edge infomation related to the deleted triangle.

                //Now we update the nodes_in_upperhem and edges_in_upperhem data structures.
                //This ensures that the newly created edges will have the correct associated spring constant.
//std::cout<<"ERROR HERE?"<<std::endl;
generalParams.edges_in_upperhem.erase(generalParams.edges_in_upperhem.begin() + iedge);
for (int i = 0; i < edges_in_upperhem.size(); i++){
    if (edges_in_upperhem[i] == iedge){
        edges_in_upperhem.erase(edges_in_upperhem.begin() + i);
        break;
    }
}

                
if (elem1 > elem2){
    generalParams.triangles_in_upperhem.erase(generalParams.triangles_in_upperhem.begin() + elem1);
    generalParams.triangles_in_upperhem.erase(generalParams.triangles_in_upperhem.begin() + elem2);
}
else if (elem2 > elem1){
    generalParams.triangles_in_upperhem.erase(generalParams.triangles_in_upperhem.begin() + elem2);
    generalParams.triangles_in_upperhem.erase(generalParams.triangles_in_upperhem.begin() + elem1);
}

                for (int q = 0; q < 4; q++){
                    int nodeP = coordInfoVecs.triangles2Nodes_1[generalParams.num_of_triangles - (4-q)]; 
                    int nodeQ = coordInfoVecs.triangles2Nodes_2[generalParams.num_of_triangles - (4-q)];
                    int nodeR = coordInfoVecs.triangles2Nodes_3[generalParams.num_of_triangles - (4-q)];

                    if (generalParams.nodes_in_upperhem[nodeP]==1 && generalParams.nodes_in_upperhem[nodeQ] ==1 && generalParams.nodes_in_upperhem[nodeR] ==1){
                        generalParams.triangles_in_upperhem.push_back(1);
                    }
                    else{
                        generalParams.triangles_in_upperhem.push_back(-1);
                    }
                }
                for (int q = 0; q < 4; q++){
                    int elem_1 = coordInfoVecs.edges2Triangles_1[generalParams.num_of_edges-(4 - q)];
                    //std::cout<<nodeP<<std::endl;
                    //std::cout<<generalParams.nodes_in_upperhem[nodeP]<<std::endl;
                    int elem_2 = coordInfoVecs.edges2Triangles_2[generalParams.num_of_edges-(4 - q)];
                    //std::cout<<nodeQ<<std::endl;
                    //std::cout<<generalParams.nodes_in_upperhem[nodeQ]<<std::endl;
                    
                    if (generalParams.triangles_in_upperhem[elem_1] == 1 && generalParams.triangles_in_upperhem[elem_2] == 1){
                        
                        generalParams.edges_in_upperhem.push_back(1);
                        edges_in_upperhem.push_back(generalParams.num_of_edges - (4 - q));
                    }
                    else{
                        
                        generalParams.edges_in_upperhem.push_back(-1);
                    }
                    
                    
                }
                
                //This should completes the dreadful data structure update associated with cell (membrane) growth.
                //Have fun modifying it if you need more function!
                //if (triggered == true){
                //	break;
                //}