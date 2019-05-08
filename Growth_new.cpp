#include "System.h"
#include <random>
#include "Growth.h"
#include <math.h>
#include <vector>

void Growth::growth(
    int iedge, 
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs) {


////////////////STRATEGY////////////////////////////////////////////////////////////
/*          n2
           /  \
          /    \
         /  t1  \
        /        \
     n3 ---------- n1
        \        /
         \  t2  /
          \    /
           \  /
            n4
*/
//The key concept is to establish n1, n2, n3, n4 in the above orientation based on the information extracted
//from t1 and t2. This will help the data structure build tremendously.
/////////////////////////////////////////////////////////////////////////////////////


int elem1 = coordInfoVecs.edges2Triangles_1[iedge];
int elem2 = coordInfoVecs.edges2Triangles_2[iedge];

int t1e1, t1e2, t1e3, t2e1, t2e2, t2e3;

if (coordinfoVecs.triangles2Edges_1[elem1] == iedge){
    t1e1 = coordInfoVecs.triangles2Edges_2[elem1];
    t1e2 = coordInfoVecs.triangles2Edges_3[elem1];
    t1e3 = coordInfoVecs.triangles2Edges_1[elem1];
}
else if (coordinfoVecs.triangles2Edges_2[elem1] == iedge){
    t1e1 = coordInfoVecs.triangles2Edges_3[elem1];
    t1e2 = coordInfoVecs.triangles2Edges_1[elem1];
    t1e3 = coordInfoVecs.triangles2Edges_2[elem1];
} 
else if (coordinfoVecs.triangles2Edges_3[elem1] == iedge){
    t1e1 = coordInfoVecs.triangles2Edges_1[elem1];
    t1e2 = coordInfoVecs.triangles2Edges_2[elem1];
    t1e3 = coordInfoVecs.triangles2Edges_3[elem1];
}

if (coordinfoVecs.triangles2Edges_1[elem2] == iedge){
    t2e1 = coordInfoVecs.triangles2Edges_2[elem2];
    t2e2 = coordInfoVecs.triangles2Edges_3[elem2];
    t2e3 = coordInfoVecs.triangles2Edges_1[elem2];
}
else if (coordinfoVecs.triangles2Edges_2[elem2] == iedge){
    t2e1 = coordInfoVecs.triangles2Edges_3[elem2];
    t2e2 = coordInfoVecs.triangles2Edges_1[elem2];
    t2e3 = coordInfoVecs.triangles2Edges_2[elem2];
} 
else if (coordinfoVecs.triangles2Edges_3[elem2] == iedge){
    t2e1 = coordInfoVecs.triangles2Edges_1[elem2];
    t2e2 = coordInfoVecs.triangles2Edges_2[elem2];
    t2e3 = coordInfoVecs.triangles2Edges_3[elem2];
}


int n1, n2, n3, n4;
if (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]){
    n1 = coordInfoVecs.edges2Nodes_1[t1e1];
    n2 = coordInfoVecs.edges2Nodes_2[t1e1];
    n3 = coordInfoVecs.edges2Nodes_2[iedge];
}
else if (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]){
    n1 = coordInfoVecs.edges2Nodes_1[t1e1];
    n2 = coordInfoVecs.edges2Nodes_2[t1e1];
    n3 = coordInfoVecs.edges2Nodes_1[iedge];
}
else if (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]){
    n1 = coordInfoVecs.edges2Nodes_2[t1e1];
    n2 = coordInfoVecs.edges2Nodes_1[t1e1];
    n3 = coordInfoVecs.edges2Nodes_1[iedge];
}
else if (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]){
    n1 = coordInfoVecs.edges2Nodes_2[t1e1];
    n2 = coordInfoVecs.edges2Nodes_1[t1e1];
    n3 = coordInfoVecs.edges2Nodes_1[iedge];
}

if (coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
    n4 = coordInfoVecs.edges2Nodes_2[t2e1];
}
else if (coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
    n4 = coordInfoVecs.edges2Nodes_1[t2e1];
}
//These extract the indices of vertices of the selected triangles "elem1" and "elem2". Now we have n1, n2, n3, n4 in the correct orientation (supposedly).



double newx = (coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
coordInfoVecs.nodeLocX.push_back(newx);
double newy = (coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
coordInfoVecs.nodeLocY.push_back(newy);
double newz = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
coordInfoVecs.nodeLocZ.push_back(newz);
//These are the coordinate of the new vertex. Its index is "coordInfoVecs.nodeLocX.size()"
generalParams.maxNodeCount += 1;


int first_deletion;
int second_deletion;
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

//Finally, we will delete the edge chosen for growth (expansion)
coordInfoVecs.edges2Nodes_1.erase(coordInfoVecs.edges2Nodes_1.begin() + iedge);
coordInfoVecs.edges2Nodes_2.erase(coordInfoVecs.edges2Nodes_2.begin() + iedge);
coordInfoVecs.edges2Nodes_3.erase(coordInfoVecs.edges2Nodes_3.begin() + iedge);
for (int i = 0; i < generalParams.num_of_triangles; i++){
    if (coordInfoVecs.triangles2Edges_1[i] > iedge){
        coordInfoVecs.triangles2Edges_1[i] = coordInfoVecs.triangles2Edges_1[i] - 1;
    }
    if (coordInfoVecs.triangles2Edges_2[i] > iedge){
        coordInfoVecs.triangles2Edges_2[i] = coordInfoVecs.triangles2Edges_2[i] - 1;
    }
}
coordInfoVecs.edges2Triangles_1.erase(coordInfoVecs.edges2Triangles_1.begin() + iedge);
coordInfoVecs.edges2Triangles_2.erase(coordInfoVecs.edges2Triangles_2.begin() + iedge);
coordInfoVecs.edges2Triangles_3.erase(coordInfoVecs.edges2Triangles_3.begin() + iedge);
generalParams.num_of_edges -= 1;

//This should completes the dreadful data structure update associated with cell (membrane) growth. Have fun modifying it if you need more function!
};




//copy configuration from device to host
void Growth::transferDtoH(CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end(),hostSetInfoVecs.nodeLocX.begin());
    thrust::copy(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end(),hostSetInfoVecs.nodeLocY.begin());
    thrust::copy(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end(),hostSetInfoVecs.nodeLocZ.begin());

    
    thrust::copy(coordInfoVecs.triangles2Nodes_1.begin(),coordInfoVecs.triangles2Nodes_1.end(),hostSetInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(coordInfoVecs.triangles2Nodes_2.begin(),coordInfoVecs.triangles2Nodes_2.end(),hostSetInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(coordInfoVecs.triangles2Nodes_3.begin(),coordInfoVecs.triangles2Nodes_3.end(),hostSetInfoVecs.triangles2Nodes_3.begin());
    
    thrust::copy(coordInfoVecs.edges2Nodes_1.begin(),coordInfoVecs.edges2Nodes_1.end(),hostSetInfoVecs.edges2Nodes_1.begin());
    thrust::copy(coordInfoVecs.edges2Nodes_2.begin(),coordInfoVecs.edges2Nodes_2.end(),hostSetInfoVecs.edges2Nodes_2.begin());
    
    thrust::copy(coordInfoVecs.edges2Triangles_1.begin(),coordInfoVecs.edges2Triangles_1.end(),hostSetInfoVecs.edges2Triangles_1.begin());
    thrust::copy(coordInfoVecs.edges2Triangles_2.begin(),coordInfoVecs.edges2Triangles_2.end(),hostSetInfoVecs.edges2Triangles_2.begin());
    
    thrust::copy(coordInfoVecs.triangles2Edges_1.begin(),coordInfoVecs.triangles2Edges_1.end(),hostSetInfoVecs.triangles2Edges_1.begin());
    thrust::copy(coordInfoVecs.triangles2Edges_2.begin(),coordInfoVecs.triangles2Edges_2.end(),hostSetInfoVecs.triangles2Edges_2.begin());
    thrust::copy(coordInfoVecs.triangles2Edges_3.begin(),coordInfoVecs.triangles2Edges_3.end(),hostSetInfoVecs.triangles2Edges_3.begin());

};

//copy configuration from host to device
void Growth::transferHtoD(CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(hostSetInfoVecs.nodeLocX.begin(),hostSetInfoVecs.nodeLocX.end(),coordInfoVecs.nodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(),hostSetInfoVecs.nodeLocY.end(),coordInfoVecs.nodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(),hostSetInfoVecs.nodeLocZ.end(),coordInfoVecs.nodeLocZ.begin());

    
    thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(),hostSetInfoVecs.triangles2Nodes_1.end(),coordInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(),hostSetInfoVecs.triangles2Nodes_2.end(),coordInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(),hostSetInfoVecs.triangles2Nodes_3.end(),coordInfoVecs.triangles2Nodes_3.begin());
    
    thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(),hostSetInfoVecs.edges2Nodes_1.end(),coordInfoVecs.edges2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(),hostSetInfoVecs.edges2Nodes_2.end(),coordInfoVecs.edges2Nodes_2.begin());
    
    thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(),hostSetInfoVecs.edges2Triangles_1.end(),coordInfoVecs.edges2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(),hostSetInfoVecs.edges2Triangles_2.end(),coordInfoVecs.edges2Triangles_2.begin());
    
    thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(),hostSetInfoVecs.triangles2Edges_1.end(),coordInfoVecs.triangles2Edges_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(),hostSetInfoVecs.triangles2Edges_2.end(),coordInfoVecs.triangles2Edges_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(),hostSetInfoVecs.triangles2Edges_3.end(),coordInfoVecs.triangles2Edges_3.begin());

};
