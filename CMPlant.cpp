
//
//  CMPlant.cpp
//  CMPlant
//
//  Created by Mihails Delmans on 25/07/2013.
//  Copyright (c) 2013 Mihails Delmans. All rights reserved.
//

#include "CMPlant.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenCL/cl.h>
#include <utility>

#include "CLGeo.h"

#include "rapidxml-1.13/rapidxml.hpp"

#define N_IDX 0
#define BUF_ID 0
#define MERGE_LIM 600

#define P_AREA  s[0]
#define P_THIC  s[1]
#define P_VISC  s[2]
#define P_YMOD  s[3]

#define P_SLEN  s[4]
#define P_SRFT  s[5]

#define MAX_SOURCE_SIZE 100000

// Private Functions



cl_int2 CMPlant::PrevInCell(int cellIdx, int nodeIdx){
    
    for (int i = 1; i <= cellNodes[cellIdx][0]; i++){
        
        if (cellNodes[cellIdx][i] == nodeIdx){
            if ( i == 1 ) return {cellNodes[cellIdx][ cellNodes[cellIdx][0] ],cellNodes[cellIdx][0]}; // return the last index
            else return {cellNodes[cellIdx][i-1],i-1};
        }
    }
    return {-1,-1};
}

cl_int2 CMPlant::NextInCell(int cellIdx, int nodeIdx){
    
    for (int i = 1; i <= cellNodes[cellIdx][0]; i++){
        
        if (cellNodes[cellIdx][i] == nodeIdx){
            if ( i == cellNodes[cellIdx][0] ) return {cellNodes[cellIdx][1],1}; // return the first index
            else return {cellNodes[cellIdx][i+1],i+1};
        }
    }
    return {-1,-1};
    
}

int CMPlant::CellSeqInNode(int cellIdx, int nodeIdx){
    for (int i = 1; i <= nodeCells[nodeIdx][N_IDX].x; i++) 
        if (nodeCells[nodeIdx][i].x == cellIdx) return i;
    
    return 0;
}

float CMPlant::GetAngle(int cNode, int cellIdx){
    
    int nNode = NextInCell(cellIdx, cNode).x;
    int pNode = PrevInCell(cellIdx, cNode).x;
    
    cl_float4 prevEdge = nodesPos[pNode] - nodesPos[cNode];
    cl_float4 nextEdge = nodesPos[nNode] - nodesPos[cNode];
    float cosine =  fmin( fmax((prevEdge * nextEdge) / (SecondNorm(prevEdge) * SecondNorm(nextEdge)) , -1) , 1);
    
    return acos(cosine);
    
}

//##################################################################################################################################


float CMPlant::CellArea(int cellIdx){
    float area = 0;
    
    for (int i =1; i <= cellNodes[cellIdx][N_IDX]; i++ ){
        
        int cNode = cellNodes[cellIdx][i];
        int nNode = NextInCell(cellIdx, cNode).x;
        
        area += Cross(nodesPos[cNode], nodesPos[nNode] - nodesPos[cNode]).z;
    }
    return area;
}

void CMPlant::CreateOneCell(int n, float r)
{
    srand(time(NULL));
    
    float ang;
    float rad;
    
    cellsN = 1;
    nodesN = n;
    
    cellDivide[1] = true;
    
    cellNodes[1][N_IDX] = n;  // Store number of nodes for the first cell at 0'th index of the vector
    
    for (int i = 1; i <= n; i++) {
        
        ang = 2 * M_PI * i / n;
        rad = r * (1 + (rand() % 90 - 45) / 100.0 );
        
        nodesPos[i] = {(float)(rad * cos(ang)), (float)(rad * sin(ang)), 0 , 0};
        cellNodes[1][i] = i;
        
        nodeCells[i][N_IDX] = {1,1}; // Set the number of cells to 1
        nodeCells[i][1] = {1,i};
        
    }
    
    for (int i =1; i <= n; i++) {
        wallProperties[i][1] = newWallProperties;
        wallProperties[i][1].P_SLEN = SecondNorm(nodesPos[i] - nodesPos[ NextInCell(1, i).x ]);
    }
}


//##################################################################################################################################

void CMPlant::ReadCellsFromFile(const char *filename, int dispW, int dispH, float scale){
    
    bool usedWall[100000];
    
    rapidxml::xml_document<> doc;
    
    std::ifstream file;
    
    file.open(filename);
    
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();
    
    
    std::string content(buffer.str());
    doc.parse<0>(&content[0]);
    
    
    char *pEnd;
    
    
    int oDispW = (int)strtol(doc.first_node("cellsetdata")->first_node("metadata")->first_node("image")->first_attribute("width")->value(),&pEnd,10);
    int oDispH = (int)strtol(doc.first_node("cellsetdata")->first_node("metadata")->first_node("image")->first_attribute("height")->value(),&pEnd,10);
    
    rapidxml::xml_node<> *pCells = doc.first_node()->first_node("cells");
    rapidxml::xml_node<> *pWalls = doc.first_node()->first_node("walls");
    
    cellsN = 0;
    nodesN = 0;
    
    
    for (rapidxml::xml_node<> *pWall = pWalls->first_node("wall"); pWall; pWall = pWall->next_sibling()) {
        
        for (rapidxml::xml_node<> *pPoint = pWall->first_node("points")->first_node("point"); pPoint ; pPoint = pPoint->next_sibling()) {
            
            bool nodeFound = false;
            cl_int2 pointPos;
            pointPos.x = (int)strtol( pPoint->first_attribute("x")->value(),&pEnd,10);
            pointPos.y = (int)strtol( pPoint->first_attribute("y")->value(),&pEnd,10);
            
            for(int i = 1; i <= nodesN; i++){
                if ( (int)nodesPos[i].x == pointPos.x && (int)nodesPos[i].y == pointPos.y ){
                    nodeFound = true;
                    break;
                }
            }
            
            if (!nodeFound) {
                nodesN++;
                nodesPos[nodesN] = {(float)pointPos.x, (float)pointPos.y};
    
            }
                
        }
        
    }
    
    
    for (rapidxml::xml_node<> *pCell = pCells->first_node("cell"); pCell; pCell = pCell->next_sibling()){
        
        int cCell = ++cellsN;
        
        int cGroup = (int)strtol(pCell->first_attribute("group")->value(),&pEnd,10);
        
        switch (cGroup) {
            case 0:
                cellDivide[cCell] = false;
                break;
            case 1:
                cellDivide[cCell] = true;
            default:
                break;
        }
        
        
        for (rapidxml::xml_node<> *pCellWall = pCell->first_node("walls")->first_node("wall"); pCellWall; pCellWall=pCellWall->next_sibling()) {
            
            int tWallIdx = (int)strtol(pCellWall->first_attribute("id")->value(), &pEnd, 10);
            
            for (rapidxml::xml_node<> *pWall = pWalls->first_node("wall"); pWall; pWall = pWall->next_sibling() ) {
                
                int cWallIdx = (int)strtol(pWall->first_attribute("id")->value(), &pEnd, 10);
                
                if (tWallIdx == cWallIdx)
                    if (!usedWall[tWallIdx]){
                        for (rapidxml::xml_node<> *pPoint = pWall->first_node("points")->first_node("point")->next_sibling(); pPoint; pPoint = pPoint->next_sibling()) {
                            
                            cl_int2 pointPos;
                            pointPos.x = (int)strtol( pPoint->first_attribute("x")->value(),&pEnd,10);
                            pointPos.y = (int)strtol( pPoint->first_attribute("y")->value(),&pEnd,10);
                            
                            int cNodeIdx = 0;
                            
                            for (int i = 1; i <= nodesN; i++)
                                if ( (int)nodesPos[i].x == pointPos.x && (int)nodesPos[i].y == pointPos.y){
                                    cNodeIdx = i;
                                    break;
                                }
                            
                            cellNodes[cCell][N_IDX]++;
                            cellNodes[cCell][cellNodes[cCell][N_IDX]] = cNodeIdx;
                            
                            nodeCells[cNodeIdx][N_IDX].x++;
                            nodeCells[cNodeIdx][ nodeCells[cNodeIdx][N_IDX].x ] = {cCell, cellNodes[cCell][N_IDX]};
                        }
                        usedWall[tWallIdx] = true;
                    }
                
                    else
                        for (rapidxml::xml_node<> *pPoint = pWall->first_node("points")->last_node("point")->previous_sibling(); pPoint; pPoint = pPoint->previous_sibling()) {
                            
                            cl_int2 pointPos;
                            pointPos.x = (int)strtol( pPoint->first_attribute("x")->value(),&pEnd,10);
                            pointPos.y = (int)strtol( pPoint->first_attribute("y")->value(),&pEnd,10);
                            
                            int cNodeIdx = 0;
                            
                            for (int i = 1; i <= nodesN; i++)
                                if ( (int)nodesPos[i].x == pointPos.x && (int)nodesPos[i].y == pointPos.y) {
                                    cNodeIdx = i;
                                    break;
                                }
                            
                            cellNodes[cCell][N_IDX]++;
                            cellNodes[cCell][cellNodes[cCell][N_IDX]] = cNodeIdx;
                            
                            nodeCells[cNodeIdx][N_IDX].x++;
                            nodeCells[cNodeIdx][ nodeCells[cNodeIdx][N_IDX].x ] = {cCell, cellNodes[cCell][N_IDX]};
                        }

                    
                
            }
            
        }
        
    }
    
    for (int i = 1; i <= nodesN; i++) {
        nodesPos[i] = {nodesPos[i].x * ((float)dispW / oDispW) - dispW * 0.5f , nodesPos[i].y * ((float)dispH / oDispH) - dispH * 0.5f };
        nodesPos[i] = nodesPos[i] * scale;
    }
    
    for (int i =1; i <= nodesN; i++) {
        for (int j = 1; j<= nodeCells[i][N_IDX].x; j++) {
            
        wallProperties[i][j] = newWallProperties;
        wallProperties[i][j].P_SLEN = SecondNorm(nodesPos[i] - nodesPos[ NextInCell(nodeCells[i][j].x, i).x ]);
        }
    }
}

//##################################################################################################################################

int CMPlant::InsertNewNode(int nodeIdx1, int nodeIdx2, cl_float4 newNodePos){
    
    int insertNodeIdx = ++nodesN;
    
    nodeCells[insertNodeIdx][N_IDX] = {0};
    wallProperties[insertNodeIdx][N_IDX] = {0};
    
    for (int i = 1; i<= nodeCells[nodeIdx1][N_IDX].x; i++) { 
        
        cl_int2 cCell = nodeCells[nodeIdx1][i];
        
        if (NextInCell(cCell.x, nodeIdx1).x == nodeIdx2 || NextInCell(cCell.x, nodeIdx2).x == nodeIdx1 ){ // In all cells look for the edge nodeIdx1 - nodeIdx2
            
            int buffIdx1, buffIdx2;
            //cl_float8 buffProp1, buffProp2;
            
            bool useBuffs = false;
            
            nodeCells[insertNodeIdx][N_IDX].x ++;
    
            nodesPos[insertNodeIdx] = newNodePos;
            
            int triggerIdx, triggerPos;
            
            
            if ( NextInCell(cCell.x, nodeIdx1).x == nodeIdx2 ){
                triggerIdx = nodeIdx1;
                triggerPos = i;
            }
            else {
                triggerIdx = nodeIdx2;
                triggerPos = PrevInCell(cCell.x, nodeIdx1).y;
            }
       
            
            
            for (int j = 1; j <= cellNodes[cCell.x][N_IDX]; j++){ // For each node in the cell
                
                
                if (useBuffs) { 
                    buffIdx2 = cellNodes[cCell.x][j+1];
                    cellNodes[cCell.x][j+1] = buffIdx1;
                    
                    for(int k = 1; k <= nodeCells[buffIdx1][0].x; k++)
                        if (nodeCells[buffIdx1][k].x == cCell.x) nodeCells[buffIdx1][k].y = j+1;
                            
                    buffIdx1 = buffIdx2;
                }
                
                if (cellNodes[cCell.x][j] == triggerIdx) { // If insertNode is found insert new Node
                    
                    buffIdx1 = cellNodes[cCell.x][j+1];
                    
                    float lNext = SecondNorm(nodesPos[triggerIdx] - nodesPos[ NextInCell(cCell.x, triggerIdx).x ] );
                    
                    cellNodes[cCell.x][j+1] = insertNodeIdx;
                    
                    nodeCells[insertNodeIdx][ nodeCells[insertNodeIdx][N_IDX].x ] = {cCell.x, j+1};
                    
                    wallProperties[insertNodeIdx] [ nodeCells[insertNodeIdx][N_IDX].x ] = wallProperties[triggerIdx][CellSeqInNode(cCell.x, triggerIdx)];
                    
                    float lProp = SecondNorm(newNodePos - nodesPos[ triggerIdx ] ) / lNext;
                    
                    wallProperties[insertNodeIdx] [ nodeCells[insertNodeIdx][N_IDX].x ].P_SLEN = (1-lProp) * wallProperties[triggerIdx][CellSeqInNode(cCell.x, triggerIdx)].P_SLEN;
                    
                    wallProperties[triggerIdx][CellSeqInNode(cCell.x, triggerIdx)].P_SLEN *= (lProp);
                    
                    
                    useBuffs = true;
                }
            }
            cellNodes[cCell.x][N_IDX]++;

            
        }
        
        
    }
    return insertNodeIdx;
}

//##################################################################################################################################

cl_int2 CMPlant::InsertNewWall(int cellIdx, int nodeIdx1, int nodeIdx2, float lFromNode1, float lFromNode2, cl_float8 newWallProperties){
    int insertNodeIdx1, insertNodeIdx2; // Retrurn of a function nodes of a new wall
    
    //int nodeIdx1 = cellNodes[cellIdx][nodeSN1];
    //int nodeIdx2 = cellNodes[cellIdx][nodeSN2];
    
    
    // Setting new points or sticking them to existing ones based on lengths
    // ######################################################################
    
    
    lFromNode1 = std::min(lFromNode1, SecondNorm( nodesPos[nodeIdx1] - nodesPos[NextInCell(cellIdx, nodeIdx1).x] ) );
    lFromNode2 = std::min(lFromNode2, SecondNorm( nodesPos[nodeIdx2] - nodesPos[NextInCell(cellIdx, nodeIdx2).x] ) );
    
    float lFromNextNode1 = SecondNorm( nodesPos[nodeIdx1] - nodesPos[NextInCell(cellIdx, nodeIdx1).x] ) - lFromNode1;
    float lFromNextNode2 = SecondNorm( nodesPos[nodeIdx2] - nodesPos[NextInCell(cellIdx, nodeIdx2).x] ) - lFromNode2;
    
    //float lProp1 = lFromNode1 / SecondNorm( nodesPos[nodeIdx1] - nodesPos[NextInCell(cellIdx, nodeIdx1).x] );
    //float lProp2 = lFromNode2 / SecondNorm( nodesPos[nodeIdx2] - nodesPos[NextInCell(cellIdx, nodeIdx2).x] );
    
    
    cl_float4 uTang1 = UnitTang(nodesPos[NextInCell(cellIdx, nodeIdx1).x] - nodesPos[nodeIdx1] );
    cl_float4 uTang2 = UnitTang(nodesPos[NextInCell(cellIdx, nodeIdx2).x] - nodesPos[nodeIdx2] );
    
    cl_float4 nodePos1 = nodesPos[nodeIdx1] + uTang1 * lFromNode1;
    cl_float4 nodePos2 = nodesPos[nodeIdx2] + uTang2 * lFromNode2;
    
    newWallProperties.P_SLEN = SecondNorm(nodePos1 - nodePos2) * 0.8;
    
    //std::cout << "\nlProps:"<< lProp1 << " " << lProp2 << "\n";
    
    if ( lFromNode1 > MERGE_LIM && lFromNextNode1 > MERGE_LIM ){
        insertNodeIdx1 = InsertNewNode(nodeIdx1, NextInCell(cellIdx, nodeIdx1).x, nodePos1);
    }
    else if (lFromNode1 <= MERGE_LIM){
        insertNodeIdx1 = nodeIdx1;
        nodeIdx1 = PrevInCell(cellIdx, nodeIdx1).x;
    }
    else insertNodeIdx1 = NextInCell(cellIdx,nodeIdx1).x;
    
    if ( lFromNode2 > MERGE_LIM && lFromNextNode2 > MERGE_LIM){
        insertNodeIdx2 = InsertNewNode(nodeIdx2, NextInCell(cellIdx, nodeIdx2).x, nodePos2);
    }
    else if (lFromNode2 <= MERGE_LIM){
        insertNodeIdx2 = nodeIdx2;
        nodeIdx2 = PrevInCell(cellIdx, nodeIdx2).x;
    }
    else insertNodeIdx2 = NextInCell(cellIdx, nodeIdx2).x;
    
    
    // Create a new daughter cell
    // ###########################################################################
    
    int newCellIdx = ++cellsN;
    cellNodes[newCellIdx][N_IDX] = 0;
    cellDivide[newCellIdx] = true;
    
    bool write = false;
    
    int i = 0;
    
    do {
        
        if (i == cellNodes[cellIdx][N_IDX]) i = 1;
        else i++;
        
        int cNodeIdx = cellNodes[cellIdx][i];
        
        
        if (write == true) {
            cellNodes[newCellIdx][N_IDX] ++;
            cellNodes[newCellIdx][ cellNodes[newCellIdx] [N_IDX]  ] = cNodeIdx;
            
            if (cNodeIdx != insertNodeIdx1)
                nodeCells[cNodeIdx][CellSeqInNode(cellIdx, cNodeIdx)] = {newCellIdx, cellNodes[newCellIdx][N_IDX]};
                
                
            else {
                nodeCells[cNodeIdx][N_IDX].x++;
                        
                nodeCells[cNodeIdx][ nodeCells[cNodeIdx][N_IDX].x ] = {newCellIdx, cellNodes[newCellIdx][N_IDX]};
                
                wallProperties[cNodeIdx][ nodeCells[cNodeIdx][N_IDX].x ] = wallProperties[cNodeIdx][CellSeqInNode(cellIdx, cNodeIdx)];
                
            }
            
        }
        
        
        if (write == true && cellNodes[cellIdx][i] == nodeIdx2){
            cl_int2 nNodeIdx = NextInCell(cellIdx, cNodeIdx);
            
            cellNodes[newCellIdx][N_IDX] ++;
            cellNodes[newCellIdx][ cellNodes[newCellIdx] [N_IDX]  ] = nNodeIdx.x;
            
            nodeCells[nNodeIdx.x][N_IDX].x++;
            
            nodeCells[nNodeIdx.x][ nodeCells[nNodeIdx.x][N_IDX].x ] = {newCellIdx, cellNodes[newCellIdx][N_IDX]};
            
            wallProperties[nNodeIdx.x][ nodeCells[nNodeIdx.x][N_IDX].x ] = newWallProperties;
            
            break;
        }
        
        
        if (cellNodes[cellIdx][i] == nodeIdx1) write = true;
        
        
    } while (1);
    
    // Modifying mother cell
    // ###########################################################################
    
    write = false;
    
    i = 0;
    
    cellNodes[BUF_ID][N_IDX] = 0; // Clear Buffer Cell
    
    do {
        
        if (i == cellNodes[cellIdx][N_IDX]) i = 1;
        else i++;
        
        int cNodeIdx = cellNodes[cellIdx][i];
        
        
        if (write == true) {
            cellNodes[BUF_ID][N_IDX] ++;
            cellNodes[BUF_ID][ cellNodes[BUF_ID] [N_IDX]  ] = cNodeIdx;

            nodeCells[cNodeIdx][CellSeqInNode(cellIdx, cNodeIdx)] = {cellIdx, cellNodes[BUF_ID][N_IDX]};
            
            
        }
        
        
        if (write == true && cellNodes[cellIdx][i] == nodeIdx1){
            cl_int2 nNodeIdx = NextInCell(cellIdx, cNodeIdx);
            
            cellNodes[BUF_ID][N_IDX] ++;
            cellNodes[BUF_ID][ cellNodes[BUF_ID] [N_IDX]  ] = nNodeIdx.x;
            
            int cSeq = CellSeqInNode(cellIdx, nNodeIdx.x);
        
            nodeCells[nNodeIdx.x][cSeq] = {cellIdx, cellNodes[BUF_ID][N_IDX]};
            wallProperties[nNodeIdx.x][cSeq] = newWallProperties;

        
        
            break;
        }
        
        
        if (cellNodes[cellIdx][i] == nodeIdx2) write = true;
        
    } while (1);
    
    for (i = 0; i <= cellNodes[BUF_ID][N_IDX]; i++ ){
        cellNodes[cellIdx][i] = cellNodes[BUF_ID][i];
    }
    
    
    return {insertNodeIdx1,insertNodeIdx2};
}

//##################################################################################################################################

void CMPlant::Divide(int cellIdx){
    
    cl_float8 newWallProperties = {0,0,0,0,0,0,0,0};
    cl_float4 centroid = {0,0,0,0};
    
    for (int i = 1; i <= cellNodes[cellIdx][N_IDX]; i++) {
        centroid = centroid + nodesPos[ cellNodes[cellIdx][i] ];
    }
    
    centroid = centroid / cellNodes[cellIdx][N_IDX];
    
    std::cout << "Centroid: " << centroid << "\n";
    
    float minLength1 = 100000, minLength2 = 100000;
    int nodeSN1, nodeSN2;
    float lFromNode1, lFromNode2;
    
    for (int i =1; i <= cellNodes[cellIdx][N_IDX]; i++){
        
        int cNode = cellNodes[cellIdx][i];
        cl_float4 cNodePos = nodesPos[cNode];
        
        int nNode = NextInCell(cellIdx, cNode).x;
        cl_float4 nNodePos = nodesPos[nNode];
        
        cl_float4 uTang = UnitTang(nNodePos - cNodePos);
        
        cl_float4 buff = uTang * ( (cNodePos - centroid) * uTang);
        
        cl_float4 buff2 = (cNodePos - centroid);
        
        float lToEdge = SecondNorm( buff2 - buff );
        
        if (lToEdge < minLength1 ){
            minLength1 = lToEdge;
            nodeSN1 = i;
            lFromNode1 = abs( (cNodePos-centroid) * uTang );
            
        }
        else if (lToEdge < minLength2) {
            minLength2 = lToEdge;
            nodeSN2 = i;
            lFromNode2 =abs((cNodePos-centroid) * uTang);
        }
        
    }
    
    std::cout << "lFromNode1:" << lFromNode1 << " lFromNode2:" << lFromNode2 << "\n";
    std::cout << "nodeSN1:" << nodeSN1 << " nodeSN2:" << nodeSN2 << "\n";
    
    cl_int2 newNodesIdx;
    
    newNodesIdx =InsertNewWall(cellIdx, nodeSN1, nodeSN2, lFromNode1, lFromNode2, newWallProperties);
    
    //InsertNewNode(newNodesIdx.x, newNodesIdx.y, centroid);
    
}

//##################################################################################################################################

void CMPlant::Divide2(int cellIdx){
    int insertNodeIdx1, insertNodeIdx2;
    float lFromNewNode1, lFromNewNode2;
    cl_float4 firstEdgeVector, nP1Intersect, nP2Intersect, nP3Intersect, nP4Intersect;
    cl_float4 intersectPos;
    float finalR;
    cl_float4 finalIntersect;
    cl_float4 finalnP1, finalnP3;
    bool flag = false;
    
    // Finding Area of the Cell
    // ########################################################
    float cellArea = 0;
    
    for (int i =1; i <= cellNodes[cellIdx][N_IDX]; i++){
        int cNode = cellNodes[cellIdx][i];
        cellArea += Cross( nodesPos[cNode], nodesPos[NextInCell(cellIdx, cNode).x ] - nodesPos[cNode] ).z;
        
    }
    cellArea *= 0.5;
    
    float minLength = 1000000;
    
    //
    // ########################################################
    for (int i = 1; i < cellNodes[cellIdx][N_IDX]; i++){
        for (int j = i+1; j <= cellNodes[cellIdx][N_IDX]; j++){
            
            int nodeIdx1 = cellNodes[cellIdx][i];
            int nodeIdx2 = cellNodes[cellIdx][j];
            
            cl_float4 nP1 = nodesPos[nodeIdx1];
            cl_float4 nP2 = nodesPos[ NextInCell(cellIdx, nodeIdx1).x ];
            cl_float4 nP3 = nodesPos[nodeIdx2];
            cl_float4 nP4 = nodesPos[ NextInCell(cellIdx, nodeIdx2).x ];
            
            intersectPos.x = ( (nP1.x * nP2.y - nP1.y * nP2.x ) * (nP3.x - nP4.x) - ( nP1.x - nP2.x ) * (nP3.x * nP4.y - nP3.y * nP4.x) ) / ( (nP1.x - nP2.x) * (nP3.y - nP4.y) - (nP1.y - nP2.y) * (nP3.x - nP4.x)  );
            intersectPos.y = ( (nP1.x * nP2.y - nP1.y * nP2.x ) * (nP3.y - nP4.y) - ( nP1.y - nP2.y ) * (nP3.x * nP4.y - nP3.y * nP4.x) ) / ( (nP1.x - nP2.x) * (nP3.y - nP4.y) - (nP1.y - nP2.y) * (nP3.x - nP4.x)  );
            intersectPos.z = 0;
            
            firstEdgeVector = nP2 - nP1;
            nP1Intersect = intersectPos - nP1;
            nP2Intersect = intersectPos - nP2;
            nP3Intersect = intersectPos - nP3;
            nP4Intersect = intersectPos - nP4;
            
            int sign;
            
            if (abs( firstEdgeVector * nP1Intersect ) > abs (firstEdgeVector * nP2Intersect))
                sign = ( firstEdgeVector * nP1Intersect > 0) ? 1 : -1;
            else
                sign = ( firstEdgeVector * nP2Intersect > 0) ? 1 : -1;
            
            float mArea = 0;
            float R = 0;
            float angle;
            
            angle = acos( (firstEdgeVector * (nP3-nP4) ) / ( SecondNorm(firstEdgeVector) * SecondNorm((nP3-nP4)) ) );
            
            if ( angle >= 0.05  ){
                
                if (sign > 0){
                    //std::cout << "i,j: " << i << " " << j << "\n";
                    int cNode = NextInCell(cellIdx, nodeIdx1).x;
                    
                    mArea += Cross ( intersectPos, intersectPos - nodesPos[cNode] ).z;
                    
                    if (cNode !=nodeIdx2)
                        
                        do {
                            mArea += Cross( nodesPos[cNode], nodesPos[cNode] - nodesPos[ NextInCell(cellIdx, cNode).x]  ).z;
                            //std::cout << cNode << " " << NextInCell(cellIdx, cNode).x << "\n";
                            cNode = NextInCell(cellIdx, cNode).x;
                        } while (cNode != nodeIdx2);
                    
                    mArea += Cross (nodesPos[nodeIdx2], nodesPos[nodeIdx2] - intersectPos ).z;
                    
                    mArea *= 0.5;
                    
                    R = sqrt( (cellArea + 2 * mArea) / angle );
                    
                    //std::cout << sign << " " << R << " " << R*angle << " "<< cellArea << " " << mArea << " " << angle << "||" << SecondNorm(nP1Intersect) << " " << SecondNorm(nP4Intersect) << " " << SecondNorm(nP2Intersect) << " " << SecondNorm(nP3Intersect) << "\n";
                    if (R >= std::max(SecondNorm(nP3Intersect), SecondNorm(nP2Intersect))  && R <= std::min( SecondNorm(nP1Intersect), SecondNorm(nP4Intersect)) && (R * angle) < minLength  ){
                        minLength = R * angle;
                        
                        insertNodeIdx1 = nodeIdx1;
                        insertNodeIdx2 = nodeIdx2;
                        
                        lFromNewNode1 = SecondNorm(nP1Intersect) - R;
                        lFromNewNode2 = R - SecondNorm(nP3Intersect);
                        
                        finalR = R;
                        finalIntersect = intersectPos;
                        finalnP1 = nP1;
                        finalnP3 = nP3;
                        flag = false;
                    }
                    
                }
                else {
                    int cNode = NextInCell(cellIdx, nodeIdx2).x;
                    //std::cout << "i,j: " << i << " " << j << "\n";
                    
                    mArea += Cross ( intersectPos, intersectPos - nodesPos[cNode] ).z;
                    
                    if (cNode !=nodeIdx1)
                        
                        do {
                            mArea += Cross( nodesPos[cNode] , nodesPos[cNode] - nodesPos[ NextInCell(cellIdx, cNode).x]  ).z;
                            //std::cout << cNode << " " << NextInCell(cellIdx, cNode).x << "\n";
                            cNode = NextInCell(cellIdx, cNode).x;
                        } while (cNode != nodeIdx1);
                    
                    mArea += Cross(nodesPos[nodeIdx1], nodesPos[nodeIdx1] - intersectPos ).z;
                    
                    mArea *= 0.5;
                    
                    R = sqrt( (cellArea + 2 * mArea) / angle );
                    //std::cout << sign << " " << R << " " << R*angle << " "<< cellArea << " " << mArea << " " << angle << "||" << SecondNorm(nP1Intersect) << " " << SecondNorm(nP4Intersect) << " " << SecondNorm(nP2Intersect) << " " << SecondNorm(nP3Intersect) << "\n";
                    
                    if (R >= std::max(SecondNorm(nP1Intersect), SecondNorm(nP4Intersect))  && R <= std::min( SecondNorm(nP2Intersect), SecondNorm(nP3Intersect)) && (R * angle) < minLength  ){
                        minLength = R * angle;
                        
                        insertNodeIdx1 = nodeIdx1;
                        insertNodeIdx2 = nodeIdx2;
                        
                        lFromNewNode1 = R - SecondNorm(nP1Intersect);
                        lFromNewNode2 = SecondNorm(nP3Intersect) - R;
                        
                        finalR = R;
                        finalIntersect = intersectPos;
                        finalnP1 = nP1;
                        finalnP3 = nP3;
                        flag = false;
                    }
                }
                
            }
            
            else{
                
                
                //Finding D
                cl_float4 uNorm = UnitTang( Cross(nP1-nP2, uZ)  );

                cl_float4 uTang = UnitTang( nP2-nP1 );
                
                float d = abs(uNorm * (nP3 - nP1));
                
                if (d < minLength){
                    
                    int cNode = NextInCell(cellIdx, nodeIdx2).x;
                    
                    intersectPos = nP1 + uNorm * d;
                    
                    mArea += Cross ( intersectPos, intersectPos - nodesPos[cNode] ).z;
                    
                    if (cNode != nodeIdx1) {
                        
                        do {
                            mArea += Cross( nodesPos[cNode] , nodesPos[cNode] - nodesPos[ NextInCell(cellIdx, cNode).x]  ).z;
                            //std::cout << cNode << " " << NextInCell(cellIdx, cNode).x << "\n";
                            cNode = NextInCell(cellIdx, cNode).x;
                        } while (cNode != nodeIdx1);
                    }
                    
                    mArea += Cross(nodesPos[nodeIdx1], nodesPos[nodeIdx1] - intersectPos ).z;
                    
                    mArea *= 0.5;
                    
                    R = round(( (cellArea * 0.5) + mArea ) / d);
                    
                    
                    
                    float l1 = round(uTang * (nP4 - nP1)), l2 = round(uTang * (nP3 - nP1));
                    float l3 = round(SecondNorm(nP2-nP1));
                    
                    if ( R>=l1 && R<=l2 && R <= l3 && R>=0 ){
                        
                        minLength = d;
                        
                        insertNodeIdx1 = nodeIdx1;
                        insertNodeIdx2 = nodeIdx2;
                        
                        lFromNewNode1 = R;
                        lFromNewNode2 = l2 - R;
                        flag = true;
                        /*
                        std::cout << "cellIDX " << cellIdx << "\n";
                        std::cout << "i: " << i << "\n";
                        std::cout << "j: " << j << "\n";
                        std::cout << "uNorm: " << uNorm << "\n";
                        std::cout << "uTang: " << uTang << "\n";
                        std::cout << "d: " << d << "\n";
                        std::cout << "minLength: " << minLength << "\n";
                        std::cout << "nodeIdx1: " << nodeIdx1<< "\n";
                        std::cout << "nodeIdx2: " << nodeIdx2<< "\n";
                        std::cout << "intersectPos: " << intersectPos<< "\n";
                        std::cout << "mArea: " << mArea<< "\n";
                        std::cout << "cellArea : " << cellArea<< "\n";
                        std::cout << "R: " << R<< "\n";
                        std::cout << "l1: " << l1<< "\n";
                        std::cout << "l2: " << l2<< "\n";
                        std::cout << "l3: " << l3<< "\n";
                        */
                    }
                    
                }
                
            }
            
        }
    }
    //std::cout << "cellIdx: " << cellIdx << " insertNodeIdx1: " << insertNodeIdx1 << " insertNodeIdx2: " << insertNodeIdx2 << " lFromNode1: " << lFromNewNode1 << " lFromNode2: " << lFromNewNode2 << "\n";
    
    
    cl_int2 newNodeIdx;
    newNodeIdx = InsertNewWall(cellIdx, insertNodeIdx1, insertNodeIdx2, lFromNewNode1, lFromNewNode2, newWallProperties);
    
    
    
    /*
    if (flag){
    cl_float4 uTnP1Int = UnitTang(finalIntersect - finalnP1);
    cl_float4 uTnP3Int = UnitTang(finalIntersect - finalnP3);
     
    cl_float4 uMiddleNode = UnitTang(uTnP1Int + uTnP3Int);
     
    InsertNewNode(newNodeIdx.x, newNodeIdx.y, finalIntersect + uMiddleNode * -finalR);
    }
    */
}


void CMPlant::CalculateForces(){
    
    for (int i = 1; i<=nodesN; i++) {
        angularForce[i] = {0,0,0,0};
    }
    
    
    for (int cNode = 1; cNode <= nodesN; cNode++) {
        nodesForce[cNode] = {0,0,0,0};
        
        // Finding Pressure Forces
        for (int j = 1; j <= nodeCells[cNode][N_IDX].x; j++){

            int cCell = nodeCells[cNode][j].x;
            
            float Angle = GetAngle(cNode, cCell);
            
            int nextNodeIdx = NextInCell(cCell, cNode).x;
            int prevNodeIdx = PrevInCell(cCell, cNode).x;
            
            cl_float4 nextEdge = nodesPos[ nextNodeIdx ] - nodesPos[cNode];
            cl_float4 prevEdge = nodesPos[cNode] - nodesPos[ prevNodeIdx ];
            
            
            cl_float4 uTangNext = UnitTang( nextEdge );
            cl_float4 uTangPrev = UnitTang( prevEdge );
            
            
            cl_float4 NormNext = Cross( nextEdge, uZ );
            cl_float4 NormPrev = Cross( prevEdge, uZ );
            
            cl_float4 uNormNext = UnitTang(NormNext);
            cl_float4 uNormPrev = UnitTang(NormPrev);
            
                       
            float lNextEdge = SecondNorm(nextEdge);
            float lPrevEdge = SecondNorm(prevEdge);
            
           
            
            nodesForce[cNode] = nodesForce[cNode] + uNormNext * lNextEdge * pressure;
            

            //std::cout << nodesForce[i] << "\n";
            
            nodesForce[cNode] = nodesForce[cNode] + uNormPrev * lPrevEdge * pressure;
            
            
            //std::cout << nodesForce[i] << "\n";
            
            cl_float8 wP = wallProperties[cNode][j];
            
            nodesForce[cNode] = nodesForce[cNode] + uTangNext * wP.P_AREA * wP.P_YMOD * ( (lNextEdge - wP.P_SLEN) / wP.P_SLEN );
            
            
            //std::cout << i << " next) " << (lNextEdge - wP.P_SLEN) / wP.P_SLEN << "\n";
        
            wP = wallProperties[prevNodeIdx][ CellSeqInNode(cCell, prevNodeIdx) ];
            
            nodesForce[cNode] = nodesForce[cNode] + uTangPrev * -wP.P_AREA * wP.P_YMOD * ( (lPrevEdge - wP.P_SLEN) / wP.P_SLEN );
                
            //std::cout << i << " prev) " << (lPrevEdge - wP.P_SLEN) / wP.P_SLEN << "\n";
            
            
            nodesForce[cNode] = nodesForce[cNode] + uTangNext * (40 / tan( Angle/2 ));
            nodesForce[cNode] = nodesForce[cNode] + uTangPrev * (-40 / tan(Angle/2));
            
            }
        }
}


void CMPlant::Step(){
    
    int N = 100;
    
    for(int i = 1; i<=nodesN; i++)
        buffNodesPos[i] = nodesPos[i];
    
    for (int c=0; c<N; c++ ){
        CalculateForces();
        
        for (int cNode = 1; cNode <= nodesN; cNode++) {
            buffNodesPos[cNode] = buffNodesPos[cNode] + nodesForce[cNode] * (dt / kEnv);
        }
        
    }
    
    for (int cNode = 1; cNode <= nodesN; cNode++) {
        
        for (int i = 1; i <= nodeCells[cNode][N_IDX].x; i++){
            int cCell = nodeCells[cNode][i].x;
            
            float lNext = SecondNorm( nodesPos[cNode] - nodesPos[ NextInCell(cCell, cNode).x ] );
            
            cl_float8 wP = wallProperties[cNode][i];
            
            wallProperties[cNode][i].P_SLEN = fmax( wP.P_SLEN + (wP.P_YMOD / wP.P_VISC ) * (lNext - wP.P_SLEN) * dt , wP.P_SLEN);
    
            
            bool changeS = true;
            for (int j = 1; j <= nodeCells[cNode][N_IDX].x; j++){
                
                int jCell = nodeCells[cNode][j].x;
                
                if ( i!=j && PrevInCell(cCell, cNode).x == NextInCell(jCell, cNode).x ) {
                    changeS = false;
                    break;
                }
            }
            
            if (changeS){
                int pNode = PrevInCell(cCell, cNode).x;
                float lPrev = SecondNorm( nodesPos[cNode] - nodesPos[pNode] );
                
                int cCellPosInpNode = CellSeqInNode(cCell, pNode);
                cl_float8 wP = wallProperties[ pNode] [ cCellPosInpNode ];
                wallProperties[ pNode ][ cCellPosInpNode ].P_SLEN = fmax( wP.P_SLEN + (wP.P_YMOD / wP.P_VISC) * (lPrev - wP.P_SLEN) * dt, wP.P_SLEN );
            }
            
            
        }
        
    }
    
    for(int i = 1; i<=nodesN; i++)
        nodesPos[i] = buffNodesPos[i];
    
}


//##################################################################################################################################


void CMPlant::SetCellParams(int cellIdx, float visc, float ymod){
    
    for (int i = 1; i <= cellNodes[cellIdx][N_IDX]; i++) {
        int cNode = cellNodes[cellIdx][i];
        
        for (int j = 1; j <= nodeCells[cNode][N_IDX].x; j++){
            int cCell = nodeCells[cNode][j].x;
            
            if (cCell == cellIdx) {
                wallProperties[cNode][j].P_VISC = visc;
                wallProperties[cNode][j].P_YMOD = ymod;
            }
            
            int PrevNodeIncCell = PrevInCell(cCell, cNode).x;
            
            if (PrevNodeIncCell == NextInCell(cellIdx, cNode).x) {
                
                int cellSeq = CellSeqInNode(cCell, PrevNodeIncCell);
                
                wallProperties[PrevNodeIncCell][cellSeq].P_VISC = visc;
                wallProperties[PrevNodeIncCell][cellSeq].P_YMOD = ymod;
                
            }
            
            
        }
    }

    
}





//##################################################################################################################################

void CMPlant::GLDisplay(int disp_w, int disp_h, bool displayNodes, bool printNumbers, bool displayForces, float scale){
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_FRONT_FACE);
    
    
    glPointSize(11);
    glLineWidth(5);
    
    

    
    for (int i = 1; i <= cellsN; i++){
        glBegin(GL_POLYGON);
        glColor3f(0.34, 0.86, 0.37);
        if (i == guiCellSelected) glColor3f(0.17, 0.43, 0.19);
        for (int j = 1; j <= cellNodes[i][N_IDX]; j++ ){
            glVertex2f( disp_w/2 +(nodesPos[ cellNodes[i][j] ].x) * scale, disp_h/2 + (nodesPos[ cellNodes[i][j] ].y)*scale   );
            
        }
        
        glEnd();
    }
    
    glPointSize(1);
    glLineWidth(3);
    for (int i = 1; i <= cellsN; i++){
        glBegin(GL_LINE_LOOP);
        
        for (int j = 1; j <= cellNodes[i][N_IDX]; j++ ){
            glColor3f(0.25, 0.61, 0.27);
            glVertex2f( disp_w/2 + (nodesPos[ cellNodes[i][j] ].x) * scale, disp_h/2 + (nodesPos[ cellNodes[i][j] ].y)*scale   );
            
        }
        
        glEnd();
    }
    
    glPointSize(5);
    
    
    for (int i = 1; i <= cellsN; i++){
        
        
        for (int j = 1; j <= cellNodes[i][N_IDX]; j++ ){
            
            if (displayNodes){
                glBegin(GL_POINTS);
                glColor3f(1, 0, 0);
                glVertex2f( disp_w/2 + (nodesPos[ cellNodes[i][j] ].x) * scale, disp_h/2 + (nodesPos[ cellNodes[i][j] ].y) * scale   );
                glEnd();
            }
            
            if (printNumbers){
                char text[5];
                sprintf(text,"%d", cellNodes[i][j] );
                GLPrintText(disp_w/2 + (nodesPos[ cellNodes[i][j] ].x) * scale, disp_h/2 + (nodesPos[ cellNodes[i][j] ].y + 8) * scale, text);
            }
            
            if (displayForces) GLDisplayForces(disp_w, disp_h, scale);
        }
        
        
    }

 
    
    glutSwapBuffers();

}

//##################################################################################################################################

void CMPlant::GLDisplayForces(int disp_w, int disp_h, float scale){
    
    
    glLineWidth(1);
    for (int i = 1; i <= nodesN; i++){
        
        glBegin(GL_LINES);
        glColor3f(1,0,0);
        glVertex2f(disp_w/2 + nodesPos[i].x * scale,  disp_h/2 + nodesPos[i].y * scale);
        glVertex2f(disp_w/2 + (nodesPos[i].x + 50 * nodesForce[i].x) * scale, disp_h/2 + (nodesPos[i].y + 50 * nodesForce[i].y)*scale);
        
        glEnd();
        
        glBegin(GL_LINES);
        glColor3f(0,0,1);
        glVertex2f(disp_w/2 + nodesPos[i].x * scale,  disp_h/2 + nodesPos[i].y * scale);
        glVertex2f(disp_w/2 + (nodesPos[i].x + 50 * angularForce[i].x) * scale, disp_h/2 + (nodesPos[i].y + 50 * angularForce[i].y)*scale);
        
        glEnd();
    }
    
}


//##################################################################################################################################

void CMPlant::GLPrintText(float x, float y, char *text){
    
    glColor3f(1,0,0);
    glRasterPos2f(x, y);
    
    while (*text) {
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *text);
        text++;
    }
}


//##################################################################################################################################

void CMPlant::CLCheckErr(cl_int err, const char *name){
    
    if ( err != CL_SUCCESS){
        std::cerr << "CL_ERROR: " << name << "(" << err << ")" << "\n";
        exit(EXIT_FAILURE);
    }
}

//##################################################################################################################################

void CMPlant::CLCreateContext(){
    
    cl_int err;
    
    cl_platform_id platformID;
    
    cl_uint platformN;
    
    err = clGetPlatformIDs(1, &platformID, &platformN);
    CLCheckErr(err, "clGetPlatformID");
    
    err = clGetDeviceIDs(platformID, CL_DEVICE_TYPE_CPU, 1, &deviceID, NULL);
    CLCheckErr(err, "clGetDeviceID");
    
 
    
    cl_ulong maxMemSize;
    size_t max_size;
    cl_uint add;
    
    err = clGetDeviceInfo(deviceID, CL_DEVICE_ADDRESS_BITS, sizeof(cl_uint), &add, NULL);
    
    std::cout << "Max size: " << add;
    CLCheckErr(err, "clGetDeviceInfo");
    
    clContext = clCreateContext(NULL, 1, &deviceID, NULL, NULL, &err);
    CLCheckErr(err, "clCreateContext");
    
    clQueue = clCreateCommandQueue(clContext, deviceID, NULL, &err);
    CLCheckErr(err, "clCreateCommandQueue");
    
    cl_buffNodesPos = clCreateBuffer(clContext, CL_MEM_READ_WRITE, maxNodes * sizeof(cl_float4), NULL, &err);
    CLCheckErr(err, "clCreateBuffer(cl_buffNodesPos)");
    
    cl_cellNodes = clCreateBuffer(clContext, CL_MEM_READ_ONLY, maxNodes * maxNodesPerCell * sizeof(int), NULL, &err);
    CLCheckErr(err, "clCreateBuffer(cl_cellNodes)");
    
    cl_nodeCells = clCreateBuffer(clContext, CL_MEM_READ_ONLY, maxNodes * maxCellsPerNode * sizeof(cl_int2) , NULL, &err);
    CLCheckErr(err, "clCreateBuffer(cl_nodeCells)");
    
    cl_nodesForce = clCreateBuffer(clContext, CL_MEM_READ_WRITE, maxNodes * sizeof(cl_float4), NULL, &err);
    CLCheckErr(err, "clCreateBuffer(cl_nodeForce)");
    
    cl_wallProperties = clCreateBuffer(clContext, CL_MEM_READ_ONLY, maxNodes * maxCellsPerNode * sizeof(cl_float8), NULL, &err);
    CLCheckErr(err, "clCreateBuffer(cl_wallProperties)");
    
    cl_Setup = clCreateBuffer(clContext, CL_MEM_READ_ONLY, 7 * sizeof (float), NULL, &err);
    CLCheckErr(err, "clCreateBuffer(cl_Setup)");
    
    FILE *fp;
    char *source;
    size_t source_size;
    
    fp = fopen("/Users/mdelmans/Desktop/CMPlant/Kernels/calculateForce.cl", "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );
    
    pCalculateForces =clCreateProgramWithSource(clContext, 1, (const char **)&source, (const size_t *)&source_size, &err);
    CLCheckErr(err, "clCreateProgramWithSource");
    
    fp = fopen("/Users/mdelmans/Desktop/CMPlant/Kernels/getNewPos.cl", "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );

    pGetNewPos = clCreateProgramWithSource(clContext, 1, (const char **)&source, (const size_t *)&source_size, &err);
    CLCheckErr(err, "clCreateProgramWithSource");
    
    fp = fopen("/Users/mdelmans/Desktop/CMPlant/Kernels/getNewStruct.cl", "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );
    
    pGetNewStruct = clCreateProgramWithSource(clContext, 1, (const char **)&source, (const size_t *)&source_size, &err);
    CLCheckErr(err, "clCreateProgramWithSource");

    err = clBuildProgram(pCalculateForces, 1, &deviceID, NULL, NULL, NULL);
    CLCheckErr(err, "clBuild(pCalculateForces)");
    
    err = clBuildProgram(pGetNewPos, 1, &deviceID, NULL, NULL, NULL);
    CLCheckErr(err, "clBuild(pGetNewPos)");
    
    err = clBuildProgram(pGetNewStruct, 1, &deviceID, NULL, NULL, NULL);
    //CLCheckErr(err, "clBuild(pGetNewStruct)");
    
    char buildLog[1000];
    
    
    
    err = clGetProgramBuildInfo(pCalculateForces, deviceID, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
    CLCheckErr(err, "clGetProgramBuildInfo");
    
    std::cout << "pCalculateForces Log:\n" << buildLog << "\n";
    
    err = clGetProgramBuildInfo(pGetNewPos, deviceID, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
    CLCheckErr(err, "clGetProgramBuildInfo");
    
    std::cout << "pGetNewPos Log:\n" << buildLog << "\n";
    
    err = clGetProgramBuildInfo(pGetNewStruct, deviceID, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
    CLCheckErr(err, "clGetProgramBuildInfo");
    
    std::cout << "pGetNewPos Log:\n" << buildLog << "\n";
    
    float setup[] = {pressure, (float)maxNodesPerCell, (float)maxCellsPerNode, (float)(nodesN), (float)cellsN, kEnv, dt};
    err = clEnqueueWriteBuffer(clQueue, cl_Setup, CL_TRUE, 0, 7 * sizeof(float), setup, 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(setup)");
    
    kCalculateForces = clCreateKernel(pCalculateForces, "calculate_forces", &err);
    CLCheckErr(err, "clCreateKernel");
    
    kGetNewPos = clCreateKernel(pGetNewPos, "get_new_pos", &err);
    CLCheckErr(err, "clCreateKernel");
    
    kGetNewStruct = clCreateKernel(pGetNewStruct, "get_new_struct", &err);
    CLCheckErr(err, "clCreateKernel");
    
    
    err = clEnqueueWriteBuffer(clQueue, cl_buffNodesPos, CL_TRUE, 0, (nodesN+1) * sizeof(cl_float4), &nodesPos[0], 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(buffNodesPos)");
    
    err = clEnqueueWriteBuffer(clQueue, cl_cellNodes, CL_TRUE, 0, (cellsN+1) * maxNodesPerCell * sizeof(int), &cellNodes[0][0], 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(cellNodes)");
    
    err = clEnqueueWriteBuffer(clQueue, cl_nodeCells, CL_TRUE, 0, (nodesN+1) * maxCellsPerNode * sizeof(cl_int2), &nodeCells[0][0], 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(nodeCells)");
    
    err = clEnqueueWriteBuffer(clQueue, cl_wallProperties, CL_TRUE, 0, (nodesN+1) * maxCellsPerNode * sizeof(cl_float8), &wallProperties[0][0], 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(wallProperties)");

}

//##################################################################################################################################

void CMPlant::CLCalculateForces(int coreSize){
    cl_int err;
    
    err = clSetKernelArg(kCalculateForces, 0, sizeof(cl_mem), (void *)&cl_nodesForce);
    CLCheckErr(err, "clSetKernelArg0");
    
    err = clSetKernelArg(kCalculateForces, 1, sizeof(cl_mem), (void *)&cl_buffNodesPos);
    CLCheckErr(err, "clSetKernelArg1");
    
    err = clSetKernelArg(kCalculateForces, 2, sizeof(cl_mem), (void *)&cl_cellNodes);
    CLCheckErr(err, "clSetKernelArg2");
    
    err = clSetKernelArg(kCalculateForces, 3, sizeof(cl_mem), (void *)&cl_nodeCells);
    CLCheckErr(err, "clSetKernelArg3");
    
    err = clSetKernelArg(kCalculateForces, 4, sizeof(cl_mem), (void *)&cl_wallProperties);
    CLCheckErr(err, "clSetKernelArg4");
    
    err = clSetKernelArg(kCalculateForces, 5, sizeof(cl_mem), (void *)&cl_Setup);
    CLCheckErr(err, "clSetKernelArg5");
    
    size_t localItemSize = coreSize;
    size_t globalItemSize = (int)ceil((nodesN+1) / (float)localItemSize) * localItemSize;
    //size_t globalItemSize = 1024;
    
    err = clEnqueueNDRangeKernel(clQueue, kCalculateForces, 1, NULL, &globalItemSize, &localItemSize, 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueNDRangeKernel(forces)");
    clFinish(clQueue);
    
}

void CMPlant::CLStep(){
    
    int N = 100;
    
    for (int c=0; c<N; c++ ){
        CLCalculateForces(32);
        CLGetNewPos();
    }
     
    CLGetNewStruct();
}

void CMPlant::CLWriteStructure(){
    
    cl_int err;
    
    err = clEnqueueWriteBuffer(clQueue, cl_cellNodes, CL_TRUE, 0, (cellsN+1) * maxNodesPerCell * sizeof(int), &cellNodes[0][0], 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(cellNodes)");
    
    err = clEnqueueWriteBuffer(clQueue, cl_nodeCells, CL_TRUE, 0, (nodesN+1) * maxCellsPerNode * sizeof(cl_int2), &nodeCells[0][0], 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(nodeCells)");
    
    err = clEnqueueWriteBuffer(clQueue, cl_wallProperties, CL_TRUE, 0, (nodesN+1) * maxCellsPerNode * sizeof(cl_float8), &wallProperties[0][0], 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueWriteBuffer(wallProperties)");

}


void CMPlant::CLGetNewPos(){
    cl_int err;
    
    err = clSetKernelArg(kGetNewPos, 0, sizeof(cl_mem), (void *)&cl_nodesForce);
    CLCheckErr(err, "clSetKernelArg0");
    
    err = clSetKernelArg(kGetNewPos, 1, sizeof(cl_mem), (void *)&cl_buffNodesPos);
    CLCheckErr(err, "clSetKernelArg1");
    
    err = clSetKernelArg(kGetNewPos, 2, sizeof(cl_mem), (void *)&cl_Setup);
    CLCheckErr(err, "clSetKernelArg2");
    
    size_t localItemSize = 32;
    size_t globalItemSize = (int)ceil((nodesN+1) / (float)localItemSize) * localItemSize;
    err = clEnqueueNDRangeKernel(clQueue, kGetNewPos, 1, NULL, &globalItemSize, &localItemSize, 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueNDRangeKernel(GetNewPos)");
}

void CMPlant::CLReadPos(){
    
    cl_int err;
    
    err = clEnqueueReadBuffer(clQueue, cl_buffNodesPos, CL_TRUE, 0, (nodesN+1) * sizeof(cl_float4), nodesPos, 0, NULL, NULL);
    
    /*
    std::cout << "nodesPos:\n";
    for (int i = 1; i<=nodesN; i++) {
        std::cout << nodesPos[i] << "\n";
    }
    std::cout << "\n\n";
    */
}

void CMPlant::CLReadForce(){
    cl_int err;
    
    err = clEnqueueReadBuffer(clQueue, cl_nodesForce, CL_TRUE, 0, (nodesN+1) * sizeof(cl_float4), nodesForce, 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueReadBuffer(forces)");
    
    
    std::cout << "nodesForce:\n";
    for (int i = 1; i<=nodesN; i++) {
        std::cout << nodesForce[i] << "\n";
    }
    std::cout << "\n\n";
    

    
    
}

void CMPlant::CLGetNewStruct(){
    cl_int err;
    
    err = clSetKernelArg(kGetNewStruct, 0, sizeof(cl_mem), (void *)&cl_nodesForce);
    CLCheckErr(err, "clSetKernelArg0");
    
    err = clSetKernelArg(kGetNewStruct, 1, sizeof(cl_mem), (void *)&cl_buffNodesPos);
    CLCheckErr(err, "clSetKernelArg1");
    
    err = clSetKernelArg(kGetNewStruct, 2, sizeof(cl_mem), (void *)&cl_cellNodes);
    CLCheckErr(err, "clSetKernelArg2");
    
    err = clSetKernelArg(kGetNewStruct, 3, sizeof(cl_mem), (void *)&cl_nodeCells);
    CLCheckErr(err, "clSetKernelArg3");
    
    err = clSetKernelArg(kGetNewStruct, 4, sizeof(cl_mem), (void *)&cl_wallProperties);
    CLCheckErr(err, "clSetKernelArg4");
    
    err = clSetKernelArg(kGetNewStruct, 5, sizeof(cl_mem), (void *)&cl_Setup);
    CLCheckErr(err, "clSetKernelArg5");
    
    size_t localItemSize = 32;
    size_t globalItemSize = (int)ceil((nodesN+1) / (float)localItemSize) * localItemSize;
    //size_t globalItemSize = 1024;
    
    err = clEnqueueNDRangeKernel(clQueue, kGetNewStruct, 1, NULL, &globalItemSize, &localItemSize, 0, NULL, NULL);
    CLCheckErr(err, "clEnqueueNDRangeKernel(forces)");

}

void CMPlant::CLStop(){
    cl_int err;
    err = clReleaseKernel(kCalculateForces);
    err = clReleaseKernel(kGetNewPos);
    err = clReleaseKernel(kGetNewStruct);
    
        
    err = clReleaseProgram(pGetNewPos);
    err = clReleaseProgram(pCalculateForces);
    err = clReleaseProgram(pGetNewStruct);
    
    err = clReleaseMemObject(cl_buffNodesPos);
    err = clReleaseMemObject(cl_nodeCells);
    err = clReleaseMemObject(cl_cellNodes);
    err = clReleaseMemObject(cl_nodesForce);
    err = clReleaseMemObject(cl_nodesPos);
    err = clReleaseMemObject(cl_Setup);
    
    err = clReleaseCommandQueue(clQueue);
    err = clReleaseContext(clContext);
    err = clReleaseDevice(deviceID);

}

//##################################################################################################################################


int CMPlant::GUICursorInCell(int x, int y, float scale){
    
    float rX = x * scale;
    float rY = y * scale;
    bool inCell = false;
    
    for (int cCell = 1; cCell <= cellsN; cCell++){
        
        for (int i = 1 ; i <= cellNodes[cCell][N_IDX]; i++) {
            int cNode = cellNodes[cCell][i];
            float X1 = nodesPos[cNode].x;
            float Y1 = nodesPos[cNode].y;
            
            int nNode = NextInCell(cCell, cNode).x;
            
            float X2 = nodesPos[nNode].x;
            float Y2 = nodesPos[nNode].y;
            
            if (( (Y1 < rY && Y2 >=rY) || ( Y2 < rY && Y1>=rY )  ) && ( X1 + (rY - Y1) / (Y2-Y1) * (X2-X1) < rX ) )
                inCell = !inCell;
            
        }
        
        if (inCell){
            guiCellSelected = cCell;
            return cCell;
        }
        
    }
    guiCellSelected = -1;
    return -1;
    
}

void CMPlant::PrintAll(){
    int i,j;
    
    std::cout << "cellsN: " << cellsN << "\n\n";
    
    std::cout << "nodesPos:\n\n";
    
    for (i = 1; i <= nodesN; i++){
        std::cout << i << ")" << nodesPos[i].x << "," << nodesPos[i].y << "," << nodesPos[i].z << "\n";
    }
    
    std::cout << "\n\ncellNodes\n\n";
    
    for (i = 1;  i<= cellsN; i++ ){
        std::cout << i << ") ";
        for (j = 1; j <= cellNodes[i][N_IDX]; j++){
            std::cout << cellNodes[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "\n\nnodeCells\n\n";
    
    for (i = 1; i <= nodesN; i++) {
        
        std::cout << i << ")";
        for (j = 1; j <= nodeCells[i][0].x; j++) {
            std::cout << nodeCells[i][j].x << "(" << nodeCells[i][j].y << ")";
        }
        std::cout << "\n";
    }
    
    std::cout << "\n\nwallProperties\n\n";
    
    for (i = 1; i <= nodesN; i++) {
        
        std::cout << i << ")";
        for (j = 1; j <= nodeCells[i][0].x; j++) {
            std::cout << wallProperties[i][j] << " " ;
        }
        std::cout << "\n";
    }

    
}


