//
//  CMPlant.h
//  CMPlant
//
//  Created by Mihails Delmans on 25/07/2013.
//  Copyright (c) 2013 Mihails Delmans. All rights reserved.
//


#ifndef __CMPlant__CMPCell__
#define __CMPlant__CMPCell__

#include <iostream>
#include <OpenCL/opencl.h>

#define MAX_NODES 5000
#define MAX_CELLS 1000
#define MAX_NODES_PER_CELL  50
#define MAX_CELLS_PER_NODE  10

class CMPlant
{
private:
    int maxNodes = MAX_NODES;
    int maxCells = MAX_CELLS;
    
    int maxNodesPerCell = MAX_NODES_PER_CELL;
    int maxCellsPerNode = MAX_CELLS_PER_NODE;
    cl_float4 uZ = {0,0,1,0};
    
    cl_int2 PrevInCell(int cellIdx, int nodeIdx);
    cl_int2 NextInCell(int cellIdx, int nodeIdx);
    
    int CellSeqInNode (int cellIdx, int nodeIdx);
    
    float kEnv = 1;
    float dt = 1;
    
    cl_float8 newWallProperties;
    
    cl_device_id deviceID;
    cl_context clContext;
    cl_command_queue clQueue;
    
    cl_mem cl_nodesPos;
    cl_mem cl_buffNodesPos;
    cl_mem cl_cellNodes;
    cl_mem cl_nodeCells;
    cl_mem cl_nodesForce;
    cl_mem cl_wallProperties;
    cl_mem cl_Setup;
    
    cl_program pCalculateForces;
    cl_program pGetNewPos;
    cl_program pGetNewStruct;
    
    cl_kernel kCalculateForces;
    cl_kernel kGetNewPos;
    cl_kernel kGetNewStruct;

public:
    float GetAngle(int cNode, int cellIdx);
    int cellsN = 0;
    int nodesN = 0;
    float pressure = 1;

    cl_float4 nodesPos[MAX_NODES]; // Stores position vector of each node (x,y,z,0)
    
    cl_float4 buffNodesPos[MAX_NODES];
    
    int cellNodes[MAX_CELLS][MAX_NODES_PER_CELL]; // Stores indexes of all nodes that belong to the cell[i], cellNodes[i][0] is the count of nodes belonging to cell[i]. cellNodes[i] = [N, NIdx1, NIdx2, ...  ]
    
    cl_int2 nodeCells[MAX_NODES][MAX_CELLS_PER_NODE]; // Stores (cell index, index inside the cellNodes). nodeCells[i][0] stores the number of entries for each node. nodeCells[i] = {N, CIdx, CNIdx, CIdx2, CNIdx, ... }
    
    cl_float8 wallProperties[MAX_NODES][MAX_CELLS_PER_NODE];
    
    cl_float4 nodesForce[MAX_NODES];
    cl_float4 angularForce[MAX_NODES];
    
    bool *cellDivide = new bool[maxCells];
    
    int guiCellSelected = -1;
    
    CMPlant(float _pressure, float _kEnv, float _dt, cl_float8 _newWallProperties){
        
        
        pressure = _pressure;
        kEnv = _kEnv;
        newWallProperties = _newWallProperties;
        dt = _dt;
        
    }
    
    void CreateOneCell(int n, float r);
    
    void ReadCellsFromFile(const char *filename, int dispW, int dispH, float scale);
    
    int InsertNewNode (int nodeIdx1, int nodeIdx2, cl_float4 newNodePos);
    cl_int2 InsertNewWall (int cellIdx, int nodePos1, int nodePos2, float lFromNode1, float lFromNode2, cl_float8 newWallProperties);
    
    float CellArea(int cellIdx);
    
    void Divide(int cellIdx);
    void Divide2(int cellIdx);
    
    void CalculateForces();
    void Step();
    
    void SetCellParams(int cellIdx, float visc, float ymod);
    
    void GLDisplay(int disp_w, int disp_h, bool displayNodes, bool printNumbers, bool displayForces, float scale);
    void GLPrintText (float x, float y, char *text);
    void GLDisplayForces(int disp_w, int disp_h, float scale);
    
    inline void CLCheckErr (cl_int err, const char *name);
    void CLCreateContext();
    void CLCalculateForces(int coreSize);
    void CLStep();
    void CLWriteStructure();
    void CLGetNewPos();
    void CLReadPos();
    void CLReadForce();
    void CLGetNewStruct();
    void CLStop();
    
    int GUICursorInCell(int x, int y, float scale);
    
    void PrintAll();
};


#endif /* defined(__CMPlant__CMPlant__) */


