#define PRESS setup[0]
#define MNPCL (int)setup[1]
#define MCLPN (int)setup[2]
#define NODEN (int)setup[3]
#define CELLN (int)setup[4]
#define K_ENV setup[5]
#define DT setup[6]

#define P_AREA  s0
#define P_THIC  s1
#define P_VISC  s2
#define P_YMOD  s3

#define P_SLEN  s4

__kernel void get_new_struct (__global       float4   *nodesForce,
                              __global       float4   *nodesPos,
                              __global const int      *cellNodes,
                              __global const int2     *nodeCells,
                              __global const float8   *wallProperties,
                              __global const float    *setup){

    
    int cNode = get_global_id(0);
    
    for (int i = 1; i <= nodeCells[cNode * MCLPN].x; i++){
        int cCell = nodeCells[cNode * MCLPN + i].x;
        
        //float lNext = SecondNorm( nodesPos[cNode] - nodesPos[ NextInCell(cCell, cNode).x ] );
        
        
        int nextIncCell;
        
        for (int i = 1; i <= cellNodes[cCell * MNPCL]; i++){
            
            if (cellNodes[cCell * MNPCL + i] == cNode){
                if ( i == cellNodes[cCell * MNPCL] ) nextIncCell = cellNodes[cCell * MNPCL +1 ]; // return the first index
                else nextIncCell = cellNodes[cCell * MNPCL + i+1];
                break;
            }
        }
        
        float lNext = fast_distance(nodesPos[cNode], nodesPos[nextIncCell]  );
        
        float8 wP = wallProperties[cNode * MCLPN + i];
        
        wallProperties[cNode *MCLPN + i].P_SLEN = wP.P_SLEN + (wP.P_YMOD / wP.P_VISC ) * (lNext - wP.P_SLEN) * DT;
        
        bool changeS = true;
        
        int j = 0;
        
        do{
            j++;
            int jCell = nodeCells[cNode * MCLPN + j].x;
            
            int prevNodeIdx;
            int nextNodeIdx;
            
            for (int i = 1; i <= cellNodes[jCell * MNPCL]; i++){
                
                if (cellNodes[jCell * MNPCL + i] == cNode){
                    if ( i == cellNodes[jCell * MNPCL] ) nextNodeIdx = cellNodes[jCell * MNPCL +1 ]; // return the first index
                    else nextNodeIdx = cellNodes[jCell * MNPCL + i+1];
                    break;
                }
            }
            
            
            for (int i = 1; i <= cellNodes[cCell * MNPCL]; i++){
                
                if (cellNodes[cCell * MNPCL + i] == cNode){
                    if ( i == 1 ) prevNodeIdx = cellNodes[cCell * MNPCL + cellNodes[cCell * MNPCL] ];
                    else prevNodeIdx = cellNodes[cCell * MNPCL + i-1];
                    break;
                }
            }

            
            
            if ( i!=j && prevNodeIdx == nextNodeIdx ) 
                changeS = false;
            
            
        }while(changeS == true && j <= nodeCells[cNode * MCLPN].x);
        
        if (changeS){
            //int pNode = PrevInCell(cCell, cNode).x;
            
            int pNode;
            
            for (int i = 1; i <= cellNodes[cCell * MNPCL]; i++){
                
                if (cellNodes[cCell * MNPCL + i] == cNode){
                    if ( i == 1 ) pNode = cellNodes[cCell * MNPCL + cellNodes[cCell * MNPCL] ];
                    else pNode = cellNodes[cCell * MNPCL + i-1];
                    break;
                }
            }

            
            //float lPrev = SecondNorm( nodesPos[cNode] - nodesPos[pNode] );
            
            float lPrev = fast_distance ( nodesPos[cNode], nodesPos[pNode] );
            
            //int cCellPosInpNode = CellSeqInNode(cCell, pNode);
            
            
            int cCellPosInpNode;
            
            for (int i = 1; i <= nodeCells[pNode * MCLPN].x; i++)
                if (nodeCells[pNode * MCLPN + i].x == cCell)  cCellPosInpNode = i;
            
            
            float8 wP = wallProperties[ pNode * MCLPN + cCellPosInpNode ];
            wallProperties[ pNode * MCLPN + cCellPosInpNode ].P_SLEN = wP.P_SLEN + (wP.P_YMOD / wP.P_VISC) * (lPrev - wP.P_SLEN) * DT;
            
        }
        
    
    }
    
}