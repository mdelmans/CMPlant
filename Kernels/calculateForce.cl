
#define PRESS setup[0]
#define MNPCL (int)setup[1]
#define MCLPN (int)setup[2]
#define NODEN (int)setup[3]
#define CELLN (int)setup[4]
#define K_ENV setup[5]

#define P_AREA  s0
#define P_THIC  s1
#define P_VISC  s2
#define P_YMOD  s3

#define P_SLEN  s4

__kernel void calculate_forces (__global       float4   *nodesForce,
                                __global       float4   *nodesPos,
                                __global const int      *cellNodes,
                                __global const int2     *nodeCells,
                                __global const float8   *wallProperties,
                                __global const float    *setup){
    
    int cNode = get_global_id(0);
    
    float4 zero = (float4)(0,0,0,0);
    
    nodesForce[cNode] = zero;
    
    float4 uZ = (float4)(0,0,1,0);
    
    
    for (int j = 1; j <= nodeCells[cNode * MCLPN].x; j++){
        
        int cCell = nodeCells[cNode * MCLPN + j].x;
        
        int nextNodeIdx;
        int prevNodeIdx;
        
        for (int i = 1; i <= cellNodes[cCell * MNPCL]; i++){
            
            if (cellNodes[cCell * MNPCL + i] == cNode){
                if ( i == cellNodes[cCell * MNPCL] ) nextNodeIdx = cellNodes[cCell * MNPCL +1 ]; // return the first index
                else nextNodeIdx = cellNodes[cCell * MNPCL + i+1];
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
        
        
        
        float4 nextEdge = nodesPos[ nextNodeIdx ] - nodesPos[cNode];
        
        float4 prevEdge = nodesPos[cNode] - nodesPos[ prevNodeIdx ];
        
        float4 NormNext = (float4)(nextEdge.y, -nextEdge.x,0,0);
        float4 NormPrev = (float4)(prevEdge.y, -prevEdge.x,0,0);
        
        float4 uTangNext = normalize(nextEdge);
        float4 uTangPrev = normalize(prevEdge);
        
        float4 uNormNext = normalize(NormNext);
        float4 uNormPrev = normalize(NormPrev);
        
        float lNextEdge = length(nextEdge);
        float lPrevEdge = length(prevEdge);

        
        nodesForce[cNode] = nodesForce[cNode] + uNormNext * lNextEdge * PRESS;
        
        nodesForce[cNode] = nodesForce[cNode] + uNormPrev * lPrevEdge * PRESS;
        
        
        float8 wP = wallProperties[cNode * MCLPN + j];
        
        nodesForce[cNode] += uTangNext * wP.P_AREA * wP.P_YMOD * ( (lNextEdge - wP.P_SLEN) / wP.P_SLEN );
        
       // bool addForce = true;
        /*
        for (int k = 1; k <= nodeCells[cNode * MCLPN].x; k++) {
            int tCell = nodeCells[cNode * MCLPN +k].x;
            
            int nextNodeInTCell;
            
            for (int i = 1; i <= cellNodes[cCell * MNPCL]; i++){
                
                if (cellNodes[tCell * MNPCL + i] == cNode){
                    if ( i == cellNodes[tCell * MNPCL] ) nextNodeInTCell = cellNodes[tCell * MNPCL +1 ]; // return the first index
                    else nextNodeInTCell = cellNodes[tCell * MNPCL + i+1];
                }
            }
            
            
            
            if (k!=j && prevNodeIdx == nextNodeInTCell ) {
                addForce = false;
                break;
            }
        }
        */
        //if  (addForce){
            
            
            int cellSeqInNode;
            
            for (int i = 1; i <= nodeCells[prevNodeIdx * MCLPN].x; i++)
                if (nodeCells[prevNodeIdx * MCLPN + i].x == cCell)  cellSeqInNode = i;
            
            
            
            wP = wallProperties[prevNodeIdx * MCLPN + cellSeqInNode];
            
            nodesForce[cNode] +=  uTangPrev * -wP.P_AREA * wP.P_YMOD * ( (lPrevEdge - wP.P_SLEN) / wP.P_SLEN );
            
        //}
    }
    
}