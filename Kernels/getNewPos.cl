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

__kernel void get_new_pos (__global       float4   *nodesForce,
                           __global       float4   *buffNodesPos,
                           __global       float    *setup){
    
    int cNode = get_global_id(0);
        
    buffNodesPos[cNode] = buffNodesPos[cNode]+ nodesForce[cNode] * (DT / K_ENV);
   
}