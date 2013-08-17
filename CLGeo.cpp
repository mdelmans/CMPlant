//
//  CLGeo.cpp
//  CMPlant
//
//  Created by Mihails Delmans on 29/07/2013.
//  Copyright (c) 2013 Mihails Delmans. All rights reserved.
//

#include "CLGeo.h"
#include <math.h>

cl_float4 out;
float fout;

cl_float4 & operator+ (cl_float4 &V1, cl_float4 &V2){
    extern cl_float4 out;
    out = {V1.x+V2.x, V1.y+V2.y, V1.z + V2.z, V1.w + V2.w};
    return out;
}

cl_float4 & operator- (cl_float4 &V1, cl_float4 &V2){
    extern cl_float4 out;
    out = {V1.x - V2.x, V1.y - V2.y, V1.z - V2.z, V1.w - V2.w};
    return out;
}

cl_float4 & operator/ (cl_float4 &V1, float k){
    extern cl_float4 out;
    out = {V1.x / k, V1.y / k, V1.z / k, V1.w / k};
    return out;
}

cl_float4 & operator* (cl_float4 &V1, float k){
    extern cl_float4 out;
    out = {V1.x * k, V1.y * k, V1.z * k, V1.w * k};
    return out;
}

std::ostream & operator<< (std::ostream &out, cl_float4 &V){
    out << "(" << V.x << "," << V.y << "," <<V.z << "," << V.w << ")";
    return out;
}

std::ostream & operator<< (std::ostream &out, cl_float8 &V){
    out << "(" << V.s[0] << ","<< V.s[1] << ","<< V.s[2] << ","<< V.s[3] << ","<< V.s[4] << ","<< V.s[5] << ","<< V.s[6] << "," << V.s[7] << ")";
    return out;
}

float SecondNorm (cl_float4 V1){
    return sqrt( pow(V1.x,2) + pow(V1.y,2) );
}

cl_float4 UnitTang (cl_float4 V){
    
    return V / SecondNorm(V);
}


cl_float4 Cross(cl_float4 V1, cl_float4 V2){
    return {V1.y * V2.z - V1.z * V2.y, V1.z * V2.x - V1.x * V2.z, V1.x * V2.y - V1.y * V2.x};
}

float & operator* (cl_float4 &V1, cl_float4 &V2){
    extern float fout;
    
    fout = V1.x * V2.x + V1.y * V2.y;
    
    return fout;
}
