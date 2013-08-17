//
//  CLGeo.h
//  CMPlant
//
//  Created by Mihails Delmans on 29/07/2013.
//  Copyright (c) 2013 Mihails Delmans. All rights reserved.
//

#ifndef __CMPlant__CLGeo__
#define __CMPlant__CLGeo__

#include <iostream>
#include <OpenCL/OpenCL.h>

cl_float4 & operator+ (cl_float4 &V1, cl_float4 &V2);
cl_float4 & operator- (cl_float4 &V1, cl_float4 &V2);

cl_float4 & operator/ (cl_float4 &V1, float k);
cl_float4 & operator* (cl_float4 &V1, float k);

float & operator* (cl_float4 &V1, cl_float4 &V2);

cl_float4 Cross (cl_float4 V1, cl_float4 V2);

std::ostream& operator<< (std::ostream &out, cl_float4 &V);
std::ostream& operator<< (std::ostream &out, cl_float8 &V);

float SecondNorm (cl_float4 V);
cl_float4 UnitTang (cl_float4 V);

#endif /* defined(__CMPlant__CLGeo__) */
