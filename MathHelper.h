#ifndef MATHHELPER_H
#define MATHHELPER_H

#include "antmath.h"
#include<stdlib.h>
class MathHelper
{
public:
	// Returns random float in [0, 1).
	static float RandF()
	{
		return (float)(rand()) / (float)RAND_MAX;
	}

	// Returns random float in [a, b).
	static float RandF(float a, float b)
	{
		return a + RandF()*(b-a);
	}

	template<typename T>
	static T Min(const T& a, const T& b)
	{
		return a < b ? a : b;
	}

	template<typename T>
	static T Max(const T& a, const T& b)
	{
		return a > b ? a : b;
	}
	 
	template<typename T>
	static T Lerp(const T& a, const T& b, float t)
	{
		return a + (b-a)*t;
	}

	template<typename T>
	static T Clamp(const T& x, const T& low, const T& high)
	{
		return x < low ? low : (x > high ? high : x); 
	}

	// Returns the polar angle of the point (x,y) in [0, 2*PI).
	static float AngleFromXY(float x, float y);

	static AV4X4FLOAT InverseTranspose(const AV4X4FLOAT M)
	{
		// Inverse-transpose is just applied to normals.  So zero out 
		// translation row so that it doesn't get into our inverse-transpose
		// calculation--we don't want the inverse-transpose of the translation.
		AV4FLOAT X(0.0f, 0.0f, 0.0f, 1.0f);
		AV4X4FLOAT A = M;
		A.setRow(3,X);
		
		AV4X4FLOAT Ainv = A.inverse();
		return Ainv.transpose();
	}

	//static AV3FLOAT RandUnitVec3();
	//static AV3FLOAT RandHemisphereUnitVec3(AV3FLOAT n);

	static const float Infinity;
	static const float Pi;


};

#endif // MATHHELPER_H
