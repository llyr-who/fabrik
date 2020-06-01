#pragma once
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <X11/X.h>
#include <X11/Xlib.h>
#include <math.h>
#include <vector>
#include "antmath.h"

class Cloth {
private:
    UINT numRows;
    UINT numCols;

    UINT vertexCount;
    UINT triangleCount;

    float dt;  // time step
    float dx;  // spatial step
    float mass;
    // all we want is the magnitude in the y-direction.
    // as gravity always points in the neg y direction.
    float gravity;

    // some ability to modify the spring constant parameters
    // should be added at some point.
    float shortDamp;
    float shortSpring;
    float longDamp;
    float longSpring;
    AV3FLOAT* prevPos;
    AV3FLOAT* currPos;
    AV3FLOAT* velocity;
    AV3FLOAT* normals;
    AV3FLOAT* tangents;
    AV3FLOAT* bitangents;
    AV3FLOAT* force;

public:
    Cloth();
    UINT RowCount() const;
    UINT ColumnCount() const;
    UINT TriangleCount() const;
    UINT VertexCount() const;
    float Width() const;
    float Depth() const;

    // this returns the solution at the ith grid point
    const AV3FLOAT& operator[](int i) const { return currPos[i]; }

    // Returns the solution normal at the ith grid point.
    const AV3FLOAT& Normal(int i) const { return normals[i]; }
    // Returns the unit tangent vector at the ith grid point
    const AV3FLOAT& Tangent(int i) const { return tangents[i]; }
    // Returns the unit bitangent vector at the ith grid point
    const AV3FLOAT& Bitangent(int i) const { return bitangents[i]; }

    void Init(UINT m, UINT n, float ddx, float ddt, float spring1,
              float spring2, float damp1, float damp2, float M);
    void Update(float dt, float windX, float windY, float windZ);
};
