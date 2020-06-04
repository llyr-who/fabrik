#pragma once
#include <math.h>
#include <cmath>
#include <iostream>
#include <vector>

#include "solvant/core/vector.hpp"

// TODO: use boost for constants
#define AL_SQRT2 0.70710678118

using UINT = unsigned int;

class cloth {
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
    solvant::vector<float, 3>* prevPos;
    solvant::vector<float, 3>* currPos;
    solvant::vector<float, 3>* velocity;
    solvant::vector<float, 3>* normals;
    solvant::vector<float, 3>* tangents;
    solvant::vector<float, 3>* bitangents;
    solvant::vector<float, 3>* force;

public:
    cloth();
    UINT RowCount() const;
    UINT ColumnCount() const;
    UINT TriangleCount() const;
    UINT VertexCount() const;
    float Width() const;
    float Depth() const;

    // this returns the solution at the ith grid point
    const solvant::vector<float, 3>& operator[](int i) const {
        return currPos[i];
    }

    // Returns the solution normal at the ith grid point.
    const solvant::vector<float, 3>& Normal(int i) const { return normals[i]; }
    // Returns the unit tangent vector at the ith grid point
    const solvant::vector<float, 3>& Tangent(int i) const {
        return tangents[i];
    }
    // Returns the unit bitangent vector at the ith grid point
    const solvant::vector<float, 3>& Bitangent(int i) const {
        return bitangents[i];
    }

    void Init(UINT m, UINT n, float ddx, float ddt, float spring1,
              float spring2, float damp1, float damp2, float M);
    void Update(float dt, float windX, float windY, float windZ);
};
