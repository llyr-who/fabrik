#include "cloth.h"
#include <iostream>

Cloth::Cloth()
    : numRows(0),
      numCols(0),
      vertexCount(0),
      triangleCount(0),
      dt(0.0f),
      dx(0.0f),
      mass(0.0f),
      gravity(0),
      shortDamp(0.0f),
      shortSpring(0.0f),
      longDamp(0.0f),
      longSpring(0),
      prevPos(0),
      currPos(0),
      velocity(0),
      normals(0),
      tangents(0),
      bitangents(0),
      force(0) {}

UINT Cloth::RowCount() const { return numRows; }

UINT Cloth::ColumnCount() const { return numCols; }

UINT Cloth::VertexCount() const { return vertexCount; }

UINT Cloth::TriangleCount() const { return triangleCount; }

float Cloth::Width() const { return numCols * dx; }

float Cloth::Depth() const { return numRows * dx; }

void Cloth::Init(UINT m, UINT n, float ddx, float ddt, float spring1,
                 float spring2, float damp1, float damp2, float M) {
    mass = M;
    numRows = m;
    numCols = n;

    vertexCount = m * n;
    triangleCount = (m - 1) * (n - 1) * 2;

    dt = ddt;
    dx = ddx;

    shortDamp = damp1;
    longDamp = damp2;
    shortSpring = spring1;
    longSpring = spring2;

    // In case Init() called again.
    delete[] prevPos;
    delete[] currPos;
    delete[] velocity;
    delete[] normals;

    prevPos = new AV3FLOAT[m * n];
    currPos = new AV3FLOAT[m * n];
    velocity = new AV3FLOAT[m * n];
    force = new AV3FLOAT[m * n];
    normals = new AV3FLOAT[m * n];
    tangents = new AV3FLOAT[m * n];
    bitangents = new AV3FLOAT[m * n];

    // Generate grid vertices in system memory.

    float halfWidth = (n - 1) * dx * 0.5f;
    float halfDepth = (m - 1) * dx * 0.5f;
    for (UINT i = 0; i < m; ++i) {
        float z = halfDepth - i * dx;
        for (UINT j = 0; j < n; ++j) {
            float x = -halfWidth + j * dx;

            prevPos[i * n + j] = AV3FLOAT(x, 0.1 * sin(x * z), z);
            currPos[i * n + j] = AV3FLOAT(x, 0.1 * sin(x * z), z);
            velocity[i * n + j] = AV3FLOAT(0, 0, 0);
            normals[i * n + j] = AV3FLOAT(0, 1, 0);
        }
    }
}

void Cloth::Update(float ddt, float windX, float windY, float windZ) {
    static float t = 0;
    t += ddt;
    UINT n = numCols;
    UINT m = numRows;
    if (t >= dt) {
        for (UINT i = 0; i < n; ++i) {
            for (UINT j = 0; j < m; ++j) {
                force[j * n + i].x = 0;
                force[j * n + i].y = 0;
                force[j * n + i].z = 0;
            }
        }

        for (UINT i = 0; i < n; ++i) {
            for (UINT j = 0; j < m; ++j) {
                float WFx, WFy, WFz;
                WFx = normals[j * n + i].x * (windX - velocity[j * n + i].x);
                WFy = normals[j * n + i].y * (windY - velocity[j * n + i].y);
                WFz = normals[j * n + i].z * (windZ - velocity[j * n + i].z);
                float Wind = ((WFx + WFy + WFz) > 0 ? WFx + WFy + WFz
                                                    : -(WFx + WFy + WFz));
                force[j * n + i].x += 0.3 * Wind * windX;
                force[j * n + i].y += mass * gravity + 0.3 * Wind * windY;
                force[j * n + i].z += 0.3 * Wind * windZ;
            }
        }
        // we do the double links first

        for (UINT i = 0; i < n - 2; ++i) {
            for (UINT j = 0; j < m - 2; ++j) {
                float F1x, F1y, F1z;
                float F2x, F2y, F2z;
                float diffx;
                float diffy;
                float diffz;
                float velDiffx;
                float velDiffy;
                float velDiffz;
                float diffNorm;
                float k;

                diffx = currPos[(j + 2) * n + i].x - currPos[j * n + i].x;
                diffy = currPos[(j + 2) * n + i].y - currPos[j * n + i].y;
                diffz = currPos[(j + 2) * n + i].z - currPos[j * n + i].z;
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[(j + 2) * n + i].x - velocity[j * n + i].x;
                velDiffy = velocity[(j + 2) * n + i].y - velocity[j * n + i].y;
                velDiffz = velocity[(j + 2) * n + i].z - velocity[j * n + i].z;

                F1x = k * diffx + longDamp * velDiffx;
                F1y = k * diffy + longDamp * velDiffy;
                F1z = k * diffz + longDamp * velDiffz;

                diffx = currPos[j * n + i + 2].x - currPos[j * n + i].x;
                diffy = currPos[j * n + i + 2].y - currPos[j * n + i].y;
                diffz = currPos[j * n + i + 2].z - currPos[j * n + i].z;
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[j * n + i + 2].x - velocity[j * n + i].x;
                velDiffy = velocity[j * n + i + 2].y - velocity[j * n + i].y;
                velDiffz = velocity[j * n + i + 2].z - velocity[j * n + i].z;

                F2x = k * diffx + longDamp * velDiffx;
                F2y = k * diffy + longDamp * velDiffy;
                F2z = k * diffz + longDamp * velDiffz;

                force[(j + 2) * n + i].x -= F1x;
                force[(j + 2) * n + i].y -= F1y;
                force[(j + 2) * n + i].z -= F1z;

                force[j * n + i].x += F1x;
                force[j * n + i].y += F1y;
                force[j * n + i].z += F1z;

                force[j * n + i + 2].x -= F2x;
                force[j * n + i + 2].y -= F2y;
                force[j * n + i + 2].z -= F2z;

                force[j * n + i].x += F2x;
                force[j * n + i].y += F2y;
                force[j * n + i].z += F2z;
            }
        }

        for (UINT i = 0; i < n - 2; ++i) {
            float F2x, F2y, F2z;
            float diffx;
            float diffy;
            float diffz;
            float velDiffx;
            float velDiffy;
            float velDiffz;
            float diffNorm;
            float k;

            diffx = currPos[(m - 1) * n + i + 2].x - currPos[(m - 1) * n + i].x;
            diffy = currPos[(m - 1) * n + i + 2].y - currPos[(m - 1) * n + i].y;
            diffz = currPos[(m - 1) * n + i + 2].z - currPos[(m - 1) * n + i].z;
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(m - 1) * n + i + 2].x - velocity[(m - 1) * n + i].x;
            velDiffy =
                velocity[(m - 1) * n + i + 2].y - velocity[(m - 1) * n + i].y;
            velDiffz =
                velocity[(m - 1) * n + i + 2].z - velocity[(m - 1) * n + i].z;

            F2x = k * diffx + longDamp * velDiffx;
            F2y = k * diffy + longDamp * velDiffy;
            F2z = k * diffz + longDamp * velDiffz;

            force[(m - 1) * n + i + 2].x -= F2x;
            force[(m - 1) * n + i + 2].y -= F2y;
            force[(m - 1) * n + i + 2].z -= F2z;

            force[(m - 1) * n + i].x += F2x;
            force[(m - 1) * n + i].y += F2y;
            force[(m - 1) * n + i].z += F2z;
        }
        for (UINT j = 0; j < m - 2; ++j) {
            float F1x, F1y, F1z;

            float diffx;
            float diffy;
            float diffz;
            float velDiffx;
            float velDiffy;
            float velDiffz;
            float diffNorm;
            float k;

            diffx = currPos[(j + 2) * n + n - 1].x - currPos[j * n + n - 1].x;
            diffy = currPos[(j + 2) * n + n - 1].y - currPos[j * n + n - 1].y;
            diffz = currPos[(j + 2) * n + n - 1].z - currPos[j * n + n - 1].z;
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(j + 2) * n + n - 1].x - velocity[j * n + n - 1].x;
            velDiffy =
                velocity[(j + 2) * n + n - 1].y - velocity[j * n + n - 1].y;
            velDiffz =
                velocity[(j + 2) * n + n - 1].z - velocity[j * n + n - 1].z;

            F1x = k * diffx + longDamp * velDiffx;
            F1y = k * diffy + longDamp * velDiffy;
            F1z = k * diffz + longDamp * velDiffz;

            force[(j + 2) * n + n - 1].x -= F1x;
            force[(j + 2) * n + n - 1].y -= F1y;
            force[(j + 2) * n + n - 1].z -= F1z;

            force[j * n + n - 1].x += F1x;
            force[j * n + n - 1].y += F1y;
            force[j * n + n - 1].z += F1z;
        }

        for (UINT i = 0; i < n - 1; ++i) {
            for (UINT j = 0; j < m - 1; ++j) {
                float F1x, F1y, F1z;
                float F2x, F2y, F2z;
                float F3x, F3y, F3z;
                float F4x, F4y, F4z;
                float diffx;
                float diffy;
                float diffz;
                float velDiffx;
                float velDiffy;
                float velDiffz;
                float diffNorm;
                float k;
                // the "from left to right" diagonal connection in our
                // computational molecule force due to spring first
                // NOTE THE SQUARE ROOT IN DIAG CONNECTIONS

                diffx = currPos[(j + 1) * n + (i + 1)].x - currPos[j * n + i].x;
                diffy = currPos[(j + 1) * n + (i + 1)].y - currPos[j * n + i].y;
                diffz = currPos[(j + 1) * n + (i + 1)].z - currPos[j * n + i].z;
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - AL_SQRT2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx =
                    velocity[(j + 1) * n + (i + 1)].x - velocity[j * n + i].x;
                velDiffy =
                    velocity[(j + 1) * n + (i + 1)].y - velocity[j * n + i].y;
                velDiffz =
                    velocity[(j + 1) * n + (i + 1)].z - velocity[j * n + i].z;

                F1x = k * diffx + shortDamp * velDiffx;
                F1y = k * diffy + shortDamp * velDiffy;
                F1z = k * diffz + shortDamp * velDiffz;

                // the "from right to left" diagonal connection in our
                // computational molecule force due to spring first
                // NOTE THE SQUARE ROOT IN DIAG CONNECTIONS

                diffx = currPos[(j + 1) * n + i].x - currPos[j * n + (i + 1)].x;
                diffy = currPos[(j + 1) * n + i].y - currPos[j * n + (i + 1)].y;
                diffz = currPos[(j + 1) * n + i].z - currPos[j * n + (i + 1)].z;
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - AL_SQRT2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx =
                    velocity[(j + 1) * n + i].x - velocity[j * n + (i + 1)].x;
                velDiffy =
                    velocity[(j + 1) * n + i].y - velocity[j * n + (i + 1)].y;
                velDiffz =
                    velocity[(j + 1) * n + i].z - velocity[j * n + (i + 1)].z;

                F2x = k * diffx + shortDamp * velDiffx;
                F2y = k * diffy + shortDamp * velDiffy;
                F2z = k * diffz + shortDamp * velDiffz;

                // non-diag connection "up to down" horizontal component

                diffx = currPos[(j + 1) * n + i].x - currPos[j * n + i].x;
                diffy = currPos[(j + 1) * n + i].y - currPos[j * n + i].y;
                diffz = currPos[(j + 1) * n + i].z - currPos[j * n + i].z;
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[(j + 1) * n + i].x - velocity[j * n + i].x;
                velDiffy = velocity[(j + 1) * n + i].y - velocity[j * n + i].y;
                velDiffz = velocity[(j + 1) * n + i].z - velocity[j * n + i].z;

                F4x = k * diffx + shortDamp * velDiffx;
                F4y = k * diffy + shortDamp * velDiffy;
                F4z = k * diffz + shortDamp * velDiffz;

                // this is the "from dleft to right" horizontal component

                diffx = currPos[j * n + (i + 1)].x - currPos[j * n + i].x;
                diffy = currPos[j * n + (i + 1)].y - currPos[j * n + i].y;
                diffz = currPos[j * n + (i + 1)].z - currPos[j * n + i].z;
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[j * n + (i + 1)].x - velocity[j * n + i].x;
                velDiffy = velocity[j * n + (i + 1)].y - velocity[j * n + i].y;
                velDiffz = velocity[j * n + (i + 1)].z - velocity[j * n + i].z;

                F3x = k * diffx + shortDamp * velDiffx;
                F3y = k * diffy + shortDamp * velDiffy;
                F3z = k * diffz + shortDamp * velDiffz;
                /*
                 *here is where we add the force to the force array
                 */

                force[j * n + i].x += F1x;
                force[j * n + i].y += F1y;
                force[j * n + i].z += F1z;

                force[j * n + i].x += F3x;
                force[j * n + i].y += F3y;
                force[j * n + i].z += F3z;

                force[j * n + i].x += F4x;
                force[j * n + i].y += F4y;
                force[j * n + i].z += F4z;

                force[(j + 1) * n + (i + 1)].x -= F1x;
                force[(j + 1) * n + (i + 1)].y -= F1y;
                force[(j + 1) * n + (i + 1)].z -= F1z;

                force[j * n + (i + 1)].x -= F3x;
                force[j * n + (i + 1)].y -= F3y;
                force[j * n + (i + 1)].z -= F3z;

                force[j * n + (i + 1)].x += F2x;
                force[j * n + (i + 1)].y += F2y;
                force[j * n + (i + 1)].z += F2z;

                force[(j + 1) * n + i].x -= F4x;
                force[(j + 1) * n + i].y -= F4y;
                force[(j + 1) * n + i].z -= F4z;

                force[(j + 1) * n + i].x -= F2x;
                force[(j + 1) * n + i].y -= F2y;
                force[(j + 1) * n + i].z -= F2z;
            }
        }
        for (UINT i = 0; i < n - 1; ++i) {
            float F4x, F4y, F4z;
            float diffx;
            float diffy;
            float diffz;
            float velDiffx;
            float velDiffy;
            float velDiffz;
            float diffNorm;
            float k;

            // non-diag connection "up to down" horizontal component

            diffx = currPos[(m - 1) * n + i + 1].x - currPos[(m - 1) * n + i].x;
            diffy = currPos[(m - 1) * n + i + 1].y - currPos[(m - 1) * n + i].y;
            diffz = currPos[(m - 1) * n + i + 1].z - currPos[(m - 1) * n + i].z;
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(m - 1) * n + i + 1].x - velocity[(m - 1) * n + i].x;
            velDiffy =
                velocity[(m - 1) * n + i + 1].y - velocity[(m - 1) * n + i].y;
            velDiffz =
                velocity[(m - 1) * n + i + 1].z - velocity[(m - 1) * n + i].z;

            F4x = k * diffx + shortDamp * velDiffx;
            F4y = k * diffy + shortDamp * velDiffy;
            F4z = k * diffz + shortDamp * velDiffz;

            force[(m - 1) * n + i + 1].x -= F4x;
            force[(m - 1) * n + i + 1].y -= F4y;
            force[(m - 1) * n + i + 1].z -= F4z;

            force[(m - 1) * n + i].x += F4x;
            force[(m - 1) * n + i].y += F4y;
            force[(m - 1) * n + i].z += F4z;
        }
        for (UINT j = 0; j < m - 1; ++j) {
            float F4x, F4y, F4z;
            float diffx;
            float diffy;
            float diffz;
            float velDiffx;
            float velDiffy;
            float velDiffz;
            float diffNorm;
            float k;

            // non-diag connection "up to down" horizontal component

            diffx = currPos[(j + 1) * n + n - 1].x - currPos[j * n + n - 1].x;
            diffy = currPos[(j + 1) * n + n - 1].y - currPos[j * n + n - 1].y;
            diffz = currPos[(j + 1) * n + n - 1].z - currPos[j * n + n - 1].z;
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(j + 1) * n + n - 1].x - velocity[j * n + n - 1].x;
            velDiffy =
                velocity[(j + 1) * n + n - 1].y - velocity[j * n + n - 1].y;
            velDiffz =
                velocity[(j + 1) * n + n - 1].z - velocity[j * n + n - 1].z;

            F4x = k * diffx + shortDamp * velDiffx;
            F4y = k * diffy + shortDamp * velDiffy;
            F4z = k * diffz + shortDamp * velDiffz;

            force[j * n + n - 1].x += F4x;
            force[j * n + n - 1].y += F4y;
            force[j * n + n - 1].z += F4z;

            force[(j + 1) * n + n - 1].x -= F4x;
            force[(j + 1) * n + n - 1].y -= F4y;
            force[(j + 1) * n + n - 1].z -= F4z;
        }

        // update the position's of the elements and velocities
        for (UINT i = 1; i < n; ++i) {
            for (UINT j = 0; j < m; ++j) {
                prevPos[j * numRows + i].x =
                    currPos[j * n + i].x + velocity[j * n + i].x * dt +
                    force[j * n + i].x * 0.5 * (1 / mass) * dt * dt;

                prevPos[j * numRows + i].y =
                    currPos[j * n + i].y + velocity[j * n + i].y * dt +
                    force[j * n + i].y * 0.5 * (1 / mass) * dt * dt;

                prevPos[j * numRows + i].z =
                    currPos[j * n + i].z + velocity[j * n + i].z * dt +
                    force[j * n + i].z * 0.5 * (1 / mass) * dt * dt;

                velocity[j * n + i].x =
                    (prevPos[j * n + i].x - currPos[j * n + i].x) * (1 / dt);
                velocity[j * n + i].y =
                    (prevPos[j * n + i].y - currPos[j * n + i].y) * (1 / dt);
                velocity[j * n + i].z =
                    (prevPos[j * n + i].z - currPos[j * n + i].z) * (1 / dt);
            }
        }

        for (UINT i = 0; i < n - 1; ++i) {
            for (UINT j = 0; j < m - 1; ++j) {
                tangents[j * n + i].x =
                    currPos[(j + 1) * n + i].x - currPos[j * n + i].x;
                tangents[j * n + i].y =
                    currPos[(j + 1) * n + i].y - currPos[j * n + i].y;
                tangents[j * n + i].z =
                    currPos[(j + 1) * n + i].z - currPos[j * n + i].z;
                tangents[j * n + i].normalize();

                bitangents[j * n + i].x =
                    currPos[j * n + i + 1].x - currPos[j * n + i].x;
                bitangents[j * n + i].y =
                    currPos[j * n + i + 1].y - currPos[j * n + i].y;
                bitangents[j * n + i].z =
                    currPos[j * n + i + 1].z - currPos[j * n + i].z;
                bitangents[j * n + i].normalize();

                crossProduct3(bitangents[j * n + i], tangents[j * n + i],
                              normals[j * n + i]);
            }
        }

        for (UINT j = 0; j < m; ++j) {
            normals[j * n + n - 1] = normals[j * n + n - 2];
        }
        for (UINT i = 0; i < n; ++i) {
            normals[(m - 1) * n + i] = normals[(m - 2) * n + i];
        }

        std::swap(prevPos, currPos);

        t = 0.0f;
    }
}
