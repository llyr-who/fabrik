#include "cloth.h"
#include <iostream>

cloth::cloth()
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

UINT cloth::RowCount() const { return numRows; }

UINT cloth::ColumnCount() const { return numCols; }

UINT cloth::VertexCount() const { return vertexCount; }

UINT cloth::TriangleCount() const { return triangleCount; }

float cloth::Width() const { return numCols * dx; }

float cloth::Depth() const { return numRows * dx; }

void cloth::Init(UINT m, UINT n, float ddx, float ddt, float spring1,
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

    prevPos = new solvant::vector<float,3>[m * n];
    currPos = new solvant::vector<float,3>[m * n];
    velocity = new solvant::vector<float,3>[m * n];
    force = new solvant::vector<float,3>[m * n];
    normals = new solvant::vector<float,3>[m * n];
    tangents = new solvant::vector<float,3>[m * n];
    bitangents = new solvant::vector<float,3>[m * n];

    // Generate grid vertices in system memory.

    float halfWidth = (n - 1) * dx * 0.5f;
    float halfDepth = (m - 1) * dx * 0.5f;
    for (UINT i = 0; i < m; ++i) {
        float z = halfDepth - i * dx;
        for (UINT j = 0; j < n; ++j) {
            float x = -halfWidth + j * dx;

            prevPos[i * n + j] = solvant::vector<float,3>(x, 0.1 * sin(x * z), z);
            currPos[i * n + j] = solvant::vector<float,3>(x, 0.1 * sin(x * z), z);
            velocity[i * n + j] = solvant::vector<float,3>(0, 0, 0);
            normals[i * n + j] = solvant::vector<float,3>(0, 1, 0);
        }
    }
}

void cloth::Update(float ddt, float windX, float windY, float windZ) {
    static float t = 0;
    t += ddt;
    UINT n = numCols;
    UINT m = numRows;
    if (t >= dt) {
        for (UINT i = 0; i < n; ++i) {
            for (UINT j = 0; j < m; ++j) {
                force[j * n + i][0] = 0;
                force[j * n + i][1] = 0;
                force[j * n + i][2] = 0;
            }
        }

        for (UINT i = 0; i < n; ++i) {
            for (UINT j = 0; j < m; ++j) {
                float WFx, WFy, WFz;
                WFx = normals[j * n + i][0] * (windX - velocity[j * n + i][0]);
                WFy = normals[j * n + i][1] * (windY - velocity[j * n + i][1]);
                WFz = normals[j * n + i][2] * (windZ - velocity[j * n + i][2]);
                float Wind = ((WFx + WFy + WFz) > 0 ? WFx + WFy + WFz
                                                    : -(WFx + WFy + WFz));
                force[j * n + i][0] += 0.3 * Wind * windX;
                force[j * n + i][1] += mass * gravity + 0.3 * Wind * windY;
                force[j * n + i][2] += 0.3 * Wind * windZ;
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

                diffx = currPos[(j + 2) * n + i][0] - currPos[j * n + i][0];
                diffy = currPos[(j + 2) * n + i][1] - currPos[j * n + i][1];
                diffz = currPos[(j + 2) * n + i][2] - currPos[j * n + i][2];
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[(j + 2) * n + i][0] - velocity[j * n + i][0];
                velDiffy = velocity[(j + 2) * n + i][1] - velocity[j * n + i][1];
                velDiffz = velocity[(j + 2) * n + i][2] - velocity[j * n + i][2];

                F1x = k * diffx + longDamp * velDiffx;
                F1y = k * diffy + longDamp * velDiffy;
                F1z = k * diffz + longDamp * velDiffz;

                diffx = currPos[j * n + i + 2][0] - currPos[j * n + i][0];
                diffy = currPos[j * n + i + 2][1] - currPos[j * n + i][1];
                diffz = currPos[j * n + i + 2][2] - currPos[j * n + i][2];
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[j * n + i + 2][0] - velocity[j * n + i][0];
                velDiffy = velocity[j * n + i + 2][1] - velocity[j * n + i][1];
                velDiffz = velocity[j * n + i + 2][2] - velocity[j * n + i][2];

                F2x = k * diffx + longDamp * velDiffx;
                F2y = k * diffy + longDamp * velDiffy;
                F2z = k * diffz + longDamp * velDiffz;

                force[(j + 2) * n + i][0] -= F1x;
                force[(j + 2) * n + i][1] -= F1y;
                force[(j + 2) * n + i][2] -= F1z;

                force[j * n + i][0] += F1x;
                force[j * n + i][1] += F1y;
                force[j * n + i][2] += F1z;

                force[j * n + i + 2][0] -= F2x;
                force[j * n + i + 2][1] -= F2y;
                force[j * n + i + 2][2] -= F2z;

                force[j * n + i][0] += F2x;
                force[j * n + i][1] += F2y;
                force[j * n + i][2] += F2z;
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

            diffx = currPos[(m - 1) * n + i + 2][0] - currPos[(m - 1) * n + i][0];
            diffy = currPos[(m - 1) * n + i + 2][1] - currPos[(m - 1) * n + i][1];
            diffz = currPos[(m - 1) * n + i + 2][2] - currPos[(m - 1) * n + i][2];
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(m - 1) * n + i + 2][0] - velocity[(m - 1) * n + i][0];
            velDiffy =
                velocity[(m - 1) * n + i + 2][1] - velocity[(m - 1) * n + i][1];
            velDiffz =
                velocity[(m - 1) * n + i + 2][2] - velocity[(m - 1) * n + i][2];

            F2x = k * diffx + longDamp * velDiffx;
            F2y = k * diffy + longDamp * velDiffy;
            F2z = k * diffz + longDamp * velDiffz;

            force[(m - 1) * n + i + 2][0] -= F2x;
            force[(m - 1) * n + i + 2][1] -= F2y;
            force[(m - 1) * n + i + 2][2] -= F2z;

            force[(m - 1) * n + i][0] += F2x;
            force[(m - 1) * n + i][1] += F2y;
            force[(m - 1) * n + i][2] += F2z;
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

            diffx = currPos[(j + 2) * n + n - 1][0] - currPos[j * n + n - 1][0];
            diffy = currPos[(j + 2) * n + n - 1][1] - currPos[j * n + n - 1][1];
            diffz = currPos[(j + 2) * n + n - 1][2] - currPos[j * n + n - 1][2];
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = longSpring * (diffNorm - 2 * dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(j + 2) * n + n - 1][0] - velocity[j * n + n - 1][0];
            velDiffy =
                velocity[(j + 2) * n + n - 1][1] - velocity[j * n + n - 1][1];
            velDiffz =
                velocity[(j + 2) * n + n - 1][2] - velocity[j * n + n - 1][2];

            F1x = k * diffx + longDamp * velDiffx;
            F1y = k * diffy + longDamp * velDiffy;
            F1z = k * diffz + longDamp * velDiffz;

            force[(j + 2) * n + n - 1][0] -= F1x;
            force[(j + 2) * n + n - 1][1] -= F1y;
            force[(j + 2) * n + n - 1][2] -= F1z;

            force[j * n + n - 1][0] += F1x;
            force[j * n + n - 1][1] += F1y;
            force[j * n + n - 1][2] += F1z;
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

                diffx = currPos[(j + 1) * n + (i + 1)][0] - currPos[j * n + i][0];
                diffy = currPos[(j + 1) * n + (i + 1)][1] - currPos[j * n + i][1];
                diffz = currPos[(j + 1) * n + (i + 1)][2] - currPos[j * n + i][2];
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - AL_SQRT2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx =
                    velocity[(j + 1) * n + (i + 1)][0] - velocity[j * n + i][0];
                velDiffy =
                    velocity[(j + 1) * n + (i + 1)][1] - velocity[j * n + i][1];
                velDiffz =
                    velocity[(j + 1) * n + (i + 1)][2] - velocity[j * n + i][2];

                F1x = k * diffx + shortDamp * velDiffx;
                F1y = k * diffy + shortDamp * velDiffy;
                F1z = k * diffz + shortDamp * velDiffz;

                // the "from right to left" diagonal connection in our
                // computational molecule force due to spring first
                // NOTE THE SQUARE ROOT IN DIAG CONNECTIONS

                diffx = currPos[(j + 1) * n + i][0] - currPos[j * n + (i + 1)][0];
                diffy = currPos[(j + 1) * n + i][1] - currPos[j * n + (i + 1)][1];
                diffz = currPos[(j + 1) * n + i][2] - currPos[j * n + (i + 1)][2];
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - AL_SQRT2 * dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx =
                    velocity[(j + 1) * n + i][0] - velocity[j * n + (i + 1)][0];
                velDiffy =
                    velocity[(j + 1) * n + i][1] - velocity[j * n + (i + 1)][1];
                velDiffz =
                    velocity[(j + 1) * n + i][2] - velocity[j * n + (i + 1)][2];

                F2x = k * diffx + shortDamp * velDiffx;
                F2y = k * diffy + shortDamp * velDiffy;
                F2z = k * diffz + shortDamp * velDiffz;

                // non-diag connection "up to down" horizontal component

                diffx = currPos[(j + 1) * n + i][0] - currPos[j * n + i][0];
                diffy = currPos[(j + 1) * n + i][1] - currPos[j * n + i][1];
                diffz = currPos[(j + 1) * n + i][2] - currPos[j * n + i][2];
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[(j + 1) * n + i][0] - velocity[j * n + i][0];
                velDiffy = velocity[(j + 1) * n + i][1] - velocity[j * n + i][1];
                velDiffz = velocity[(j + 1) * n + i][2] - velocity[j * n + i][2];

                F4x = k * diffx + shortDamp * velDiffx;
                F4y = k * diffy + shortDamp * velDiffy;
                F4z = k * diffz + shortDamp * velDiffz;

                // this is the "from dleft to right" horizontal component

                diffx = currPos[j * n + (i + 1)][0] - currPos[j * n + i][0];
                diffy = currPos[j * n + (i + 1)][1] - currPos[j * n + i][1];
                diffz = currPos[j * n + (i + 1)][2] - currPos[j * n + i][2];
                diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
                // now force due to damper
                velDiffx = velocity[j * n + (i + 1)][0] - velocity[j * n + i][0];
                velDiffy = velocity[j * n + (i + 1)][1] - velocity[j * n + i][1];
                velDiffz = velocity[j * n + (i + 1)][2] - velocity[j * n + i][2];

                F3x = k * diffx + shortDamp * velDiffx;
                F3y = k * diffy + shortDamp * velDiffy;
                F3z = k * diffz + shortDamp * velDiffz;
                /*
                 *here is where we add the force to the force array
                 */

                force[j * n + i][0] += F1x;
                force[j * n + i][1] += F1y;
                force[j * n + i][2] += F1z;

                force[j * n + i][0] += F3x;
                force[j * n + i][1] += F3y;
                force[j * n + i][2] += F3z;

                force[j * n + i][0] += F4x;
                force[j * n + i][1] += F4y;
                force[j * n + i][2] += F4z;

                force[(j + 1) * n + (i + 1)][0] -= F1x;
                force[(j + 1) * n + (i + 1)][1] -= F1y;
                force[(j + 1) * n + (i + 1)][2] -= F1z;

                force[j * n + (i + 1)][0] -= F3x;
                force[j * n + (i + 1)][1] -= F3y;
                force[j * n + (i + 1)][2] -= F3z;

                force[j * n + (i + 1)][0] += F2x;
                force[j * n + (i + 1)][1] += F2y;
                force[j * n + (i + 1)][2] += F2z;

                force[(j + 1) * n + i][0] -= F4x;
                force[(j + 1) * n + i][1] -= F4y;
                force[(j + 1) * n + i][2] -= F4z;

                force[(j + 1) * n + i][0] -= F2x;
                force[(j + 1) * n + i][1] -= F2y;
                force[(j + 1) * n + i][2] -= F2z;
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

            diffx = currPos[(m - 1) * n + i + 1][0] - currPos[(m - 1) * n + i][0];
            diffy = currPos[(m - 1) * n + i + 1][1] - currPos[(m - 1) * n + i][1];
            diffz = currPos[(m - 1) * n + i + 1][2] - currPos[(m - 1) * n + i][2];
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(m - 1) * n + i + 1][0] - velocity[(m - 1) * n + i][0];
            velDiffy =
                velocity[(m - 1) * n + i + 1][1] - velocity[(m - 1) * n + i][1];
            velDiffz =
                velocity[(m - 1) * n + i + 1][2] - velocity[(m - 1) * n + i][2];

            F4x = k * diffx + shortDamp * velDiffx;
            F4y = k * diffy + shortDamp * velDiffy;
            F4z = k * diffz + shortDamp * velDiffz;

            force[(m - 1) * n + i + 1][0] -= F4x;
            force[(m - 1) * n + i + 1][1] -= F4y;
            force[(m - 1) * n + i + 1][2] -= F4z;

            force[(m - 1) * n + i][0] += F4x;
            force[(m - 1) * n + i][1] += F4y;
            force[(m - 1) * n + i][2] += F4z;
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

            diffx = currPos[(j + 1) * n + n - 1][0] - currPos[j * n + n - 1][0];
            diffy = currPos[(j + 1) * n + n - 1][1] - currPos[j * n + n - 1][1];
            diffz = currPos[(j + 1) * n + n - 1][2] - currPos[j * n + n - 1][2];
            diffNorm = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
            k = shortSpring * (diffNorm - dx) * (1 / diffNorm);
            // now force due to damper
            velDiffx =
                velocity[(j + 1) * n + n - 1][0] - velocity[j * n + n - 1][0];
            velDiffy =
                velocity[(j + 1) * n + n - 1][1] - velocity[j * n + n - 1][1];
            velDiffz =
                velocity[(j + 1) * n + n - 1][2] - velocity[j * n + n - 1][2];

            F4x = k * diffx + shortDamp * velDiffx;
            F4y = k * diffy + shortDamp * velDiffy;
            F4z = k * diffz + shortDamp * velDiffz;

            force[j * n + n - 1][0] += F4x;
            force[j * n + n - 1][1] += F4y;
            force[j * n + n - 1][2] += F4z;

            force[(j + 1) * n + n - 1][0] -= F4x;
            force[(j + 1) * n + n - 1][1] -= F4y;
            force[(j + 1) * n + n - 1][2] -= F4z;
        }

        // update the position's of the elements and velocities
        for (UINT i = 1; i < n; ++i) {
            for (UINT j = 0; j < m; ++j) {
                prevPos[j * numRows + i][0] =
                    currPos[j * n + i][0] + velocity[j * n + i][0] * dt +
                    force[j * n + i][0] * 0.5 * (1 / mass) * dt * dt;

                prevPos[j * numRows + i][1] =
                    currPos[j * n + i][1] + velocity[j * n + i][1] * dt +
                    force[j * n + i][1] * 0.5 * (1 / mass) * dt * dt;

                prevPos[j * numRows + i][2] =
                    currPos[j * n + i][2] + velocity[j * n + i][2] * dt +
                    force[j * n + i][2] * 0.5 * (1 / mass) * dt * dt;

                velocity[j * n + i][0] =
                    (prevPos[j * n + i][0] - currPos[j * n + i][0]) * (1 / dt);
                velocity[j * n + i][1] =
                    (prevPos[j * n + i][1] - currPos[j * n + i][1]) * (1 / dt);
                velocity[j * n + i][2] =
                    (prevPos[j * n + i][2] - currPos[j * n + i][2]) * (1 / dt);
            }
        }

        for (UINT i = 0; i < n - 1; ++i) {
            for (UINT j = 0; j < m - 1; ++j) {
                tangents[j * n + i][0] =
                    currPos[(j + 1) * n + i][0] - currPos[j * n + i][0];
                tangents[j * n + i][1] =
                    currPos[(j + 1) * n + i][1] - currPos[j * n + i][1];
                tangents[j * n + i][2] =
                    currPos[(j + 1) * n + i][2] - currPos[j * n + i][2];
                tangents[j * n + i].normalize();

                bitangents[j * n + i][0] =
                    currPos[j * n + i + 1][0] - currPos[j * n + i][0];
                bitangents[j * n + i][1] =
                    currPos[j * n + i + 1][1] - currPos[j * n + i][1];
                bitangents[j * n + i][2] =
                    currPos[j * n + i + 1][2] - currPos[j * n + i][2];
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
