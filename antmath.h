#pragma once
#include <math.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#define AL_PI 3.141592654f
#define AL_2PI 6.283185307f
#define AL_1DIVPI 0.318309886f
#define AL_1DIV2PI 0.159154943f
#define AL_PIDIV2 1.570796327f
#define AL_PIDIV4 0.785398163f
#define AL_SQRT2 1.414213562f

typedef unsigned int UINT;

typedef struct _AV2FLOAT {
    float x, y;
    _AV2FLOAT(){};
    _AV2FLOAT(float X, float Y) : x(X), y(Y){};

    _AV2FLOAT& operator=(const _AV2FLOAT& V) {
        x = V.x;
        y = V.y;
        return *this;
    }
    bool operator==(const _AV2FLOAT& v) const {
        return (fabs(x - v.x) < 0.001 && fabs(y - v.y) < 0.001);
    }
    bool operator!=(const _AV2FLOAT& v) { return !(*this == v); }

} AV2FLOAT;

typedef struct _AV3FLOAT {
    float x, y, z;
    _AV3FLOAT(){};
    _AV3FLOAT(float X, float Y, float Z) : x(X), y(Y), z(Z){};

    _AV3FLOAT& operator=(const _AV3FLOAT& V) {
        x = V.x;
        y = V.y;
        z = V.z;
        return *this;
    };

    void normalize() {
        float mag = sqrt(x * x + y * y + z * z);
        x = x / mag;
        y = y / mag;
        z = z / mag;
    }
} AV3FLOAT;

typedef struct _AV4FLOAT {
    float x, y, z, w;
    _AV4FLOAT(){};
    _AV4FLOAT(float X, float Y, float Z, float W) : x(X), y(Y), z(Z), w(W){};

    _AV4FLOAT& operator=(const _AV4FLOAT& V) {
        x = V.x;
        y = V.y;
        z = V.z;
        w = V.w;
        return *this;
    }
    void normalize() {
        float mag = sqrt(x * x + y * y + z * z + w * w);
        if (mag == 0) {
        } else {
            x = x / mag;
            y = y / mag;
            z = z / mag;
            w = w / mag;
        }
    }

} AV4FLOAT;

typedef struct _AV3INT {
    int x, y, z;
    _AV3INT(){};
    _AV3INT(int X, int Y, int Z) : x(X), y(Y), z(Z){};
    int operator[](int index) {
        if (index > 2) {
            std::cout << "We should be throwing errors now" << std::endl;
            return -1;
        } else if (index == 0) {
            return x;
        } else if (index == 1) {
            return y;
        } else {
            return z;
        }
    }
    void setVertex(int index, int value) {
        if (index > 2) {
            std::cout << "We should be throwing errors now" << std::endl;

        } else if (index == 0) {
            x = value;
        } else if (index == 1) {
            y = value;
        } else {
            z = value;
        }
    }
    _AV3INT& operator=(const _AV3INT& threeint) {
        x = threeint.x;
        y = threeint.y;
        z = threeint.z;
        return *this;
    }
} AV3INT;

typedef struct _AV4INT {
    int x, y, z, w;
    _AV4INT(){};
    _AV4INT(int X, int Y, int Z, int W) : x(X), y(Y), z(Z), w(W){};
} AV4INT;

typedef struct _Vertex {
    UINT VertexNumber;
    AV2FLOAT Position;
    bool Ear;
    int Boundary; /*0 is interior node, 1 Dirichlet, 2 Neumann*/
    _Vertex() {
        Boundary = 0;
        VertexNumber = 0;
        Position.x = 0;
        Position.y = 0;
        Ear = 0;
    }
    _Vertex(const UINT vnum, const AV2FLOAT& p, const bool ear,
            const int bdy = 0)
        : VertexNumber(vnum), Position(p), Ear(ear), Boundary(bdy) {}
    _Vertex& operator=(const _Vertex& v) {
        Boundary = v.Boundary;
        VertexNumber = v.VertexNumber;
        Position.x = v.Position.x;
        Position.y = v.Position.y;
        Ear = v.Ear;
        return *this;
    }
    bool operator==(const _Vertex& v) { return VertexNumber == v.VertexNumber; }
    bool operator!=(const _Vertex& v) { return !(*this == v); }
    float distance(const _Vertex& v) {
        return sqrt(
            ((Position.x - v.Position.x) * (Position.x - v.Position.x)) +
            (((Position.y - v.Position.y) * (Position.y - v.Position.y))));
    }

} Vertex;

struct Triangle {
    // vertex numbers of the vertices
    int x, y, z;
    // edge numbers of the edges
    int e1, e2, e3;
    // marked to be swapped
    int swap;
    Triangle() { swap = 0; }
    Triangle(int xt, int yt, int zt) {
        x = xt;
        y = yt;
        z = zt;
        swap = 0;
    }
    void setVertices(int xt, int yt, int zt) {
        x = xt;
        y = yt;
        z = zt;
        swap = 0;
    }
};

// Edge structure
struct Edge {
    // store two vertex numbers
    int vnum1, vnum2;
    // boundary tag
    int boundary;
    // triangle numbers shared
    //
    Edge() {
        vnum1 = 0;
        vnum2 = 0;
        boundary = 0;
    }
    Edge(int v1, int v2) {
        vnum1 = v1;
        vnum2 = v2;
        boundary = 0;
    }

    Edge(int v1, int v2, int bdy) {
        vnum1 = v1;
        vnum2 = v2;
        boundary = bdy;
    }
    bool operator==(const Edge d) {
        if (vnum1 == d.vnum1) {
            return vnum2 == d.vnum2;
        }
        if (vnum2 == d.vnum1) {
            return vnum1 == d.vnum2;
        }
        return 0;
    }
    void setEndPoints(int v1, int v2) {
        vnum1 = v1;
        vnum2 = v2;
    }
};

class AV4X4FLOAT {
public:
    // consider reserving memory here.
    float m[16];
    AV4X4FLOAT();
    AV4X4FLOAT(float m00, float m01, float m02, float m03, float m10, float m11,
               float m12, float m13, float m20, float m21, float m22, float m23,
               float m30, float m31, float m32, float m33);
    AV4X4FLOAT(const float* a);
    void diag(AV4FLOAT d);
    // what the fuk is this Anthony?!
    float operator()(UINT r, UINT c) const;
    float& operator()(UINT r, UINT c);
    AV4X4FLOAT operator*(const AV4X4FLOAT M);
    AV4X4FLOAT inverse();
    float det();
    AV4FLOAT getRow(UINT i);
    void setRow(UINT i, const AV4FLOAT r);
    void setCol(UINT i, const AV4FLOAT r);
    AV4X4FLOAT transpose();
    void display();
};

extern AV4X4FLOAT formViewModelMatrix(AV4FLOAT pos, AV4FLOAT target,
                                      AV4FLOAT up);
extern AV4X4FLOAT formProjMatrix(float FOVangle, float aspect, float nearz,
                                 float farz);
extern float ANTMATHConvertToRadians(float degree);

// operations on structs
// for AV4FLOAT
extern float dotProduct4(AV4FLOAT u, AV4FLOAT v);
extern void disp4(AV4FLOAT u);
// for AV3FLOAT
extern void add3(AV3FLOAT& u, AV3FLOAT& v, AV3FLOAT& result);
extern void addEquals3(AV3FLOAT& result, AV3FLOAT& v);
extern void subtract3(AV3FLOAT& u, AV3FLOAT& v, AV3FLOAT& result);
extern void disp3(AV3FLOAT& u);
extern void crossProduct3(AV3FLOAT& u, AV3FLOAT& v, AV3FLOAT& result);
extern float dotProduct3(AV3FLOAT u, AV3FLOAT v);

/*
 *Overides ouput/input operator to allow easy output/input of AV2FLOATS
 */

inline std::ostream& operator<<(std::ostream& os, const AV2FLOAT& v) {
    std::ostream wrap(os.rdbuf());
    wrap.imbue(std::locale("C"));

    wrap << '(' << std::setprecision(4) << v.x << ", " << v.y << ')';
    os.setstate(wrap.rdstate());
    return os;
}

inline std::istream& operator>>(std::istream& is, AV2FLOAT& v) {
    std::istream wrap(is.rdbuf());
    wrap.imbue(std::locale("C"));

    auto expect = [&](char const expected) {
        char actual;
        if (!(wrap >> actual)) {
            return false;
        }
        if (actual != expected) {
            wrap.setstate(wrap.rdstate() | std::ios::failbit);
        }
        return actual == expected;
    };

    expect('(') && (wrap >> v.x) && expect(',') && (wrap >> v.y) && expect(')');

    is.setstate(wrap.rdstate());
    return is;
}

