#include "antmath.h"
#include <iostream>

AV4X4FLOAT::AV4X4FLOAT(const float* a) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) m[j * 4 + i] = a[j * 4 + i];
};

AV4X4FLOAT::AV4X4FLOAT() {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) m[j * 4 + i] = 0;
};
AV4X4FLOAT AV4X4FLOAT::operator*(const AV4X4FLOAT M) {
    float data[16];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            float d = 0;
            for (int k = 0; k < 4; k++) {
                d += m[k * 4 + i] * M(k, j);
            }
            data[j * 4 + i] = d;
        }
    }
    AV4X4FLOAT C(data);
    return C;
};

float AV4X4FLOAT::operator()(UINT r, UINT c) const { return m[c * 4 + r]; }
float& AV4X4FLOAT::operator()(UINT r, UINT c) { return m[c * 4 + r]; }

float AV4X4FLOAT::det() {
    // will we ever need this??!?!!?!?
    return 0;
};

AV4FLOAT AV4X4FLOAT::getRow(UINT i) {
    AV4FLOAT r(m[0 * 4 + i], m[1 * 4 + i], m[2 * 4 + i], m[3 * 4 + i]);
    return r;
};

void AV4X4FLOAT::setRow(UINT i, const AV4FLOAT r) {
    m[0 * 4 + i] = r.x;
    m[1 * 4 + i] = r.y;
    m[2 * 4 + i] = r.z;
    m[3 * 4 + i] = r.w;
};
void AV4X4FLOAT::setCol(UINT i, const AV4FLOAT r) {
    m[i * 4 + 0] = r.x;
    m[i * 4 + 1] = r.y;
    m[i * 4 + 2] = r.z;
    m[i * 4 + 3] = r.w;
}

AV4X4FLOAT AV4X4FLOAT::inverse() {
    // hard coding an inverse, we can test at some later point
    // it seems as though matrix inversion is not needed often.

    float inv[16], det;
    int i;

    inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] +
             m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

    inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] +
             m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] +
             m[12] * m[7] * m[10];

    inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] +
             m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

    inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] +
              m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] +
              m[12] * m[6] * m[9];

    inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] +
             m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] +
             m[13] * m[3] * m[10];

    inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] +
             m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

    inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] -
             m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

    inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] +
              m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

    inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] +
             m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];

    inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] -
             m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

    inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] +
              m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

    inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] -
              m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] -
             m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] -
              m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
    det = 1.0 / det;
    for (i = 0; i < 16; ++i) inv[i] = inv[i] * det;

    AV4X4FLOAT C(inv);
    return C;
}

AV4X4FLOAT AV4X4FLOAT::transpose() {
    float t[16];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            t[j * 4 + i] = m[i * 4 + j];
        }
    }
    AV4X4FLOAT C(t);
    return C;
}

void AV4X4FLOAT::diag(AV4FLOAT d) {
    m[0] = d.x;
    m[5] = d.y;
    m[10] = d.z;
    m[15] = d.w;
}

void AV4X4FLOAT::display() {
    std::cout << std::endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cout << m[j * 4 + i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

float dotProduct4(AV4FLOAT u, AV4FLOAT v) {
    return u.x * v.x + u.y * v.y + u.w * v.w + u.z * v.z;
}
float dotProduct3(AV3FLOAT u, AV3FLOAT v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

AV4X4FLOAT formViewModelMatrix(AV4FLOAT pos, AV4FLOAT target, AV4FLOAT up) {
    AV4FLOAT mz;
    mz.x = pos.x - target.x;
    mz.y = pos.y - target.y;
    mz.z = pos.z - target.z;
    mz.w = 0.0f;
    mz.normalize();

    AV4FLOAT my;
    my.x = up.x;
    my.y = up.y;
    my.z = up.z;
    my.w = 0.0f;

    AV4FLOAT mx;
    mx.x = my.y * mz.z - my.z * mz.y;
    mx.y = my.z * mz.x - my.x * mz.z;
    mx.z = my.x * mz.y - my.y * mz.x;
    mx.w = 0.0f;
    mx.normalize();

    my.x = mz.y * mx.z - mz.z * mx.y;
    my.y = mz.z * mx.x - mz.x * mx.z;
    my.z = mz.x * mx.y - mz.y * mx.x;
    my.w = 0.0f;

    AV4FLOAT t;
    t.x = mx.x * pos.x + mx.y * pos.y + mx.z * pos.z;
    t.y = my.x * pos.x + my.y * pos.y + my.z * pos.z;
    t.z = -(mz.x * pos.x + mz.y * pos.y + mz.z * pos.z);

    AV4X4FLOAT m;
    m.m[0] = mx.x;
    m.m[1] = my.x;
    m.m[2] = mz.x;
    m.m[3] = 0.0f;
    m.m[4] = mx.y;
    m.m[5] = my.y;
    m.m[6] = mz.y;
    m.m[7] = 0.0f;
    m.m[8] = mx.z;
    m.m[9] = my.z;
    m.m[10] = mz.z;
    m.m[11] = 0.0f;
    m.m[12] = t.x;
    m.m[13] = t.y;
    m.m[14] = t.z;
    m.m[15] = 1.0f;

    return m;
}

AV4X4FLOAT formProjMatrix(float FOVangle, float aspect, float nearz,
                          float farz) {
    AV4X4FLOAT A;

    A.m[0] = 1 / (aspect * tanf(FOVangle / 2));
    A.m[5] = 1 / tanf(FOVangle / 2);
    A.m[10] = (nearz + farz) / (farz - nearz);
    A.m[11] = -2.0 * nearz * farz / (farz - nearz);
    A.m[14] = -1.0;

    return A;
}

float ANTMATHConvertToRadians(float degree) {
    return degree * AL_PI * (1 / 180.0);
}

void disp4(AV4FLOAT u) {
    std::cout << u.x << " " << u.y << " " << u.z << " " << u.w << std::endl;
}

void disp3(AV3FLOAT& u) {
    std::cout << u.x << " " << u.y << " " << u.z << std::endl;
}

void add3(AV3FLOAT& u, AV3FLOAT& v, AV3FLOAT& result) {
    result.x = u.x + v.x;
    result.y = u.y + v.y;
    result.z = u.z + v.z;
}

void addEquals3(AV3FLOAT& result, AV3FLOAT& v) {
    result.x += v.x;
    result.y += v.y;
    result.z += v.z;
}
void subtract3(AV3FLOAT& u, AV3FLOAT& v, AV3FLOAT& result) {
    result.x = u.x - v.x;
    result.y = u.y - v.y;
    result.z = u.z - v.z;
}

void crossProduct3(AV3FLOAT& u, AV3FLOAT& v, AV3FLOAT& result) {
    result.x = u.y * v.z - v.y * u.z;
    result.y = -u.x * v.z + v.x * u.z;
    result.z = u.x * v.y - v.x * u.y;
}
