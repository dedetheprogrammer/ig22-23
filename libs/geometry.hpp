#pragma once
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <iostream>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define N 4
#define EPSILON_ERROR 0.000001

enum Rotation { X_ROT = 0x1, Y_ROT = 0x2, Z_ROT = 0x4 };

//===============================================================//
// Vector3: 3d vector
//===============================================================//
/**
 * @brief Three dimension vector. Four components:
 *  - Coordinate x.
 *  - Coordinate y.
 *  - Coordinate z.
 *  - Homogeneous coordinate h. Used for matrix transformators:
 *      - 0 if Direction.
 *      - 1 if Point.
 */
class Vector3 {
private:

    // Para limpiar números feos que se generen por redondeo o por 
    // operaciones matemáticas. Me da muchisimo toc ver un 0 negativo
    // o un número elevado a la -12.
    void clean() {
        if (std::abs(x) < EPSILON_ERROR || x == -0) x = 0;
        if (std::abs(y) < EPSILON_ERROR || y == -0) y = 0;
        if (std::abs(z) < EPSILON_ERROR || z == -0) z = 0;
    }

public:

    double x, y, z, h;

    Vector3 (double x = 0, double y = 0, double z = 0)
        : x(x), y(y), z(z), h(0) { clean(); }
    Vector3 (double x, double y, double z, int h)
        : x(x), y(y), z(z), h(h) { clean(); }
    Vector3 (Vector3 v, int h)
        : x(v.x), y(v.y), z(v.z), h(h) { clean(); }
    Vector3 (double m[4])
        : x(m[0]), y(m[1]), z(m[2]), h(m[3]) { clean(); }

    // Vector module
    double mod() const {
        return sqrt((x * x) + (y * y) + (z * z));
    }

    // Vector asignation.
    void operator=(Vector3 v) {
        x = v.x;
        y = v.y;
        z = v.z;
        clean();
    }

    // Vector add.
    void operator+=(Vector3 v) {
        x += v.x;
        y += v.y; 
        z += v.z;
    }

    // Vector substract.
    void operator-=(Vector3 v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    void operator*=(double d) {
        x *= d;
        y *= d;
        z *= d; 
    }

};

// Vector add.
Vector3 operator+(Vector3 v, Vector3 w) {
    return Vector3(v.x+w.x, v.y+w.y, v.z+w.z);
}

// Vector substract.
Vector3 operator-(Vector3 v, Vector3 w) {
    return Vector3(v.x-w.x, v.y-w.y, v.z-w.z);
}

// Vector sign change.
Vector3 operator-(Vector3 v) {
    return Vector3(-v.x, -v.y, -v.z);
}

// Vector dot product.
double operator*(Vector3 v, Vector3 w) {
    return (v.x * w.x) + (v.y * w.y) + (v.z * w.z) + (v.h * w.h);
}

// Vector scalar product.
Vector3 operator*(Vector3 v, double r) {
    return Vector3(v.x * r, v.y * r, v.z * r);
}

// Vector scalar product.
Vector3 operator*(double r, Vector3 v) {
    return v * r;
}

// Vector scalar division.
Vector3 operator/(Vector3 v, double r) {
    return Vector3(v.x/r, v.y/r, v.z/r);
}

// Vector division (be careful, this doesn't make sense).
Vector3 operator/(Vector3 v, Vector3 w) {
    return Vector3(v.x / w.x, v.y / w.y, v.z / w.z);
}

// Vector absolute.
inline Vector3 abs(Vector3 v) {
    if (v.x < 0) v.x *= -1;
    if (v.y < 0) v.y *= -1;
    if (v.z < 0) v.z *= -1;
    return v;
}

// Angle between two vectors in radians.
double rad(Vector3 v, Vector3 w) {
    return acos((v * w)/(v.mod() * w.mod()));
}

// Angle between two vectors in grades.
double grd(Vector3 v, Vector3 w) {
    return rad(v, w) * (180/M_PI);
}

// Vector equality.
bool operator==(Vector3 v, Vector3 w) {
    return (v.x - w.x < EPSILON_ERROR) 
        && (v.y - w.y < EPSILON_ERROR)
        && (v.z - w.z < EPSILON_ERROR);
}

// Vector inequality.
bool operator!=(Vector3 v, Vector3 w) {
    return !(v == w);
}

// Direction normalization. By default, it normalize to the unitary value.
Vector3 nor(Vector3 v, double mod = 1.0) {
    return v * mod/v.mod();
}

// Vector cross product.
Vector3 crs(Vector3 v, Vector3 w) {
    return Vector3(
        (v.y * w.z) - (v.z * w.y),
        (v.z * w.x) - (v.x * w.z),
        (v.x * w.y) - (v.y * w.x)
    );
}

// How to rotate a vector in reference to other vector. Rotating a respect of b.
// https://math.stackexchange.com/questions/511370/how-to-rotate-one-vector-about-another
Vector3 rot(Vector3 a, Vector3 b, double r = M_PI/2) {

    Vector3 abb = ((a*b)/(b*b))*b; // a component in the direction of b.
    Vector3 abp = a - abb;         // a component in the direction orthogonal to b.
    Vector3 w = crs(b, abp);       // w component, orthogonal to a and b.

    double abp_m = abp.mod();
    double x1 = cos(r)/abp_m;
    double x2 = sin(r)/abp_m;
    return (abp_m * (x1*abp + x2*w)) + abb;
}

// Orthonormal basis
// https://graphics.pixar.com/library/OrthonormalB/paper.pdf
std::vector<Vector3> orthonormal_basis(const Vector3& n) {

    double sign = copysignf(1.0f, n.z);
    const double a = -1.0f / (sign + n.z);
    const double b = n.x * n.y * a;
    return {
        Vector3(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x),
        Vector3(b, sign + n.y * n.y *a, -n.y)
    };

}

// Print vector out.
std::ostream& operator<<(std::ostream& os, const Vector3& v) {
    return os << "(" << v.x << "," << v.y << "," << v.z << ")";
}

// Initialize vector.
void operator>>(std::istream& in, Vector3& v) {
    in >> v.x >> v.y >> v.z;
}

/*
void operator>>(std::istream& is, Vector3& v) {
    std::cout << "\n    px: "; is >> v.x;
    std::cout <<   "    py: "; is >> v.y;
    std::cout <<   "    pz: "; is >> v.z; std::cout << "  )";
    return is;
}
*/

//===============================================================//
// Ray: 3d ray.
//===============================================================//

class Ray {
private:
    // ...
public:
    Vector3 p, d;
    Ray () : p(Vector3()) {}
    Ray (Vector3 p, Vector3 d) : p(p), d(d) {}
};

std::ostream& operator<<(std::ostream& os, const Ray& r) {
    return os << "RAY {" << "ORIGIN: " << r.p << ", DIRECTION: " << r.d << "}";
}

//===============================================================//
// Matrix3: 3d matrix.
//===============================================================//
class Matrix3 {
private:

    // Coeficient matrix.
    void cof3(double A[N][N], double temp[N][N], int row, int col, int n) {
        int i = 0, j = 0;
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                if (x != row && y != col) {
                    temp[i][j++] = A[x][y];
                    if (j == n-1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }

    // Matrix determinant.
    double det3(double A[N][N], int n) {
        if (n == 1) return A[0][0];
        
        double d = 0, temp[N][N];
        for (int f = 0; f < n; f++) {
            if (A[0][f]) {
                cof3(A, temp, 0, f, n);
                d += (f % 2 == 0 ? 1 : -1) * A[0][f] * det3(temp, n - 1);
            }
        }
        return d;
    }

    // Adjunct matrix.
    void adj3(double A[N][N], double adj[N][N]) {
        double temp[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                cof3(A, temp, i, j, N);
                adj[j][i] = (det3(temp, N - 1));
                if (adj[j][i]) adj[j][i] *= (((i + j) % 2 == 0) ? 1 : -1);
            }
        }
    }

public:

    double m[4][4]; // Matrix values.

    Matrix3(bool identity = true) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) m[i][j] = 0;
        }
        if (identity) { 
            for (int i = 0; i < 4; i++) m[i][i] = 1;
        }
    }
    Matrix3(double n[4][4]) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m[i][j] = n[i][j];
            }
        }
    }
    
    virtual Matrix3 invert() {
        double det = det3(m, N);
        if (!det) return *this;

        double aux[N][N];
        adj3(m, aux);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                aux[i][j] = aux[i][j] / det;
            }
        }
        return Matrix3(aux);
    }
};

// Matrix product.
Matrix3 operator* (Matrix3 t1, Matrix3 t2) {
    Matrix3 tf(false);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            tf.m[i][j] = Vector3(t1.m[i]) * 
                Vector3(t2.m[0][j], t2.m[1][j], t2.m[2][j], t2.m[3][j]);
        }
    }
    return tf;
}

// Matrix and vector product.
Vector3 operator* (Matrix3 t, Vector3 v) {
    return Vector3(
        (Vector3(t.m[0]) * v),
        (Vector3(t.m[1]) * v),
        (Vector3(t.m[2]) * v),
        (Vector3(t.m[3]) * v)
    );
}

// Print matrix out.
std::ostream& operator<< (std::ostream& os, const Matrix3& t) {
    os << std::endl;
    for (int i = 0; i < 4; i++) {
        os << "[ ";
        for (int j = 0; j < 4; j++) os << t.m[i][j] << " ";
        os << "]" << std::endl;
    }
    return os;
}

//==============================//
// 3d matrix translation.
//==============================//
class Matrix3Translation : public Matrix3 {
private:
    // ...
public:

    Matrix3Translation(double translation) {
        m[0][3] = m[1][3] = m[2][3] = translation;
    }

    Matrix3Translation(double tx, double ty, double tz) {
        m[0][3] = tx;
        m[1][3] = ty;
        m[2][3] = tz;
    }

    Matrix3 invert() override {
        return Matrix3Translation(m[0][3] * -1, m[1][3] * -1, m[2][3] * -1);
    }
};

//==============================//
// 3d matrix scaling.
//==============================//
class Matrix3Scale : public Matrix3 {
private:
    // ...
public:
    Matrix3Scale(double scale) {
        m[0][0] = m[1][1] = m[2][2] = scale;
    }

    Matrix3Scale(double sx, double sy, double sz) {
        m[0][0] = sx;
        m[1][1] = sy;
        m[2][2] = sz;
    }

    Matrix3 invert() override {
        return Matrix3Scale(1/m[0][0], 1/m[1][1], 1/m[2][2]);
    }
};

//==============================//
// 3d matrix rotation.
//==============================//
class Matrix3Rotation : public Matrix3 {
private:
    //double rotation;
public:

    Matrix3Rotation(int flags, double degrees) {
        double r = (degrees * M_PI)/180;
        if (flags & X_ROT) {
            m[1][1] =  cos(r); if (std::abs(m[1][1]) < EPSILON_ERROR) m[1][1] = 0;
            m[1][2] = -sin(r); if (std::abs(m[1][2]) < EPSILON_ERROR) m[1][2] = 0;
            m[2][1] =  sin(r); if (std::abs(m[2][1]) < EPSILON_ERROR) m[2][1] = 0;
            m[2][2] =  cos(r); if (std::abs(m[2][2]) < EPSILON_ERROR) m[2][2] = 0;   
        } else if (flags & Y_ROT) {
            m[0][0] =  cos(r); if (std::abs(m[0][0]) < EPSILON_ERROR) m[0][0] = 0;
            m[0][2] =  sin(r); if (std::abs(m[0][2]) < EPSILON_ERROR) m[0][2] = 0;
            m[2][0] = -sin(r); if (std::abs(m[2][0]) < EPSILON_ERROR) m[2][0] = 0;
            m[2][2] =  cos(r); if (std::abs(m[2][2]) < EPSILON_ERROR) m[2][2] = 0;
        } else if (flags & Z_ROT) {
            m[0][0] =  cos(r); if (std::abs(m[0][0]) < EPSILON_ERROR) m[0][0] = 0;
            m[0][1] = -sin(r); if (std::abs(m[0][1]) < EPSILON_ERROR) m[0][1] = 0;
            m[1][0] =  sin(r); if (std::abs(m[1][0]) < EPSILON_ERROR) m[1][0] = 0;
            m[1][1] =  cos(r); if (std::abs(m[1][1]) < EPSILON_ERROR) m[1][1] = 0;
        }
    }

};

// Ortogonal vector.
Vector3 orto(Vector3 v, int flags, int sense) {
    if (flags & X_ROT) v = Matrix3Rotation(X_ROT, sense * 90) * v;
    if (flags & Y_ROT) v = Matrix3Rotation(Y_ROT, sense * 90) * v;
    if (flags & Z_ROT) v = Matrix3Rotation(Z_ROT, sense * 90) * v;
    return v;
}

// Vector rotation.
Vector3 rot(Vector3 v, int flags, std::vector<double> r) {
    if (flags & X_ROT) {
        v = Matrix3Rotation(X_ROT, r[0]) * v;
        r.erase(r.begin());
    }
    if (flags & Y_ROT) {
        v = Matrix3Rotation(Y_ROT, r[0]) * v;
        r.erase(r.begin());
    }
    if (flags & Z_ROT) {
        v = Matrix3Rotation(Z_ROT, r[0]) * v;
        r.erase(r.begin());
    }
    return v;
}

//==============================//
// 3d matrix base change.
//==============================//
class Matrix3BaseChange : public Matrix3 {
private:
    // ...
public:
    Matrix3BaseChange(Vector3 u, Vector3 v, Vector3 w, Vector3 o) {
        m[0][0] = u.x; m[0][1] = v.x; m[0][2] = w.x; m[0][3] = o.x;
        m[1][0] = u.y; m[1][1] = v.y; m[1][2] = w.y; m[1][3] = o.y;
        m[2][0] = u.z; m[2][1] = v.z; m[2][2] = w.z; m[2][3] = o.z;
        m[3][0] = u.h; m[3][1] = v.h; m[3][2] = w.h; m[3][3] = o.h;
    }
};

#endif