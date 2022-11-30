#pragma once
#include <algorithm>
#include "geometry.hpp"
#include "image.hpp"
#include "utils.hpp"

//===============================================================//
// Textures
//===============================================================//
struct Texture_ref {
    
    // Texture start reference.
    Vector3 p;
    bool point_ref;

    // Texture dimension.
    int r_flags;
    std::vector<float> r_vals;

    Texture_ref(int r_flags, std::vector<float> r_vals)
        : point_ref(false), r_flags(r_flags), r_vals(r_vals) {}

    Texture_ref(Vector3 p, int r_flags, std::vector<float> r_vals)
        : p(p), point_ref(true), r_flags(r_flags), r_vals(r_vals) {}

};

class Texture {
private:

    struct Quad {
    private:
        //...
    public:
        bool alpha;
        RGB  color;
    };

public:
    int width;
    int height;
    std::vector<std::vector<Quad>> quads;

    Texture() {}
    Texture(Image i) {
        this->width  = i.width;
        this->height = i.height;
        quads = std::vector<std::vector<Quad>>(height, std::vector<Quad>(width));
        for (int h = 0; h < i.height; h++) {
            for (int w = 0; w < i.width; w++) {
                quads[h][w] = Quad{i.has_color_key && i.color_key == i.pixels[h][w],
                    i.pixels[h][w]};
            }
        }
    }

};

//===============================================================//
// Directional light.
//===============================================================//
class Light {
private:
    // ...
public:
    Vector3 c;
    RGB pow;
    Light(Vector3 c = Vector3(), RGB pow = RGB()) : c(c), pow(pow) {}
};

//===============================================================//
// Object
//===============================================================//

class Object {
private:
    virtual std::ostream& print(std::ostream& os) = 0; // He tenido que cargarme un 
                                                       // poquito el overload para que 
                                                       // se pudiese printear cualquier
                                                       // objeto.
public:
    // Object info.
    RGB kd;

    // Texture things.
    Texture t; // Object texture.
    bool has_texture;

    Object(Texture t) : t(t), has_texture(true) {}
    Object(RGB kd = RGB(0,0,0)) : kd(kd), has_texture(false) {}
    virtual ~Object() = default; // Por alguna razón, el compilador de C++ necesita un
                                 // destructor virtual, no hace absolutamente nada, pero
                                 // no toca los huevos.
                                 
    virtual RGB get_kd() { return kd; }
    RGB fr() const { return kd/M_PI; }

    virtual std::vector<float> intersects(const Ray& ray) = 0;
    virtual Vector3 normal(Vector3 p, Vector3 wi) = 0;


    friend std::ostream& operator<<(std::ostream& os, Object& p) {
        return p.print(os);
    }
};

//===============================================================//
// Plane
//===============================================================//

class Plane : public Object {
private:
    // Print plane attributes.
    std::ostream& print(std::ostream& os) override {
        return os << "PLANE {"
            << "\n  normal: "   << n 
            << "\n  distance: " << D 
            << "\n  kd: " << kd
            << "\n  finite plane bounds: " << b
            << "\n  texture reference: " << t_ref
            << "\n}";
    }

    // Texture things.
    std::vector<Vector3> t_ref; // Texture start reference and orientation.
                                //  - 0: reference point.
                                //  - 1: texture geometrical height;
                                //  - 2: texture geometrical width;
    //int qi, qj; // Index of the texture quad that has been intersected.
    float qw, qh; // Width and height of a single quad.
    float tw, th; // Width and height of the texture.
    RGB q_color;  // Color of the texture quad that has been intersected.

public:
    // Geometrical things
    float D;                 // Implicit equation A*x+B*y+C*z+D
                             //     (= 0 if the point is in the plane).
    Vector3 n;               // Normal of the plane = (A,B,C).
    std::vector<Vector3> b;  // Vertex bounds of the plane (if finite).

    Plane() {}

    // ==========================
    // Plane constructors
    // ==========================

    // Solid color plane defined by a normal and its distance to the origin.
    Plane(float D, Vector3 n, RGB kd = RGB(185,185,185)) : Object(kd), D(D), n(nor(n)) {}

    // Plane defined by a normal and a plane contained point with solid color.
    Plane(Vector3 p, Vector3 n, RGB kd = RGB(185,185,185)) : Object(kd) {
        this->n = nor(n);
        this->D = -this->n*p;
    }

    // Texturized plane defined by a normal and the distance of the plane to the origin.
    //  - t is the texture.
    //  - p is the point reference where the texture will start to be drawn.
    //  - tw and th are the texture geometrical width and height.
    //  - r is the rotation of the texture.
    //  - o is the orientation (positive or negative).
    Plane(float D, Vector3 n, Texture t, Texture_ref r, float tw, float th)
        : Object(t), tw(tw), th(th), D(D)
    {

        // Plane things.
        this->n = nor(n);

        // Texture reference.
        t_ref.push_back((r.point_ref) ? r.p : -D*this->n);
        t_ref.push_back(nor(rot(this->n, r.r_flags, r.r_vals), th));
        t_ref.push_back(nor(crs(t_ref[1], this->n), tw));
        qw = tw/t.width;
        qh = th/t.height;
    }

    // ==========================
    // Finite plane constructors
    // ==========================

    Plane(float D, std::vector<Vector3> b, RGB kd = RGB(185,185,185)) : Object(kd), D(D), b(b) {
        // Geometrical things.
        this->n = nor(crs(b[1]-b[0], b[3]-b[0]));
    }

    Plane(std::vector<Vector3> b, RGB kd = RGB(185,185,185)) : Object(kd), b(b) {
        // Geometrical things.
        this->n = nor(crs(b[1]-b[0], b[3]-b[0]));
        this->D = -n*b[0];
    }

    Plane(std::vector<Vector3> b, Texture t, Texture_ref r) : Object(t), b(b) {
        // Geometrical things.
        n = nor(crs(b[1]-b[0], b[3]-b[0]));
        D = -n*b[0];

        // Texture reference.
        tw = (b[3]-b[0]).mod();
        th = (b[1]-b[0]).mod();
        t_ref.push_back((r.point_ref) ? r.p : -D*n);
        t_ref.push_back(nor(rot(n, r.r_flags, r.r_vals), th));
        t_ref.push_back(nor(crs(t_ref[1], n), tw));
        qw = tw/t.width;
        qh = th/t.height;
    }

    RGB get_kd() override {
        if (!has_texture) return kd;
        else return q_color;
    }

    std::vector<float> intersects(const Ray& r) override {

        if (n*r.d == 0.f) return {};     // Si la división es 0 no hay corte.
        float ts = -(n*r.p + D)/(n*r.d); // Calcula la distancia desde el
                                         //   centro del rayo hasta el punto de corte.
        
        Vector3 p = (r.p + r.d*ts);
        if (b.size()) {
            if( n * crs(b[1]-b[0], p-b[0]) < 0 || // Point p inside edge 1 (v1 to v2).
                n * crs(b[2]-b[1], p-b[1]) < 0 || // Point p inside edge 2 (v2 to v3).
                n * crs(b[3]-b[2], p-b[2]) < 0 || // Point p inside edge 3 (v3 to v4).
                n * crs(b[0]-b[3], p-b[3]) < 0)   // Point p inside edge 4 (v4 to v1).
            {
                return {};
            }        
        }
        
        // Texture things.
        if (has_texture) {
            p -= t_ref[0];
            // Obtaining proportional width and height.
            float hs = t_ref[1] * (p/th);
            float ws = t_ref[2] * (p/tw);

            // Obtaining correspondent texture indexes.
            int qi = std::abs(std::fmod(hs, th)/ qh);
            int qj = std::abs(std::fmod(ws, tw)/ qw);

            // If the proportional distance is negative, we have to
            // flip the index to avoid fliped texture tiles.
            if (hs < 0) qi = (t.height - 1) - qi;
            if (ws < 0) qj = (t.width  - 1) - qj;

            // If the texture quad is transparent, it "doesn't intersects".
            if (t.quads[qi][qj].alpha) return {};

            q_color = t.quads[qi][qj].color;
        }
        return {ts}; // Devolver la distancia al punto de corte.
    }

    Vector3 normal(Vector3 p, Vector3 wi) override {
        if (n * wi <= 0) return n;
        else return -n;
    }

};

//===============================================================//
// Sphere
//===============================================================//

class Sphere : public Object {
private:
    std::ostream& print(std::ostream& os) override {
        return os << "SPHERE {"
            << "\n  center: "   << c
            << "\n  radius: "   << r 
            << "\n  kd: " << kd
            << "\n}";
    }
public:
    Vector3 c;
    float r;

    Sphere(Vector3 c = Vector3(), float r = 1.0, RGB kd = RGB()) : Object(kd), c(c), r(r) {}

    std::vector<float> intersects(const Ray& ray) override {
        float a, b, c_, q, x0, x1; Vector3 L;
        L = ray.p-c;
        a = ray.d*ray.d;
        b = 2*ray.d*L;
        c_ = L*L - r*r;
        // Solve quadratic formula.
        float discr = b*b - 4*a*c_; 
        if (discr<0) return {};
        else if (discr == 0) return {-0.5f*b/a};
        else {
            q = (b > 0) ? -0.5f*(b+sqrt(discr)) : -0.5f*(b-sqrt(discr));
            x0 = q/a; x1 = c_/q;
            if (x0 > x1) swap(x0,x1);
            return {x0, x1};
        }
    }

    Vector3 normal(Vector3 p, Vector3 wi) override {
        return c-p;
    }
};

//===============================================================//
// Box
//===============================================================//

class Box : public Object {
private:
    std::ostream& print(std::ostream& os) override {
        return os << "BOX {"
            << "\n  bounds: "   << bounds
            << "\n  kd: " << kd
            << "\n}";
    }
public:
    std::vector<Vector3> bounds;
    Vector3 center;

    Box() {}
    Box(std::vector<Vector3> bounds)
        : bounds(bounds), center((bounds[1] - bounds[0])/2) {}
    Box(Vector3 min, Vector3 max, RGB kd = RGB(185,185,185))
        : Object(kd), bounds({ min, max }), center( (max-min)/2) {}

    std::vector<float> intersects(const Ray& ray) override {

        Vector3 inv_d(1.f/ray.d.x, 1.f/ray.d.y, 1.f/ray.d.z);
        bool sign_dir_x = inv_d.x < 0;
        bool sign_dir_y = inv_d.y < 0;
        bool sign_dir_z = inv_d.z < 0;

        float p0 = (bounds[sign_dir_x].x - ray.p.x) * inv_d.x;
        float p1 = (bounds[1 - sign_dir_x].x - ray.p.x) * inv_d.x;
        float b0 = (bounds[sign_dir_y].y - ray.p.y) * inv_d.y;
        float b1 = (bounds[1 - sign_dir_y].y - ray.p.y) * inv_d.y;
        if ((p0 > b1) || (b0 > p1)) return {};
        if (b0 > p0) p0 = b0;
        if (b1 < p1) p1 = b1;

        b0 = (bounds[sign_dir_z].z - ray.p.z) * inv_d.z;
        b1 = (bounds[1 - sign_dir_z].z - ray.p.z) * inv_d.z;
        if ((p0 > b1) || (b0 > p1)) return {};
        if (b0 > p0) p0 = b0;
        if (b1 < p1) p1 = b1;

        if (p0 < 0) {
            if (p1 < 0) return {};
            else return { p1 };
        } else return { p0 };
    }

    // Translate one of the bounds along half the distance of the direction
    // created from that bound to the other to obtain the center. Then, n = p-c.
    Vector3 normal(Vector3 p, Vector3 wi) override { 
        return p-(bounds[0] + 0.5f * (bounds[1] - bounds[0]));
    }
};

//===============================================================//
// Meshes and its faces: triangles.
//===============================================================//

class Triangle : public Plane {
private:
    std::ostream& print(std::ostream& os) override {
        return os << "TRIANGLE {" 
            << "\n vertex: "   << v 
            << "\n kd: " << kd
            << "\n}";
    }
public:
    std::vector<Vector3> v; // Triangle vertex.

    Triangle() {}
    Triangle(std::vector<Vector3> v, RGB kd = RGB(185,185,185)) 
        : Plane(v[0], crs(v[1]-v[0], v[2]-v[0]), kd), v(v) {
    }

    std::vector<float> intersects(const Ray& r) override {
        std::vector<float> t;
        if ((t = Plane::intersects(r)).empty()) return {};

        Vector3 p = r.p + r.d*t[0];
        if (n * crs(v[1]-v[0], p-v[0]) < 0 ||          // Point p inside edge 1 (v1 to v2).
            n * crs(v[2]-v[1], p-v[1]) < 0 ||          // Point p inside edge 2 (v2 to v3).
            n * crs(v[0]-v[2], p-v[2]) < 0) return {}; // Point p inside edge 3 (v3 to v1).
        return t;
    }

};

class Mesh : public Object {
private:

    std::ostream& print(std::ostream& os) override {
        return os << "MESH {"
            << "\n  faces: "        << faces
            << "\n  kd: "     << kd
            << "\n  box collider: " << collider
            << "\n}";
    }

    // Ply metadata
    // nothing for the moment.
    Box collider;

public:

    Triangle q_dot;
    std::vector<Triangle> faces;

    Mesh(std::vector<Triangle> faces) : Object(RGB(185,185,185)), collider(Vector3(INFINITY, INFINITY, INFINITY), 
        Vector3(-INFINITY, -INFINITY, -INFINITY)), faces(faces)
    {   
        for (auto& f : faces) {
            for (auto& v : f.v) {
                collider.bounds[0].x = std::min(v.x, collider.bounds[0].x);
                collider.bounds[0].y = std::min(v.y, collider.bounds[0].y);
                collider.bounds[0].z = std::min(v.z, collider.bounds[0].z);
                collider.bounds[1].x = std::max(v.x, collider.bounds[1].x);
                collider.bounds[1].y = std::max(v.y, collider.bounds[1].y);
                collider.bounds[1].z = std::max(v.z, collider.bounds[1].z);
            }
        }
    }

    Mesh(std::string ply_file) : Object(RGB(185,185,185)), 
        collider(Vector3(INFINITY, INFINITY, INFINITY), Vector3(-INFINITY, -INFINITY, -INFINITY))
    {
        std::string s("");
        std::ifstream in(ply_file);
        assert(in.is_open() && "file not found, check it out.");
        std::vector<Vector3> vertex;

        while (s.compare("end_header")) {
            s = get_line(in);
            if (s.find("element vertex") != std::string::npos) {
                vertex = std::vector<Vector3>(
                    std::stoi(replace(s, int_d))
                );
            } else if (s.find("element face") != std::string::npos) {
                faces  = std::vector<Triangle>(
                    std::stoi(replace(s, int_d))
                );
            }
        }
        assert(vertex.size() && faces.size() && "no vertex..? no faces..?"); // no bitches..? Sorry

        // Reading the vertex list.
        for (auto& v : vertex) {
            in >> v;
            collider.bounds[0].x = std::min(v.x, collider.bounds[0].x);
            collider.bounds[0].y = std::min(v.y, collider.bounds[0].y);
            collider.bounds[0].z = std::min(v.z, collider.bounds[0].z);
            collider.bounds[1].x = std::max(v.x, collider.bounds[1].x);
            collider.bounds[1].y = std::max(v.y, collider.bounds[1].y);
            collider.bounds[1].z = std::max(v.z, collider.bounds[1].z);
        }

        // Reading the faces list.
        int n, v1, v2, v3;
        for (auto& f : faces) {
            in >> n >> v1 >> v2 >> v3;
            f = Triangle({vertex[v1], vertex[v2], vertex[v3]}, kd);
        }
    }

    RGB get_kd() override {
        return q_dot.kd;
    }

    RGB get_reflectance() const {
        return q_dot.kd/M_PI; 
    }

    std::vector<float> intersects(const Ray& r) override {
        if (collider.intersects(r).empty()) return {};
        
        std::vector<float> ts;
        for (auto& f : faces) {
            std::vector<float> tp;
            if ((tp = f.intersects(r)).size() && tp[0] > 0 && !insert(ts, tp[0])) {
                q_dot = f;
            }
        }
        return ts;
    }

    // Calcular en que triangulo queda el punto.
    Vector3 normal(Vector3 p, Vector3 wi) override {
        return q_dot.n;
    }

    friend Mesh operator*(Mesh m, Matrix3 transform);
};

Mesh operator*(Mesh m, Matrix3 transform) {
    std::vector<Vector3> b{
        Vector3(INFINITY, INFINITY, INFINITY),
        Vector3(-INFINITY, -INFINITY, -INFINITY)
    };
    for (auto& f : m.faces) {
        for (auto& v : f.v) {
            v.h = 1;
            v = transform * v;
            b[0].x = std::min(v.x, b[0].x);
            b[0].y = std::min(v.y, b[0].y);
            b[0].z = std::min(v.z, b[0].z);
            b[1].x = std::max(v.x, b[1].x);
            b[1].y = std::max(v.y, b[1].y);
            b[1].z = std::max(v.z, b[1].z);
        }
        f = Triangle(f.v, f.kd);
    }
    m.collider.bounds = b;
    return m;
};