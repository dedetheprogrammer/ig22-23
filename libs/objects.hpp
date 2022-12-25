#pragma once
#include <algorithm>
#include <random>
#include "geometry.hpp"
#include "image.hpp"
#include "utils.hpp"

//===============================================================//
// Textures
//===============================================================//
/*struct Texture_ref {
    
    // Texture start reference.
    Vector3 p;
    bool point_ref;

    // Texture dimension.
    int r_flags;
    std::vector<double> r_vals;

    Texture_ref(int r_flags, std::vector<double> r_vals)
        : point_ref(false), r_flags(r_flags), r_vals(r_vals) {}

    Texture_ref(Vector3 p, int r_flags, std::vector<double> r_vals)
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

};*/

//===============================================================//
// Materials
//===============================================================//

struct Sample {

    Vector3 wi;   // Resultant direction.
    RGB fr;       // Color factor.

    // In some materials getting the direct light doesn't make sense with the
    // way the scattering is simulated. For example:
    //
    // With specularity, we redirect the direction based in the Snell's theorem,
    // the probability that the resultant direction intersects with a point 
    // light is practically 0 as we redirect in purpose the direction (yes,
    // could be a light right on the direction trayectory, but the chance is 
    // too small to take in account).
    //
    // Then we indicate if the d_light has to be calculated or not.
    bool is_delta;

    Sample()
        : wi(Vector3()), fr(RGB()), is_delta(true) {}
    Sample(Vector3 wi, RGB fr, bool is_delta)
        : wi(wi), fr(fr), is_delta(is_delta) {}
};

class Material {
private:

    // First version before Fresnell.
    void coefficient_correction() {

        pd = max(kd);
        ps = max(ks);
        pt = max(kt);

        double coeff = pd + ps + pt;
        if (coeff > 1) {
            pd /= coeff;
            ps /= coeff;
            pt /= coeff;
        }
    }

    // etat = ref_index_i, etai = ref_index_o.
    double fresnel_coefficients(Vector3& n, Vector3 wo, double ref_index_o, double ref_index_i) {
        
        if ((n * wo) > 0) {
            n *= -1;
            std::swap(ref_index_o, ref_index_i);
        }
        wo = nor(wo);
        double ref_coef = ref_index_o/ref_index_i;
        double cos_i = n * wo;
        double cos_t2 = 1.0 - ref_coef * ref_coef * (1 - cos_i * cos_i);
        if (cos_t2 < 0) {
            ps = 1;
            pt = 0;
        } else {
            double cos_t = sqrt(cos_t2);
            double Rs = ((ref_index_i * cos_i) - (ref_index_o * cos_t))
                / ((ref_index_i * cos_i) * (ref_index_o * cos_t));
            double Rp = ((ref_index_o * cos_i) - (ref_index_i * cos_t))
                / ((ref_index_o * cos_i) * (ref_index_i * cos_t));
            ps = (Rs * Rs + Rp * Rp)/2;
            pt = 1 - ps;
        }

        return ref_coef;
    }

public:

    // Lambertian diffuse parameters.
    RGB kd;   // Lambertian diffuse coefficient.
    double pd; // Lambertian diffuse probability.

    // Perfect specular reflectance parameters.
    RGB ks;   // Perfect specular reflectante coefficient.
    double ps; // Perfect specular reflectance probability.

    // Perfect refrectation parameters.
    RGB kt;   // Perfect refrectation coefficient.
    double pt; // Perfect refrectation probability.

    // Material emission.
    RGB ke;

    // Material refractance index.
    double ref_index_i;

    Material(RGB kd = RGB(185), RGB ks = RGB(), RGB kt = RGB(), RGB ke = RGB(), double ref_index_i = 0) {

        // Lambertian diffuse parameters.
        this->kd = kd;
        // Perfect specular reflectance parameters.
        this->ks = ks;
        // Perfect refrectation parameters.
        this->kt = kt;
        // Material emission.
        this->ke = ke;

        // Material refraction index.
        this->ref_index_i = ref_index_i;

        // Coefficients correction.
        coefficient_correction();

    }

    Sample scattering(Vector3 n, Vector3 wo = Vector3(), double ref_index_o = 1) {

        // Russian roulette event generator.
        double rr_event = E(e2);

        // Fresnel coefficients evaluation:
        // double ref_coef = fresnel_coefficients(n, wo, ref_index_o, ref_index_i);

        // Lambertian diffuse event:
        if (pd > 0 && rr_event < pd) {
            double lat = acos(sqrt(1 - E(e2))); // SE GENERAN COMO RADIANES
            double azi = 2*M_PI*E(e2);          // LO HE COMPROBADO.
            // Orthonormal basis:
            std::vector<Vector3> b = orthonormal_basis(n);
            // New direction sampling:
            Vector3 dir = Matrix3BaseChange(b[0],b[1], n, Vector3())
                * Vector3(sin(lat)*cos(azi), sin(lat)*sin(azi), cos(lat));
            return Sample(dir, kd/pd, false);
        }
        // Perfect specular reflectance event:
        else if (ps > 0 && rr_event < (pd + ps)) {
            return Sample(wo - ((2*n)*(wo*n)), ks/ps, true);
        }
        // Perfect refrectation event:
        // - https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
        // - https://stackoverflow.com/a/58676386
        else if (pt > 0 && rr_event < (pd + ps + pt)) {

            double ref_coef = ref_index_o/ref_index_i;
            if ((n * wo) >  0) {
                n *= -1;
                ref_coef = ref_index_i/ref_index_o;
            }
            wo = nor(wo);
            double cos_i  = n * wo;
            double cos_t2 = 1.0 - ref_coef * ref_coef * (1 - cos_i * cos_i);
            if (cos_t2 < 0) {
                return Sample(wo - ((2*n)*(wo*n)), kt/pt, true);
            } else {
                return Sample(ref_coef * (wo - n * cos_i) - n * sqrtf(cos_t2), kt/pt, true);
            }
        }
        // Ray death event:
        else {
            return Sample();
        }
    }

};

std::ostream& operator<<(std::ostream& os, const Material& m) {
    return os << "[ kd: " << m.kd << ", ks: " << m.ks << ", kt: " << m.kt << "]";
}

//=================================================================//
// Light
//=================================================================//

//===============================================================//
// Light photon
//===============================================================//
class Photon {
private:
    // ... 
public:
    // Photon position.
    Vector3 pos;
    // Photon flux.
    RGB flux;
    // Photon next direction.
    Vector3 wi;
    Photon (Vector3 pos, RGB flux, Vector3 wi) : pos(pos), flux(flux), wi(wi) {}
};

std::ostream& operator<<(std::ostream& os, const Photon& p) {
    return os << "Photon { " << p.pos << ", " << p.flux << ", " << p.wi << " }";
}

struct PhotonAxisPosition {
    double operator()(const Photon& p, std::size_t i) const {
        return p.pos[i];
    }
};

//===============================================================//
// Point light
//===============================================================//
class Light {
private:
    // ...
public:
    // Point light center.
    Vector3 c;
    // Point light power.
    RGB pow;

    Light(Vector3 c = Vector3(), RGB pow = RGB()) : c(c), pow(pow) {}
};

//=================================================================//
// Objects
//=================================================================//

class Object;

struct Collision {
    std::shared_ptr<Object> obj; // Collisioned object.
    Vector3 normal; // Collisioned object normal.
    Vector3 point;  // Collision point.
    double dist;    // Collision distance.
};

class Object {
private:
    virtual std::ostream& print(std::ostream& os) = 0; // He tenido que cargarme un 
                                                       // poquito el overload para que 
                                                       // se pudiese printear cualquier
                                                       // objeto.
public:

    // Object info.
    Material m;

    Object(Material m = Material()) : m(m) {}

    // Texture things.
    // Texture t; // Object texture.
    // bool has_texture;
    //Object() : m(Material(RGB(185), RGB(0), RGB(0))), has_texture(false) {}
    // Object(Texture t) : t(t), has_texture(true) {}
    // Object(Material m) : m(m), has_texture(false) {}
    virtual ~Object() = default; // Por alguna razón, el compilador de C++ necesita un
                                 // destructor virtual, no hace absolutamente nada, pero
                                 // no toca los huevos.
    
    // Ahora el fr ya no pertenece al objeto sino al material y se calcula a lo
    // ruleta rusa, entonces, aqui ya no hace nada.
    virtual RGB emission() const {
        return m.kd/M_PI;
    }

    virtual double intersects(const Ray& ray) = 0;
    virtual Vector3 normal(Vector3 wo = Vector3(), Vector3 p = Vector3()) = 0;

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
            << "\n  finite plane bounds: " << b
            << "\n  material: " << m
            //<< "\n  texture reference: " << t_ref
            << "\n}";
    }

    //// Texture things.
    //std::vector<Vector3> t_ref; // Texture start reference and orientation.
    //                            //  - 0: reference point.
    //                            //  - 1: texture geometrical height;
    //                            //  - 2: texture geometrical width;
    ////int qi, qj; // Index of the texture quad that has been intersected.
    //double qw, qh; // Width and height of a single quad.
    //double tw, th; // Width and height of the texture.
    //RGB q_color;  // Color of the texture quad that has been intersected.

public:
    // Geometrical things
    double D;                 // Implicit equation A*x+B*y+C*z+D
                             //     (= 0 if the point is in the plane).
    Vector3 n;               // Normal of the plane = (A,B,C).
    std::vector<Vector3> b;  // Vertex bounds of the plane (if finite).

    Plane() {}

    // ==========================
    // Plane constructors
    // ==========================

    // Solid color plane defined by a normal and its distance to the origin.
    Plane(double D, Vector3 n, Material m = Material()) : Object(m), D(D), n(nor(n)) {}

    // Plane defined by a normal and a plane contained point with solid color.
    Plane(Vector3 p, Vector3 n, Material m = Material()) : Object(m) {
        this->n = nor(n);
        this->D = -this->n*p;
    }

    // Texturized plane defined by a normal and the distance of the plane to the origin.
    //  - t is the texture.
    //  - p is the point reference where the texture will start to be drawn.
    //  - tw and th are the texture geometrical width and height.
    //  - r is the rotation of the texture.
    //  - o is the orientation (positive or negative).
    /*
    Plane(double D, Vector3 n, Texture t, Texture_ref r, double tw, double th)
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
    */

    // ==========================
    // Finite plane constructors
    // ==========================
    Plane(double D, std::vector<Vector3> b, Material m = Material())
        : Object(m), D(D), n(nor(crs(b[1]-b[0], b[3]-b[0]))), b(b) {}

    Plane(std::vector<Vector3> b, Material m = Material()) : Object(m), b(b) {
        this->n = nor(crs(b[1]-b[0], b[3]-b[0]));
        this->D = -n*b[0];
    }

    /*
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
    */

    double intersects(const Ray& r) override {

        if (n*r.d == 0.f) return {};     // Si la division es 0 no hay corte.
        double t = -(n*r.p + D)/(n*r.d); // Calcula la distancia desde el
                                         //   centro del rayo hasta el punto de corte.
        
        Vector3 x = (r.p + r.d*t);
        if (b.size()) {
            if( n * crs(b[1]-b[0], x-b[0]) < 0 || // Point p inside edge 1 (v1 to v2).
                n * crs(b[2]-b[1], x-b[1]) < 0 || // Point p inside edge 2 (v2 to v3).
                n * crs(b[3]-b[2], x-b[2]) < 0 || // Point p inside edge 3 (v3 to v4).
                n * crs(b[0]-b[3], x-b[3]) < 0)   // Point p inside edge 4 (v4 to v1).
            {
                return -1;
            }  
        }

        // Texture things.
        /*
        if (has_texture) {
            p -= t_ref[0];
            // Obtaining proportional width and height.
            double hs = t_ref[1] * (p/th);
            double ws = t_ref[2] * (p/tw);

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
        }*/

        // Devolver la distancia al punto de corte.
        return t; 
    }

    Vector3 normal(Vector3 wo = Vector3(), Vector3 p = Vector3()) override {
        return (n * wo <= 0) ? n : -n;
    }

};

//===============================================================//
// Sphere
//===============================================================//

class Sphere : public Object {
private:

    bool solve_quatratic(double a, double b, double c, double& x0, double& x1) {
        double discr = b * b - 4 * a * c;
        if (discr < 0) return false;
        else if (discr == 0) x0 = x1 = - 0.5 * b / a;
        else {
            double q = (b > 0) ?
                -0.5 * (b + sqrt(discr)) :
                -0.5 * (b - sqrt(discr));
            x0 = q / a;
            x1 = c / q;
        }
        if (x0 > x1) std::swap(x0, x1);

        return true;
    }

    std::ostream& print(std::ostream& os) override {
        return os << "SPHERE {"
            << "\n  center: "   << center
            << "\n  radius: "   << r 
            << "\n  material: " << m
            << "\n}";
    }
public:
    Vector3 center;
    double r;

    Sphere(Vector3 center, double r, Material m = Material()) 
        : Object(m), center(center), r(r) {}

    double intersects(const Ray& ray) override {

        Vector3 L = ray.p - center;
        double a = ray.d * ray.d;
        double b = 2 * ray.d * L;
        double c = L * L - r * r;
        double delta = b*b - 4*a*c;

        if (delta < EPSILON_ERROR) return -1;
        double t0 = (-b - sqrt(delta)) / (2*a);
        double t1 = (-b + sqrt(delta)) / (2*a);
        if (t0 <= EPSILON_ERROR && t1 <= EPSILON_ERROR) return -1;
        return t0 > EPSILON_ERROR ? t0 : t1;

    }

    Vector3 normal(Vector3 wo = Vector3(), Vector3 p = Vector3()) override {
        return p-center;
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
            << "\n  material: " << m
            << "\n}";
    }
public:
    std::vector<Vector3> bounds;
    Vector3 center;

    Box() {}
    Box(std::vector<Vector3> bounds, Material m = Material())
        : bounds(bounds), center((bounds[1] - bounds[0])/2) {}

    Box(Vector3 min, Vector3 max, Material m = Material())
        : Object(m), bounds({ min, max }), center( (max-min)/2) {}

    double intersects(const Ray& ray) override {
        return -1;
        /*
        Vector3 inv_d(1.f/ray.d.x, 1.f/ray.d.y, 1.f/ray.d.z);
        bool sign_dir_x = inv_d.x < 0;
        bool sign_dir_y = inv_d.y < 0;
        bool sign_dir_z = inv_d.z < 0;

        double p0 = (bounds[sign_dir_x].x - ray.p.x) * inv_d.x;
        double p1 = (bounds[1 - sign_dir_x].x - ray.p.x) * inv_d.x;
        double b0 = (bounds[sign_dir_y].y - ray.p.y) * inv_d.y;
        double b1 = (bounds[1 - sign_dir_y].y - ray.p.y) * inv_d.y;
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
        } else return { p0 };*/
    }

    Vector3 normal(Vector3 wo = Vector3(), Vector3 p = Vector3()) override { 
        // TODO: corregirlo.
        return Vector3();
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
            << "\n material: " << m
            << "\n}";
    }
public:
    std::vector<Vector3> v; // Triangle vertex.

    Triangle() {}
    Triangle(std::vector<Vector3> v, Material m = Material()) 
        : Plane(v[0], crs(v[1]-v[0], v[2]-v[0]), m), v(v) {}

    double intersects(const Ray& r) override {
        double t;
        if ((t = Plane::intersects(r)) < 0) return t;

        Vector3 x = r.d*t + r.p;
        if (n * crs(v[1]-v[0], x-v[0]) < 0 ||          // Point p inside edge 1 (v1 to v2).
            n * crs(v[2]-v[1], x-v[1]) < 0 ||          // Point p inside edge 2 (v2 to v3).
            n * crs(v[0]-v[2], x-v[2]) < 0) return -1; // Point p inside edge 3 (v3 to v1).
        return t;
    }

};

class Mesh : public Object {
private:

    std::ostream& print(std::ostream& os) override {
        return os << "MESH {"
            << "\n  faces: "        << faces
            << "\n  material: "     << m
            << "\n  box collider: " << collider
            << "\n}";
    }

    // Ply metadata
    // nothing for the moment.
    Box collider;

public:

    Triangle q_dot;
    std::vector<Triangle> faces;

    Mesh(std::vector<Triangle> faces) : collider(Vector3(INFINITY, INFINITY, INFINITY), 
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

    Mesh(std::string ply_file, Material m = Material()) : Object(m), 
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
            f = Triangle({vertex[v1], vertex[v2], vertex[v3]}, m);
        }
    }

    double intersects(const Ray& r) override {
        /*
        if (collider.intersects(r).empty()) return {};
        
        std::vector<double> ts;
        for (auto& f : faces) {
            std::vector<double> tp;
            if ((tp = f.intersects(r)).size() && tp[0] > 0 && !insert(ts, tp[0])) {
                q_dot = f;
            }
        }
        return ts;*/
        return -1;
    }

    // Calcular en que triangulo queda el punto.
    Vector3 normal(Vector3 wo = Vector3(), Vector3 p = Vector3()) override {
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
        f = Triangle(f.v, f.m);
    }
    m.collider.bounds = b;
    return m;
};