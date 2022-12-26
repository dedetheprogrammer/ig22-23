#pragma once
#include <algorithm>
#include <random>
#include "geometry.hpp"
#include "image.hpp"
#include "utils.hpp"

//===============================================================//
// Properties
//===============================================================//
class Box_collider {
private:
    //...
public:
    std::vector<Vector3> b;
    Vector3 center;

    Box_collider() : b({Vector3(INFINITY, INFINITY, INFINITY), Vector3(-INFINITY, -INFINITY,-INFINITY)}) {}
    Box_collider(Vector3 min, Vector3 max) : b({ min, max }), center((max-min)/2) {}

    bool intersects(const Ray& r) {

        Vector3 inv_d(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z);
        bool sign_dir_x = inv_d.x < 0;
        bool sign_dir_y = inv_d.y < 0;
        bool sign_dir_z = inv_d.z < 0;

        double p0 = (b[sign_dir_x].x - r.p.x) * inv_d.x;
        double p1 = (b[1 - sign_dir_x].x - r.p.x) * inv_d.x;
        double b0 = (b[sign_dir_y].y - r.p.y) * inv_d.y;
        double b1 = (b[1 - sign_dir_y].y - r.p.y) * inv_d.y;
        if ((p0 > b1) || (b0 > p1)) return false;
        if (b0 > p0) p0 = b0;
        if (b1 < p1) p1 = b1;

        b0 = (b[sign_dir_z].z - r.p.z) * inv_d.z;
        b1 = (b[1 - sign_dir_z].z - r.p.z) * inv_d.z;
        if ((p0 > b1) || (b0 > p1)) return false;
        if (b0 > p0) p0 = b0;
        if (b1 < p1) p1 = b1;

        return (p0 < 0 && p1 < 0) ? false : true;

    }

};

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
    Vector3 normal;              // Collisioned object normal.
    Vector3 point;               // Collision point.
    double dist;                 // Collision distance.

    Collision()
        : normal(Vector3()), point(Vector3()), dist(INFINITY) {}
    Collision(double d)
        : normal(Vector3()), point(Vector3()), dist(d) {}
    Collision(Vector3 n, Vector3 p, double d)
        : normal(n), point(p), dist(d) {}

};

//--------------------------------
class Object {
private:
    virtual std::ostream& print(std::ostream& os) = 0; // He tenido que cargarme un 
                                                       // poquito el overload para que 
                                                       // se pudiese printear cualquier
                                                       // objeto.
public:

    // Object info.
    Material m; // Object material.
    // Box_collider collider; // Object box collider.

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

    virtual Collision intersects(const Ray& ray) = 0;
    
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

    // Default constructor.
    Plane() {}

    // ==========================
    // Plane constructors
    // ==========================

    // Solid color plane defined by a normal and its distance to the origin.
    Plane(double D, Vector3 n, Material m = Material()) : Object(m), D(D), n(nor(n)) {}

    // Plane defined by a normal and a plane contained point with solid color.
    Plane(Vector3 p, Vector3 n, Material m = Material()) : Object(m) {
        this->n = nor(n);
        this->D = -(this->n)*p;
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
    //Plane(double D, std::vector<Vector3> b, Material m = Material())
    //    : Object(m), D(D), n(nor(crs(b[1]-b[0], b[3]-b[0]))), b(b) {}

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

    Collision intersects(const Ray& r) override {

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
                return Collision(-1);
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
        return Collision(nor((n * r.d <= 0) ? n : -n), x, t);
    }

};

//===============================================================//
// Sphere
//===============================================================//

class Sphere : public Object {
private:

    std::ostream& print(std::ostream& os) override {
        return os << "SPHERE {"
            << "\n  center: "   << center
            << "\n  radius: "   << radius
            << "\n  material: " << m
            << "\n}";
    }

public:

    Vector3 center; // Sphere center.
    double radius;  // Sphere radius.

    Sphere(Vector3 center, double radius, Material m = Material()) 
        : Object(m), center(center), radius(radius) {}

    Collision intersects(const Ray& r) override {

        Vector3 L = r.p - center;
        double a = r.d * r.d;
        double b = 2 * r.d * L;
        double c = L * L - radius * radius;
        double delta = b*b - 4*a*c;

        if (delta < EPSILON_ERROR) return Collision(-1);
        double t0 = (-b - sqrt(delta)) / (2*a);
        double t1 = (-b + sqrt(delta)) / (2*a);
        if (t0 <= EPSILON_ERROR && t1 <= EPSILON_ERROR) return Collision(-1);
        if (t0 > EPSILON_ERROR) {
            Vector3 x = r.d * t0 + r.p;
            return Collision(nor(x-center), x, t0);
        } else {
            Vector3 x = r.d * t1 + r.p;
            return Collision(nor(x-center), x, t1);
        }

    }

};

//===============================================================//
// Cube
//===============================================================//

class Cube : public Object {
private:
    std::ostream& print(std::ostream& os) override {
        return os << "CUBE {"
            //<< "\n  bounds: "   << bounds
            << "\n  material: " << m
            << "\n}";
    }
public:
    // Cube vertex representation:
    //        b[5]             max (b[6])
    //            +-----------+
    //           /.          /|
    //     b[1] +-----------+<--- b[2]
    //          | .         | |
    //     b[4]-| ..FRONT...|.+ b[7]
    //          |.          |/
    //          +-----------+
    //     min (b[0])         b[3]
    std::vector<Vector3> v;   // Cube vertex.
    std::vector<Plane> faces; // Cube faces.

    // Default constructor.
    Cube() {}
    // Cube constructor with the minimum and maximum bounds.
    Cube(Vector3 min, Vector3 max, Material m = Material()) : Object(m) {

        // Vertex order:
        v.push_back(min);                          // v[0]: -x, -y, -z.
        v.push_back(Vector3(min.x, max.y, min.z)); // v[1]: -x, +y, -z.
        v.push_back(Vector3(max.x, max.y, min.z)); // v[2]: +x, +y, -z.
        v.push_back(Vector3(max.x, min.y, min.z)); // v[3]: +x, -y, -z.
        v.push_back(Vector3(min.x, min.y, max.z)); // v[4]: -x, -y, +z.
        v.push_back(Vector3(min.x, max.y, max.z)); // v[5]: -x, +y, +z.
        v.push_back(max);                          // v[6]: +x, +y, +z.
        v.push_back(Vector3(max.x, min.y, max.z)); // v[7]: +x, -y, +z.

        // Faces order:
        faces.push_back(Plane({v[0],v[1],v[2],v[3]})); // Front face.
        faces.push_back(Plane({v[3],v[2],v[6],v[7]})); // Right face.
        faces.push_back(Plane({v[7],v[6],v[5],v[4]})); // Back face.
        faces.push_back(Plane({v[4],v[5],v[1],v[0]})); // Left face.
        faces.push_back(Plane({v[1],v[5],v[6],v[2]})); // Top face.
        faces.push_back(Plane({v[0],v[4],v[7],v[4]})); // Bottom face.

    }
    // Cube constructor with each vertice.
    Cube(std::vector<Vector3> v, Material m = Material()) : Object(m) {

        // Cube vertices:
        this->v = v;
        // Cube faces:
        faces.push_back(Plane({v[0],v[1],v[2],v[3]})); // Front face.
        faces.push_back(Plane({v[3],v[2],v[6],v[7]})); // Right face.
        faces.push_back(Plane({v[7],v[6],v[5],v[4]})); // Back face.
        faces.push_back(Plane({v[4],v[5],v[1],v[0]})); // Left face.
        faces.push_back(Plane({v[1],v[5],v[6],v[2]})); // Top face.
        faces.push_back(Plane({v[0],v[4],v[7],v[4]})); // Bottom face.

    }
    // Cube constructor with each plane.
    Cube(std::vector<Plane> faces, Material m = Material())
        : Object(m), faces(faces) {}

    Collision intersects(const Ray& ray) override {

        Collision c;
        for (auto& f : faces) {
            auto t = f.intersects(ray);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
            }
        }
        return (c.dist != INFINITY) ? c : Collision(-1);
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

    Collision intersects(const Ray& r) override {

        Collision c;
        if ((c = Plane::intersects(r)).dist < 0) return c;

        Vector3 x = r.d*c.dist + r.p;
        if (n * crs(v[1]-v[0], x-v[0]) < 0 || // Point p inside edge 1 (v1 to v2).
            n * crs(v[2]-v[1], x-v[1]) < 0 || // Point p inside edge 2 (v2 to v3).
            n * crs(v[0]-v[2], x-v[2]) < 0)   // Point p inside edge 3 (v3 to v1).
        {
            return Collision(-1);
        }
        return Collision(nor(n), x, c.dist);
    }

};

class Mesh : public Object {
private:

    std::ostream& print(std::ostream& os) override {
        return os << "MESH {"
            << "\n  faces: "        << faces
            << "\n  material: "     << m
            //<< "\n  box collider: " << c
            << "\n}";
    }

    // Ply metadata
    // nothing for the moment.

public:

    Box_collider collider;       // Mesh box collider.
    std::vector<Triangle> faces; // Mesh faces.

    Mesh(std::vector<Triangle> faces, Material m = Material()) 
        : Object(m), faces(faces)
    {   
        for (auto& f : faces) {
            for (auto& v : f.v) {
                collider.b[0] = min(collider.b[0], v);
                collider.b[1] = max(collider.b[1], v);
            }
        }
    }

    Mesh(std::string ply_file, Material m = Material()) : Object(m) {

        // Reading the PLY file:
        std::string s("");
        std::ifstream in(ply_file);
        assert(in.is_open() && "file not found, check it out.");
        std::vector<Vector3> vertex;
        // Reading the PLY header:
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

        // Reading the PLY vertex list:
        for (auto& v : vertex) {
            in >> v;
            // Defining the mesh collider.
            collider.b[0] = min(collider.b[0], v);
            collider.b[1] = max(collider.b[1], v);
        }

        // Reading the PLY faces list.
        int n, v1, v2, v3;
        for (auto& f : faces) {
            in >> n >> v1 >> v2 >> v3;
            f = Triangle({vertex[v1], vertex[v2], vertex[v3]});
        }
    }

    Collision intersects(const Ray& r) override {
        // If the ray doesn't intersects the mesh collider, then it won't
        // intersect anything, just return no collision.
        if (!collider.intersects(r)) return Collision(-1);

        // Mesh internal collision.
        Collision c; 
        // Calculating possible intersections with the mesh faces:
        for (auto& f : faces) {
            auto t = f.intersects(r);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
            }
        }
        return (c.dist != INFINITY) ? c : Collision(-1);
    }

};

Mesh operator*(Mesh m, Matrix3 transform) {

    m.collider = Box_collider();
    for (auto& f : m.faces) {
        for (auto& v : f.v) {
            v.h = 1;
            v = transform * v;
            m.collider.b[0] = min(m.collider.b[0], v);
            m.collider.b[1] = max(m.collider.b[1], v);
        }
        f = Triangle(f.v, f.m);
    }
    return m;

};