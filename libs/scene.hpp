#pragma once
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <mutex>
#include <random>
#include <thread>
#include "geometry.hpp"
#include "kd_tree.h"
#include "objects.hpp"

//===============================================================//
// Definitions
//===============================================================//

using Photon_map = nn::KDTree<Photon, 3, PhotonAxisPosition>;
using Object_ptr = std::shared_ptr<Object>; 
using Objects    = std::vector<Object_ptr>;
using Lights     = std::vector<Light>;

/**
 * @brief Camera render mode.
 *  1. RAY_TRACING: just rasterization, without any light.
 *  2. PATH_TRACING: ray path + direct & indirect light.
 *  3. PHOTON_MAPPING: photon path + direct & indirect light.
 * 
 */
enum Render_mode { RAY_TRACING, PATH_TRACING, PHOTON_MAPPING };

/**
 * @brief Render parameters, created for some coding magic.
 * 
 */
struct Render_params {
public:
    // Renderable Objects.
    Objects objects;

    // Renderable lights.
    Lights lights;

    // Renderable photons.
    Photon_map pmap;

    // Incident ray (wo).
    Ray ray;

    // Object that was previously intersected. Null if no object was.
    Object_ptr p_obj;

    // How many bounces has our ray bounced already.
    int depth;

    Render_params(Objects o, Lights l, Photon_map pmap)
        : objects(o), lights(l), pmap(pmap), p_obj(nullptr), depth(0) {}
    Render_params(Objects o, Lights l, Ray r, Object_ptr p_obj, int depth)
        : objects(o), lights(l), ray(r), p_obj(p_obj), depth(depth) {}
    Render_params(Objects o, Lights l, Photon_map pmap, Ray r, Object_ptr p_obj)
        : objects(o), lights(l), pmap(pmap), ray(r), p_obj(p_obj) {}
};

std::ostream& operator<<(std::ostream& os, const Render_params& r) {
    return os << r.ray << ", " << r.depth << std::endl;
}

//===============================================================//
// Camera
//===============================================================//

// Intersection offset to avoid Shadow Acne and self intersection:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/ligth-and-shadows
#define SHADOW_BIAS 0.00001

// Camera settings
struct Camera_settings {

    // General settings:
    // -- Camera geometrical base:
    Vector3 center;         // Camera position.
    Vector3 left;           // Camera width projection.
    Vector3 up;             // Camera height projection.
    Vector3 front;          // Camera depth projection.
    // -- Camera pixel base:
    int ppp;                // Point/paths per pixel.
    // -- Camera renderization base:
    Render_mode m;          // Camera render mode.

    // Path tracing settings:
    int RAY_MAX_DEPTH;      // Path tracing max depth.

    // Photon mapping settings:
    unsigned long nphotons; // Max number of photons to look for.
    double radius;          // Photon search radius.

    // Default constructor.
    Camera_settings() {}
    // Ray tracing camera settings.
    Camera_settings(std::vector<Vector3> base, int ppp) {

        // Camera geometrical base:
        center = base[0];
        left   = base[1];
        up     = base[2];
        front  = base[3];

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;

        // Ray tracing settings:
        m = RAY_TRACING;
    }
    // Path tracing camera settings.
    Camera_settings(std::vector<Vector3> base, int ppp, int RAY_MAX_DEPTH) {
        
        // Camera geometrical base:
        center = base[0];
        left   = base[1];
        up     = base[2];
        front  = base[3];

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;

        // Path tracing settings:
        m = PATH_TRACING;
        this->RAY_MAX_DEPTH = RAY_MAX_DEPTH;

    }
    // Photon mapping camera settings.
    Camera_settings(std::vector<Vector3> base, int ppp, unsigned long nphotons, double radius) {

        // Camera geometrical base:
        center = base[0];
        left   = base[1];
        up     = base[2];
        front  = base[3];

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;

        // Photon mapping settings:
        m = PHOTON_MAPPING;
        this->nphotons = nphotons;
        this->radius   = radius;

    }
    // Free camera settings.
    Camera_settings(std::vector<Vector3> base, int ppp, Render_mode m, 
        int RAY_MAX_DEPTH, int nphotons, double radius)
    {
        // Camera geometrical base:
        center = base[0];
        left   = base[1];
        up     = base[2];
        front  = base[3];

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;

        // Render settings:
        this->m = m;
        this->RAY_MAX_DEPTH = RAY_MAX_DEPTH;
        this->nphotons = nphotons;
        this->radius   = radius;

        

    }

};

std::ostream& operator<<(std::ostream& os, const Camera_settings& s) {
    return os << "Camera coords: " << s.center << ", " << s.left << ", " << s.up << ", " << s.front << "\n"
        << "PPP: " << s.ppp << "\n"
        << "RAY_MAX_DEPTH: " << s.RAY_MAX_DEPTH << "\n" 
        << "NPHOTONS: " << s.nphotons << "\n" 
        << "RADIUS SEARCH: " << s.radius << "\n";
}   

class Camera {
private:
public:

    /**
     * @brief Generates an antialiasing pixel point.
     * 
     * @param ref Pixel center reference.
     * @return Vector3 
     */
    inline Vector3 antialiasing_sample(Vector3 ref) {
        std::uniform_real_distribution<> x( ref.x - p_left.x, ref.x + p_left.x);
        std::uniform_real_distribution<> y( ref.y - p_up.y, ref.y + p_up.y);
        return Vector3(x(e2), y(e2), ref.z);
    }

    /**
     * @brief Ray tracing.
     * 
     * @param s Set of parameters (objects, lights, etc) needed by the render.
     * @return RGB 
     */
    inline RGB ray_tracing(Render_params s) {

        // Closest collision.
        Collision c;

        // Calculating possible intersections in the Scene.
        for (auto& o : s.objects) {
            auto t = o->intersects(s.ray);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
                c.obj = o;
            }
        }
        return c.obj->m.kd;

    }

    /** 
     * @brief Path tracing.
     * 
     * @param s Set of parameters (objects, lights, etc) needed by the render.
     * @return RGB
     */
    inline RGB path_tracing(Render_params s) {

        if (s.depth > ini.RAY_MAX_DEPTH) return RGB();

        // Light contributions.
        RGB direct_light_contrib, indirect_light_contrib;
        // Closest collision.
        Collision c;
        // Calculating possible intersection with the scene:
        for (auto& o : s.objects) {
            auto t = o->intersects(s.ray);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
                c.obj = o;
            }
        }
        // If the ray doesn't intersect with any object, just return black color.
        if (c.obj == nullptr) return RGB();
        // If the ray intersects with an object with emission, return its light 
        // emission coefficient.
        if (c.obj->m.ke  > 0) return c.obj->m.ke;

        // Montecarlo sample:
        Sample samp = c.obj->m.scattering(c.normal, s.ray.d, 1.0 /*
            ((s.p_obj == nullptr || s.p_obj == c.obj || c.dist > SHADOW_BIAS) ? 
                1.0 : 
                s.p_obj->m.ref_index_i
            )*/
        );

        // Direct light computation. Delta explanation in Material class
        // (objects.hpp):
        if (!samp.is_delta) {
            for (auto& l : s.lights) {
                // Shadow ray. See shadow bias explanation in the variable
                // declaration:
                Ray lr(c.point + c.normal * SHADOW_BIAS, l.c - c.point);
                double l_mod = lr.d.mod();
                if (l_mod == 0.0) continue;

                bool collides = false;
                for (auto& o : s.objects) {
                    auto t = o->intersects(lr);
                    if (t.dist > 0 && t.dist < 1 && (t.normal * lr.d) < 0) {
                        collides = true;
                        break;
                    }
                }
                if (collides) continue;

                // Render equation:
                direct_light_contrib +=
                    (l.pow / (l_mod * l_mod)) *
                    c.obj->emission() *
                    std::abs(c.normal * (lr.d / l_mod));
            }
        }

        // Indirect light computation.
        // Bouncing ray. Now the direction depends on the material properties:
        //  - Uniform cousine sample (Total difussion).
        //  - Reflection.
        //  - Refraction.
        //  - More..?
        indirect_light_contrib += path_tracing(
            Render_params(s.objects, s.lights, Ray(c.point, samp.wi), c.obj, s.depth+1)
        );

        return direct_light_contrib + indirect_light_contrib * samp.fr;
    }

    /**
     * @brief Photon mapping 
     * 
     * @param s Set of parameters (objects, lights, etc) needed by the render.
     * @return RGB
     */
    RGB photon_mapping(Render_params s) {
        return RGB();

        // Light contribution.
        RGB photon_light_contrib;
        // Closest collision.
        Collision c;
        // Calculating possible intersections in the Scene.
        for (auto& o : s.objects) {
            auto t = o->intersects(s.ray);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
                c.obj = o;
            }
        }
        // If the ray doesn't intersect with any object, just return black color.
        if (c.obj == nullptr) return RGB();
        // If the ray intersects with an object with emission, return its light 
        // emission coefficient.
        if (c.obj->m.ke  > 0) return c.obj->m.ke;

        // Object scattering and sampling:
        Sample samp = c.obj->m.scattering(c.normal, s.ray.d, 1.0 /*/
            ((s.p_obj == nullptr || s.p_obj == c_obj || c_dist > SHADOW_BIAS) ? 
                1.0 : 
                s.p_obj->m.ref_index_i
            )*/
        );
        if (samp.fr == 0) return RGB();
        // If the object material is a delta material, we have to follow the ray 
        // path until it intersects with a diffuse object or dies.
        if (samp.is_delta) {
            return photon_mapping(
                Render_params(s.objects, s.lights, s.pmap, Ray(c.point, samp.wi), c.obj)
            ) * samp.fr;
        }

        // Maximum distance to look for photons:
        double real_rad = 0;

        // Nearest photons to the collision point:
        auto neighbors = s.pmap.nearest_neighbors(c.point, ini.nphotons, ini.radius);

        // Compute real radius:
        for (auto& p : neighbors) {
            double dist = (p->pos - c.point).mod();
            if (dist > real_rad) real_rad = dist;
        }

        // Compute light contribution:
        for (auto& p : neighbors) {
            // Radial basis function kernel
            // https://towardsdatascience.com/gaussian-process-kernels-96bafb4dd63e
            double dist = (p->pos - c.point).mod();
            double RBFK = exp(-(dist*dist)/(2*real_rad*real_rad));

            photon_light_contrib +=
                // (c_obj->emission() * (p->flux / (M_PI * rad * rad)));
                (c.obj->emission() * p->flux * RBFK);
        }

        return photon_light_contrib;
    }

    // Camera pixel base.
    Vector3 p_left;     // Pixel width projection.
    Vector3 p_up;       // Pixel height projection.
    Vector3 p_ref;      // Pixel dot reference.

public:

    Camera_settings ini; // Camera settings.
    Image view;          // Camera image view.
    Progress_bar bar;    // Render progress bar.

    Camera(Camera_settings ini, int width, int height) {

        // Camera settings.
        this->ini = ini;

        // Camera pixel base.
        p_left = ini.left/width;
        p_up   = ini.up/height;
        p_ref  = ini.center + (ini.up - p_up) + (ini.left - p_left) + ini.front;

        // Camera image view.
        view = Image(width, height, "Camera view");

        // Render progress.
        bar  = Progress_bar("RENDERING SCENE", 80, width * height, STYLE1, 100);

    }

    void render(Render_params s) {

        // Render function to execute.
        std::function<RGB(Camera*, Render_params)> render_function;
        // Evaluating which render function has to be executed.
        if (ini.m == RAY_TRACING) {
            render_function = ray_tracing;
        } else if (ini.m == PATH_TRACING) {
            render_function = path_tracing;
        } else if (ini.m == PHOTON_MAPPING) {
            if (s.pmap.empty()) {
                std::cerr << "warning: your scene doesn't have any preloaded photon map.\n";
                return;
            }
            render_function = photon_mapping;
        }

        // First pixel center as reference.
        auto ref = p_ref;
        // Iterating over camera view pixels:

        for (auto& h : view.pixels) {
            for (auto& w : h) {
                // For ETA purposes.
                auto t = now();
                // First point will be the pixel center.
                s.ray = Ray(ini.center, ref - ini.center);
                w += render_function(this, s);
                // If more points needed, these will be randomly calculated:
                for (int a = 1; a < ini.ppp; a++) {
                    s.ray = Ray(ini.center, antialiasing_sample(ref) - ini.center);
                    w += render_function(this, s);
                }
                // Dividing total pixel color with the ppp.
                w /= ini.ppp;
                // Estimating new memory color resolution value.
                if (w.R > view.maxval) view.maxval = view.memval = w.R;
                if (w.G > view.maxval) view.maxval = view.memval = w.G;
                if (w.B > view.maxval) view.maxval = view.memval = w.B;
                // Move pixel center reference horizontally.
                ref -= 2*p_left;
                // Updating bar status.
                bar.update(std::cout, now()-t);
            }
            ref += 2*ini.left - 2*p_up;
        }
        flush_stream(std::cout);
    }

    void export_render(std::string render_name = "./new_scene.ppm", int res = 10e8) {
        view.colres = view.memres = (view.maxval * res);
        export_image(view, render_name);
    }

};

std::ostream& operator<<(std::ostream& os, const Camera& c) {
    return os << c.ini
        << "Pixel ref: " << c.p_ref  << "\n"
        << "Pixel lft: " << c.p_left << "\n"
        << "Pixel up:  " << c.p_up   << "\n";
}

//===============================================================//
// Scene: a set of cameras, objects and lights.
//===============================================================//
class Scene {
private:

    Progress_bar bar;

    inline void generate_photon_map() {

        RGB total_pow;
        for (auto& l : lights) {
            total_pow += l.pow;
        }

        // Scene photons.
        std::vector<Photon> photons;
        for (auto& l : lights) {

            // Light proportional photons.
            int S = INITIAL_PHOTONS * l.pow.rad() / total_pow.rad();

            // Casting photons
            for (int p = 0; p < S; p++) {
                // Uniform spherical sampling : Angular sampling.
                double lat = acos(2*E(e2) - 1); //acos(C(e2));
                double azi = 2*M_PI*E(e2);
                Ray r(l.c, Vector3(sin(lat)*cos(azi), sin(lat)*sin(azi), cos(lat)));
                RGB flux = (4*M_PI*l.pow)/S;

                auto t = now();
                for (int depth = 0; depth < MAX_PHOTON_DEPTH; depth++) {
                    // If the photon flux is near to 0, just stop iterating.
                    if (flux.rad() < EPSILON_ERROR) break;
                    // Closest collision.
                    Collision c;
                    // Calculating possible intersections with the scene:
                    for (auto& o : objects) {
                        auto t = o->intersects(r);
                        if (t.dist > 0 && t.dist < c.dist) {
                            c = t;
                            c.obj = o;
                        }
                    }
                    // If the photon doesn't intersect with any object of the scene
                    // just stop iterating.
                    if (c.obj == nullptr) break;
                    // If the photon intersects with an object with light emission, 
                    // stop iterating.
                    if (c.obj->m.ke  > 0) break;

                    // Colision point.
                    Sample samp = c.obj->m.scattering(c.normal, r.d);

                    // Changing to the new direction.
                    r = Ray(c.point, samp.wi);

                    // Adding a new photon. If is delta material, doesn't make 
                    // sense to store the photon:
                    if (!samp.is_delta) {
                        photons.push_back(Photon(c.point, flux, samp.wi));
                    }

                    // New flux resulting of the bounce.
                    flux *= (samp.fr);
                } 
                bar.update(std::cout, now() - t);
            }   
        }
    }

public:

    // Scene cameras.
    std::vector<Camera> cameras; 
    
    // Scene objects.
    Objects objects;

    // Scene lights.        
    Lights lights;

    // Scene photon map
    Photon_map pmap;
    // How many photons to shoot in the first iteration between every light source.
    int INITIAL_PHOTONS;
    // Max number of total photons.
    int MAX_PHOTONS;
    // Max bounce photon depth.
    int MAX_PHOTON_DEPTH;

    Scene(std::vector<Camera> cameras, Objects objects)
        : cameras(cameras), objects(objects) {}
    Scene(std::vector<Camera> cameras, Objects objects, Lights lights)
        : cameras(cameras), objects(objects), lights(lights) {}
    Scene(std::vector<Camera> cameras, Objects objects, Lights lights,
        int INITIAL_PHOTONS/*, int MAX_PHOTONS*/, int MAX_PHOTON_DEPTH)
        : cameras(cameras), objects(objects), lights(lights)
    {
        this->INITIAL_PHOTONS = INITIAL_PHOTONS;
        //this->MAX_PHOTONS = MAX_PHOTONS;
        this->MAX_PHOTON_DEPTH = MAX_PHOTON_DEPTH;
        bar = Progress_bar("CREATING PHOTON MAP", 60, INITIAL_PHOTONS, STYLE1, 50);
        generate_photon_map();
    }

    void render(int i = 0) {
        cameras[i].render(Render_params(objects, lights, pmap));
    }

    void export_render(int i = 0, std::string render_name = "./new_scene.ppm", int res = 10e8) {
        cameras[i].export_render(render_name, res);
    }

};