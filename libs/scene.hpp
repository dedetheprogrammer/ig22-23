#pragma once
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <mutex>
#include <random>
#include <thread>
#include "geometry.hpp"
#include "objects.hpp"

//===============================================================//
// Definitions
//===============================================================//

using Object_ptr = std::shared_ptr<Object>; 
using Objects    = std::vector<Object_ptr>;
using Lights     = std::vector<Light>;

struct RenderObject {
public:
    Objects objects;
    Lights  lights;

    RenderObject(Objects o, Lights l) : objects(o), lights(l) {}
};

//===============================================================//
// Camera
//===============================================================//

// Intersection offset to avoid Shadow Acne and self intersection:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/ligth-and-shadows
#define SHADOW_BIAS 0.00001f

enum Mode { RAY_TRACING, PATH_TRACING, PHOTON_MAPPING };

class Camera {
private:

    // Returns an antialiasing point.
    std::vector<Vector3> antialiasing_samples(Vector3 ref) {
        if (ppp == 1) {
            return { ref };
        } else {
            std::vector<Vector3> samples(ppp);
            std::uniform_real_distribution<> x( ref.x - p_left.x, ref.x + p_left.x);
            std::uniform_real_distribution<> y( ref.y - p_up.y, ref.y + p_up.y);
            for (auto& s : samples) {
                s = Vector3(x(e2), y(e2), ref.z);
            }
            return samples;
        }
    }

    // Casting a ray.
    /**
     * @brief Casting a ray.
     * 
     * @param s     Set of renderable objects (Geometrical objects and lights).
     * @param r     Incident ray.
     * @param p_obj The object that the previous cast intersected with.
     * @param depth How many bounces do our ray to bounce.
     * @return RGB 
     */
    RGB cast_ray(RenderObject& s, Ray r, Object_ptr p_obj, int depth) {

        if (depth > MAX_DEPTH) return RGB();
        // Light contributios.
        RGB direct_light_contrib, indirect_light_contrib;
        // Distance of the closest intersection point.
        double c_dist = INFINITY;
        // Closest object.
        Object_ptr c_obj;
        // Calculating possible intersections in the Scene.
        for (auto& o : s.objects) {
            auto t = o->intersects(r);
            if (t > 0 && t < c_dist) {
                c_dist = t;
                c_obj  = o;
            }
        }
        if (c_obj == nullptr) return RGB();
        if (c_obj->m.ke  > 0) return c_obj->m.ke;
        
        // Aqui toca meter si intersecta con un objeto que emite luz devuelva
        // la luz o yo que se, hay que mirarlo. También pienso que hay que
        // hacer algo en la direct_light del tipo, hemos recorrido las luces
        // pero tenemos que recorrer los objetos también ya que algunos podrían
        // emitir luz y es luz que también influye.
        //if (c_obj.emissor()) return 

        // Montecarlo sample:
        Vector3 x = r.d * c_dist + r.p, n = nor(c_obj->normal(r.d, x));
        Sample samp = c_obj->m.scattering(n, r.d, 1.0 /*
            ((p_obj == nullptr || p_obj == c_obj || c_dist > BIAS) ? 
                1.0 : 
                p_obj->m.ref_index_i
            )*/
        );

        // Compute direct light: explanation in Material class.
        if (samp.wd_light) {
            for (auto& l : s.lights) {

                // Shadow ray:
                Ray lr(x + n * SHADOW_BIAS, l.c - x);
                double l_mod = lr.d.mod();
                if (l_mod == 0.0f) continue;

                bool collides = false;
                for (auto& o : s.objects) {
                    auto t = o->intersects(lr);
                    Vector3 ln = nor(o->normal(lr.d, lr.d * t + lr.p));
                    if (t > 0 && t < 1 && (ln * lr.d) < 0) {
                        collides = true;
                        break;
                    }
                }
                if (collides) continue;
                
                direct_light_contrib += 
                    (l.pow / (l_mod * l_mod)) *
                    c_obj->emission() *
                    std::abs(n * (lr.d / l_mod));
            }
        }

        // Compute indirect light.
        // Bouncing ray. Now the direction depends on the material properties:
        //  - Uniform cousine sample (Total difussion).
        //  - Reflection.
        //  - Refraction.
        //  - More..?        
        indirect_light_contrib += cast_ray(s, Ray(x, samp.wi), c_obj, depth+1);

        return direct_light_contrib + indirect_light_contrib * samp.fr;
    }

    // Path tracing.
    inline RGB path_tracing(RenderObject& s, Ray r) {
        return cast_ray(s, r, nullptr, 0);
    }

    // Ray tracing.
    inline RGB ray_tracing(RenderObject& s, Ray r) {

        // Distance of the closest intersection point.
        double c_dist = INFINITY;
        // Closest object.
        Object_ptr c_obj;

        // Calculating possible intersections in the Scene.
        for (auto& o : s.objects) {
            auto t = o->intersects(r);
            if (t > 0 && t < c_dist) {
                c_dist = t;
                c_obj  = o;
            }
        }
        return c_obj->m.kd;

    }

public:

    // Camera geometrical base.
    Vector3 center;     // Camera position.
    Vector3 left;       // Camera width projection.
    Vector3 up;         // Camera height projection.
    Vector3 front;      // Camera depth projection.

    // Camera pixel base.
    Vector3 p_left;     // Pixel width projection.
    Vector3 p_up;       // Pixel height projection.
    Vector3 p_ref;      // Pixel dot reference.
    int ppp;            // Point/Paths per pixel.

    // Camera view.
    Image view;         // Camera image view.

    
    Mode mode;          // Render mode.
    int MAX_DEPTH;      // Path tracing max depth.

    // Other.
    progress_bar bar;   // Render progress bar.

    Camera(std::vector<Vector3> base, int width, int height, int ppp = 1,
        Mode mode = PATH_TRACING, int MAX_DEPTH = 1) {

        // Camera geometrical base.
        center = base[0];
        left   = base[1];
        up     = base[2];
        front  = base[3];

        // Camera pixel base.
        p_left = left/width;
        p_up   = up/height;
        p_ref  = center + (up - p_up) + (left - p_left) + front;
        this->ppp  = ppp;

        this->mode      = mode;
        this->MAX_DEPTH = MAX_DEPTH;

        // Camera image view.
        view = Image(width, height, "Camera View");

        // Render progress.
        bar  = progress_bar(80, width * height, STYLE1, 100);

    }

    void render(RenderObject s) {
        auto ref = p_ref;
        for (auto& h : view.pixels) {
            for (auto& w : h) {
                auto p = now();
                for (auto &as : antialiasing_samples(ref)) {
                    if (mode == RAY_TRACING) {
                        w += ray_tracing(s, Ray(center, as - center)); 
                    } else if (mode == PATH_TRACING) {
                        w += path_tracing(s, Ray(center, as - center));
                    }
                    
                }
                w /= ppp;
                if (w.R > view.maxval) view.maxval = view.memval = w.R;
                if (w.G > view.maxval) view.maxval = view.memval = w.G;
                if (w.B > view.maxval) view.maxval = view.memval = w.B;
                ref -= 2*p_left;
                bar.update(std::cout, now()-p);
            }
            ref += 2*left - 2*p_up; 
        }
        flush_stream(std::cout);
    
    }

    void export_render(std::string render_name = "./new_scene.ppm", int res = 10e8) {
        view.colres = view.memres = (view.maxval * res);
        export_image(view, render_name);
    }

};

//===============================================================//
// Scene: a set of cameras, objects and lights.
//===============================================================//

class Scene {
private:
    // ...
public:

    // Scene cameras.
    std::vector<Camera> cameras; 
    
    // Scene objects.
    Objects objects;

    // Scene lights.        
    Lights lights;

    Scene(std::vector<Camera> cameras, Objects objects, Lights lights)
        : cameras(cameras), objects(objects), lights(lights) {}

    void render(int i = 0) {
        cameras[i].render(RenderObject(objects, lights));
    }

    void export_render(int i = 0, std::string render_name = "./new_scene.ppm", int res = 10e8) {
        cameras[i].export_render(render_name, res);
    }

};