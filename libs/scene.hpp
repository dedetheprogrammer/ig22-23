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
// Scene: a set of objects and lights.
//===============================================================//

using Object_ptr = std::shared_ptr<Object>; 
using Light_ptr  = std::shared_ptr<Light>;
using Objects    = std::vector<std::shared_ptr<Object>>;
using Lights     = std::vector<std::shared_ptr<Light>>;

class Scene {
private:
    // ...
public:
    // Scene cameras.
    //std::vector<Camera> cameras; 
    
    // Scene objects.
    Objects objects;

    // Scene lights.        
    Lights  lights;

    Scene() {}
    Scene(Objects objects, Lights lights) : objects(objects), lights(lights) {}

};

//===============================================================//
// Pixel: the elemental of the grid.
//===============================================================//

#define BIAS 0.001f
#define MAX_DEPTH 1

class Pixel {
private:

    RGB cast_ray(Scene& s, Ray r, int depth) {

        if (depth > MAX_DEPTH) return RGB();

        RGB direct_light_contrib, indirect_light_contrib;
        // Distance of the closest intersection point.
        float c_dist = INFINITY;
        // Closest object.
        Object_ptr c_obj;
        // Calculating possible intersections in the Scene.
        for (auto& o : s.objects) {
            auto vt = o->intersects(r);
            if (vt.size() > 0 && vt[0] > 0 && (vt[0] < c_dist)) {
                c_dist = vt[0];
                c_obj  = o;
            }
        }
        if (c_obj == nullptr) return RGB();

        // Compute direct light.
        // Colision point.
        Vector3 x = r.d * c_dist + r.p, n = nor(c_obj->normal(x, r.d));
        for (auto& l : s.lights) {

            Ray lr(x, l->c - x);
            float l_mod = lr.d.mod();
            if (std::abs(l_mod) < EPSILON_ERROR) continue;

            bool collides = false;
            for (auto& o : s.objects) {
                for (auto& t : o->intersects(lr)) {
                    if (t > BIAS && t < 1) {
                        collides = true;
                        break;
                    }
                }
            }
            if (collides) continue;

            direct_light_contrib += 
                (l->pow / (l_mod * l_mod)) *
                /*c_obj->fr() **/
                std::abs(n * (lr.d / l_mod));
        }

        // Compute indirect light.
        // Bouncing ray. Now the direction depends on the material properties:
        //  - Uniform cousine sample.
        //  - Reflection.
        //  - Refraction.
        Vector3 new_dir = c_obj->m.scattering(n); 
        indirect_light_contrib += cast_ray(s, Ray(x, new_dir), depth+1);

        return (direct_light_contrib + indirect_light_contrib) * c_obj->fr(Vector3(), Vector3(), Vector3()) * M_PI;
    }

public:
    // Pixel center.
    Vector3 c;

    // Pixel points.              
    std::vector<Vector3> dots;

    // Pixel generated color.
    RGB color;
    
    Pixel() {}
    Pixel(Vector3 c, std::vector<Vector3> dots) : c(c), dots(dots), color(RGB()) {}

    inline void ray_tracing(Scene s, Vector3 cc) {
        for (auto& d : dots) {
            // Distance of the closest intersection point.
            float c_dist = INFINITY;

            // Closest object.
            Object_ptr c_obj;

            // Calculating possible intersections in the Scene.
            for (auto& o : s.objects) {
                auto vt = o->intersects(Ray(cc, d - cc));
                if (vt.size() > 0 && vt[0] > 0 && vt[0] < c_dist) {
                    c_dist = vt[0];
                    c_obj = o;
                }
            }
            color += c_obj->m.kd;
        }
        color /= dots.size();
    }

    inline void path_tracing(Scene s, Vector3 cc) {
        for (auto& d : dots) {
            color += cast_ray(s, Ray(cc, d - cc), 0);
        }
        color /= dots.size();
    }

};

//===============================================================//
// Camera
//===============================================================//

using CameraGridRow = std::vector<Pixel>;
using CameraGrid    = std::vector<CameraGridRow>;
enum Render { RAY_TRACING, PATH_TRACING };

class Camera {
private:
    // ...
public:
    Vector3 c, l, u, f; // Center. Left + Up + Front dimensions.
    Vector3 pl, pu;     // Pixel Left + Up dimensions.
                        // - u = center->up border; size = down border->up border = 2u
                        // - l = center->left; size = right border->left border = 2l
    int w, h, ppp;      // Height, width and points per pixel.
    CameraGrid grid;    // Camera grid.
    progress_bar bar;   // Progress bar.

    Camera(Vector3 c, Vector3 l, Vector3 u, Vector3 f, int w, int h, int ppp = 1) 
        : c(c), l(l), u(u), f(f), pl(2*l/w), pu(2*u/h), w(w), h(h), ppp(ppp)
    {
        grid = CameraGrid(h, CameraGridRow(w));
        bar  = progress_bar(70, w * h, STYLE1, 250);

        Vector3 pc = c + u-pu/2 + l-pl/2 + f; // Pixel center.
        // Randomice.
        for (auto& hi : grid) {
            for (auto& wi : hi) {
                std::vector<Vector3> dots(ppp);
                if (ppp == 1) dots[0] = pc;
                else {
                    std::uniform_real_distribution<> x( pc.x - pl.x/2, pc.x + pl.x/2);
                    std::uniform_real_distribution<> y( pc.y - pu.y/2, pc.y + pu.y/2);
                    for (auto& d : dots) {
                        d = Vector3(x(e2), y(e2), pc.z);
                    }
                }
                wi = Pixel(pc, dots);
                pc -= pl;
            }
            pc += (2*l - pu);
        }
    }

    void render(Scene scene, Render type) {
        for (auto& h: grid) {
            for (auto& w : h) {
                auto p = now();
                w.path_tracing(scene, c);
                bar.update(std::cout, now()-p);
            }
        }
        flush_stream(std::cout);
    }

    void export_render(std::string render_name = "./new_scene.ppm", int res = 10e8) {
        Image r(w, h);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                auto pixel = grid[i][j].color;
                r.pixels[i][j] = pixel;
                if (pixel.R > r.maxval) r.maxval = r.memval = pixel.R;
                if (pixel.G > r.maxval) r.maxval = r.memval = pixel.G;
                if (pixel.B > r.maxval) r.maxval = r.memval = pixel.B; 
            }
        }
        r.colres = r.memres = (r.maxval * res);
        export_image(r, render_name);
    }
};