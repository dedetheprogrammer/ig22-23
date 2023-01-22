#pragma once
#ifndef RENDER_H
#define RENDER_H

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <mutex>
#include <random>
#include "camera.hpp"
#include "geometry.hpp"
#include "kd_tree.hpp"
#include "objects.hpp"

//===============================================================//
// Scene: a set of cameras, objects and lights.
//===============================================================//
using Photon_map = nn::KDTree<Photon, 3, PhotonAxisPosition>;
using Camera_ptr = std::shared_ptr<Camera>;
using Object_ptr = std::shared_ptr<Object>; 

class Scene {
private:
    //...
public:

    // Scene cameras.
    std::vector<Camera_ptr> cameras; 
    
    // Scene objects.
    std::vector<Object_ptr> objects;

    // Scene lights.        
    std::vector<Light> lights;

    // Scene photon map
    Photon_map pmap;
    // How many photons to shoot in the first iteration between every light source.
    int INITIAL_PHOTONS;
    // Max number of total photons.
    int MAX_PHOTONS;
    // Max bounce photon depth.
    int MAX_PHOTON_DEPTH;

    Scene(std::vector<Camera_ptr> cameras, std::vector<Object_ptr> objects)
        : cameras(cameras), objects(objects) {}
    Scene(std::vector<Camera_ptr> cameras, std::vector<Object_ptr> objects, std::vector<Light> lights)
        : cameras(cameras), objects(objects), lights(lights) {}
    Scene(std::vector<Camera_ptr> cameras, std::vector<Object_ptr> objects, std::vector<Light> lights,
        int INITIAL_PHOTONS, int MAX_PHOTONS, int MAX_PHOTON_DEPTH)
        : cameras(cameras), objects(objects), lights(lights)
    {
        this->INITIAL_PHOTONS = INITIAL_PHOTONS;
        this->MAX_PHOTONS = MAX_PHOTONS; // No se usa de momento.
        this->MAX_PHOTON_DEPTH = MAX_PHOTON_DEPTH;
    }

};

//===============================================================//
// Render: the engine in charge of doing the magic.
//===============================================================//

// Intersection offset to avoid Shadow Acne and self intersection:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/ligth-and-shadows
#define SHADOW_BIAS 0.00001

class Render {
private:

    /**
     * @brief Ray tracing.
     * 
     * @param s Set of parameters (objects, lights, etc) needed by the render.
     * @return RGB 
     */
    inline RGB ray_tracing(Scene& s, Ray wo, Object_ptr p_obj, int depth) {

        // Closest collision.
        Collision c;

        // Calculating possible intersections in the Scene.
        for (auto& o : s.objects) {
            auto t = o->intersects(wo);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
                if (c.obj == nullptr) c.obj = o;
            }
        }
        return c.obj->m.kd(c.obj->b, c.normal, c.point, Vector3(), -wo.d);
    }

    /** 
     * @brief Path tracing.
     * 
     * @param s Set of parameters (objects, lights, etc) needed by the render.
     * @return RGB
     */
    
    inline RGB path_tracing(Scene& s, Ray wo, Object_ptr p_obj, int depth) {

        if (depth > RAY_MAX_DEPTH) return RGB();
        if (wo.d.mod() == 0) return RGB();

        // Light contributions.
        RGB direct_light_contrib, indirect_light_contrib;
        // Closest collision.
        Collision c; 
        // Calculating possible intersection with the scene:
        for (auto& o : s.objects) {
            auto t = o->intersects(wo);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
                if (c.obj == nullptr) c.obj = o;
            }
        }
        // If the ray doesn't intersect with any object, just return black color.
        if (c.obj == nullptr) return RGB();
        // If the ray intersects with an object with emission, return its light 
        // emission coefficient.
        if (c.obj->m.ke  > 0) return c.obj->m.ke;

        // Montecarlo sample:
        Sample samp = c.obj->m.scattering(c.obj->b, c.normal, c.point, wo.d,
            ((p_obj == nullptr) || (p_obj == c.obj) || (c.dist > SHADOW_BIAS) ?
                1.0 :
                p_obj->m.ref_index_i
            )
        );

        // Direct light computation. Delta explanation in Material class
        // (objects.hpp):
        if (!samp.is_delta) {
            for (auto& l : s.lights) {
                // Shadow ray. See shadow bias explanation in the variable
                // declaration:
                auto fix = c.point + (c.normal * SHADOW_BIAS);
                Ray lr(fix, l.c - fix);
                double l_mod = lr.d.mod();
                if (l_mod == 0.0) continue;

                bool collides = false;
                for (auto& o : s.objects) {
                    auto t = o->intersects(lr);
                    if (t.dist > 0 && t.dist < 1 /*&& (t.normal * lr.d) < 0*/) {
                        collides = true;
                        break;
                    }
                }
                if (collides) continue;

                // Render equation:
                direct_light_contrib +=
                    (l.pow / (l_mod * l_mod)) *
                    c.obj->m.kd(c.obj->b, c.normal, c.point, lr.d, -wo.d) / M_PI *
                    std::abs(c.normal * (lr.d / l_mod));
            }
        }

        // Indirect light computation.
        // Bouncing ray. Now the direction depends on the material properties:
        //  - Uniform cousine sample (Total difussion).
        //  - Reflection.
        //  - Refraction.
        //  - More..?
        indirect_light_contrib += path_tracing(s, Ray(c.point, samp.wi), c.obj, depth+1);
        return direct_light_contrib
            + indirect_light_contrib * samp.fr;
    }


    // ===========================
    // Photon mapping.
    // ===========================
    // Photon mapping parameters.
    int SEARCH_PHOTONS; // Number of photons to search.
    double RADIUS;      // Radius of the area where the photons will be searched.

    inline void import_pmap(Scene& s, std::string pmap_file) {
        // Opening the PMAP file.
        std::ifstream in(pmap_file);
        debug(!in.is_open(), "file '" + pmap_file + "' not found, check it out!", 1);

        int saved_photons;
        int dummy;
        in >> saved_photons;
        // Esto seria para leer la dimension de Vector posicion, del RGB flujo y
        // el vector Wi. Pero no nos hace falta, lo pongo por que queda más fancy.
        in >> dummy >> dummy >> dummy;

        bar.init(std::cout, "PHOTON MAP IMPORTATION", saved_photons);
        std::vector<Photon> photons(saved_photons);
        for (auto& p : photons) {
            auto t = now();
            in >> p;
            bar.update(std::cout, now()-t);
        }
        bar.kill(std::cout);
        s.pmap = Photon_map(photons);

    }

    // For casting a photon.
    inline void photon_cast_th(std::vector<Photon>& photons, const Scene& s,
        const Light& l, const int S, const int min, const int max)
    {
        // Casting photons
        Object_ptr p_obj;
        for (int p = min; p < max; p++) {
            // Uniform spherical sampling : Angular sampling.
            double lat = acos(2*E(e2) - 1); //acos(C(e2));
            double azi = 2*M_PI*E(e2);
            Ray r(l.c, Vector3(sin(lat)*cos(azi), sin(lat)*sin(azi), cos(lat)));
            RGB flux = (4*M_PI*l.pow)/S;

            auto t = now();
            for (int depth = 0; depth < s.MAX_PHOTON_DEPTH; depth++) {
                // If the photon flux is near to 0, just stop iterating.
                if (flux.rad() < EPSILON_ERROR) break;
                // Closest collision.
                Collision c;
                // Calculating possible intersections with the scene:
                for (auto& o : s.objects) {
                    auto t = o->intersects(r);
                    if (t.dist > 0 && t.dist < c.dist) {
                        c = t;
                        if (c.obj == nullptr) c.obj = o;
                    }
                }
                // If the photon doesn't intersect with any object of the scene
                // just stop iterating.
                if (c.obj == nullptr) break;
                // If the photon intersects with an object with light emission, 
                // stop iterating.
                if (c.obj->m.ke  > 0) break;

                // Colision point.
                Sample samp = c.obj->m.scattering(c.obj->b, c.normal, c.point, r.d,
                    ((p_obj == nullptr) || (p_obj == c.obj) || (c.dist > SHADOW_BIAS) ?
                        1.0 :
                        p_obj->m.ref_index_i
                    )
                );

                // Changing to the new direction.
                r = Ray(c.point, samp.wi);

                // Adding a new photon. If is delta material, doesn't make 
                // sense to store the photon:
                if (!samp.is_delta) {
                    m.lock();
                    photons.push_back(Photon(c.point, flux, samp.wi));
                    m.unlock();
                }

                // New flux resulting of the bounce.
                flux *= (samp.fr);
                p_obj = c.obj;
            }
            m.lock();
            bar.update(std::cout, now() - t);
            m.unlock();
        }
    }

    // For generating the scene photon map.
    void generate_photon_map(Scene& s) {
        // Total scene power.
        RGB total_pow;
        for (auto& l : s.lights) {
            total_pow += l.pow;
        }

        // Scene photons.
        std::vector<Photon> photons;
        for (auto& l : s.lights) {

            // Light proportional number of photons.
            int S = s.INITIAL_PHOTONS * l.pow.rad() / total_pow.rad();
            // Threading things.
            std::vector<std::thread> threads(N_THREADS);
            int task_range = S/N_THREADS;
            int extra_task = S%N_THREADS;
            for (int i = 0; i < 8; i++) {
                int ini = i*task_range + (i<extra_task)*i + (i>=extra_task)*extra_task;
                int end = (i+1)*task_range + (i<extra_task)*(i+1) + (i>=extra_task)*extra_task;
                threads[i] = std::thread(&Render::photon_cast_th, this,
                    std::ref(photons), std::ref(s), std::ref(l), S, ini, end+1);
            }
            for (auto& th : threads) {
                th.join();
            }
        }
        s.pmap = Photon_map(photons, PhotonAxisPosition());
    }

    /**
     * @brief Photon mapping 
     * 
     * @param s Set of parameters (objects, lights, etc) needed by the render.
     * @return RGB
     */
    RGB photon_mapping(Scene& s, Ray wo, Object_ptr p_obj, int depth) {

        if (depth > RAY_MAX_DEPTH) return RGB();
        // Light contribution.
        RGB photon_light_contrib;
        // Closest collision.
        Collision c;
        // Calculating possible intersections in the Scene.
        for (auto& o : s.objects) {
            auto t = o->intersects(wo);
            if (t.dist > 0 && t.dist < c.dist) {
                c = t;
                if (c.obj == nullptr) c.obj = o;
            }
        }
        // If the ray doesn't intersect with any object, just return black color.
        if (c.obj == nullptr) return RGB();
        // If the ray intersects with an object with emission, return its light 
        // emission coefficient.
        if (c.obj->m.ke  > 0) return c.obj->m.ke;

        // Object scattering and sampling:
        Sample samp = c.obj->m.scattering(c.obj->b, c.normal, c.point, wo.d,
            ((p_obj == nullptr || p_obj == c.obj || c.dist > SHADOW_BIAS) ? 
                1.0 : 
                p_obj->m.ref_index_i
            )
        );
        if (samp.fr == 0) return RGB();
        // If the object material is a delta material, we have to follow the ray 
        // path until it intersects with a diffuse object or dies.
        if (samp.is_delta) {
            return photon_mapping(s, Ray(c.point, samp.wi), c.obj, depth+1) * samp.fr;
        }

        // Maximum distance to look for photons:
        double real_rad = 0;

        // Nearest photons to the collision point:
        auto neighbors = (RADIUS > 0.0) ?
                    s.pmap.nearest_neighbors(c.point, SEARCH_PHOTONS, RADIUS) :
                    s.pmap.nearest_neighbors(c.point, SEARCH_PHOTONS);

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
                (c.obj->m.kd(c.obj->b, c.normal, c.point, p->wi, wo.d)/M_PI * p->flux * RBFK);
        }

        return photon_light_contrib;
    }

    // ====================================
    // Threading.
    // ====================================
    int N_THREADS; // Number of threads.
    std::mutex m;  // Mutex for concurrency.
    void render_th(Scene& s, Camera_ptr c, int min, int max) {
        for (int x = min; x < max; x++) {
            for (int y = 0; y < c->image_width; y++) {
                // For render ETA purposes.
                auto t = now();
                // First point will be the pixel center.
                c->view[x][y] += render_function(s, c->cast_ray(x,y),nullptr,0);
                for (int a = 1; a < c->ppp; a++) {
                    c->view[x][y] += render_function(s,c->cast_antialiasing_ray(x,y),nullptr,0);
                }
                m.lock();
                c->update_view(x,y);
                // Updating bar status.
                bar.update(std::cout, now()-t);
                m.unlock();
            }
        }
    }

    // ====================================
    // Render mode.
    // ====================================
    enum Mode { RAY_TRACING, PATH_TRACING, PHOTON_MAPPING };
    friend std::ostream& operator<<(std::ostream& os, const Mode& m) {
        if (m == RAY_TRACING) return os << "Ray tracing";
        else if (m == PATH_TRACING) return os << "Path tracing";
        else if (m == PHOTON_MAPPING) return os << "Photon mapping";
    }
    Mode mode; // Render mode.
    // Render function.
    std::function<RGB(Scene&, Ray, Object_ptr, int)> render_function;

public:

    // Path tracing parameters.
    int RAY_MAX_DEPTH;

    Progress_bar bar = Progress_bar(80, STYLE2, 100); // Render progress bar.

    /**
     * @brief Construct a new Ray tracing renderer.
     * 
     */
    Render(int N_THREADS) {
        this->N_THREADS = (N_THREADS == 0) ? 1 : N_THREADS;
        mode = RAY_TRACING;
        render_function = std::bind(&Render::ray_tracing, this, std::placeholders::_1,
            std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
    }

    /**
     * @brief Construct a new Path tracing renderer.
     * 
     * @param RAY_MAX_DEPTH Max bounces per ray.
     */
    Render(int RAY_MAX_DEPTH, int N_THREADS) {
        this->N_THREADS = (N_THREADS == 0) ? 1 : N_THREADS;
        this->RAY_MAX_DEPTH = RAY_MAX_DEPTH;
        mode = PATH_TRACING;
        render_function = std::bind(&Render::path_tracing, this, std::placeholders::_1,
            std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
    }

    /**
     * @brief Construct a new Photon mapping renderer.
     * 
     * @param SEARCH_PHOTONS Number of photons to search in each intersection.
     * @param RADIUS Area radius to search photons into.
     */
    Render(int RAY_MAX_DEPTH, int SEARCH_PHOTONS, double RADIUS, int N_THREADS) {
        this->N_THREADS = (N_THREADS == 0) ? 1 : N_THREADS;
        this->RAY_MAX_DEPTH  = RAY_MAX_DEPTH;
        this->SEARCH_PHOTONS = SEARCH_PHOTONS;
        this->RADIUS = RADIUS;
        mode = PHOTON_MAPPING;
        render_function = std::bind(&Render::photon_mapping, this, std::placeholders::_1,
            std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
        
    }

    void render(Scene& s, int camera_index = 0, std::string pmap_file = "") {

        Camera_ptr c = s.cameras[camera_index];
        if (mode == PHOTON_MAPPING && s.pmap.empty()) {
            if (pmap_file.empty()) {
                bar.init(std::cout, "PHOTON MAPPING GENERATION", s.INITIAL_PHOTONS);
                std::cout << "Photon mappingn render: "
                    << "\n> Search photons number: " << SEARCH_PHOTONS
                    << "\n> Search radius: " << RADIUS
                    << "\n> Scene specificactions: "
                    << "\n   > INITIAL PHOTONS: " << s.INITIAL_PHOTONS
                    << "\n   > MAX PHOTONS: " << s.MAX_PHOTONS 
                    << "\n   > MAX PHOTONS DEPTH: " << s.MAX_PHOTON_DEPTH
                    << "\n"; 
                generate_photon_map(s);
                bar.kill(std::cout);
            } else {
                import_pmap(s, pmap_file);
            }
        }

        bar.init(std::cout, "SCENE RENDERIZATION", c->image_width * c->image_height);
        std::cout << *c << std::endl;

        // The thread thing.
        int task_range = c->image_height / N_THREADS;
        int extra_task = c->image_height % N_THREADS;
        std::vector<std::thread> threads(N_THREADS);
        for (int i = 0; i < N_THREADS; i++) {
            int ini = i*task_range + (i<extra_task)*i + (i>=extra_task)*extra_task;
            int end = (i+1)*task_range + (i<extra_task)*(i+1) + (i>=extra_task)*extra_task;
            threads[i] = std::thread(&Render::render_th, this, std::ref(s), c, ini, end);
        }
        for (auto& th : threads) {
            th.join();
        }

        /*
        std::vector<std::thread> threads;
        int thread_range = c->image_width / N_THREADS;
        int leftover_ran = c->image_width % N_THREADS;
        for (int x = 0; x < c->image_height; x++) {
            for (int i = 0; i < N_THREADS; i++) {
                int ini = i * thread_range;
                int end = (i+1) * thread_range - 1;
                if (i < leftover_ran) {
                    ini += i;
                    end += i+1;
                } else {
                    ini += leftover_ran;
                    end += leftover_ran;
                }
                threads.push_back(std::thread(&Render::render_th, this, std::ref(s), c, x, ini, end+1));
            }
        }
        for (auto& th : threads) {
            th.join();
        }
        */

        /*
        for (int x = 0; x < c->image_height; x++) {
            for (int y = 0; y < c->image_width; y++) {
                // For render ETA purposes.
                auto t = now();
                // First point will be the pixel center.
                c->view[x][y] += render_function(s, c->cast_ray(x,y),nullptr,0);
                for (int a = 1; a < c->ppp; a++) {
                    c->view[x][y] += render_function(s,c->cast_antialiasing_ray(x,y),nullptr,0);
                }
                c->update_view(x,y);
                // Updating bar status.
                bar.update(std::cout, now()-t);
            }
        }*/
        bar.kill(std::cout);

    }

    void export_render(Scene& s, int camera_index = 0, std::string render_name = "./new_scene.ppm", int res = 10e7) {
        s.cameras[camera_index]->export_view(render_name, res);
    }

    void export_pmap(Scene& s, std::string pmap_file = "./new_scene.pmap") {
        std::ofstream os(pmap_file);
        os << s.pmap.size() << "\n";
        os << "3 3 3\n";
        for (auto& p : s.pmap.get_elements()) {
            os << p.pos.x  << " " << p.pos.y  << " " << p.pos.z  << " "
               << p.flux.R << " " << p.flux.G << " " << p.flux.B << " "
               << p.wi.x   << " " << p.wi.y   << " " << p.wi.z   << "\n";
        }
        os.close();
    }
};

#endif