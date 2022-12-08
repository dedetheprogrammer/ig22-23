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

// TODO:
// - RGB luminance.
// - Import indirect light.
// - Start looking into photon mapping.
// - Threading.
// - Implement kd-tree algorythm.
// - 3d meshes textures.
// - Normal mapping.

//===============================================================//
// Scene: a set of objects and lights.
//===============================================================//

using Object_ptr = std::shared_ptr<Object>; 
using Light_ptr  = std::shared_ptr<Light>;
using Objects    = std::vector<std::shared_ptr<Object>>;
using Lights     = std::vector<std::shared_ptr<Light>>;

std::mutex m;
std::random_device rd;
std::mt19937 e2(rd());

class Scene {
private:
    // ...
public:
    Objects objects; // Scene objects.
    Lights  lights;  // Scene lights.

    Scene() {}
    Scene(Objects objects, Lights lights) : objects(objects), lights(lights) {}

};

//===============================================================//
// Pixel: the elemental of the grid.
//===============================================================//

bool once = true;
#define BIAS 0.001f
#define MAX_DEPTH 2

class Pixel {
private:

    // Render equation.
    /*inline RGB RenderEquation(Object_ptr c_obj, Vector3 casted_dir, RGB pow, Ray light_ray, float l_mod) {
        Vector3 n = c_obj->normal(light_ray.p, casted_dir);
        //     Incoming light:         BDRF:         Object geometry:
        return pow / (l_mod * l_mod) * c_obj->fr() * std::abs(n * (light_ray.d / l_mod));
    }

    inline RGB find_path(Scene s, Ray r, int n_bounce = 0) {
        RGB ret;
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
        if (c_obj != nullptr) {
            Vector3 x = r.d * c_dist + r.p;
            RGB lDir = direct_light(s, r, c_obj, x);
            ret += lDir;
            if (n_bounce < MAX_BOUNCES) {
                ret += indirect_light(s, r, c_obj, x, n_bounce);
            }
            ret *= c_obj->fr();
        }
        return ret;
    }

    // Direct light renderization.
    inline RGB direct_light(const Scene& s, Ray r, Object_ptr c_obj, Vector3 x) {
        
        RGB dl_color;
        for (auto& l : s.lights) {
        
            Ray lr(x, l->c - x);
            float l_mod = lr.d.mod();
            if (l_mod == 0) continue;

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
            dl_color += RenderEquation(c_obj, r.d, l->pow, lr, l_mod);
        }
        return dl_color;
    }

    // Indirect light renderization.
    inline RGB indirect_light(Scene& s, Ray r, Object_ptr c_obj, Vector3 x, int n_bounce) {

        std::uniform_real_distribution<> E(0, 1);
        float lat = acos(sqrt(1 - E(e2)));
        float azi = 2*M_PI*E(e2);

        std::vector<Vector3> b = orthonormal_basis(nor(c_obj->normal(x, r.d)));
        Vector3 dir = Matrix3BaseChange(b[1], b[0], b[2], x) * 
            Vector3(sin(lat) * cos(azi), sin(lat) * sin(azi), cos(lat));
 
        // Bouncing ray
        return find_path(s, Ray(x, dir), n_bounce+1) * M_PI;
    }*/

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
                /*(c_obj->fr()) **/
                std::abs(n * (lr.d / l_mod));
        }

        // Compute indirect light.
        std::uniform_real_distribution<> E(0, 1);
        float lat = acos(sqrt(1 - E(e2))); // SE GENERAN COMO RADIANES
        float azi = 2*M_PI*E(e2);          // LO HE COMPROBADO.
        // EXPLICACIÓN DEL NUEVO RAYO GENERADO: 
        // 1. Utilizando una matriz de cambio de base y/o las propiedades 
        //    trigonometricas para pasar de polares a cartesianas daría malos
        //    resultados porque estas propiedades se basan suponiendo que 
        //    nuestro vector esta alineado a los ejes cartesianos.
        // 2. Coge por ejemplo sin(lat) * cos(azi), sabemos que el seno calcula
        //    el cateto opuesto y el coseno el contiguo, si tu vector es paralelo
        //    al eje Y, la componente x del nuevo vector si es sin(lat) * cos(azi).
        //    Pero si tu vector es paralelo al eje x, es la componente y quien debería
        //    tener el valor sin(lat) * cos(azi). Esto hace que en todos los puntos de
        //    intersección la hemiesfera este alineada al eje Y, por eso hay rayos
        //    que se generan detras de la pared al rotarlos.
        // 3. Usar esto nos hace totalmente dependientes de los ejes, por lo que 
        //    he hecho un truquito de magia para que la rotación sea correcta.
        //
        // GENERAMOS UNA BASE ORTONORMAL RESPECTO DE LA NORMAL OBTENIDA.
        // Se devuelven los dos vectores restantes que conforman la base ortonormal.
        std::vector<Vector3> b = orthonormal_basis(n);
        // Ahora, lo que hacemos es rotar uno de estos dos vectores respecto del 
        // otro no la latitud sino la colatitud, ya que los dos estan contenidos
        // en el mismo plano. Además, debe ser negativa la rotación ya que estamos
        // rotando este vector en el sentido de las agujas del reloj.
        Vector3 dir = rot(b[0], b[1], -(lat));
        // Lo último que queda es rotar el azimuth. Para ello, cogemos el vector
        // rotado colatitud radianes y lo rotamos respecto de la normal azimuth 
        // radianes. Da igual el sentido, sea uno u otro, como rotas hasta 2π,
        // pues te da igual.
            dir = nor(rot(dir, n, azi));
        
        // Bouncing ray
        indirect_light_contrib += cast_ray(s, Ray(x, dir), depth+1);

        return (direct_light_contrib + indirect_light_contrib) * c_obj->fr() * M_PI;
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
            color += c_obj->get_kd();
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