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

// TODO: investigate KD-trees.
// TODO: restructuration of the code.

// Obtiene el número de nucleos del sistema. Parece ser que si usas todos los núcleos (
// thread por nucleo) el rendimiento empeora, de ahi que se use dos menos.
int THREADS = std::thread::hardware_concurrency() - 2;
std::mutex m;

using Scene = std::vector<std::shared_ptr<Object>>;

class Pixel {
private:

    inline RGB generate_pixel(Scene figures, Vector3 cc, Vector3 dot) {

        // Closest color.
        RGB dColor;
        // Distance of the closest intersection point.
        float closest = INFINITY;
        // Generated ray.
        Ray r = Ray(cc, dot-cc);

        // Calculating possible intersections in the Scene.
        for (auto& f : figures) {
            auto vt = f->intersects(r);
            if (vt.size() > 0 && vt[0] > 0 && (vt[0] < closest)) {
                closest = vt[0];
                dColor = f->get_emission();
            }
        }
        return dColor;
    }

    // Función para el thread.
    void calculate_th(Scene figures, Vector3 cc, int min, int max) {
        for (int i = min; i < max; i++) {
            RGB dColor = generate_pixel(figures, cc, dots[i]);
            m.lock();
            color += dColor;
            m.unlock();
        }
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

    void calculate_wout_threads(Scene figures, Vector3 cc) {
        for (auto& dot : dots) {
            color += generate_pixel(figures, cc, dot);
        }
        color /= dots.size();
    }

    void calculate_with_threads(Scene figures, Vector3 cc) {

        // Dots sobrantes si el número de dots no es exacto al número de núcleos, y dots por
        // thread.
        int first_take = dots.size()%THREADS, dots_th = (dots.size()-first_take)/THREADS;
        // Vector de threads.
        std::vector<std::thread> ths(THREADS);

        // Es posible que al pasarle una serie de muestras no sean exactas al número de núcleos
        // por lo que sobraran tantas como modulo de número de núcleos, por lo que estos primeros
        // threads recibiran un dot más.
        for (int i = 0; i < first_take; i++) {
            ths[i] = std::thread(&Pixel::calculate_th, std::ref(*this), figures, cc, 
                i * (dots_th+1), (i+1) * (dots_th+1));
        }
        // Una vez se repartan los dots sobrantes, el resto de threads recibe el número que
        // recibiría si la division fuese exacta. Si es exacta desde el principio, el primer for
        // se lo suda.
        for (int i = first_take; i < THREADS; i++) {
            ths[i] = std::thread(&Pixel::calculate_th, std::ref(*this), figures, cc, 
                i * dots_th, (i+1) * dots_th);
        }

        // Refresco de PSCD, join de los threads para que el proceso padre espere a que terminen.
        for (auto& th : ths) {
            th.join();
        }

        // Color medio del pixel.
        color /= dots.size();
    }
};

//===============================================================//
// Camera
//===============================================================//

using CameraGridRow = std::vector<Pixel>;
using CameraGrid    = std::vector<CameraGridRow>;

class Camera {
private:
    // ...
public:
    Vector3 c, l, u, f; // Center. Left + Up + Front dimensions.
    Vector3 pl, pu;     // Pixel Left + Up dimensions.
                        // - u = center->up border; size = down border->up border = 2u
                        // - l = center->left; size = right border->left border = 2l
    int h, w, ppp;      // Height, width and points per pixel.
    CameraGrid pixels;  // Camera grid.

    Camera(Vector3 c, Vector3 l, Vector3 u, Vector3 f, int h, int w, int ppp = 1) 
        : c(c), l(l), u(u), f(f), pl(2*l/w), pu(2*u/h), h(h), w(w), ppp(ppp)
    {
        // std::cout << "Cores: " << THREADS << std::endl;
        pixels = CameraGrid(h, CameraGridRow(w));

        Vector3 pc = c + u-pu/2 + l-pl/2 + f; // Pixel center.
        // Randomice.
        std::random_device rd;
        std::mt19937 e2(rd());

        // srand (static_cast <unsigned> (time(0))); // Generate seed for the pixels' dots' positions
        for (auto& hi : pixels) {
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
                pc = pc - pl;
            }
            pc = pc + 2*l - pu;
        }
    }


    // Without using threads (8): 
    //  ppp(1024): 108025ms
    //  ppp(2048): 216056ms
    // Using threads (8):
    //  ppp(1024): 21947ms
    //  ppp(2048): 39171ms, 32130ms
    void rayTrace(Scene scene, int m) {

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        
        // Al parecer, si hay pocos ppp, va mucho peor que haciendolo secuencial.
        if (!m/*|| ppp < THREADS*5*/) {
            for (auto& row : pixels) {
                for (auto& pixel : row) {
                    pixel.calculate_wout_threads(scene, c);
                }
            }
        } else {
            for (auto& row : pixels) {
                for (auto& pixel : row) {
                    pixel.calculate_with_threads(scene, c);
                }
            }
        }

        std::cout << "Time cost = " <<
           std::chrono::duration_cast<std::chrono::milliseconds>
           (std::chrono::steady_clock::now() - begin).count() << "[ms]" << std::endl;
    }

    Channels getRGB(int& colres) {
        float max = 0;
        Channels colors(h, std::vector<RGB>(w));
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                colors[i][j] = pixels[i][j].color;
                if (colors[i][j].R > max) max = colors[i][j].R;
                if (colors[i][j].G > max) max = colors[i][j].G;
                if (colors[i][j].B > max) max = colors[i][j].B; 
            }
        }
        colres = max * 255.0;
        return colors;
    }
};