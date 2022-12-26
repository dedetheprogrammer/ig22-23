#include <iostream>
#include <functional>
#include "scene.hpp"

int main(int argc, char* argv[]) {

    Scene s({
        // Camaras.
        Camera(Camera_settings({Vector3(0,0,-3.5), Vector3(-1,0,0), Vector3(0,1,0), Vector3(0,0,3)}, 64, 100, 0.05), 256, 256)
    }, {
        // Objetos.
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(210,0,0)))),    //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(0,210,0)))),   //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(240,240,240)))),    //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240)))),   //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))),     //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(220,145,220), RGB(80,80,80)))),  //Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, Material(RGB(125,200,200), RGB(30,30,30), RGB(255,255,255), RGB(), 1.5))),  //Esfera de la derecha
    }, {
        // Luces puntuales.
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0))
    }, 10000000, 5);

    s.render();
    s.export_render();

}