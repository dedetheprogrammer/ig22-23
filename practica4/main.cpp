#include <iostream>
#include "scene.hpp"
#include "objects.hpp"

// TODO: Hacer modificadores poco a poco (decidir cuales si y no).
// - OK ... Other geometries
// - OK ... Import indirect light.
// - OK ... Refactor Scene.
// - OK ... Russian Roulette.
// - OK ... BXDFs (Scattering).
// - OK ...... Diffussion.
// - OK ...... Reflectance.
// - OK ...... Refractance.
// - ?? ...... More?
// - IP ... Fresnell effects.  
// - OK ... Light areas.
// - OK ... Start looking into photon mapping.
// - OK ... Make the photon mapping.
// - OK ... Refactor scene: Ray tracing + Path tracing + Photon mapping options.
// -    ... Debug photon mapping.
// - OK ... Finite planes
// - OK ... Boxes
// - IP ... Meshes
// -    ... Acceleration structures.
// -    ... Parallelization.
// -    ... Debug and log system.
// -    ... Textures.
// -    ... Normal mapping.
// -    ... Spectral rendering.
// -    ... Motion blur.
// -    ... Depth of field.
// -    ... Importance sampling next event estimation. (I think we have done this already...).
// -    ... Bidirectional path tracing.
// -    ... Participating media.
// -    ... RGB luminance.
int main (int argc, char* argv[]) {

    Scene s({
        // Camaras.
        Camera(Camera_settings({Vector3(0,0,-4), Vector3(-1,0,0), Vector3(0,1,0), Vector3(0,0,3)}, 256, 5), 512, 512)
    }, {
        // Objetos.
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(210,0,0)))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(0,210,0)))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(240)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240), RGB(), RGB(), RGB(/*1.0, 1.0, 1.0*/)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240)))), //Plano de atrás
        std::make_shared<Plane>(Plane({Vector3(-0.5,0.998,-0.5), Vector3(-0.5,0.998,0.5), Vector3(0.5,0.998,0.5), Vector3(0.5,0.998,-0.5)}, Material(RGB(180,180,180), RGB(), RGB(), RGB(1.0, 1.0, 1.0)))),
        //std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(220,145,220), RGB(80,80,80)))), // Esfera de la izquierda
        //std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, Material(RGB(/*125,200,200*/), RGB(30,30,30), RGB(255,255,255), RGB(), 1.5))), // Esfera de la derecha
        //std::make_shared<Sphere>(Sphere(Vector3(0,0.9,0), 0.3, Material(RGB(/*125,200,200*/), RGB(30,30,30), RGB(255,255,255), RGB(0,0,80), 1.5))),  //Esfera de la derecha
        //std::make_shared<Cube>(Cube(Vector3(-0.25,-1,-0.25), Vector3(0.25,-0.5,0.25), Material(RGB(240,80,0)))),
        //std::make_shared<Cube>(Cube({
        //        Vector3(0,-1,-0.5), Vector3(0,-0.5,-0.5), Vector3(0.5,-0.5,0), Vector3(0.5,-1,0),
        //        Vector3(-0.5,-1,0), Vector3(-0.5,-0.5,0), Vector3(0,-0.5,0.5), Vector3(0,-1,0.5)
        //    }, 
        //    Material(RGB(0,180,20), RGB(), RGB(), RGB(/*0.5, 0.75, 0.5*/)))
        //),
        // std::make_shared<Sphere>(Sphere(Vector3(-0.25,-0.7,-1), 0.3, Material(RGB(8,140,30)))),  //Esfera de la derecha
        // std::make_shared<Sphere>(Sphere(Vector3(0.7,-0.7,-0.5), 0.3, Material(RGB(140,80,10)))),  //Esfera de la derecha
        // std::make_shared<Sphere>(Sphere(Vector3(-0.9,-0.7,0), 0.3, Material(RGB(70,120,20)))),  //Esfera de la derecha
        // std::make_shared<Sphere>(Sphere(Vector3(0.1,-0.7,1), 0.3, Material(RGB(56,209,210)))),  //Esfera de la derecha
        // std::make_shared<Sphere>(Sphere(Vector3(0.9,-0.7,0.9), 0.3, Material(RGB(203,0,0)))),  //Esfera de la derecha
        // std::make_shared<Sphere>(Sphere(Vector3(0,-0.7,0), 0.3, Material(RGB(9,207,0)))),  //Esfera de la derecha
        //std::make_shared<Mesh>(Mesh("../src/ply/cow.ply") * Matrix3Scale(0.4) * Matrix3Translation(0, -0.6, 0) * Matrix3Rotation(Y_ROT, 45)),
        std::make_shared<Mesh>(Mesh("../src/ply/mario.ply") * Matrix3Scale(0.07) * Matrix3Rotation(Y_ROT, -135) * Matrix3Translation(0,-0.1,0)),
        //std::make_shared<Mesh>(Mesh ({
        //    Triangle({Vector3(-0.5,0,0.5), Vector3(0,0.5,0), Vector3(0,-0.5,0)},   RGB(205, 170, 240)),
        //    Triangle({Vector3(0,-0.5,0), Vector3(0,0.5,0), Vector3(0.5,0,0.5)},    RGB(205, 170, 240)),
        //    Triangle({Vector3(-0.5,0,0.5), Vector3(0,0.5,0), Vector3(0.5,0,0.5)},  RGB(205, 170, 240)),
        //    Triangle({Vector3(-0.5,0,0.5), Vector3(0,-0.5,0), Vector3(0.5,0,0.5)}, RGB(205, 170, 240)),
        //    },
        //    Material(RGB("#FCD0B4"))) * Matrix3Rotation(X_ROT, -45)
        //),
    }, {
        // Luces puntuales.
        //Light(Vector3(0,0.5,0), RGB(1.0,1.0,1.0))
    });

    s.render();
    s.export_render();

}