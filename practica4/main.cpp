#include <iostream>
#include "render.hpp"

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
// - OK ... Fresnell effects.  
// - OK ... Light areas.
// - OK ... Start looking into photon mapping.
// - OK ... Make the photon mapping.
// - OK ... Refactor scene: Ray tracing + Path tracing + Photon mapping options.
// - OK ... Debug photon mapping.
// - OK ... Finite planes
// - OK ... Boxes
// - OK ... Meshes
// - IP ... Debug and log system.
// - OK ...... Objects and materials.
// - OK ...... Camera and render.
// - ?? ... Acceleration structures.
// - OK ... Parallelization.
// - OK ... Textures.
// - ?? ... Normal mapping.
// - NO ... Spectral rendering.
// - ?? ... Motion blur.
// - ?? ... Depth of field.
// - ?? ... Importance sampling next event estimation. (I think we have done this already...).
// - NO ... Bidirectional path tracing.
// -    ... Participating media.
// -    ... RGB luminance.
int main (int argc, char* argv[]) {

    Scene s({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,-0.9, -2), -22, 0, -4, 120, Camera::HORIZONTAL, 3.1, 128, 1280, 720)),
    },
    {
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(), RGB(255,255,255)))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(), RGB(255,255,255)))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(30,30,30), RGB(40,40,40)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB()))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(), RGB(255,255,255)))), //Plano de atrás
        std::make_shared<Plane>(Plane(1, Vector3(0,0,1), Material(RGB(), RGB(255,255,255)))), //Plano de atrás
        std::make_shared<Cube>(Cube(Vector3(-0.6,0.5,-0.6), Vector3(-0.55,0.995,-0.55), Material(RGB(45, 45, 45), RGB(), RGB(), RGB(210,210,210)))),
        std::make_shared<Cube>(Cube(Vector3(0.6,0.25,0.6), Vector3(0.65,0.995,0.65), Material(RGB(45, 45, 45), RGB(), RGB(), RGB(210,210,210)))),
        std::make_shared<Cube>(Cube(Vector3(-0.45,0.7,0.45), Vector3(-0.4,0.995,0.4), Material(RGB(45, 45, 45), RGB(), RGB(), RGB(210,210,210)))),
        std::make_shared<Cube>(Cube(Vector3(0.65,0.4,-0.45), Vector3(0.7,0.995,-0.4), Material(RGB(45, 45, 45), RGB(), RGB(), RGB(210,210,210)))),
        std::make_shared<Cube>(Cube(Vector3(-0.7,0.4,-0.05), Vector3(-0.75,0.995,0), Material(RGB(45, 45, 45), RGB(), RGB(), RGB(210,210,210)))),

        //std::make_shared<Plane>(Plane({Vector3(-0.5,-0.995,-0.5), Vector3(-0.5,-0.995,0.5), Vector3(0.5,-0.995,0.5), Vector3(0.5,-0.995,-0.5) }, {}, {}, Material(RGB(), RGB(), RGB(), RGB(1.0,1.0,1.0)))),
        //std::make_shared<Mesh>(Mesh("../src/ply/cherry-neon.ply", Material(RGB(), RGB(), RGB(), RGB(204, 214, 21))) * Matrix3Rotation(Y_ROT, -90) * Matrix3Scale(2.5) * Matrix3Translation(0, 0.2, 0.999)),
        //std::make_shared<Mesh>(Mesh("../src/ply/cherry-neon.ply", Material(RGB(), RGB(), RGB(), RGB(186, 56, 209))) * Matrix3Rotation(Y_ROT, -90) * Matrix3Scale(2.5) * Matrix3Translation(0.6, 0.2, 0.999)),
        //std::make_shared<Mesh>(Mesh("../src/ply/cherry-neon.ply", Material(RGB(), RGB(), RGB(), RGB(29, 29, 171))) * Matrix3Rotation(Y_ROT, -90) * Matrix3Scale(2.5) * Matrix3Translation(-0.6, 0.2, 0.999)),
        std::make_shared<Cube>(Cube(Vector3(-0.2,-0.1,-0.2), Vector3(0.2,0.3,0.2), Material(RGB(), RGB(), RGB(), RGB(15, 15, 15)))),
        std::make_shared<Mesh>(Mesh("../src/ply/pedestal.ply", Material(RGB(107, 107, 107), RGB(40,40,40))) * Matrix3Scale(7) * Matrix3Translation(0, -1, 0)),
        std::make_shared<Mesh>(Mesh("../src/ply/portal-out1.ply", Material(RGB(255, 255, 255))) * Matrix3Scale(0.1) * Matrix3Translation(0, 0.1, 0)),
        std::make_shared<Mesh>(Mesh("../src/ply/portal-out2.ply", Material(RGB(143, 203, 227), RGB(), RGB(), RGB(143, 203, 227))) * Matrix3Scale(0.1) * Matrix3Translation(0, 0.1, 0)),
        //std::make_shared<Mesh>(Mesh("../src/ply/water.ply", Material(RGB(), RGB(20,20,20), RGB(235,235,235), RGB(), 1.33)) * Matrix3Scale(0.35) * Matrix3Translation(0, -0.95, 0)),
        std::make_shared<Sphere>(Sphere(Vector3(-0.4, -1.18, -0.3), 0.2, Material(RGB(), RGB(), RGB(), RGB(1.0,1.0,1.0)))),
        std::make_shared<Sphere>(Sphere(Vector3(0.4, -1.18, -0.3), 0.2, Material(RGB(), RGB(), RGB(), RGB(1.0,1.0,1.0)))),
        std::make_shared<Sphere>(Sphere(Vector3(-0.4, -1.18, -0.8), 0.2, Material(RGB(), RGB(), RGB(), RGB(1.0,1.0,1.0)))),
        std::make_shared<Sphere>(Sphere(Vector3(0.4, -1.18, -0.8), 0.2, Material(RGB(), RGB(), RGB(), RGB(1.0,1.0,1.0)))),
    
    }, {

    });
    //s.objects[7]->collider.print(std::cout, 0);

    // Cornell box scene with point light.
    Scene s0({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 130, Camera::HORIZONTAL, 8, 512, 1800, 1350)),
        std::make_shared<Orthographic_camera>(Orthographic_camera(Vector3(0,0,-3.5), Vector3(-1,0,0), Vector3(0,1,0), Vector3(0,0,3), 256, 512, 512)),
    }, {
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(210,0,0), RGB(), RGB()))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(0,210,0), RGB(), RGB()))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(240,240,240)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(240, 110, 203), RGB(80,80,80)))), // Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, Material(RGB(/*125,200,200*/), RGB(25,25,25), RGB(230,230,230), RGB(), 1.5))), // Esfera de la derecha
    }, {
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0)),
    });

    // Cornell box scene with area light.
    Scene s1({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 130, Camera::HORIZONTAL, 8, 512, 1800, 1350)),
    }, {
        // Objects
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(210,0,0), RGB(), RGB()))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(0,210,0), RGB(), RGB()))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(240,240,240)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240), RGB(), RGB(), RGB(1.0, 1.0, 1.0)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(240, 110, 203), RGB(80,80,80)))),                  // Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, Material(RGB(), RGB(25,25,25), RGB(230,230,230), RGB(), 1.5))), // Esfera de la derecha
    }, {
        // Lights.
    });


    Scene s2({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 100, Camera::HORIZONTAL, 3.5, 512, 1000, 1000)),
    }, {
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(214, 200, 163), RGB(), RGB()))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(214, 200, 163), RGB(), RGB()))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(150, 93, 44)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(214, 200, 163)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(214, 200, 163)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(157, 209, 94), RGB(80,80,80)))),
        std::make_shared<Sphere>(Sphere(Vector3(0,-0.9,0), 0.1, Material(RGB(94, 117, 209)))), 
        std::make_shared<Sphere>(Sphere(Vector3(-0.3,-0.6,0.8), 0.4, Material(RGB(222, 118, 49)))),
        std::make_shared<Sphere>(Sphere(Vector3(0.4,-0.7,-0.4), 0.3, Material(RGB(), RGB(25,25,25), RGB(230,230,230), RGB(), 1.5))),
    }, {
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0)),
    });


    Scene s3({
        // Camaras.
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,-0.8, -3.5), -10, 0, 0, 130, Camera::HORIZONTAL, 8, 512, 1800, 1350)),
    }, {
        // Objetos.
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(210,0,0)))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(0,210,0)))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(240,240,240)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(-0.25,-0.7,1), 0.3, Material(RGB(8,140,30)))), 
        std::make_shared<Sphere>(Sphere(Vector3(0.7,-0.7,-0.5), 0.3, Material(RGB(140,80,10)))),
        std::make_shared<Sphere>(Sphere(Vector3(-0.8,0.6,0.9), 0.3, Material(RGB(), RGB(255,255,255)))), // Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(-0.9,-0.7,0), 0.3, Material(RGB(70,120,20)))), 
        std::make_shared<Sphere>(Sphere(Vector3(0,-0.7,1), 0.3, Material(RGB(56,209,210), RGB(100, 100, 100)))),
        std::make_shared<Sphere>(Sphere(Vector3(0.9,-0.7,0.9), 0.3, Material(RGB(203,0,0)))),
        std::make_shared<Sphere>(Sphere(Vector3(0,-0.7,0), 0.3, Material(RGB(9,207,0), RGB(), RGB(), RGB(200, 10, 40)))), // Esfera rosa inferior.
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.85,-0.7), 0.15, Material(RGB(), RGB(255,255,255)))), // Esfera metalica izquierda mas cercana.
        std::make_shared<Sphere>(Sphere(Vector3(0.0,0.65,0.25), 0.35, Material(RGB(), RGB(), RGB(), RGB(20,235,20)))), // Esfera verde superior.
        std::make_shared<Sphere>(Sphere(Vector3(-0.95,0.1,0.1), 0.15, Material(RGB(), RGB(), RGB(), RGB(227,115,30)))), // Esfera verde superior.
        std::make_shared<Sphere>(Sphere(Vector3(0.75,-0.3,0.75), 0.5, Material(RGB(220,90,76), RGB(10, 10, 10)))), // Esfera de plastico a la derecha.
        std::make_shared<Sphere>(Sphere(Vector3(0.3,-0.8,-0.8), 0.2, Material(RGB(), RGB(10, 10, 10), RGB(245,245,245), RGB(), 1.2))), // Esfera de plastico a la derecha.
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,0,0.9), 0.3, Material(RGB(207, 137, 39)))), // Esfera de plastico a la derecha.
        std::make_shared<Sphere>(Sphere(Vector3(0.75,0.6,0.9), 0.25, Material(RGB(43, 72, 237)))), // Esfera de plastico a la derecha.
        std::make_shared<Sphere>(Sphere(Vector3(0.2,0.35,1), 0.4, Material(RGB(43, 237, 185), RGB(25,25,25)))), // Esfera de plastico a la derecha.
        std::make_shared<Sphere>(Sphere(Vector3(0.8,-0.9,-0.9), 0.05, Material(RGB(), RGB(), RGB(), RGB(52, 103, 234)))), // Esfera de plastico a la derecha.
    }, {});

    // Color bleeding
    Scene s4({
        // Camaras.
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 100, Camera::HORIZONTAL, 3.5, 256, 1000, 1000)),
    }, {
        // Objetos.
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(240,240,240)))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(240,240,240)))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(210,0,0)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))), //Plano de atrás
        std::make_shared<Cube>(Cube(Vector3(-0.9,-1,-0.1), Vector3(-0.2,-0.4,0.5), Material(RGB(14, 230, 50)))),
        std::make_shared<Sphere>(Sphere(Vector3(0.7,-0.7,0.7), 0.3, Material(RGB(178, 61, 224)))),
        std::make_shared<Mesh>(Mesh("../src/ply/pawn.ply", Material(RGB(205,185,121), RGB(50,50,50))) * Matrix3Scale(0.2) * Matrix3Translation(1.5, -1.305, 0.6)),

    }, {
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0)),
    });

    // Hard shadows.
    Scene s5({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 90, Camera::HORIZONTAL, 3, 512, 1000, 1000))
    }, {
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(210,0,0), RGB(), RGB()))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(13, 99, 219), RGB(), RGB()))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(240,240,240)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240), RGB(), RGB(), RGB(/*1.0, 1.0, 1.0*/)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(240, 110, 203)))), // Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, Material(RGB(125,200,200)))), // Esfera de la derecha
    }, {
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0)),
    });

    // Soft shadows.
    Scene s6({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 90, Camera::HORIZONTAL, 3, 512, 1000, 1000))
    }, {
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(210,0,0), RGB(), RGB()))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(13, 99, 219), RGB(), RGB()))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(240,240,240)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240), RGB(), RGB(), RGB(/*1.0, 1.0, 1.0*/)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(240, 110, 203)))), // Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, Material(RGB(125,200,200)))), // Esfera de la derecha
    }, {
        Light(Vector3(-0.5,0.5,-0.5), RGB(1.0, 1.0, 1.0)),
        Light(Vector3(-0.5,0.5,0.5), RGB(1.0, 1.0, 1.0)),
        Light(Vector3(0.5,0.5,0.5), RGB(1.0, 1.0, 1.0)),
        Light(Vector3(0.5,0.5,-0.5), RGB(1.0, 1.0, 1.0)),
    });

    // Escena causticas:
    Scene s7({
        // Camaras.
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 100, Camera::HORIZONTAL, 3.5, 512, 1000, 1000)),
    }, {
        // Objetos.
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(13, 99, 219)))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(0,210,0)))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(210,0,0)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(140, 81, 38)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(0,-0.5,-0.5), 0.5, Material(RGB(), RGB(20,20,20), RGB(235,235,235), RGB(), 1.3))),
    }, {
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0)),
    });

    Scene s8({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0.75,-0.8, -0.3), -20, -80, -7, 130, Camera::HORIZONTAL, 4, 512, 1600, 900)),
    }, {
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(41, 154, 230), RGB(), RGB()))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(41, 154, 230), RGB(), RGB()))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(227, 174, 59)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240), RGB(), RGB(), RGB(/*1.0, 1.0, 1.0*/)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(41, 154, 230)))), //Plano de atrás
        std::make_shared<Sphere>(Sphere(Vector3(0, -0.7, 0), 0.3, Material(RGB(125,200,200)))), // Esfera de la derecha
    }, {
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0)),
    });

    // Materials scene.
    Scene s9({
        std::make_shared<Pinhole_camera>(Pinhole_camera(Vector3(0,0, -3.5), 100, Camera::HORIZONTAL, 3.6, 512, 1000, 1000)),
    }, {
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , Material(RGB(240,240,240), RGB(), RGB()))), //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(240,240,240), RGB(), RGB()))), //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , Material(RGB(107, 99, 67)))), //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240,240,240), RGB(), RGB(), RGB(/*1.0, 1.0, 1.0*/)))), //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(240,240,240)))), //Plano de atrás
        std::make_shared<Cube>(Cube(Vector3(-1,-1,-0.1), Vector3(1,-0.7,1), Material(RGB(107, 99, 67)))),
        std::make_shared<Sphere>(Sphere(Vector3(-0.65,-0.4,0.3), 0.3, Material(RGB(207, 219, 31)))), // Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(0,-0.4,0.3), 0.3, Material(RGB(200, 14, 28), RGB(55,55,55)))), // Esfera de la derecha
        std::make_shared<Sphere>(Sphere(Vector3(0.65,-0.4,0.3), 0.3, Material(RGB(), RGB(), RGB(), RGB(18, 196, 24)))), // Esfera de la derecha
        std::make_shared<Sphere>(Sphere(Vector3(-0.35, -0.7, -0.6), 0.3, Material(RGB(), RGB(255,255,255)))), // Esfera de la derecha
        std::make_shared<Sphere>(Sphere(Vector3(0.35, -0.7, -0.6), 0.3, Material(RGB(), RGB(45,45,45), RGB(210,210,210), RGB(), 1.35))), // Esfera de la derecha
    }, {
        Light(Vector3(0,0.5,0), RGB(1.0, 1.0, 1.0)),
    });


    int N_BOUNCES = 15;
    int N_THREADS = 8;
    int using_camera  = 0;
    Scene using_scene = s;
    std::string export_name = "./new_scene.ppm";
    Render path_tracer(N_BOUNCES, N_THREADS);
    path_tracer.render(using_scene, using_camera);
    path_tracer.export_render(using_scene, using_camera, export_name );
}