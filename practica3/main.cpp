#include <iostream>
#include "scene.hpp"

int main (int argc, char* argv[]) {

    //Texture_ref r(b[0], X_ROT | Z_ROT, {-90, grd(nor(crs(b[1]-b[0], b[3]-b[0])), Vector3(0,0,-1))});
    Scene s{{
        // Camaras.
        Camera(Camera_settings({Vector3(0,0,-3.5), Vector3(-1,0,0), Vector3(0,1,0), Vector3(0,0,3)}, 512), 512, 512),
    }, {
        // Objetos.
        // Cornell's box:
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0),  Material(RGB(210,0,0)))), // Plano de la izquierda.
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), Material(RGB(0,210,0)))), // Plano de la derecha.
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0),  Material(RGB(240)))), // Plano de abajo.
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), Material(RGB(240)))), // Plano de arriba.
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), Material(RGB(200)))), // Plano de atrás.
        
        // Prinitives test:
        std::make_shared<Triangle>(Triangle({Vector3(-0.95, -0.5, -0.5), Vector3(-0.95, -0.5, 0.5), Vector3(-0.95, 0.5, 0)}, Material(RGB(200,150,0)))),
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,0.3,-0.6),  0.3, Material(RGB(80,67,145)))), // Una primera esfera.
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, Material(RGB(220,145,220)))), // Esfera de la izquierda.
        std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, Material(RGB(125,200,200)))), // Esfera de la derecha.
        std::make_shared<Cube>(Cube(Vector3(0.45,0.45,0.45), Vector3(0.95,0.95,0.95), Material(RGB(200,0,60)))),

        // 3d meshes test:
        //std::make_shared<Mesh>(Mesh({
        //    Triangle({Vector3(-0.5,0,0.5), Vector3(0,0.5,0),  Vector3(0,-0.5,0)}),
        //    Triangle({Vector3(0,-0.5,0),   Vector3(0,0.5,0),  Vector3(0.5,0,0.5)}),
        //    Triangle({Vector3(-0.5,0,0.5), Vector3(0,0.5,0),  Vector3(0.5,0,0.5)}),
        //    Triangle({Vector3(-0.5,0,0.5), Vector3(0,-0.5,0), Vector3(0.5,0,0.5)}),
        //}) * Matrix3Rotation(Y_ROT, -90) * Matrix3Rotation(X_ROT, -45)),
        //std::make_shared<Mesh>(Mesh("../src/cow.ply") * Matrix3Scale(0.4) * Matrix3Translation(0, -0.6, 0) * Matrix3Rotation(Y_ROT, 45)),

        // Finite planes and textures test:
        //std::make_shared<Plane>(Plane({
        //    Vector3(-0.5,-1,0.5),
        //    Vector3(-0.5,0.4,0.5),
        //    Vector3(0.5,0.4,-0.5),
        //    Vector3(0.5,-1,-0.5)
        //})),
        // std::make_shared<Plane>(Plane(1, Vector3( 1, 0, 0), Texture(Image("doom_wall-3.ppm")) , 0.5, 0.5, Z_ROT)),    //Plano de la izquierda
        // std::make_shared<Plane>(Plane(1, Vector3(-1, 0, 0), Texture(Image("doom_wall-3.ppm")) , 0.5, 0.5, Z_ROT, 1)),   //Plano de la derecha
        // std::make_shared<Plane>(Plane(1, Vector3( 0, 1, 0), Texture(Image("doom_floor-3.ppm")), 0.5, 0.5, X_ROT)),    //Plano de abajo
        // std::make_shared<Plane>(Plane(1, Vector3( 0,-1, 0), Texture(Image("doom_sky-2.ppm"))  , Vector3(2,1,3), 4, 6, X_ROT)),   //Plano de arriba
        // std::make_shared<Plane>(Plane(3, Vector3( 0, 0,-1), Texture(Image("doom_wall-3.ppm")) , 0.5, 0.5, X_ROT)),     //Plano de atrás
        //std::make_shared<Plane>(Plane(b, Texture(Image("big_demon_sprite.ppm")), r)),

    }};

    s.render();
    s.export_render();

}
