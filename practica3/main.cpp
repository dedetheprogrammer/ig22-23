#include <iostream>
#include "camera.hpp"
#include "objects.hpp"

int main (int argc, char* argv[]) {
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    /*
    std::vector<Vector3> b {
        Vector3(-0.5,-1,0.5),
        Vector3(-0.5,0.4,0.5),
        Vector3(0.5,0.4,-0.5),
        Vector3(0.5,-1,-0.5)
    };

    std::vector<Triangle> m {
        Triangle({Vector3(-0.5,0,0.5), Vector3(0,0.5,0), Vector3(0,-0.5,0)},  RGB(50,50,50)),
        Triangle({Vector3(0,-0.5,0), Vector3(0,0.5,0), Vector3(0.5,0,0.5)},   RGB(100,100,100)),
        Triangle({Vector3(-0.5,0,0.5), Vector3(0,0.5,0), Vector3(0.5,0,0.5)}, RGB(170,170,170)),
        Triangle({Vector3(-0.5,0,0.5), Vector3(0,-0.5,0), Vector3(0.5,0,0.5)}, RGB(230,230,230)),
    };

    //Mesh a(m);

    //Texture_ref r(b[0], X_ROT | Z_ROT, {-90, grd(nor(crs(b[1]-b[0], b[3]-b[0])), Vector3(0,0,-1))});
    
    // Ray R(Vector3(0,0,0), Vector3(1, 2, 3));
    // Scene O {
    //     std::make_shared<Plane>(Plane(Vector3(1,0,0), 1)),
    //     std::make_shared<Triangle>(Triangle(Vector3(-1, -0.5, -0.5), Vector3(-1, -0.5, 0.5), Vector3(-1, 0.5, 0))),
    //     std::make_shared<Sphere>(Sphere(Vector3(1, 0.4, 0), 0.5)),
    //     std::make_shared<Box>(Box(Vector3(3, -1, -1), Vector3(6, 1, 1)))
    // };

    // for (auto& o : scene) {
    //     std::cout << *o << std::endl;
    // }

    // for (auto& o : O) {
    //     std::vector<Vector3> points;
    //     std::cout << *o.get() << "\n";
    //     std::cout << "\tR: " << R << " intersecta en ";
    //     for (auto& p : o.get()->intersects(R)) {
    //         std::cout << p << " ";
    //     }
    //     std::cout << "\n";
    // }

    */
   
    Camera C(Vector3(0,0,-3.5), Vector3(-1,0,0), Vector3(0,1,0), Vector3(0,0,3), 256, 256, 1);
    Scene scene {
        // El main escenario.
        std::make_shared<Plane>(Plane(1, Vector3(1,0,0) , RGB(200,0,0))),    //Plano de la izquierda
        std::make_shared<Plane>(Plane(1, Vector3(-1,0,0), RGB(0,200,0))),   //Plano de la derecha
        std::make_shared<Plane>(Plane(1, Vector3(0,1,0) , RGB(204,200,200))),    //Plano de abajo
        std::make_shared<Plane>(Plane(1, Vector3(0,-1,0), RGB(200,204,200))),   //Plano de arriba
        std::make_shared<Plane>(Plane(1, Vector3(0,0,-1), RGB(200,200,179))),     //Plano de atrás
        //std::make_shared<Mesh>(Mesh("./obj/cow.ply") * Matrix3Scale(0.4) * Matrix3Translation(0, -0.6, 0) * Matrix3Rotation(Y_ROT, 45)),

        //std::make_shared<Mesh>(a * Matrix3Rotation(Y_ROT, -90) * Matrix3Rotation(X_ROT, -45)),
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,0.3,-0.6), 0.3, RGB(255,0,200))),  //Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(-0.5,-0.7,0.25), 0.3, RGB(254,125,254))),  //Esfera de la izquierda
        std::make_shared<Sphere>(Sphere(Vector3(0.5,-0.7,-0.25), 0.3, RGB(125,254,254)))  //Esfera de la derecha
        // std::make_shared<Plane>(Plane(b)),
        
        // Pruebas texturas:
        // std::make_shared<Plane>(Plane(1, Vector3( 1, 0, 0), Texture(Image("doom_wall-3.ppm")) , 0.5, 0.5, Z_ROT)),    //Plano de la izquierda
        // std::make_shared<Plane>(Plane(1, Vector3(-1, 0, 0), Texture(Image("doom_wall-3.ppm")) , 0.5, 0.5, Z_ROT, 1)),   //Plano de la derecha
        // std::make_shared<Plane>(Plane(1, Vector3( 0, 1, 0), Texture(Image("doom_floor-3.ppm")), 0.5, 0.5, X_ROT)),    //Plano de abajo
        // std::make_shared<Plane>(Plane(1, Vector3( 0,-1, 0), Texture(Image("doom_sky-2.ppm"))  , Vector3(2,1,3), 4, 6, X_ROT)),   //Plano de arriba
        // std::make_shared<Plane>(Plane(3, Vector3( 0, 0,-1), Texture(Image("doom_wall-3.ppm")) , 0.5, 0.5, X_ROT)),     //Plano de atrás
        //std::make_shared<Plane>(Plane(b, Texture(Image("big_demon_sprite.ppm")), r)),
    };

    C.rayTrace(scene);
    float colres  = 0;
    auto pixels = C.getRGB(colres);
    Image i(1.0, colres * 255, pixels);
    export_image(i, "./new_scene.ppm");

    std::cout << "Time cost = " <<
        std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::steady_clock::now() - begin).count() << "[ms]" << std::endl;

}