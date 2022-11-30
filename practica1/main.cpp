#include "geometry.hpp"
#include "planet.hpp"

void geometry_test(bool verbose = false) {
    // Vector3 tests.
    std::cout << "== PRUEBA DE VECTOR3 =================" << std::endl;
    float s = 0.5;
    Vector3 a(1.0,1.0,1.0), b(5.0,4.0,4.5), d(2.0, 2.0, 2.0), v(2.0, -2.0, 0.0), w = crs(d,v);

    if (verbose) { std::cout << "Params: " << a << ", " << b << ", " << s << std::endl; }
    if (verbose) { std::cout << "Add: " << a+b      << std::endl; } assert((a+b) == Vector3(6,5,5.5));
    if (verbose) { std::cout << "Sub: " << a-b      << std::endl; } assert((a-b) == Vector3(-4,-3,-3.5));
    if (verbose) { std::cout << "Mul: " << a*s      << std::endl; } assert((a*s) == Vector3(0.5,0.5,0.5));
    if (verbose) { std::cout << "Div: " << a/s      << std::endl; } assert((a/s) == Vector3(2,2,2));
    if (verbose) { std::cout << "Mod: " << a.mod()  << std::endl; } assert(a.mod() - sqrt(3) < EPSILON_ERROR);
    if (verbose) { std::cout << "Nor: " << nor(a)   << std::endl; } assert(nor(a) == (a/a.mod()));
    if (verbose) { std::cout << "Grd: " << a/b      << std::endl; }
    if (verbose) { std::cout << "Dot: " << a*b      << std::endl; } assert(a*b == 13.5);
    if (verbose) { std::cout << "Crs: " << crs(a,b) << std::endl; } assert(crs(a,b) == Vector3(0.5,0.5,-1));
    std::cout << "------> OK\n";

    // Matrices tests.
    std::cout << "== PRUEBA DE MATRICES =======================" << std::endl;
    float A[4][4] = {{3,4,7,8},{3,2,8,7},{0,2,6,7},{4,2,2,1}};
    float B[4][4] = {{8,7,6,5},{2,4,2,1},{1,1,3,4},{7,8,9,0}};

    Matrix3 mA(A);
    Matrix3 mB(B);
    Matrix3 mC = mA * mB;
    Matrix3 mD = mC.invert();
    Matrix3Translation mT(3,5,7);
    Matrix3Scale mS(2,1,2);
    Matrix3Rotation mR(90, X_ROT);
    Matrix3BaseChange mBC(d, v, w, a);
    std::cout << "Transform Matrix A: " << mA << std::endl;
    std::cout << "Transform Matrix B: " << mB << std::endl;
    std::cout << "Transform Matrix C = A * B: " << mC << std::endl;
    std::cout << "Inverse Matrix mD: " << mD << std::endl;
    std::cout << "Translation Matrix mT: " << mT << std::endl;
    std::cout << "Scaling Matrix mS: " << mS << std::endl;
    std::cout << "Rotation Matrix mR: " << mR << std::endl;
    std::cout << "Base change Matrix mBC: " << mBC << std::endl;

    // Transformations tests.
    Vector3 hp(b, 1);
    Vector3 hd(d, 0);
    std::cout << "== PRUEBA DE TRANSFORMACIONES ===============" << std::endl;
    std::cout << "== TRASLACION" << std::endl;
    std::cout << "Point traslation: " << mT * hp << std::endl;
    std::cout << "Vector traslation: " << mT * hd << std::endl;

    std::cout << "== ESCALADO" << std::endl;
    std::cout << "Point scaling: " << mS * hp << std::endl;
    std::cout << "Vector scaling: " << mS * hd << std::endl;
    
    std::cout << "==ROTATION" << std::endl;
    std::cout << "Point rotation: " << mR * hp << std::endl;
    std::cout << "Vector rotation: " << mR * hd << std::endl;

    std::cout << "== CAMBIO DE BASE" << std::endl;
    std::cout << "Point base change: " << mBC * hp << std::endl;
    std::cout << "(No sense) Vector base change: " << mBC * hd << std::endl;

    Vector3 hc = mT * mS * hp;
    std::cout << "Point traslation + rotation: " << hc << std::endl;
    std::cout << "Recover original point: " << (mT * mS).invert() * hc << std::endl;
}

void fatal_death(Planet planeta1, int inclination1, int azymuth1, Planet planeta2, int inclination2, int azymuth2) {
    Vector3 estacion1 = planeta1.getGlobalPoint(inclination1,azymuth1);
    Vector3 estacion2 = planeta2.getGlobalPoint(inclination2,azymuth2);
    std::cout << estacion1 << " " << estacion2 << std::endl;
    Vector3 AB(estacion2 - estacion1);
    Matrix3BaseChange p1Coor(planeta1.axis, planeta1.equator, planeta1.thirdAxis, planeta1.center),
        p2Coor(planeta2.axis, planeta2.equator, planeta2.thirdAxis, planeta2.center);

    Vector3 ABp1 = p1Coor*AB;
    Vector3 ABp2 = p2Coor*AB;

    // Si la dirección de lanzamiento transformada al sistema de coordenadas del planeta 1 tiene componente hy negativa,
    // significa geométricamente que la dirección de lanzamiento es hacia dentro del planeta, por lo que hay choque.
    if (ABp1.y < 0) {
        std::cout << "Ha habido choque con el planeta 1" << std::endl;
    }
    // Si la dirección de lanzamiento transformada al sistema de coordenadas del planeta 2 tiene componente hy positiva,
    // significa geométricamente que la dirección de lanzamiento proviene desde dentro del planeta hacia fuera, atravesandolo, por lo que hay choque.
    if (ABp2.y > 0) {
        std::cout << "Ha habido choque con el planeta 2" << std::endl;
    }
    // Si la dirección de lanzamiento tiene componente hy positiva vista desde el planeta 1, y negativa vista desde el planeta 2,
    // Significa que el lanzamiento se expulsa hacia fuera del planeta 1 y llega desde fuera del planeta 2, por lo que no hay choque.
    if (ABp1.y >= 0 && ABp2.y <= 0) {
        std::cout << "No ha habido choque" << std::endl;
    }
}

int main (int argc, char* argv[]) {
    // Hay error si metes Latitud y Azimuth 90º en ambos planetas. Por? Ni idea.
    // Tengo que probar el get global point pero con matriz de transformación.
    // Si hay sugerencia de alternativas de nombre de Matrix3, soy todo oidos.
    Planet A(Vector3(1,1,1), Vector3(0,0,0), Vector3(2,2,2)), 
        B(Vector3(1,1,1), Vector3(5,5,5), Vector3(7,7,7));
    float latA, aziA, latB, aziB;

    // PLANET A:
    std::cout << "Welcome to the FTL Enterprises Interplanetary Simulator" << std::endl;
    //std::cout << "Introduce Planet A {";
    //std::cin >> A; std::cout << "} " << std::endl;
    std::cout << "Introduce Station A latitude: "; std::cin >> latA;
    std::cout << "Introduce Station A azimuth: ";  std::cin >> aziA;
    std::cout << std::endl;

    // PLANET B:
    //std::cout << "Introduce Planet B {";
    //std::cin >> B; std::cout << "} " << std::endl;
    std::cout << "Introduce Station B latitude: "; std::cin >> latB;
    std::cout << "Introduce Station B azimuth: ";  std::cin >> aziB;
    std::cout << std::endl;

    // Trayectory Evaluation:
    fatal_death(A, latA, aziA, B, latB, aziB);
}