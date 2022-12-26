#include "geometry.hpp"
#include "planet.hpp"

int main (int argc, char* argv[]) {

    std::cout << "Welcome to the FTL Enterprises Interplanetary Simulator\n";
    std::cout << "First, some test to check the system...\n";

    std::cout << "= GEOMETRY TESTS =======================\n";
    // Vector3 evaluation.
    std::cout << "== VECTOR3 TESTS =======================\n";
    double s = 0.5;
    Vector3 a(1.0,1.0,1.0), b(5.0,4.0,4.5), d(2.0, 2.0, 2.0), v(2.0, -2.0, 0.0), w = crs(d,v);
    std::cout << "Params: " << a << ", " << b << " and the scalar: " << s << "\n";
    std::cout << "Add: " << a+b      << std::endl;
    std::cout << "Sub: " << a-b      << std::endl;
    std::cout << "Mul: " << a*s      << std::endl;
    std::cout << "Div: " << a/s      << std::endl;
    std::cout << "Mod: " << a.mod()  << std::endl;
    std::cout << "Nor: " << nor(a)   << std::endl;
    std::cout << "Grd: " << a/b      << std::endl;
    std::cout << "Dot: " << a*b      << std::endl;
    std::cout << "Crs: " << crs(a,b) << std::endl;

    // Matrix3 evaluation.
    std::cout << "== MATRIX3 TESTS =======================\n";
    double ma[4][4] = {{3,4,7,8},{3,2,8,7},{0,2,6,7},{4,2,2,1}};
    double mb[4][4] = {{8,7,6,5},{2,4,2,1},{1,1,3,4},{7,8,9,0}};
    Matrix3 mA(ma);
    Matrix3 mB(mb);
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

    // Matrix3 transforms evaluation.
    std::cout << "== MATRIX3 TRANSFORMATIONS TESTS =======\n";
    Vector3 hp(b, 1);
    Vector3 hd(d, 0);
    std::cout << "== TRASLATION" << std::endl;
    std::cout << "Point traslation: " << mT * hp << std::endl;
    std::cout << "Vector traslation: " << mT * hd << std::endl;

    std::cout << "== SCALING" << std::endl;
    std::cout << "Point scaling: " << mS * hp << std::endl;
    std::cout << "Vector scaling: " << mS * hd << std::endl;
    
    std::cout << "== ROTATION" << std::endl;
    std::cout << "Point rotation: " << mR * hp << std::endl;
    std::cout << "Vector rotation: " << mR * hd << std::endl;

    std::cout << "== BASE CHANGE" << std::endl;
    std::cout << "Point base change: " << mBC * hp << std::endl;
    std::cout << "(No sense) Vector base change: " << mBC * hd << std::endl;

    Vector3 hc = mT * mS * hp;
    std::cout << "Point traslation + rotation: " << hc << std::endl;
    std::cout << "Recover original point: " << (mT * mS).invert() * hc << std::endl;

    // Planet geometry evaluation:
    std::cout << "== PLANET GEOMETRY TESTS ===============\n";
    Planet pA(Vector3(5.0, 2.0, 0.0), Vector3(1.0,1.0,1.0), Vector3(2.0, 3.0, 1.0));
    double a_lat = 45, a_azi = a_lat;
    Vector3 lpa = pA.get_local_point(a_lat, a_azi), cpa = pA.get_global_point(45, 45);
    std::cout << pA
        << "\nLatitude: " << a_lat << "º, Azimuth: " << a_azi << "º"
        << "\nPolar Coords: " << lpa
        << "\nCarts Coords: " << cpa << "\n\n";
    
    Planet pB(Vector3(1,2,2)*2, Vector3(0,0,0), Vector3(2.0,3.0,1.0));
    Vector3 o(1,2,2), lp = pB.get_local_point(o), cp = pB.get_global_point(lp.y, lp.z);
    std::cout << pB 
        << "\nOrign Coords: " << o
        << "\nPolar Coords: " << lp
        << "\nCarts Coords: " << cp << "\n\n";

    // Trajectory Evaluation:
    std::cout << "== PLANET COLLISION ====================\n";
    // Latitudes and azimthus in degrees:
    double lat_A = 80;  // Planet A latitude.
    double azi_A = 160; // Planet A azimuth.
    double lat_B = 45;  // Planet B latitude.
    double azi_B = 35;  // Planet B azimuth.
    // Planets:
    Planet A(Vector3(1,1,1), Vector3(0,0,0), Vector3(2,2,2)); // Planet A.
    Planet B(Vector3(1,1,1), Vector3(5,5,5), Vector3(7,7,7)); // Planet B.
    std::cout << "Parameters:\n"
        << A << "\n"
        << "Planet A latitude: " << lat_A << " and azimuth: " << azi_A << "\n"
        << B << "\n"
        << "Planet B latitude: " << lat_B << " and azimuth: " << azi_B << "\n";

    Vector3 station_A = A.get_global_point(lat_A,azi_A);
    Vector3 station_B = B.get_global_point(lat_B,azi_B);
    Vector3 traject(station_B - station_A);
    std::cout << "Station A: " << station_A 
        << ", station B: " << station_B 
        << ", trajectory: " << traject << "\n";
    
    Matrix3BaseChange coor_A(A.axis, A.equator, A.cardinal, A.center);
    Matrix3BaseChange coor_B(B.axis, B.equator, B.cardinal, B.center);
    Vector3 local_A = coor_A*traject;
    Vector3 local_B = coor_B*traject;

    // If the y component of the Planet A local coordinate system is negative,
    // that means that the launch direction is towards the planet, then boom.
    if (local_A.y < 0) {
        std::cout << "Ha habido choque con el planeta 1" << std::endl;
    }
    // If the y component of the Planet B local coordinates system is positive,
    // that means that the launch direction comes from inside the planet, then boom.
    if (local_B.y > 0) {
        std::cout << "Ha habido choque con el planeta 2" << std::endl;
    }
    // If both of above are false, then it means that the launch was succesful
    // and there won't be any boom.
    if (local_A.y >= 0 && local_B.y <= 0) {
        std::cout << "No ha habido choque" << std::endl;
    }
}