#pragma once
#include "geometry.hpp"
#include <iostream>

class Planet {
private:

    Vector3 surface_point(const Vector3 v) {
        return (abs(axis.mod()/2 - v.mod()) < EPSILON_ERROR)
            ? v : v * (axis/2).mod()/v.mod();
    }

public:
    Vector3 center, reference;
    Vector3 axis, equator, thirdAxis;

    Planet() : center(Vector3()), reference(Vector3()), axis(Vector3()){}

    Planet(Vector3 axis, Vector3 center, Vector3 reference) : center(center), reference(reference), axis(axis) {
        Vector3 v = surface_point(reference - center);

        equator = crs(v, crs(axis/2, v));
        equator = equator * (axis/2).mod() / equator.mod();

        thirdAxis = crs(axis/2, equator);
        thirdAxis = thirdAxis * (axis/2).mod() / thirdAxis.mod();
    }

    Vector3 getLocalPoint(Vector3 cart_coord) {
        Vector3 v = surface_point(cart_coord - center);

        float lat = acos(v.z/v.mod()) * 180/M_PI;
        float azi = acos(v.x / (sqrt(v.x * v.x + v.y * v.y))) * 180/M_PI;
        return Vector3(v.mod(), lat, azi);
    }

    Vector3 getLocalPoint(float inclination, float azimuth) {
        return Vector3(equator.mod(), inclination, azimuth);
    }

    Vector3 getGlobalPoint(float latitude, float azimuth) {
        float lat_r = latitude * M_PI/180, azi_r = azimuth * M_PI/180, r = (axis/2).mod();
        //return MatrixBaseChange3D(axis, equator, crs(axis, equator), center) * Vector3();
        return Vector3(r * sin(lat_r) * cos(azi_r), r * sin(lat_r) * sin(azi_r), r * cos(lat_r)) + center;
    }
};

std::ostream& operator<<(std::ostream& os, const Planet& s) {
    return os << "SPHERE {"<< std::endl
        << " CENTER:    " << s.center << std::endl
        << " AXIS:      " << s.axis << std::endl
        << " RADIUS:    " << (s.axis/2).mod() << std::endl
        << " GREENWICH: " << s.reference << std::endl
        << " EQUATOR:   " << s.equator << std::endl
        << "}";
}

std::istream& operator>> (std::istream& is, Planet& s) {
    Vector3 a, c, r;

    std::cout << "\n  axis (";
    is >> a; std::cout << ", ";    
    std::cout << "center (";
    is >> c; std::cout << ", ";
    std::cout << "reference (";
    is >> r; std::cout << std::endl;
    s = Planet(a, c, r);
    return is;
}

void sphere_test() {
    Planet A(Vector3(5.0, 2.0, 0.0), Vector3(1.0,1.0,1.0), Vector3(2.0, 3.0, 1.0));
    float a_lat = 45, a_azi = a_lat;
    Vector3 lpa = A.getLocalPoint(a_lat, a_azi), cpa = A.getGlobalPoint(45, 45);
    std::cout << A << std::endl;
    std::cout << "Latitude: " << a_lat << "º, Azimuth: " << a_azi << "º" << std::endl;
    std::cout << "Polar Coords: " << lpa << std::endl;
    std::cout << "Carts Coords: " << cpa << std::endl << std::endl;
    
    Planet B(Vector3(1,2,2)*2, Vector3(0,0,0), Vector3(2.0,3.0,1.0));
    std::cout << B << std::endl;
    Vector3 o(1,2,2), lp = B.getLocalPoint(o), cp = B.getGlobalPoint(lp.y, lp.z);
    std::cout << "Orign Coords: " << o << std::endl;
    std::cout << "Polar Coords: " << lp << std::endl;
    std::cout << "Carts Coords: " << cp << std::endl << std::endl;
}