#pragma once
#include <iostream>
#include "geometry.hpp"

class Planet {
private:

    Vector3 surface_point(const Vector3 v) {
        return (abs(axis.mod()/2 - v.mod()) < EPSILON_ERROR)
            ? v : v * (axis/2).mod()/v.mod();
    }

public:
    Vector3 center;    // Planet center.
    Vector3 reference; // Planet point-city reference.
    Vector3 axis;      // Planet axis.
    Vector3 equator;   // Planet equator.
    Vector3 cardinal;  // Planet third axis 

    // Default constructor.
    Planet() {}
    // Assign center, city reference and axis of the planet.
    Planet(Vector3 axis, Vector3 center, Vector3 reference) {

        // Geometrical things:
        this->center    = center;
        this->reference = reference;
        this->axis      = axis;

        // Planet equator:
        Vector3 s = surface_point(reference - center);
        equator = crs(s, crs(axis/2, s));
        equator = nor(equator, axis.mod()/2);

        // Planet third axis:
        cardinal = crs(axis/2, equator);
        cardinal = nor(cardinal, axis.mod()/2);
    }

    Vector3 get_local_point(Vector3 c_coord) {

        // Fixed surface coordinate.
        Vector3 s = surface_point(c_coord - center);
        // Latitude.
        double lat = acos(s.z/s.mod()) * 180/M_PI;
        // Azimuth.
        double azi = acos(s.x / (sqrt(s.x * s.x + s.y * s.y))) * 180/M_PI;
        return Vector3(s.mod(), lat, azi);
    }

    Vector3 get_local_point(double inclination, double azimuth) {
        return Vector3(equator.mod(), inclination, azimuth);
    }

    Vector3 get_global_point(double latitude, double azimuth) {
        double lat_r = latitude * M_PI/180, azi_r = azimuth * M_PI/180, r = (axis/2).mod();
        //return MatrixBaseChange3D(axis, equator, crs(axis, equator), center) * Vector3();
        return Vector3(r * sin(lat_r) * cos(azi_r), r * sin(lat_r) * sin(azi_r), r * cos(lat_r)) + center;
    }
};

std::ostream& operator<<(std::ostream& os, const Planet& s) {
    return os << "PLANET {"<< std::endl
        << "  Center:    " << s.center << std::endl
        << "  Axis:      " << s.axis << std::endl
        << "  Radius:    " << (s.axis/2).mod() << std::endl
        << "  Greenwich: " << s.reference << std::endl
        << "  Equator:   " << s.equator << std::endl
        << "}";
}
