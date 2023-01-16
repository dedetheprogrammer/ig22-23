#pragma once
#ifndef CAMERA_H
#define CAMERA_H

#include "geometry.hpp"
#include "image.hpp"

//===============================================================//
// Generic camera
//===============================================================//
class Camera {
private:
    // Print camera attributes.
    virtual std::ostream& print(std::ostream& os) = 0;
public:

    // Indicates the orientation of the parameter. Used in things such as the 
    // FOV in the constructor.
    enum PARAMETER_ORIENTATION { HORIZONTAL, VERTICAL };

    // Canvas/Screen window. These are the coordinates of the "canvas" in the
    // image plane. These coordinates are computed from the canvas size and 
    // the film gate aspect ratio.
    Vector3 origin; // Camera position.
    Vector3 left;   // Camera width projection.
    Vector3 up;     // Camera height projection.
    Vector3 front;  // Camera depth projection.

    // Camera pixel base:
    int ppp;        // Points/paths per pixel.
    Vector3 p_ref;  // Pixel center reference.
    Vector3 p_left; // Pixel width.
    Vector3 p_up;   // Pixel height.

    // Defines the width in pixels of the output image. The image size also
    // defines the resolution gate aspect ratio.
    int image_width;     // Image width.
    // Defines the height in pixels of the output image. The image size also
    // defines the resolution gate aspect ratio.
    int image_height;    // Image height.
    // The ratio between the image width and its height (in pixels).
    double aspect_ratio; // Image aspect ratio.
    Image view;          // Camera image view.

    Camera() {}
    virtual ~Camera() = default; // Misma razon que Object, a estas horas no me
                                 // apetece discutir con el compilador.

    virtual inline Ray cast_ray(double x, double y) = 0;
    virtual inline Ray cast_antialiasing_ray(double x, double y) = 0;
    void export_view(std::string image_name = "./new_scene.ppm", int res = 10e8) {
        view.colres = view.memres = (view.maxval * res);
        export_image(view, image_name);
    }
    void update_view(int x, int y) {
        view[x][y] /= (double)ppp;
        if (view[x][y].R > view.maxval) view.maxval = view.memval = view[x][y].R;
        if (view[x][y].G > view.maxval) view.maxval = view.memval = view[x][y].G;
        if (view[x][y].B > view.maxval) view.maxval = view.memval = view[x][y].B;
    }

    friend std::ostream& operator<<(std::ostream& os, Camera& c) {
        return c.print(os);
    }
};

//===============================================================//
// Fisheye camera
//===============================================================//
class Fisheye_camera : public Camera {
private:
    std::ostream& print(std::ostream& os) override {
        return os << "Fisheye camera:"
            << "\n> View dimensions: " << image_width << "x" << image_height
            << "\n> Camera position: " << origin
            << "\n> Camera orientation: " << front
            << "\n> Film dimensions. Width: " << left << ", Height: " << up 
            << "\n> Points/paths per pixel: " << ppp;
    }
public:

    enum PROJECTION_TYPE { EQUIDISTANT, EQUISINUSOIDAL, EQUICUBIC };

    double FOV;
    PROJECTION_TYPE projection_type;

    Fisheye_camera(Vector3 origin, double FOV, int image_width, int image_height, double focal_length)
    {

        // Camera fisheye parameters:
        this->FOV = FOV * M_PI / 180;
        this->projection_type = projection_type;

        // Camera geometric base:
        this->origin = origin;
        //this->up = ...
        //this->left = ...

        // Camera pixel base:

        // Camera image view:
        this->image_width  = image_width;
        this->image_height = image_height;
        this->aspect_ratio = (double) image_width / (double) image_height;

    }

    inline Ray cast_ray(double x, double y) override {
        Vector3 ref = p_ref - (2*y*p_left) - (2*x*p_up);
        double theta = 2 * std::asin(ref.mod())/FOV;
        double phi   = std::atan2((ref - origin).mod(), ref.z);
        return Ray(origin, Vector3(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)));
    }

    inline Ray cast_antialiasing_ray(double x, double y) override {
        Vector3 ref = (p_ref - (2*y*p_left) - (2*x*p_up));
        std::uniform_real_distribution<> px(ref.x + p_left.x, ref.x - p_left.x);
        std::uniform_real_distribution<> py(ref.y - p_up.y, ref.y + p_up.y);
        ref = Vector3(px(e2), py(e2), ref.z);
        double theta = 2 * std::asin(ref.mod())/FOV;
        double phi   = std::atan2((ref - origin).mod(), ref.z);
        return Ray(origin, Vector3(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)));
    }


};

//===============================================================//
// Orthographic camera
//===============================================================//
class Orthographic_camera : public Camera {
private:
    std::ostream& print(std::ostream& os) override {
        return os << "Orthographic camera:"
            << "\n> View dimensions: " << image_width << "x" << image_height
            << "\n> Camera position: " << origin
            << "\n> Camera orientation: " << front
            << "\n> Film dimensions. Width: " << left << ", Height: " << up 
            << "\n> Points/paths per pixel: " << ppp;
    }
public:
    Orthographic_camera(Vector3 origin, Vector3 left, Vector3 up, Vector3 front,
        int ppp, int image_width, int image_height)
    {
        // Camera geometrical base:
        this->origin = origin;
        this->left   = left;
        this->up     = up;
        this->front  = front;

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;
        p_left = left/image_width;
        p_up   = up/image_height;
        p_ref  = origin + (up - p_up) + (left - p_left) + front;

        // Camera image view:
        this->image_width  = image_width;
        this->image_height = image_height;
        this->view = Image(image_width, image_height, "Camera view");

    }

    inline Ray cast_ray(double x, double y) override {
        return Ray(p_ref - (2*y*p_left) - (2*x*p_up), nor(front));
    }

    inline Ray cast_antialiasing_ray(double x, double y) override {
        Vector3 ref = (p_ref - (2*y*p_left) - (2*x*p_up));
        std::uniform_real_distribution<> px(ref.x + p_left.x, ref.x - p_left.x);
        std::uniform_real_distribution<> py(ref.y - p_up.y, ref.y + p_up.y);
        return Ray(Vector3(px(e2), py(e2), ref.z), nor(front));
    }
};

//===============================================================//
// Pinhole camera
//===============================================================//
class Pinhole_camera : public Camera {
private:
    std::ostream& print(std::ostream& os) override {
        return os << "Pinhole camera:"
            << "\n> View dimensions: " << image_width << "x" << image_height
            << "\n> Aspect ratio: " << aspect_ratio
            << "\n> Camera position: " << origin
            << "\n> Camera orientation: " << front
            << "\n> Film dimensions. Width: " << left << ", Height: " << up
            << "\n> Focal length: " << focal_length
            << "\n> Points/paths per pixel: " << ppp;
    }
public:

    // Defines the distance between the eye (the camera position), and the image
    // plane in a pinhole camera. This parameter is used to compute the angle of
    // view.
    double focal_length;

    // Define the width dimension of the film that would be used in a real
    // camera. The angle of view depends on this value. It also defines the 
    // film gate aspect ratio.
    //double film_aperture_width;
    // Define the height dimension of the film that would be used in a real
    // camera. The angle of view depends on this value. It also defines the 
    // film gate aspect ratio.
    //double film_aperture_height;

    // Near clipping plane is an imaginary plane located at a particular distance
    // from the camera along the its sight line. Only objects after the near 
    // clipping plane are rendered in the camera's view.
    //double near_clipping_plane;
    // Far clipping plane is an imaginary plane located at a particular distance
    // from the camera along the its sight line. Only objects before the far
    // clipping plane are rendered in the camera's view.
    //double far_clipping_plane;

    // The camera to world transformation matrix. Defines the camera position and
    // orientation.
    // Matrix3BaseChange camera_to_world;

    // The angle of view is computed from the focal length and the film size
    // parameters.
    // double FOVh;
    // double FOVv;

    /**
     * @brief Construct a new Pinhole_camera object. Harder (maybe inefficient)
     * way to configure the camera, you will have to have in count the image
     * ratio manually.
     * 
     * @param origin 
     * @param up 
     * @param left 
     * @param front 
     * @param ppp 
     * @param image_width 
     * @param image_height 
     */
    Pinhole_camera(Vector3 origin, Vector3 up, Vector3 left, Vector3 front, int ppp, 
        int image_width, int image_height)
    {
        // Camera geometrical base:
        this->origin = origin;
        this->left   = left;
        this->up     = up;
        this->front  = front;

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;
        p_left = left/image_width;
        p_up   = up/image_height;
        p_ref  = origin + (up - p_up) + (left - p_left) + front;

        // Camera image view:
        this->image_width  = image_width;
        this->image_height = image_height;
        this->view = Image(image_width, image_height, "Camera view");
    }

    /**
     * @brief Construct a new Pinhole_camera object. With this constructor, your
     * camera will be aligned with the Z-axis. 
     * 
     * @param origin 
     * @param FOV 
     * @param FOV_ORIENTATION 
     * @param focal_length 
     * @param image_width 
     * @param image_height 
     */
    Pinhole_camera(Vector3 origin, double FOV, PARAMETER_ORIENTATION FOV_ORIENTATION,
        double focal_length, int ppp, int image_width, int image_height)
    {

        this->focal_length = focal_length;
        this->image_width  = image_width;
        this->image_height = image_height;
        this->aspect_ratio = (double) image_width / (double) image_height;
        double half_FOV_len = tan((FOV * M_PI / 180.0)/2);
        
        this->origin = origin;
        if (FOV_ORIENTATION == HORIZONTAL) {
            this->left = Vector3(-half_FOV_len, 0, 0);
            this->up   = Vector3(0, tan(atan(half_FOV_len / aspect_ratio)), 0);
        } else if (FOV_ORIENTATION == VERTICAL) {
            this->up   = Vector3(0, half_FOV_len, 0);
            this->left = Vector3(-tan(atan(half_FOV_len / aspect_ratio)), 0, 0); 
        }
        this->front = Vector3(0, 0, focal_length/half_FOV_len);

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;
        p_left = left/image_width;
        p_up   = up/image_height;
        p_ref  = origin + (up - p_up) + (left - p_left) + front;

        // Camera image view:
        view = Image(image_width, image_height, "Camera view");

    }

    /**
     * @brief Construct a new Pinhole_camera object. With this constructor, your
     * camera will be aligned with the direction vector.
     * 
     * @param origin 
     * @param direction
     * @param FOV 
     * @param FOV_ORIENTATION 
     * @param focal_length 
     * @param image_width 
     * @param image_height 
     */
    Pinhole_camera(Vector3 origin, Vector3 orientation, double FOV, PARAMETER_ORIENTATION FOV_ORIENTATION,
        double focal_length, int ppp, int image_width, int image_height)
    {

        this->focal_length = focal_length;
        this->image_width  = image_width;
        this->image_height = image_height;
        this->aspect_ratio = (double) image_width / (double) image_height;
        double half_FOV_len = tan((FOV * M_PI / 180.0)/2);
        
        this->origin = origin;
        if (FOV_ORIENTATION == HORIZONTAL) {
            this->left = Vector3(-half_FOV_len, 0, 0);
            this->up   = Vector3(0, tan(atan(half_FOV_len / aspect_ratio)), 0);
        } else if (FOV_ORIENTATION == VERTICAL) {
            this->up   = Vector3(0, half_FOV_len, 0);
            this->left = Vector3(-tan(atan(half_FOV_len / aspect_ratio)), 0, 0); 
        }
        this->front = Vector3(0, 0, focal_length/half_FOV_len);

        // Appling the orientation into the camera dimensions.
        // These need to be normalized so the camera vectors won't deform.
        // New camera up dimension in world coordinates.
        Vector3 world_up    = nor(crs(this->left, orientation));
        // New camera front dimension in camera coordinates.
        Vector3 local_front = nor(orientation);
        // New camera left dimension in camera coordinates.
        Vector3 local_left  = nor(crs(world_up, local_front));
        // New camera up dimension in camera coordinates.
        Vector3 local_up    = nor(crs(local_front, local_left));
        // Matrix camera rotation.
        Matrix3BaseChange rotation(local_left, local_up, local_front, origin);
        
        // Updating camera dimension vectors.
        this->left  = rotation * this->left;
        this->up    = rotation * this->up;
        this->front = rotation * this->front;

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;
        p_left = left/image_width;
        p_up   = up/image_height;
        p_ref  = origin + (up - p_up) + (left - p_left) + front;

        // Camera image view:
        view = Image(image_width, image_height, "Camera view");

    }

    /**
     * @brief Construct a new Pinhole_camera object. With this constructor, your
     * camera will be rotated with the given angles (in degrees).
     * https://stackoverflow.com/a/36604847
     * 
     * @param origin 
     * @param pitch Rotation around the x-axis.
     * @param yaw   Rotation around the y-axis.
     * @param roll  Rotation around the z-axis.
     * @param FOV 
     * @param FOV_ORIENTATION 
     * @param focal_length 
     * @param image_width 
     * @param image_height 
     */
    Pinhole_camera(Vector3 origin, double pitch, double yaw, double roll,
        double FOV, PARAMETER_ORIENTATION FOV_ORIENTATION,
        double focal_length, int ppp, int image_width, int image_height)
    {

        this->focal_length = focal_length;
        this->image_width  = image_width;
        this->image_height = image_height;
        this->aspect_ratio = (double) image_width / (double) image_height;
        double half_FOV_len = tan((FOV * M_PI / 180.0)/2);
        
        this->origin = origin;
        if (FOV_ORIENTATION == HORIZONTAL) {
            this->left = Vector3(-half_FOV_len, 0, 0);
            this->up   = Vector3(0, tan(atan(half_FOV_len / aspect_ratio)), 0);
        } else if (FOV_ORIENTATION == VERTICAL) {
            this->up   = Vector3(0, half_FOV_len, 0);
            this->left = Vector3(-tan(atan(half_FOV_len / aspect_ratio)), 0, 0); 
        }
        this->front = Vector3(0, 0, focal_length/half_FOV_len);

        // Appling the orientation into the camera dimensions.
        // X-axis matrix rotation (pitch rotation).
        Matrix3Rotation pitch_rotation(X_ROT, pitch);
        // Y-axis matrix rotation (yaw rotation).
        Matrix3Rotation yaw_rotation(Y_ROT, yaw);
        // Z-axis matrix rotation (roll rotation).
        Matrix3Rotation roll_rotation(Z_ROT, roll);
        
        // Updating camera dimension vectors.
        this->left  = pitch_rotation * yaw_rotation * roll_rotation * this->left;
        this->up    = pitch_rotation * yaw_rotation * roll_rotation * this->up;
        this->front = pitch_rotation * yaw_rotation * roll_rotation * this->front;

        // Camera pixel base:
        if (ppp < 0) { ppp = 1; } this->ppp = ppp;
        p_left = left/image_width;
        p_up   = up/image_height;
        p_ref  = origin + (up - p_up) + (left - p_left) + front;

        // Camera image view:
        view = Image(image_width, image_height, "Camera view");

    }

    inline Ray cast_ray(double x, double y) override {
        return Ray(origin, nor((p_ref - (2*y*p_left) - (2*x*p_up)) - origin));
    }

    inline Ray cast_antialiasing_ray(double x, double y) override {
        Vector3 ref = (p_ref - (2*y*p_left) - (2*x*p_up));
        std::uniform_real_distribution<> px( ref.x + p_left.x, ref.x - p_left.x);
        std::uniform_real_distribution<> py( ref.y - p_up.y, ref.y + p_up.y);
        return Ray(origin, nor(Vector3(px(e2), py(e2), ref.z) - origin));
    }

};

/*
std::ostream& operator<<(std::ostream& os, const Camera& c) {
    return os << "> View dimensions: " << c.view.width << "x" << c.view.height
        << "\n> Aspect ratio: " << c.aspect_ratio
        << "\n> Camera position: " << c.center
        << "\n> Camera orientation: " << c.front
        << "\n> Film dimensions. Width: " << c.left << ", Height: " << c.up 
        << "\n> Points/paths per pixel: " << c.ppp;
}
*/

#endif