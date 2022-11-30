#pragma once
#include <cassert>
#include <cmath>    
#include <sstream>
#include "utils.hpp"

#define clamp 0x01
#define equalize 0x02
#define gamma 0x04

//===============================================================//
// RGB
//===============================================================//

class RGB {
private:
    // ...
public:
    float R, G, B;
    RGB(int R, int G, int B) : R(R/255.0), G(G/255.0), B(B/255.0) {}
    RGB(float R = 0.0, float G = 0.0, float B = 0.0) : R(R), G(G), B(B) {}

    void operator+=(const RGB& pixel) {
        R += pixel.R;
        G += pixel.G;
        B += pixel.B;
    }

    void operator/=(const float d) {
        R /= d;
        G /= d;
        B /= d;
    }
};

RGB operator+(const RGB& fst, const RGB& snd) {
    return RGB(fst.R + snd.R, fst.G + snd.G, fst.B + snd.B);
}

RGB operator*(const RGB& pixel, const float d) {
    return RGB(pixel.R * d, pixel.G * d, pixel.B * d);
}

RGB operator*(const RGB& fst, const RGB& snd) {
    return RGB(fst.R * snd.R, fst.G * snd.G, fst.B * snd.B);
}

RGB operator/(const RGB& pixel, const float d) {
    return RGB(pixel.R / d, pixel.G / d, pixel.B / d);
}

bool operator==(const RGB& a, const RGB& b) {
    return a.R == b.R && a.G == b.G && a.B == b.B;
}

void operator>>(std::istream& is, RGB& pixel) {
    is >> pixel.R >> pixel.G >> pixel.B;
}

std::ostream& operator<<(std::ostream& os, const RGB& pixel) {
    return os << pixel.R << " " << pixel.G << " " << pixel.B;
}

//===============================================================//
// Image
//===============================================================//

typedef std::vector<std::vector<RGB>> Channels;

class Image {
private:
    
    void read_header(std::ifstream& in) {

        // Reading the netpbm file format. This should be P3.
        format = get_line(in);
        assert(!format.compare("P3") && "unrecognized format, expected P3 format...\n");

        // Reading the maximum value and the color key, if there's any. Initially the 
        // memory value is equal to the maximum.
        std::string s; getline(in, s);
        std::vector<std::string> ts = tokenize(s, ",");
        maxval = memval = std::stof(
            std::regex_replace(ts[0], std::regex("[^0-9.]"),"").c_str()
        );
        if (ts.size() > 1) {
            ts = tokenize(ts[1], " ");
            color_key = RGB(
                std::stof(ts[0].c_str()),
                std::stof(ts[1].c_str()),
                std::stof(ts[2].c_str())
            );
            has_color_key = true;
        } else {
            has_color_key = false;
        }

        // Reading the file name.
        name = get_line(in).substr(2);

        // Reading image width, height and color resolution. Initially, the memory resolution
        // is equal to the color.
        in >> width >> height >> colres;
        memres = colres;
    }

public:

    std::string format; // Netpbm file format. The PPM format is P3 if ASCII, or P6 if binary:
                        //      https://netpbm.sourceforge.net/doc/index.html
                        //      https://netpbm.sourceforge.net/doc/ppm.html#index
    std::string name;   // Name of the file.
    float maxval;       // Maximum value of the color resolution in memory.
    float memval;       // Maximum value of the color resolution generated in memory (after
                        //      tone mapping operations).
    int   width;        // Width of the ppm file.    
    int   height;       // Height of the ppm file.
    int   colres;       // Maximum color resolution of the ppm file.
    int   memres;       // Maximum color resolution of the ppm file that was generated in memory
                        //      (after tone mapping operations).
    Channels pixels;    // Image pixels data. A channel is conformed by Red, Green and Blue values 
                        //      (0..255, 0..65535). Alpha channel not supported.
    RGB color_key;      // As the alpha channel is not supported, the color key indicates which 
                        //      RGB value will representate transparency.
    bool has_color_key; // Indicates if the ppm file has a color key, necessary because the program
                        //      could have determined RGB(0,0,0) as the color key.
    bool data_per_line;
    
    Image() {}

    Image(float maxval, int colres, Channels pixels) : format("P3"), pixels(pixels) {
        this->name   = "IMAGE";
        this->width  = pixels[0].size();
        this->height = pixels.size();
        this->maxval = memval = maxval;
        this->colres = memres = colres;
    }

    Image(std::string file) {
        // Opening the PPM file.
        std::ifstream in(file);
        assert(in.is_open() && "File not found.");

        // Reading PPM file header.
        read_header(in);

        // Reading PPM pixels data.
        pixels = std::vector<std::vector<RGB>>(height, std::vector<RGB>(width));
        for (auto& h : pixels) {
            for (auto& w : h) {
                in >> w;
                w = (w * maxval) / colres;
            }
        }
        in.close();
    }
};

std::ostream& operator<<(std::ostream& os, const Image& i) {
    return os << i.name << " {"
        << "\n  Format ............. " << i.format
        << "\n  Max value .......... " << i.memval
        << "\n  Resolution ......... " << i.width << "x" << i.height
        << "\n  Color resolution ... " << i.memres
        << "\n  Color key .......... " << i.color_key
        << "\n}";
}

//===============================================================//
// Exporting an image.
//===============================================================//

void export_image(Image& i, std::string path, bool conversion = 0) {   
    std::ofstream os(path);
    os << i.format << "\n";
    os << "#MAX=" << i.memval << "\n";
    os << "# " + get_filename(path)+ "\n";
    os << i.width << " " << i.height << "\n";
    os << static_cast<int>((conversion ? 255 : i.memres)) << "\n";

    for (auto& h : i.pixels) {
        for (auto& w : h) {
            w = w * (conversion ? 255 : i.colres/i.maxval);
            os << static_cast<int>(std::ceil(w.R)) << " " 
               << static_cast<int>(std::ceil(w.G)) << " "
               << static_cast<int>(std::ceil(w.B)) << "     ";
        }
        os << "\n";
    }
    os.close();
}

//===============================================================//
// Tone mapping
//===============================================================//

struct Tone {
    int t_flags = 0;
    float c_map, e_map, g_map;

    Tone() : t_flags(0), c_map(0), e_map(0), g_map(0) {}
    Tone(int t_flags, float c_map, float e_map, float g_map)
        : t_flags(t_flags), c_map(c_map), e_map(e_map), g_map(g_map) {}
};

void tone_mapping(Image& i, Tone p) {
    float max = i.memval;
    if ((p.t_flags & clamp) && (i.memval > p.c_map)) {
        if (p.c_map >= 1.f) i.memval = p.c_map;
        i.memres = p.c_map * i.colres/i.maxval;
        max = p.c_map;
    }

    if (p.t_flags & equalize) {
        i.memres = i.colres = p.e_map * i.colres/i.maxval;
        i.memval = i.maxval = p.e_map;
    }

    for (auto& h : i.pixels) {
        for (auto& w : h) {
            if (p.t_flags & clamp) {
                if (w.R > p.c_map) w.R = p.c_map;
                if (w.G > p.c_map) w.G = p.c_map;
                if (w.B > p.c_map) w.B = p.c_map;
            }

            if (p.t_flags & (equalize | gamma)) {
                w = w * p.e_map/max;
            }

            if (p.t_flags & gamma) {
                w.R = pow(w.R, 1/p.g_map);
                w.G = pow(w.G, 1/p.g_map);
                w.B = pow(w.B, 1/p.g_map);
            }
        }
    }
}