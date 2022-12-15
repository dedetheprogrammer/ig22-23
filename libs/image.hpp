#pragma once
#include <cassert>
#include <cmath>   
#include <iomanip>   
#include <iostream>  
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
    double R, G, B;
    RGB() : R(0.0), G(0.0), B(0.0) {}
    RGB(int C) : R(C/255.0), G(C/255.0), B(C/255.0) {}
    RGB(int R, int G, int B) : R(R/255.0), G(G/255.0), B(B/255.0) {}
    RGB(double C) : R(C), G(C), B(C) {}
    RGB(double R, double G, double B) : R(R), G(G), B(B) {}

    void operator=(const int i) {
        R = G = B = double(i)/255.0;
    }

    void operator=(const double d) {
        R = G = B = d;
    }

    void operator+=(const RGB& pixel) {
        R += pixel.R;
        G += pixel.G;
        B += pixel.B;
    }

    void operator*=(const RGB& pixel) {
        R *= pixel.R;
        G *= pixel.G;
        B *= pixel.B;
    }

    void operator*=(const double d) {
        R *= d;
        G *= d;
        B *= d;
    }

    void operator/=(const double d) {
        R /= d;
        G /= d;
        B /= d;
    }

    void operator^=(const double d) {
        R = pow(R, d);
        G = pow(G, d);
        B = pow(B, d);
    }
};

RGB operator+(const RGB& fst, const RGB& snd) {
    return RGB(fst.R + snd.R, fst.G + snd.G, fst.B + snd.B);
}

RGB operator*(const RGB& pixel, const double d) {
    return RGB(pixel.R * d, pixel.G * d, pixel.B * d);
}

RGB operator*(const RGB& fst, const RGB& snd) {
    return RGB(fst.R * snd.R, fst.G * snd.G, fst.B * snd.B);
}

RGB operator/(const RGB& pixel, const double d) {
    return RGB(pixel.R / d, pixel.G / d, pixel.B / d);
}

RGB operator^(const RGB& pixel, const double d) {
    return RGB(double(pow(pixel.R, d)), pow(pixel.G, d), pow(pixel.B, d));
}

bool operator==(const RGB& a, const RGB& b) {
    return a.R == b.R && a.G == b.G && a.B == b.B;
}

template <typename T>
bool operator>(const RGB& pixel, const T& d) {
    return pixe.R > d || pixel.G > d || pixel.B > d;
}

double max(const RGB& pixel) {
    double max = pixel.R;
    if (max < pixel.G) max = pixel.G;
    if (max < pixel.B) max = pixel.B;
    return max;
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

using Channels = std::vector<std::vector<RGB>>;

class Image {
private:
    
    void read_header(std::ifstream& in) {

        // Reading the netpbm file format. This should be P3.
        format = get_line(in);

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
    double maxval;       // Maximum value of the color resolution in memory.
    double memval;       // Maximum value of the color resolution generated in memory (after
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
    
    Image() {}

    Image(int width, int height, std::string name = "IMAGE") {
        this->format = "P3";
        this->name   = name;
        this->width  = width;
        this->height = height;
        this->pixels = Channels(height, std::vector<RGB>(width));
        this->maxval = this->memval = 0;
    }

    Image(std::string file) {
        // Opening the PPM file.
        std::ifstream in(file);
        assert(in.is_open() && "File not found.");

        // Reading PPM file header.
        has_color_key = false;
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

    std::vector<RGB>& operator[ ](const int& i) {
        return pixels[i];
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
    for (auto& h : i.pixels) {
        for (auto& w : h) {
            std::cout << w << "\t";
        }
        std::cout << std::endl;
    }
}



//===============================================================//
// Exporting an image.
//===============================================================//

void export_image(Image& i, std::string name, bool LDR = false, bool convert = false) {

    std::ofstream os(name);
    os << i.format << "\n";
    os << "#MAX=" << ((LDR || convert) ? 1.0 : i.memval) << "\n";
    os << "# " + get_filename(name)+ "\n";
    os << i.width << " " << i.height << "\n";
    os << (convert ? 255 : i.memres) << "\n";

    for (auto& h : i.pixels) {
        for (auto& w : h) {
            w = w * (convert && !LDR ? 255 : i.colres/i.maxval);
            os << static_cast<int>(w.R) << " " 
               << static_cast<int>(w.G) << " "
               << static_cast<int>(w.B) << "     ";
        }
        os << "\n";
    }
    os.close();
}

//===============================================================//
// Tone mapping
//===============================================================//

void tone_mapping(Image& i, int flags, double c_map = 1.0, double e_map = 1.0, double g_map = 2.2) {

    if (flags & clamp) {
        if (c_map > i.memval) return;
        i.memres = c_map * (double)i.colres/i.maxval;
        i.memval = c_map;
    }

    double e_max = i.memval;
    if (flags & (equalize | gamma)) {
        i.colres = i.memres = e_map * (double)i.colres/i.maxval;
        i.maxval = i.memval = e_map;
    }

    if (flags & gamma) {
        i.maxval = i.memval = std::pow(i.memval, 1/g_map);
        i.colres = i.memres = i.memval * (double)i.colres/i.memval;
    }

    for (auto& h : i.pixels) {
        for (auto& w : h) {
            if (flags & clamp) {
                if (w.R > c_map) w.R = c_map;
                if (w.G > c_map) w.G = c_map;
                if (w.B > c_map) w.B = c_map;
            }

            if (flags & (equalize | gamma)) {
                w *= (e_map/e_max);
            }

            if (flags & gamma) {
                w ^= (1/g_map);
            }
        }
    }
}