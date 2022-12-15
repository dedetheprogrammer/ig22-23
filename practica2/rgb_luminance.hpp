#pragma once
#include <iostream>
#include <math.h>
#include <string>

/* LUMINANCE RELATED FUNCTIONS */
/**
 * Conversion from RGB to Photometric/Digital ITU BT.709.
 * Source: https://en.wikipedia.org/wiki/Relative_luminance
 * 
 * @param pixel RGB tuple to get the luminance.
 * @return The RGB relative luminance.
 */
double StandardLum(RGB pixel) {
    return (0.2126 * pixel.R) + (0.7152 * pixel.G) + (0.0722 * pixel.B);
}

RGB StandardRGB (RGB pixel, double lum) {
    return RGB(
        (lum - 0.7152 * pixel.G - 0.0722 * pixel.B) / 0.2126,
        (lum - 0.2126 * pixel.R - 0.0722 * pixel.B) / 0.7152,
        (lum - 0.2126 * pixel.R - 0.7152 * pixel.G) / 0.0722
    );
}

/**
 * Conversion from RGB to Digital ITU BT.601.
 * Source: https://www.w3.org/TR/AERT/#color-contrast
 * 
 * @param pixel RGB tuple to get the luminance.
 * @return The RGB color-contrast luminance.
 */
double ContrastLum(RGB pixel) {
    return (0.299 * pixel.R) + (0.587 * pixel.G) + (0.114 * pixel.B);
}

/**
 * Conversion from RGB to HSP luminance.
 * Source: https://alienryderflex.com/hsp.html
 * 
 * @param pixel RGB tuple to get the luminance.
 * @return The RGB HSP luminance.
 */
double HSPLum(RGB pixel) {
    return sqrt(
        pow(0.299 * pixel.R, 2) + 
        pow(0.587 * pixel.G, 2) + 
        pow(0.114 * pixel.B, 2)
    );
}

#pragma once

#include <cmath>
#include "RGB.h"
#include "image.h"



void debug(Image image/*, double high*/) {

    RGB p(15989, 22769, 17304);
    double lum = StandardLum(p);
    std::cout << "Orginal: " << p;
    std::cout << "Luminance: " << lum << std::endl;
    std::cout << "Obtained: " << StandardRGB(p, lum) << std::endl;

    for (int h = 0; h < image.height; h++) {
        for (int w = 0; w < image.width; w++) {
            double lum = StandardLum(image.pixels[h][w]);
            image.pixels[h][w] = StandardRGB(image.pixels[h][w], lum);
        }
    }
}


void clamp(Image& image, int low = 0, int high = 255) {
    for (auto& h : image.pixels) {
        for (auto& w : h) {
            if (w.R < low) w.R = low;
            else if (w.R > high) w.R = high;
            if (w.G < low) w.G = low;
            else if (w.G > high) w.G = high;
            if (w.B < low) w.B = low;
            else if (w.B > high) w.B = high;
        }
    }
    if (high < image.color_res) image.color_res = high;
}

void equalization(Image& image, int max) {
    std::cout << "Hola: " << image.color_res << ", " << max << std::endl;
    for (auto& h : image.pixels) {
        for (auto& w : h) {
            w = (w * max)/image.color_res;
        }
    }
    image.color_res = max;
}

void gamma(Image& image, int gamma) {
    for (auto &w : image.pixels) {
        for (auto &h : w) {
            h.R = pow(h.R, gamma);
            h.G = pow(h.G, gamma);
            h.B = pow(h.B, gamma);
        }
    }
    image.color_res = pow(image.color_res, gamma);
}


/*
void clamp(HDRImage& hdr_img, double low, double high) {
    // Because of the MAX_VALUE parameter, our pixel values are MAX_VALUE times bigger
    // than the value stored, which matches the range [0..color_res] thus making the
    // calculus more confusing, so we are going to assume our low and high are on the
    // color range scale and we're going to adjust them to the memory range scale, this
    // will be done in the other functions aswell without this comment.
    low = low * (hdr_img.maxval / hdr_img.color_res);
    high = high * (hdr_img.maxval / hdr_img.color_res);
    std::cout << "low:" << low << "  high:"<< high << std::endl;


    for (int h = 0; h < hdr_img.height; h++) {
        for (int w = 0; w < hdr_img.width; w++) {
            double lum = RGBtolum(hdr_img.image[h][w]);
            lum = (lum < low) ? 0 : lum - low;
            lum = (lum > high) ? high : lum;
            hdr_img.image[h][w] = lumtoRGB(lum, hdr_img.image[h][w]);
        }
    }
    // To obtain the new color_res we're going to adjust the values back to the [0..color_res] scale
    hdr_img.color_res = (high - low) * (hdr_img.color_res / hdr_img.maxval);
}

void equalize(HDR_image& hdr_img, double high) {
    high = high * (hdr_img.maxval / hdr_img.color_res);

    double highest = 0;
    //Evaluate highest luminosity
    for (int h = 0; h < hdr_img.height; h++) {
        for (int w = 0; w < hdr_img.width; w++) {
            double lum = RGBtolum(hdr_img.image[h][w]);
            highest = (lum > highest) ? lum : highest;
        }
    }
    //Now that the highest point is found, we can equalize all values
    for (int h = 0; h < hdr_img.height; h++) {
        for (int w = 0; w < hdr_img.width; w++) {
            double lum = RGBtolum(hdr_img.image[h][w])/highest * high;
            hdr_img.image[h][w] = lumtoRGB(lum, hdr_img.image[h][w]);
        }
    }
    hdr_img.color_res = high * (hdr_img.color_res / hdr_img.maxval);
}

void equalizeANDclamp(HDR_image& hdr_img, double high) {
    high = high * (hdr_img.maxval / hdr_img.color_res);
    for (int h = 0; h < hdr_img.height; h++) {
        for (int w = 0; w < hdr_img.width; w++) {
            double lum = RGBtolum(hdr_img.image[h][w]) / (hdr_img.color_res * hdr_img.maxval) * high;
            hdr_img.image[h][w] = lumtoRGB(lum, hdr_img.image[h][w]);
        }
    }
    hdr_img.color_res = high * (hdr_img.color_res / hdr_img.maxval);
}

void gammaCorrection(HDR_image& hdr_img, double gamma) {
    for (int h = 0; h < hdr_img.height; h++) {
        for (int w = 0; w < hdr_img.width; w++) {
            double lum = pow((RGBtolum(hdr_img.image[h][w]) / (hdr_img.color_res * hdr_img.maxval)), (1.0F / gamma))
                            * (hdr_img.color_res * hdr_img.maxval);
            hdr_img.image[h][w] = lumtoRGB(lum, hdr_img.image[h][w]);
        }
    }
}

void gammaRestore(HDR_image& hdr_img, double gamma) {
    for (int h = 0; h < hdr_img.height; h++) {
        for (int w = 0; w < hdr_img.width; w++) {
            double lum = pow((RGBtolum(hdr_img.image[h][w]) / (hdr_img.color_res * hdr_img.maxval)), gamma)
                            * (hdr_img.color_res * hdr_img.maxval);
            hdr_img.image[h][w] = lumtoRGB(lum, hdr_img.image[h][w]);
        }
    }
}
*/

#include "tone_mapping.h"

int main(int argc, char *argv[]) {
    //HDRImage h("HDR_PPM_Files/forest_path.ppm");
    //std::cout << h << std::endl;
    //std::cout << RGB(523352.5,523532.2,732634.5);

    // prueba_debug(h, 255);

    // equalize(h, 255);
    // std::string path("HDR_PPM_Files/forest_path_rewritten.ppm");
    // h.write_ppm(path);
    // std::cout << h << std::endl;

    Image h("test.ppm");
    Image i = h, j = h;
    Image k = h, l = h, m = h;
    std::cout << h << std::endl;
    for (auto& rp : h.pixels) {
        for (auto& cp : rp) {
            std::cout << cp << "\t";
        }
        std::cout << std::endl;
    }

    clamp(h, 4, 7);
    h.write_ppm("test_clamp-1.ppm", true);

    clamp(i, 0, 8);
    i.write_ppm("test_clamp-2.ppm", true);

    clamp(j, 6);
    j.write_ppm("test_clamp-3.ppm", true);

    equalization(k, 189);
    k.write_ppm("test_equalize-1.ppm", true);

    equalization(l, 255);
    l.write_ppm("test_equalize-2.ppm", true);

    equalization(m, 19);
    m.write_ppm("test_equalize-3.ppm", true);

    // h.write_ppm("forest_path_ldr.ppm", true);
    // gammaCorrection(h, 2);
    // h.write_ppm("forest_path_ldr_gamma_2.ppm", true);
    //clamp(h, 0, 255);
    //h.write_ppm("forest_path_ldr_clamp.ppm", true);
}