#include <iostream>
#include "image.hpp"

void help(int status = 0, std::string msg = "") {
    std::cout << msg << "\nusage: tone_mapper -i <path> -o <path> [<flags>] <mapping options> \n"
        << "  Options: different behaviours for the tone mapper.\n"
        << "    -c,--conversion     Conversion between HDR to LDR.\n"
        << "    -h,--help           Help.\n"
        << "    -i,--input <path>   Where is your image.\n"
        << "    -l,--ldr            Your image is LDR\n"           
        << "    -o,--output <path>  Output file.\n"
        << "  Mapping options: which tone mapping options to apply. The mapping level has to be decimal.\n"
        << "    --clamp <n>         Clamping.\n"
        << "    --equalize <n>      Equalization.\n"
        << "    --gamma <n>         Gamma curve.\n\n";
    exit(status);
}

int main(int argc, char *argv[]) {

    int flags = 0; float c_map = 1.0, e_map = 1.0, g_map = 2.2;
    bool conversion = false, ldr = false;
    std::string input, output;
    for (int s = 1; s < argc; s++) {
        if (!strcmp(argv[s], "-c") || !strcmp(argv[s], "--conversion")) conversion = true;
        else if (!strcmp(argv[s], "-h") || !strcmp(argv[s], "--help")) help();
        else if (!strcmp(argv[s], "-i") || !strcmp(argv[s], "--input")) input = argv[++s];
        else if (!strcmp(argv[s], "-l") || !strcmp(argv[s], "--ldr")) ldr = true;
        else if (!strcmp(argv[s], "-o") || !strcmp(argv[s], "--output")) output = argv[++s];
        else if (!strcmp(argv[s], "--clamp")) {
            flags |= clamp;
            c_map = std::stof(argv[++s]);
        } 
        else if (!strcmp(argv[s], "--equalize")) {
            flags |= equalize;
            e_map = std::stof(argv[++s]);
        }
        else if (!strcmp(argv[s], "--gamma")) {
            flags |= gamma;
            g_map = std::stof(argv[++s]);
        } else help(1, "\nerror: unknown option '" + std::string(argv[s]) + "'.");
    }
    if (!input.length()) help(1, "\nerror: no input image given.");
    if (!output.length()) help(1, "\nerror: no output image given.");
    if (!flags) help(1, "\nerror: no mapping operations given.");

    Image i(input);
    tone_mapping(i, flags, c_map, e_map, g_map);
    export_image(i, output, ldr, conversion);
}
