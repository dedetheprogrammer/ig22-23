#include <iostream>
#include <functional>
#include "scene.hpp"

RGB photon_mapping_ex(Render_params s) {
    std::cout << "PHOTON MAPPING\n";
    return RGB(1,1,1);
}

RGB path_tracing_ex(Render_params s) {
    std::cout << "PATH TRACING\n";
    return RGB(1,1,1);
}

RGB ray_tracing_ex(Render_params s) {
    std::cout << "RAY TRACING\n";
    return RGB(1,1,1);
}

template <class ... Types>
void print(Types && ... inputs) {
    int i = 0;
    ([&] {
        ++i;
        std::cout << inputs << std::endl;
    } (), ...);
}

int main(int argc, char* argv[]) {

    std::function<RGB(Camera*, Render_params)> r;
    r = Camera::ray_tracing;

    print(1, 2, 3.14,
          "Pass me any "
          "number of arguments",
          "I will print\n",
          RGB(1.0,1.0,1.0));
}