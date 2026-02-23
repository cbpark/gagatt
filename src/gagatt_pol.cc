#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "gagatt_unpol: the output has been stored in " << argv[1]
              << '\n';
}
