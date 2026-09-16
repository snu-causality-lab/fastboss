#include "tie_break.h"

#include <cstdlib>
#include <iostream>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "fastboss_determinism_tests failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

}  // namespace

int main() {
    using fastboss::detail::better_indexed_score;
    using fastboss::detail::better_relocation;

    require(better_indexed_score(4.0, 2, 4.0, 7),
            "equal GST scores must choose the smaller variable ID");
    require(!better_indexed_score(4.0, 7, 4.0, 2),
            "a larger variable ID must not replace a tied smaller ID");
    require(better_indexed_score(5.0, 7, 4.0, 2),
            "a higher score must win regardless of variable ID");

    require(better_relocation(3.0, 3, 3.0, 5, 2),
            "equal relocation scores must choose the shorter move");
    require(better_relocation(3.0, 1, 3.0, 3, 2),
            "equidistant relocation ties must choose the smaller slot");
    require(!better_relocation(2.0, 2, 3.0, 5, 2),
            "a shorter relocation must not replace a higher score");

    std::cout << "fastboss determinism tests passed\n";
    return EXIT_SUCCESS;
}
