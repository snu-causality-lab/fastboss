#pragma once

#include <cstdlib>

namespace fastboss::detail {

inline bool better_indexed_score(double candidate_score,
                                 int candidate_index,
                                 double incumbent_score,
                                 int incumbent_index) noexcept {
    if (candidate_score != incumbent_score) {
        return candidate_score > incumbent_score;
    }
    return candidate_index < incumbent_index;
}

inline bool better_relocation(double candidate_score,
                              int candidate_slot,
                              double incumbent_score,
                              int incumbent_slot,
                              int old_position) noexcept {
    if (candidate_score != incumbent_score) {
        return candidate_score > incumbent_score;
    }
    const int candidate_distance = std::abs(candidate_slot - old_position);
    const int incumbent_distance = std::abs(incumbent_slot - old_position);
    if (candidate_distance != incumbent_distance) {
        return candidate_distance < incumbent_distance;
    }
    return candidate_slot < incumbent_slot;
}

}  // namespace fastboss::detail
