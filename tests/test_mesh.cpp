#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "core/mesh.hpp"
#include "core/params.hpp"

namespace {
MeshParm distributed_parameters() {
    MeshParm p{};
    p.NX = 32;
    p.NY = 32;
    p.NZ = 32;
    p.NXS = 16;
    p.NYS = 30;
    p.NZS = 16;
    p.NX2 = 48;
    p.NY2 = 48;
    p.NZ2 = 48;
    p.INDEX_Serial = 0;
    return p;
}
}

TEST(Mesh, DistributedBlocksCoverEveryIndexExactlyOnce) {
    for (const int ranks : {2, 3, 4}) {
        std::vector<int> x(17), y(33), y2(49);
        int zero_owners = 0;
        for (int rank = 0; rank < ranks; ++rank) {
            Mesh mesh(distributed_parameters(), 2.0, 3.0, rank, ranks);
            for (int i = 0; i < mesh.x_count(rank); ++i) ++x[mesh.x_begin(rank) + i];
            for (int i = 0; i < mesh.y_count(rank); ++i) ++y[mesh.y_begin(rank) + i];
            for (int i = 0; i < mesh.y2_count(rank); ++i) ++y2[mesh.y2_begin(rank) + i];
            zero_owners += mesh.owns_zero_x_mode() ? 1 : 0;
            EXPECT_EQ(mesh.zero_mode_owner(), 0);
        }
        EXPECT_EQ(x, std::vector<int>(17, 1));
        EXPECT_EQ(y, std::vector<int>(33, 1));
        EXPECT_EQ(y2, std::vector<int>(49, 1));
        EXPECT_EQ(zero_owners, 1);
    }
}

TEST(Mesh, RejectsEmptyDistributedBlocks) {
    auto p = distributed_parameters();
    p.NXS = 1;
    EXPECT_THROW((Mesh(p, 2.0, 3.0, 0, 3)), std::invalid_argument);
}
