#include <gtest/gtest.h>
#include <mpi.h>

#include <cstdint>
#include <limits>
#include <vector>

#include "core/mpi_context.hpp"

TEST(MpiContext, CommunicatorIdentityAndGatherPreserveIntegerWidthAndSign) {
    MpiContext mpi(MPI_COMM_WORLD);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    EXPECT_EQ(mpi.comm(), MPI_COMM_WORLD);
    EXPECT_EQ(mpi.rank(), rank);
    EXPECT_EQ(mpi.size(), size);
    const auto unsigned_values = mpi.allgather_u64(std::numeric_limits<std::uint64_t>::max() - rank, "gather unsigned");
    const auto signed_values = mpi.allgather_i64(std::numeric_limits<std::int64_t>::min() + rank, "gather signed");
    ASSERT_EQ(unsigned_values.size(), static_cast<std::size_t>(size));
    ASSERT_EQ(signed_values.size(), static_cast<std::size_t>(size));
    for (int r = 0; r < size; ++r) {
        EXPECT_EQ(unsigned_values[r], std::numeric_limits<std::uint64_t>::max() - r);
        EXPECT_EQ(signed_values[r], std::numeric_limits<std::int64_t>::min() + r);
    }
}

TEST(MpiContext, CollectivesSupportNonzeroRoot) {
    MpiContext mpi(MPI_COMM_WORLD);
    const auto root = mpi.size() - 1;
    double data[] = {mpi.rank() == root ? 3.25 : -1., mpi.rank() == root ? -7.5 : -1.};
    mpi.broadcast_doubles(data, 2, root, "broadcast last root");
    EXPECT_EQ(data[0], 3.25);
    EXPECT_EQ(data[1], -7.5);
    const auto sum = mpi.reduce_sum(mpi.rank() + .5, root, "sum last root");
    if (mpi.rank() == root) { EXPECT_EQ(sum, .5 * mpi.size() * mpi.size()); }
    EXPECT_EQ(mpi.allreduce_min(10 - mpi.rank(), "minimum"), 11 - mpi.size());
    mpi.barrier("barrier");
}

TEST(MpiContext, SendrecvExchangesPairedRankBuffers) {
    MpiContext mpi(MPI_COMM_WORLD);
    const auto peer = mpi.size() - 1 - mpi.rank();
    double send[] = {mpi.rank() + .25, -mpi.rank() - .5}, recv[2]{};
    mpi.sendrecv_doubles(send, 2, peer, 17, recv, 2, 17, "pair exchange");
    EXPECT_EQ(recv[0], peer + .25);
    EXPECT_EQ(recv[1], -peer - .5);
}

TEST(MpiContext, AlltoallvHandlesUnequalCountsZeroCountsAndGappedDisplacements) {
    MpiContext mpi(MPI_COMM_WORLD);
    const int n = mpi.size(), r = mpi.rank();
    std::vector<int> sc(n), sd(n), rc(n), rd(n);
    std::vector<double> send(4*n, -999.), recv(4*n, -999.);
    for (int peer = 0; peer < n; ++peer) {
        sc[peer] = rc[peer] = (r + peer) % 3;
        sd[peer] = rd[peer] = 4*peer;
        for (int i = 0; i < sc[peer]; ++i) send[sd[peer]+i] = 100*r + 10*peer + i;
    }
    mpi.alltoallv_doubles(send.data(), sc, sd, recv.data(), rc, rd, "alltoallv unequal");
    for (int peer = 0; peer < n; ++peer)
        for (int i = 0; i < 4; ++i)
            EXPECT_EQ(recv[rd[peer]+i], i < rc[peer] ? 100*peer + 10*r + i : -999.);
}

TEST(MpiContext, InvalidVectorSizesAreRejectedBeforeCollective) {
    MpiContext mpi(MPI_COMM_WORLD);
    std::vector<int> valid(mpi.size(), 0), invalid(mpi.size()+1, 0);
    double buffer = 0;
    for (int slot = 0; slot < 4; ++slot) {
        SCOPED_TRACE(slot);
        EXPECT_THROW(mpi.alltoallv_doubles(&buffer, slot == 0 ? invalid : valid,
            slot == 1 ? invalid : valid, &buffer, slot == 2 ? invalid : valid,
            slot == 3 ? invalid : valid, "invalid sizes"), std::invalid_argument);
    }
}

TEST(MpiContext, CommunicationErrorsIncludeOperationAndRank) {
    MpiContext mpi(MPI_COMM_WORLD);
    double value = 0;
    try {
        mpi.sendrecv_doubles(&value, 1, mpi.size(), 5, &value, 1, 5, "invalid peer");
        FAIL() << "invalid peer accepted";
    } catch (const std::runtime_error& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("invalid peer"), std::string::npos);
        EXPECT_NE(message.find("rank " + std::to_string(mpi.rank())), std::string::npos);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MpiContext mpi(MPI_COMM_WORLD);
    const int success = mpi.allreduce_min(result == 0 ? 1 : 0, "test status");
    MPI_Finalize();
    return success ? 0 : 1;
}
