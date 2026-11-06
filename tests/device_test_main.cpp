#include "core/cpu_device.hpp"
#include <gtest/gtest.h>
#include <string>

int main(int argc, char **argv) {
    int threads = 1;
    if (argc >= 3 && std::string(argv[1]) == "--threads") {
        threads = std::stoi(argv[2]);
        argc -= 2;
        for (int i = 1; i < argc; ++i) argv[i] = argv[i + 2];
    }
    auto &device = cpu_device(threads);
#ifndef CHANNEL_DNS_MPI_ONLY
    if (device.threads() != threads) return 1;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
