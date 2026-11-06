#include "core/cpu_device.hpp"
#include <gtest/gtest.h>
#include <cstdlib>
#include <vector>

TEST(CpuDevice, SingletonAndGridCoverRaggedBlocksExactly) {
    auto &device = cpu_device();
    EXPECT_EQ(&device, &cpu_device());
#ifdef CHANNEL_DNS_MPI_ONLY
    EXPECT_EQ(device.threads(), 1);
#else
    EXPECT_GE(device.threads(), 1);
#endif
    std::vector<int> visits(10007, 0);
    device.launch(Grid1D{visits.size(), 127}, [&](Block1D block) {
        EXPECT_EQ(block.begin, block.index * 127);
        EXPECT_LE(block.end - block.begin, 127u);
        for (auto i = block.begin; i < block.end; ++i)
            ++visits[i];
    });
    device.sync();
    for (int count : visits)
        EXPECT_EQ(count, 1);
}

TEST(CpuDevice, KernelCanPerformMultipleOperationsAndShareOrderedQueue) {
    auto &device = cpu_device();
    std::vector<int> a(10007), b(10007);
    bool ordered_task = false;
    device.enqueue([&] { ordered_task = true; });
    device.launch(Grid1D{a.size(), 64}, [&](Block1D block) {
        for (auto i = block.begin; i < block.end; ++i) {
            a[i] = static_cast<int>(i);
            a[i] *= 3;
            b[i] = a[i] + 2;
        }
    });
    device.sync();
    EXPECT_TRUE(ordered_task);
    for (int i = 0; i < 10007; ++i)
        EXPECT_EQ(b[i], 3 * i + 2);
    EXPECT_THROW(device.launch(Grid1D{1, 0}, [](Block1D) {}), std::invalid_argument);
    EXPECT_THROW(device.launch(Grid1D{std::numeric_limits<std::size_t>::max(), 1}, [](Block1D) {}),
                 std::overflow_error);
    device.launch(Grid1D{0, 1}, [](Block1D) { throw std::runtime_error("empty grid ran"); });
    EXPECT_NO_THROW(device.sync());
}

TEST(CpuDevice, CopyAndFillKeepGuardsAndSupportNestedExecution) {
    auto &device = cpu_device();
    constexpr std::size_t count = 10007;
    std::vector<int> source(count + 2, -7), output(count + 4, -9);
    for (std::size_t i = 0; i < count; ++i) source[i + 1] = static_cast<int>(i * 3);
    device.copy_n(source.data() + 1, count, output.data() + 2);
    for (std::size_t i = 0; i < count; ++i) EXPECT_EQ(output[i + 2], source[i + 1]);
    EXPECT_EQ(output[0], -9); EXPECT_EQ(output[1], -9);
    EXPECT_EQ(output[count + 2], -9); EXPECT_EQ(output[count + 3], -9);
    device.launch(Grid1D{1, 1}, [&](Block1D) {
        device.fill_n(output.data() + 2, count, 42);
        device.copy_n(output.data() + 2, count, source.data() + 1);
    });
    device.sync();
    for (std::size_t i = 0; i < count; ++i) EXPECT_EQ(source[i + 1], 42);
    EXPECT_EQ(source.front(), -7); EXPECT_EQ(source.back(), -7);
    EXPECT_NO_THROW(device.copy_n(output.data(), output.size(), output.data()));
    EXPECT_NO_THROW(device.copy_n<int>(nullptr, 0, nullptr));
    EXPECT_NO_THROW(device.fill_n<int>(nullptr, 0, 0));
}
