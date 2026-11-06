#include "core/cpu_device.hpp"
#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

TEST(CpuDeviceInline, DependentWorkCompletesOnCallerBeforeReturning) {
    auto &device = cpu_device();
    const auto caller = std::this_thread::get_id();
    std::vector<int> values(17);
    device.enqueue([token = std::make_unique<int>(7), &values, caller, &device] {
        EXPECT_EQ(std::this_thread::get_id(), caller);
        device.fill_n(values.data(), values.size(), *token);
    });
    EXPECT_EQ(values.back(), 7);
    device.launch(Grid1D{values.size(), 3}, [&](Block1D block) {
        EXPECT_EQ(std::this_thread::get_id(), caller);
        device.for_each(block.begin, block.end, [&](int i) { values[i] += i; });
    });
    for (int i = 0; i < 17; ++i) EXPECT_EQ(values[i], 7 + i);
    device.enqueue_for_each(0, 17, [&](int i) {
        EXPECT_EQ(std::this_thread::get_id(), caller);
        ++values[i];
    }, 4096);
    EXPECT_EQ(values.back(), 24);
    EXPECT_EQ(device.threads(), 1);
    EXPECT_NO_THROW(device.sync());
}

TEST(CpuDeviceInline, FailuresPropagateAtCallSiteAndDeviceRemainsUsable) {
    auto &device = cpu_device();
    EXPECT_THROW(device.enqueue([] { throw std::runtime_error("task failed"); }),
                 std::runtime_error);
    std::vector<int> values(9);
    EXPECT_THROW(device.enqueue_for_each(0, 9, [&](int i) {
        if (i == 4) throw std::runtime_error("range failed");
        values[i] = 1;
    }), std::runtime_error);
    for (int i = 0; i < 9; ++i) EXPECT_EQ(values[i], i < 4 ? 1 : 0);
    device.fill_n(values.data(), values.size(), 3);
    EXPECT_EQ(values.back(), 3);
    EXPECT_NO_THROW(device.sync());
}
