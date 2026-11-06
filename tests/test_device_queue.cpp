#include <gtest/gtest.h>

#include <atomic>
#include <condition_variable>
#include <chrono>
#include <cstdlib>
#include <future>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <thread>
#include <vector>

#include "core/cpu_device.hpp"

TEST(CpuDeviceQueue, SingletonAndExactCoverageAcrossRepeatedUnevenRanges) {
    auto &device = cpu_device();
    EXPECT_EQ(&device, &cpu_device());
    std::set<std::thread::id> participants;
    std::mutex mutex;
    for (int repeat = 0; repeat < 50; ++repeat) {
        std::vector<int> values(10007, 0);
        device.for_each(-3, 10004, [&](int i) {
            ++values[i + 3];
            if (repeat == 0) {
                std::lock_guard lock(mutex);
                participants.insert(std::this_thread::get_id());
            }
        });
        for (int value : values)
            ASSERT_EQ(value, 1);
    }
    EXPECT_GE(participants.size(), 1u);
    EXPECT_LE(participants.size(), static_cast<std::size_t>(device.threads()));
}

TEST(CpuDeviceQueue, NestedWorkFinishesWithoutAdditionalThreads) {
    std::vector<int> values(257 * 31, 0);
    cpu_device().for_each(0, 257, [&](int i) {
        const auto outer = std::this_thread::get_id();
        cpu_device().for_each(0, 31, [&](int j) {
            if (outer != std::this_thread::get_id())
                throw std::runtime_error("nested dispatch");
            values[i * 31 + j] = i + j + 1;
        });
    });
    for (int i = 0; i < 257; ++i)
        for (int j = 0; j < 31; ++j)
            EXPECT_EQ(values[i * 31 + j], i + j + 1);
}

TEST(CpuDeviceQueue, ExceptionsPropagateAndPoolCanBeReused) {
    for (int failing : {0, 73}) {
        EXPECT_THROW(cpu_device().for_each(0, 101,
                                  [&](int i) {
                                      if (i == failing)
                                          throw std::runtime_error("task failure");
                                  }),
                     std::runtime_error);
        std::atomic<int> count{0};
        cpu_device().for_each(0, 101, [&](int) { ++count; });
        EXPECT_EQ(count, 101);
    }
    cpu_device().for_each(10, 10, [](int) { throw std::runtime_error("empty range executed"); });
}

TEST(CpuDeviceQueue, SmallRangesCompleteAndForeignSubmissionsAreRejected) {
    auto &device = cpu_device();
    std::atomic<int> count{0};
    device.for_each(0, 7, [&](int) { ++count; }, 8);
    EXPECT_EQ(count, 7);
    bool rejected = false;
    std::thread foreign([&] {
        try {
            device.for_each(0, 100, [](int) {});
        } catch (const std::logic_error &) {
            rejected = true;
        }
    });
    foreign.join();
    EXPECT_TRUE(rejected);
}

TEST(CpuDeviceQueue, EnqueueIsAsynchronousAndSyncWaitsForAllSubmittedWork) {
    auto &device = cpu_device();
    std::promise<void> release, started;
    auto gate = release.get_future().share();
    auto entered = started.get_future();
    std::atomic<int> completed{0};
    device.enqueue([&] {
        started.set_value();
        gate.wait();
        ++completed;
    });
    device.enqueue([&] { ++completed; });
    // The first task is still blocked; submission itself must have returned.
    EXPECT_EQ(entered.wait_for(std::chrono::seconds(2)), std::future_status::ready);
    release.set_value();
    device.sync();
    EXPECT_EQ(completed, 2);
    device.sync(); // Repeated synchronization is safe.
}

TEST(CpuDeviceQueue, MultipleQueuedRangesHonorOrderWithOneFinalSync) {
    auto &device = cpu_device();
    std::vector<int> a(10007, 0), b(10007, 0), c(10007, 0);
    device.enqueue_for_each(0, a.size(), [&](int i) { a[i] = i; });
    device.enqueue_for_each(0, b.size(), [&](int i) { b[i] = 2 * i; });
    device.enqueue_for_each(0, c.size(), [&](int i) { c[i] = a[i] + b[i]; });
    device.sync();
    for (int i = 0; i < 10007; ++i)
        EXPECT_EQ(c[i], 3 * i);
    device.enqueue([] { throw std::runtime_error("asynchronous failure"); });
    std::atomic<bool> drained{false};
    device.enqueue([&] { drained = true; });
    EXPECT_THROW(device.sync(), std::runtime_error);
    EXPECT_TRUE(drained);
    device.enqueue([] {});
    EXPECT_NO_THROW(device.sync());
}

TEST(CpuDeviceQueue, WorkerCountMatchesConfigurationAndOwnerNeverExecutesTasks) {
    auto &device = cpu_device();
    ASSERT_GE(device.threads(), 1);
    std::mutex rendezvous_mutex;
    std::condition_variable rendezvous;
    int arrived = 0;
    std::vector<std::thread::id> ids(device.threads());
    device.enqueue_for_each(0, device.threads(), [&](int i) {
        ids[i] = std::this_thread::get_id();
        std::unique_lock<std::mutex> lock(rendezvous_mutex);
        if (++arrived == device.threads())
            rendezvous.notify_all();
        else
            rendezvous.wait(lock, [&] { return arrived == device.threads(); });
    });
    device.sync();
    const std::set<std::thread::id> participants(ids.begin(), ids.end());
    EXPECT_EQ(participants.size(), static_cast<std::size_t>(device.threads()));
    EXPECT_EQ(participants.count(std::this_thread::get_id()), 0u);
    std::thread::id single_task_thread;
    device.enqueue([&] { single_task_thread = std::this_thread::get_id(); });
    device.sync();
    EXPECT_NE(single_task_thread, std::this_thread::get_id());
    EXPECT_NE(single_task_thread, std::thread::id{});
}

TEST(CpuDeviceQueue, MutableCallableStateIsPrivateToEachChunk) {
    auto &device = cpu_device();
    const auto chunks = static_cast<std::size_t>(device.threads());
    const std::size_t count = chunks * 3 + 1;
    std::vector<int> values(count);
    device.enqueue_for_each(0, count, [&, visits = 0](int i) mutable { values[i] = ++visits; });
    device.sync();
    std::size_t first = 0;
    for (std::size_t lane = 0; lane < chunks; ++lane) {
        const auto length = count / chunks + (lane < count % chunks);
        for (std::size_t i = 0; i < length; ++i)
            EXPECT_EQ(values[first + i], static_cast<int>(i + 1));
        first += length;
    }
}

TEST(CpuDeviceQueue, QueuedBatchesOwnTemporaryCapturesUntilSync) {
    auto &device = cpu_device();
    std::promise<void> release;
    auto gate = release.get_future().share();
    std::weak_ptr<int> first, second;
    std::vector<int> a(103), b(107);
    {
        auto value = std::make_shared<int>(17);
        first = value;
        device.enqueue_for_each(0, a.size(), [&, gate, value](int i) {
            gate.wait();
            a[i] = *value + i;
        });
    }
    {
        auto value = std::make_shared<int>(29);
        second = value;
        device.enqueue_for_each(0, b.size(), [&, gate, value](int i) {
            gate.wait();
            b[i] = *value + i;
        });
    }
    EXPECT_FALSE(first.expired());
    EXPECT_FALSE(second.expired());
    release.set_value();
    device.sync();
    EXPECT_TRUE(first.expired());
    EXPECT_TRUE(second.expired());
    for (int i = 0; i < static_cast<int>(a.size()); ++i) EXPECT_EQ(a[i], 17 + i);
    for (int i = 0; i < static_cast<int>(b.size()); ++i) EXPECT_EQ(b[i], 29 + i);
}

TEST(CpuDeviceQueue, ThrowingCallableCopyPublishesNoPartialBatch) {
    struct Body {
        Body(int &copies, int fail_at, std::atomic<int> &calls)
            : copies(copies), fail_at(fail_at), calls(calls) {}
        Body(Body &&) = default;
        Body(const Body &other) : copies(other.copies), fail_at(other.fail_at), calls(other.calls) {
            if (++copies == fail_at) throw std::runtime_error("copy failed");
        }
        void operator()(int) { ++calls; }
        int &copies;
        int fail_at;
        std::atomic<int> &calls;
    };
    auto &device = cpu_device();
    int copies = 0;
    std::atomic<int> calls{0};
    std::atomic<bool> earlier_done{false};
    device.enqueue([&] { earlier_done = true; });
    EXPECT_THROW(device.enqueue_for_each(0, 100, Body(copies, std::min(device.threads(), 2), calls)),
                 std::runtime_error);
    device.sync();
    EXPECT_TRUE(earlier_done);
    EXPECT_EQ(calls, 0);
    device.for_each(0, 100, [&](int) { ++calls; });
    EXPECT_EQ(calls, 100);
}

TEST(CpuDeviceQueue, SyncIncludesDestructionOfSingleTaskCaptures) {
    auto &device = cpu_device();
    for (int repeat = 0; repeat < 50; ++repeat) {
        auto value = std::make_shared<int>(repeat);
        std::weak_ptr<int> lifetime = value;
        device.enqueue([value] {});
        value.reset();
        device.sync();
        EXPECT_TRUE(lifetime.expired());
    }
    EXPECT_THROW(device.enqueue(std::function<void()>{}), std::invalid_argument);
}

TEST(CpuDeviceQueue, LaterSubmissionWaitsForEveryChunkOfPreviousBatch) {
    auto &device = cpu_device();
    std::promise<void> release, blocked_chunk_started, later_started;
    auto gate = release.get_future().share();
    auto blocked = blocked_chunk_started.get_future();
    auto later = later_started.get_future();
    std::atomic<int> completed{0};
    device.enqueue_for_each(0, device.threads(), [&](int i) {
        if (i == device.threads() - 1) {
            blocked_chunk_started.set_value();
            gate.wait();
        }
        ++completed;
    });
    EXPECT_EQ(blocked.wait_for(std::chrono::seconds(2)), std::future_status::ready);
    device.enqueue([&] {
        EXPECT_EQ(completed.load(), device.threads());
        later_started.set_value();
    });
    EXPECT_EQ(later.wait_for(std::chrono::milliseconds(20)), std::future_status::timeout);
    release.set_value();
    device.sync();
    EXPECT_EQ(later.wait_for(std::chrono::seconds(0)), std::future_status::ready);
}

TEST(CpuDeviceQueue, SingleTasksCompleteInSubmissionOrder) {
    auto &device = cpu_device();
    std::atomic<int> sequence{0};
    for (int i = 0; i < 1000; ++i)
        device.enqueue([&, i] { EXPECT_EQ(sequence.fetch_add(1), i); });
    device.sync();
    EXPECT_EQ(sequence.load(), 1000);
}

TEST(CpuDeviceQueue, CaptureDestructionPrecedesNextSubmission) {
    auto &device = cpu_device();
    std::promise<void> release;
    auto gate = release.get_future().share();
    std::atomic<bool> retired{false};
    auto lifetime = std::shared_ptr<int>(new int, [&](int *p) { delete p; retired = true; });
    device.enqueue([lifetime, gate] { gate.wait(); });
    lifetime.reset();
    device.enqueue([&] { EXPECT_TRUE(retired.load()); });
    release.set_value();
    device.sync();
    EXPECT_TRUE(retired.load());
}
