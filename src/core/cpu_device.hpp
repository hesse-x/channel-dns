#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

#ifndef CHANNEL_DNS_PURE_MPI
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>
#endif

struct Grid1D {
    std::size_t elements;
    std::size_t block_size;
};

struct Block1D {
    std::size_t index;
    std::size_t begin;
    std::size_t end;
};

// The execution interface for host-resident computation. The default backend is
// async FIFO; CHANNEL_DNS_PURE_MPI executes inline on the calling MPI thread.
// A block is a CPU loop range, not a GPU group of concurrently executing threads.
class CpuDevice {
  public:
#ifdef CHANNEL_DNS_PURE_MPI
    ~CpuDevice() = default;
    constexpr int threads() const noexcept { return 1; }

    template <typename Function> void enqueue(Function &&task) {
        std::forward<Function>(task)();
    }
    void sync() const noexcept {}

    template <typename Function>
    void enqueue_for_each(std::ptrdiff_t first, std::ptrdiff_t last, Function &&function,
                          std::size_t = 1) {
        for (auto i = first; i < last; ++i) function(i);
    }

    template <typename Function>
    void for_each(std::ptrdiff_t first, std::ptrdiff_t last, Function &&function,
                  std::size_t grain = 1) {
        enqueue_for_each(first, last, std::forward<Function>(function), grain);
    }
#else
    ~CpuDevice();
    int threads() const noexcept { return threads_; }

    void enqueue(std::function<void()> task);
    void sync();

    // Each chunk owns its callable. Referenced objects must outlive sync().
    template <typename Function>
    void enqueue_for_each(std::ptrdiff_t first, std::ptrdiff_t last, Function &&function,
                          std::size_t grain = 1) {
        require_owner();
        if (last <= first) return;
        const auto count = static_cast<std::size_t>(last - first);
        const auto chunks = std::min<std::size_t>(threads(),
            1 + (count - 1) / std::max<std::size_t>(grain, 1));
        using Body = std::decay_t<Function>;
        static_assert(std::is_copy_constructible_v<Body>, "range callables must be copyable");
        enqueue_batch(std::make_unique<RangeBatch<Body>>(
            first, count, chunks, std::forward<Function>(function)));
    }

    // Synchronous convenience. Nested work stays on its current worker.
    template <typename Function>
    void for_each(std::ptrdiff_t first, std::ptrdiff_t last, Function &&function,
                  std::size_t grain = 1) {
        if (executing_) {
            for (auto i = first; i < last; ++i) function(i);
            return;
        }
        enqueue_for_each(first, last, std::forward<Function>(function), grain);
        sync();
    }

#endif
    CpuDevice(const CpuDevice &) = delete;
    CpuDevice &operator=(const CpuDevice &) = delete;

    template <typename Kernel> void launch(Grid1D grid, Kernel &&kernel) {
        if (!grid.block_size)
            throw std::invalid_argument("kernel block size must be positive");
        if (grid.elements > static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()))
            throw std::overflow_error("kernel grid exceeds addressable range");
        const auto blocks = grid.elements ? 1 + (grid.elements - 1) / grid.block_size : 0;
        enqueue_for_each(0, blocks,
            [grid, body = std::forward<Kernel>(kernel)](std::size_t block) mutable {
                const auto begin = block * grid.block_size;
                body(Block1D{block, begin,
                             begin + std::min(grid.block_size, grid.elements - begin)});
            });
    }

    template <typename T> void copy_n(const T *source, std::size_t count, T *destination) {
        if (!count || source == destination) return;
#ifdef CHANNEL_DNS_PURE_MPI
        std::copy_n(source, count, destination);
#else
        const auto chunks = std::min<std::size_t>(threads(), 1 + (count - 1) / 4096);
        for_each(0, chunks, [&](std::size_t lane) {
            const auto first = count / chunks * lane + std::min(count % chunks, lane);
            const auto length = count / chunks + (lane < count % chunks);
            std::copy_n(source + first, length, destination + first);
        });
#endif
    }

    template <typename T> void fill_n(T *destination, std::size_t count, T value) {
        if (!count) return;
#ifdef CHANNEL_DNS_PURE_MPI
        std::fill_n(destination, count, value);
#else
        const auto chunks = std::min<std::size_t>(threads(), 1 + (count - 1) / 4096);
        for_each(0, chunks, [&](std::size_t lane) {
            const auto first = count / chunks * lane + std::min(count % chunks, lane);
            const auto length = count / chunks + (lane < count % chunks);
            std::fill_n(destination + first, length, value);
        });
#endif
    }

  private:
    friend CpuDevice &cpu_device(int);
#ifdef CHANNEL_DNS_PURE_MPI
    explicit constexpr CpuDevice(int) {}
#else
    // Type erasure happens once per batch, not once per chunk. Scheduling state
    // and worker ownership are private to the device implementation.
    class Impl;
    struct Batch {
        explicit Batch(std::size_t chunks) : chunks(chunks), remaining(chunks) {}
        virtual ~Batch() = default;
        virtual void run(std::size_t lane) = 0;
        const std::size_t chunks;
      private:
        friend class Impl;
        std::size_t next = 0;
        std::size_t remaining;
    };

    template <typename Function> struct RangeBatch final : Batch {
        RangeBatch(std::ptrdiff_t first, std::size_t count, std::size_t chunks, Function function)
            : Batch(chunks), first(first), width(count / chunks), remainder(count % chunks),
              bodies(chunks, function) {}

        void run(std::size_t lane) override {
            const auto begin = first + static_cast<std::ptrdiff_t>(
                                           width * lane + std::min(remainder, lane));
            const auto end = begin + static_cast<std::ptrdiff_t>(width + (lane < remainder));
            auto &body = bodies[lane];
            for (auto i = begin; i < end; ++i)
                body(i);
        }

        const std::ptrdiff_t first;
        const std::size_t width, remainder;
        std::vector<Function> bodies;
    };

    explicit CpuDevice(int threads);
    void require_owner() const;
    void enqueue_batch(std::unique_ptr<Batch> batch);
    static thread_local bool executing_;
    const int threads_;
    std::unique_ptr<Impl> impl_;
#endif
};

// Pass a positive worker count on first access; zero selects hardware concurrency.
// Later explicit requests must match the initialized pool. Pure MPI always runs inline.
inline CpuDevice &cpu_device(int threads = 0) {
    if (threads < 0) throw std::invalid_argument("thread count must be nonnegative");
    static CpuDevice device(threads);
#ifndef CHANNEL_DNS_PURE_MPI
    if (threads > 0 && threads != device.threads())
        throw std::logic_error("CPU device thread count is already initialized");
#endif
    return device;
}
