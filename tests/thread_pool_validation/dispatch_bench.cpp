#include "core/thread_pool.hpp"
#include <array>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <new>

std::atomic<bool> count_new{false};
std::atomic<std::size_t> allocations{0};
void *operator new(std::size_t n) {
    void *p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    if (count_new.load(std::memory_order_relaxed)) allocations.fetch_add(1, std::memory_order_relaxed);
    return p;
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }

int main() {
    auto &pool = thread_pool();
    std::array<std::size_t, 1024> output{};
    auto submit = [&](std::size_t iteration) {
        const std::array<std::size_t, 8> state{iteration, 1, 2, 3, 4, 5, 6, 7};
        pool.for_each(0, output.size(), [&, state](std::size_t i) { output[i] = state[0] + i; });
    };
    for (int i = 0; i < 100; ++i) submit(i);
    allocations = 0;
    count_new = true;
    for (int i = 0; i < 1000; ++i) submit(i);
    count_new = false;
    const double per_dispatch = allocations.load() / 1000.0;
    constexpr int repeats = 10000;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) submit(i);
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (output.back() != repeats - 1 + output.size() - 1) return 1;
    std::cout << "threads=" << pool.size() << " allocations_per_dispatch=" << per_dispatch
              << " dispatch_us=" << elapsed * 1e6 / repeats << '\n';
}
