#include "core/cpu_device.hpp"
#include <type_traits>
static_assert(!std::is_copy_constructible_v<CpuDevice>);
void public_interface(CpuDevice &device, int *values) {
    device.enqueue([values] { values[0] = 1; });
    device.enqueue_for_each(1, 8, [values](int i) { values[i] = i; });
    device.launch(Grid1D{8, 3}, [values](Block1D block) {
        for (auto i = block.begin; i < block.end; ++i) ++values[i];
    });
    device.sync();
    device.fill_n(values, 8, 3);
}
