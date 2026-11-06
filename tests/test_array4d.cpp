// Array4D 单测（D4/D5）：布局等价性、z 偏移、语义。
#include "core/array4d.hpp"

#include <gtest/gtest.h>

TEST(Array4D, LayoutHasContiguousComponents) {
    // 填充线性序号，断言 data()[linear] 与 operator()(i,j,k,l) 一致
    Array4D<double> a(3, 4, 5, 6);
    int seq = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 5; ++k)
                for (int l = 0; l < 6; ++l, ++seq)
                    a(i, j, k, l) = static_cast<double>(seq);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 5; ++k)
                for (int l = 0; l < 6; ++l) {
                    int linear = ((i * 4 + j) * 5 + k) * 6 + l;
                    EXPECT_DOUBLE_EQ(a.data()[linear], a(i, j, k, l));
                }
    EXPECT_EQ(a.size(), 3 * 4 * 5 * 6);
}

TEST(Array4D, ZOffsetAccessor) {
    // iz ∈ [-NZS, NZS] → at_z(iz) = data[iz + NZS]
    const int NZS = 8;
    Array4D<double> a(2 * NZS + 1, 2, 2, 2);
    a.set_z_offset(NZS);
    a.fill(0.0);
    a.at_z(-NZS, 0, 0, 0) = 1.0;
    EXPECT_DOUBLE_EQ(a.data()[0], 1.0);                 // -NZS 命中首元素
    const int last = (2 * NZS + 1) * 2 * 2 * 2 - 1;
    a.at_z(NZS, 1, 1, 1) = 2.0;
    EXPECT_DOUBLE_EQ(a.data()[last], 2.0);              // +NZS 命中末元素
    a.at_z(0, 0, 0, 0) = 3.0;
    EXPECT_DOUBLE_EQ(a(NZS, 0, 0, 0), 3.0);             // iz=0 ↔ 物理 NZS
    EXPECT_EQ(a.z_offset(), NZS);
}

TEST(Array4D, DefaultCopyMoveFill) {
    Array4D<double> a;
    EXPECT_EQ(a.size(), 0);
    a.resize(2, 2, 2, 2);
    a.fill(7.0);
    Array4D<double> b = a;              // 拷贝
    EXPECT_DOUBLE_EQ(b(1, 1, 1, 1), 7.0);
    Array4D<double> c = std::move(a);   // 移动
    EXPECT_DOUBLE_EQ(c(0, 0, 0, 0), 7.0);
    EXPECT_EQ(c.dim(0), 2);
    EXPECT_EQ(c.dim(3), 2);
}

#ifndef NDEBUG
TEST(Array4D, OutOfBoundsAssertsInDebug) {
    Array4D<double> a(2, 2, 2, 2);
    EXPECT_DEATH(a(2, 0, 0, 0), "");
    EXPECT_DEATH(a(0, 0, 0, -1), "");
}
#endif
