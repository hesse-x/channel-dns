// wisdom 导入/导出（D3，step2 §3.4）。
// 阶段 2 单 rank 场景：构造时导入、plan 创建完成后导出（文件固定 channel_dns.wisdom）。
// 多 rank 使用包含 nproc/rank 的独立文件；每 rank 只更新自己的 plan wisdom。
#pragma once

#include <string>

namespace wisdom {

// 导入失败不报错（首次运行无 wisdom 文件），返回是否成功。
bool import(const std::string& file);

// 导出当前进程全部 wisdom（plan 创建完成后调用）。
void export_wisdom(const std::string& file);

}  // namespace wisdom
