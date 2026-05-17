# 知识库 — C. elegans 操作系统反编译

本目录包含从操作系统视角反编译线虫发育的全部分层文件。
全部基因验证来自NCBI Entrez公开数据，无需外部API依赖。

## 结构

| 层 | 文件 | 对应OS概念 |
|----|------|-----------|
| 架构总览 | celegans_os_v1.0.md, celegans_os_architecture.md | 整体框架 |
| Layer 0 | (内置在架构文档中) | 引导加载器 |
| Layer 1 | layer1_proc_tree_v2.py | 进程管理 |
| Layer 2 | layer2_ipc_stack.py | 信号通路/IPC |
| Layer 3 | layer3_syscalls.py | 系统调用 |
| Layer 4 | layer4_scheduler.py | 调度器 |
| Layer 5 | layer5_from_genes.py | 驱动程序 |
| Layer 6 | layer6_memory_manager.py | 内存管理 |
| Layer 7 | layer7_bug_tracker.py | Bug追踪 |
| AC/VU | ac_vu_decompile.py | 分布式选主算法反编译 |
| DNA | dna_as_programming_language.py | DNA编程语言框架 |
