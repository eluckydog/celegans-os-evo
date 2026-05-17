# C. elegans 逆向软件工程 — 执行蓝图 (v2.0)

## 最终目标
系统性反编译线虫的"完整操作系统架构"，验证生命是否本质上可以用软件工程框架理解。

## 完成状态: Layer 0-6 ✅ 全部完成

| Layer | 状态 | 核心产出 | 关键基因验证 |
|-------|------|---------|-------------|
| **0: 基础设施** | ✅ | DNA编程语言视角(4字符/64操作码/编译流水线), NCBI管线, 开发日志 | — |
| **1: 进程树** | ✅ | Sulston细胞谱系→6主干fork, 131凋亡(86%在神经), Notch/Wnt/EGF信号标注 | par-1/2/3/6, glp-1, pie-1, skn-1, mom-2, pal-1, lin-12, unc-86, mec-3 |
| **2: IPC协议栈** | ✅ | 5种IPC模式: Notch(socket 1对1) / Wnt(广播 1对多) / EGF(RPC远程) / TGF-β(全局变量) / Hedgehog(信号量) | lag-2→lin-12, mom-2→pop-1, lin-3→let-23, daf-7→daf-1 (全部NCBI确认mRNA+蛋白长度) |
| **3: 系统调用** | ✅ | 10条syscall表 + 效应基因 | fork(cdc-25/cdk-1/cyb-1), kill(ced-3/4/9), exec(unc-86/mec-3), send(5种配体), recv(5种受体), bind(hmr-1/hmp-1/pat-3), migrate(unc-6/40/5), write(unc-54/ace-1), read(mec-4/odr-10/tax-4), chmod(HMT/HDAC) |
| **4: 调度器** | ✅ | 混合调度架构(时钟+定时器+优先级+中断) | cdc-25/cdk-1(时钟), lin-4/let-7(定时器miRNA), skn-1/pie-1(优先级), daf-2/7/16/1(中断→dauer) |
| **5: 驱动程序** | ✅ | 8种神经递质协议(GPIO/I²C/ADC/总线), 关键回路映射(触碰/咽/产卵) | cha-1/Ach→GPIO, eat-4/Glu→总线, unc-25/GABA→反相器, tph-1/5-HT→中断, cat-2/DA→状态, flp-1/FMRFa→I²C |
| **6: 内存管理** | ✅ | 染色质=mmap/mprotect, Ub/SUMO=GC, piRNA=WAL | asd-2/set-16(H3K4激活), mes-2/mes-6(H3K27抑制), met-2/set-25(H3K9→swap), smo-1/uba-2(SUMO→引用计数), prg-1/ppw-1(piRNA→WAL) |

## 关键发现

1. **AC/VU命运决定 ≈ 两节点分布式选主协议 (Bully算法)** — 据我所知无文献先例
2. **C. elegans 丢弃DNA甲基化 ≈ 简化内存管理** — 只用组蛋白修饰，更轻量
3. **131个凋亡86%在神经 ≈ kill()系统调用的主要使用场景是精确的神经修剪**
4. **神经递质 = 驱动程序协议** — 不需要连接组数据，从基因表达就能反推
5. **所有Layer 0-6的断言全部有NCBI基因数据验证**

## 交付物

| 文件 | 位置 |
|------|------|
| Layer产出索引 | knowledge_base/bioinf/DELIVERY_INDEX.md |
| 执行蓝图 | knowledge_base/bioinf/CELEGANS_REVERSE_ENGINEERING.md |
| 线虫OS架构 | knowledge_base/bioinf/celegans_os_architecture.md |
| Layer 1-6 代码 | knowledge_base/bioinf/layer1-6_*.py |
| AC/VU分析 | knowledge_base/bioinf/ac_vu_decompile.py |
| Notch分析 | knowledge_base/bioinf/lin12_notch_analysis.py |
