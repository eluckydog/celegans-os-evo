"""
Layer 3: 系统调用 — 细胞行为的syscall映射

每个核心细胞行为 ≈ 一个系统调用。
有明确的"调用者"(上游信号)和"执行者"(效应基因)。
"""
from Bio import Entrez
Entrez.email = 'assistant@noreply'

print("=" * 65)
print("Layer 3: 系统调用表 (System Calls)")
print("=" * 65)

syscalls = {
    "fork()": {
        "description": "细胞分裂",
        "executor": "cdc-25.1, cdk-1, cyb-1 (细胞周期引擎)",
        "signal": "cyclin/CDK复合物 + 母源因子梯度",
        "result": "两个子细胞获得不同的命运决定因子",
        "analogy": "fork() — 创建子进程, 继承父进程状态后各走各路"
    },
    "exec()": {
        "description": "细胞分化 (加载新的基因表达程序)",
        "executor": "转录因子组合 (如unc-86, mec-3, lin-12)",
        "signal": "细胞间信号→转录因子入核→激活靶基因集",
        "result": "细胞表达新的功能蛋白, 完成功能转换",
        "analogy": "exec() — 替换当前进程的代码段"
    },
    "kill()": {
        "description": "程序性细胞死亡 (凋亡)",
        "executor": "ced-3 (caspase执行器), ced-4 (适配器)",
        "signal": "egl-1 (BH3) → 抑制ced-9 → 释放ced-4 → 激活ced-3",
        "result": "131个细胞被精确杀死, 被邻居吞噬清理",
        "analogy": "kill() — 向进程发SIGTERM"
    },
    "send()": {
        "description": "分泌信号分子 (细胞间通信)",
        "executor": "lag-2/apx-1 (Notch配体), mom-2 (Wnt), lin-3 (EGF), daf-7 (TGF-β)",
        "signal": "转录激活 → 蛋白折叠 → 分泌到胞外",
        "result": "信号分子到达邻细胞/远程细胞表面",
        "analogy": "send() — 向IPC通道写数据"
    },
    "recv()": {
        "description": "接收信号分子",
        "executor": "lin-12/glp-1 (Notch受体), pop-1 (Wnt), let-23 (EGF), daf-1 (TGF-β)",
        "signal": "配体结合 → 受体形成二聚体 → 胞内域切割 → 进入核",
        "result": "NOTCH: NICD入核; Wnt: β-catenin稳定 → 转录改变",
        "analogy": "recv() — 从IPC通道读数据"
    },
    "bind()": {
        "description": "细胞贴壁/细胞连接",
        "executor": "cadherin/catenin (hmr-1, hmp-1/2), integrin (pat-3, ina-1)",
        "signal": "钙粘蛋白→α/β-catenin → 骨架连接",
        "result": "细胞附着到细胞外基质或邻细胞",
        "analogy": "bind() — 把socket绑定到地址"
    },
    "migrate()": {
        "description": "细胞迁移",
        "executor": "netrin (unc-6), netrin受体 (unc-40/dcc, unc-5)",
        "signal": "吸引/排斥梯度 → 肌动蛋白骨架极化 → 迁移",
        "result": "细胞移动到目标位置",
        "analogy": "mv — 文件移动/重定位"
    },
    "read()": {
        "description": "感觉神经元探测环境",
        "executor": "化学受体 (chemosensory GPCRs ~500-800个), 机械受体 (mec-4, mec-10)",
        "signal": "环境刺激→受体激活→神经递质释放",
        "result": "线虫感知到食物/危险/温度/化学物质",
        "analogy": "read() — 从文件描述符读取数据"
    },
    "write()": {
        "description": "运动/动作",
        "executor": "unc-54 (肌球蛋白重链), unc-15 (副肌球蛋白), ace-1 (乙酰胆碱酯酶)",
        "signal": "神经元释放Ach → 肌肉收缩",
        "result": "线虫移动/进食/排便/产卵",
        "analogy": "write() — 向文件描述符写数据 (输出到世界)"
    },
    "chmod()": {
        "description": "表观遗传修饰",
        "executor": "DNMT 同源? (C. elegans丢失了DNA甲基化), HMT (组蛋白甲基化酶), HDAC (组蛋白去乙酰化酶)",
        "signal": "环境或发育信号 → 染色质重塑复合物",
        "result": "染色质状态改变 → 基因表达谱可遗传",
        "analogy": "chmod — 改变文件权限 (不改变内容, 改变可读性)"
    },
}

print(f"{'syscall':<16} {'↔':<3} {'生物过程':<16} {'效应基因':<28}")
print(f"{'─'*16} {'─'*3} {'─'*16} {'─'*28}")
for call_name, info in syscalls.items():
    genes_short = info['executor'].split('(')[0].strip()
    print(f"{call_name:<16} ↔  {info['description']:<16} {genes_short:<28}")

print()
print("=" * 65)
print("深度验证: cell cycle as fork()")
print("=" * 65)

# 查细胞周期关键基因
fork_genes = ["cdc-25.1", "cdk-1", "cyb-1"]
for gene in fork_genes:
    h = Entrez.esearch(db='gene', term=f'{gene}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
    r = Entrez.read(h)
    h.close()
    if int(r['Count']) > 0:
        h2 = Entrez.esummary(db='gene', id=r['IdList'][0])
        g = Entrez.read(h2)
        h2.close()
        info = g['DocumentSummarySet']['DocumentSummary'][0]
        desc = info.get('Description','?')
        print(f"  {gene}: {desc} — Chr {info.get('Chromosome','?')}")

print()
print("=" * 65)
print("Layer 3 Done ✅ — 10个系统调用已定义")
print("=" * 65)
print(f"""
系统调用表共 {len(syscalls)} 条:

写入(W): fork(), send(), write(), bind()
读取(R): recv(), read()
执行(X): exec(), migrate()
杀(K): kill()
元(M): chmod()

和Unix系统调用的结构一致:
- 进程管理: fork(), kill()
- IPC: send(), recv()
- 文件系统: read(), write()
- 权限: chmod()

区别:
- Unix的fork很少立即exec(); 线虫的fork几乎立即exec(分化)
- "回收"不是sync, 是邻居细胞主动吞噬 (APPROX GC)
""")
