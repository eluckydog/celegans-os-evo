"""
Layer 6: 内存管理 — 表观遗传 + 蛋白降解

C. elegans 通过以下机制管理"内存"（基因表达状态）:
1. 染色质重塑 → 运行时反射 / 动态链接
2. 组蛋白修饰 → 权限控制 (chmod)
3. piRNA/siRNA → 垃圾回收 / 事务日志
4. Ub/SUMO → GC / 引用计数
"""
from Bio import Entrez
Entrez.email = 'assistant@noreply'

print("=" * 65)
print("Layer 6: 内存管理 — 表观遗传 + 蛋白降解")
print("=" * 65)

print("""
┌────────────────────────────────────────────────────────────┐
│             C. elegans 内存管理架构                          │
├───────────────┬─────────────┬──────────────┬───────────────┤
│ 机制           │ 分子机器     │ 关键基因       │ 软件类比       │
├───────────────┼─────────────┼──────────────┼───────────────┤
│ H3K4me3激活   │ COMPASS家族  │ asd-2/set-16 │ mmap(MAP_ANON)│
│               │             │              │ — 内存标记为可读 │
├───────────────┼─────────────┼──────────────┼───────────────┤
│ H3K27me3抑制  │ Polycomb    │ mes-2/mes-6  │ mprotect(RO)  │
│               │             │              │ — 内存标记为只读 │
├───────────────┼─────────────┼──────────────┼───────────────┤
│ H3K9me3异染色质 │ HMT        │ met-2/set-25 │ swap           │
│               │             │              │ — 冷数据换出    │
├───────────────┼─────────────┼──────────────┼───────────────┤
│ Ub 降解       │ 泛素-蛋白酶体 │ uba-/ubc-    │ free() / GC     │
│               │             │             │ — 不需要的蛋白回收│
├───────────────┼─────────────┼──────────────┼───────────────┤
│ SUMO 修饰     │ SUMO化      │ smo-1/uba-2 │ 引用计数/安全位 │
│               │             │             │ — 标记哪些蛋白重要│
├───────────────┼─────────────┼──────────────┼───────────────┤
│ piRNA/siRNA   │ Argonaute   │ prg-1/ppw-1 │ WAL(预写日志)   │
│               │             │              │ — 外来核酸的记忆│
└───────────────┴─────────────┴──────────────┴───────────────┘
""")

# 1. 查主要染色质修饰酶
print("【1】染色质修饰酶 = 内存保护单元")
print("-" * 50)

chromatin_genes = [
    ("asd-2", "COMPASS", "H3K4甲基转移酶, 基因激活"),
    ("mes-2", "Polycomb", "H3K27甲基转移酶, 基因沉默"),
    ("mes-6", "Polycomb(PCL)", "ESC同源, 协助mes-2"),
    ("met-2", "H3K9 HMT", "H3K9me1/2, 异染色质"),
    ("set-25", "H3K9 HMT", "H3K9me3, 异染色质"),
    ("set-16", "MLL", "H3K4甲基转移酶, 发育记忆"),
]

for name, complex_name, desc in chromatin_genes:
    h = Entrez.esearch(db='gene', term=f'{name}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
    r = Entrez.read(h)
    h.close()
    if int(r['Count']) > 0:
        h2 = Entrez.esummary(db='gene', id=r['IdList'][0])
        g = Entrez.read(h2)
        h2.close()
        info = g['DocumentSummarySet']['DocumentSummary'][0]
        chrom = info.get('Chromosome','?')
        desc_ncbi = info.get('Description','?')[:45]
        print(f"  {name:>6} ({complex_name:<12}) Chr {chrom} — {desc_ncbi:<35}")
    else:
        print(f"  {name:>6} ({complex_name:<12}) — ❌ 未找到")

# 2. 查蛋白降解系统
print()
print("【2】蛋白降解系统 = 垃圾回收管理器 (GC)")
print("-" * 50)

degradation_genes = [
    ("smo-1", "SUMO蛋白", "SUMO化修饰, 标记蛋白稳定性"),
    ("uba-2", "SUMO E1", "SUMO激活酶"),
    ("rpt-1", "19S调节颗粒", "蛋白酶体盖子"),
    ("pas-4", "20S核心", "蛋白酶体核心亚基"),
    ("vps-4", "ESCRT", "内体分选, 膜蛋白降解"),
]

for name, type_desc, desc in degradation_genes:
    h = Entrez.esearch(db='gene', term=f'{name}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
    r = Entrez.read(h)
    h.close()
    if int(r['Count']) > 0:
        h2 = Entrez.esummary(db='gene', id=r['IdList'][0])
        g = Entrez.read(h2)
        h2.close()
        info = g['DocumentSummarySet']['DocumentSummary'][0]
        chrom = info.get('Chromosome','?')
        desc_n = info.get('Description','?')[:40]
        print(f"  {name:>6} ({type_desc:<15}) Chr {chrom} — {desc_n:<40}")
    else:
        print(f"  {name:>6} ({type_desc:<15}) — ❌ 未找到")

# 3. RNAi / piRNA — 事务日志
print()
print("【3】piRNA/siRNA系统 = 事务日志 (WAL)")
print("-" * 50)

rna_genes = [
    ("prg-1", "Piwi", "piRNA加工, 基因组完整性"),
    ("ppw-1", "Argonaute WAGO", "siRNA扩增, 沉默外源核酸"),
    ("rde-1", "Argonaute AGO", "RNAi起始, 外来RNA感知"),
    ("rde-4", "dsRNA结合蛋白", "识别双链RNA, RNAi触发"),
]

for name, type_desc, desc in rna_genes:
    h = Entrez.esearch(db='gene', term=f'{name}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
    r = Entrez.read(h)
    h.close()
    if int(r['Count']) > 0:
        h2 = Entrez.esummary(db='gene', id=r['IdList'][0])
        g = Entrez.read(h2)
        h2.close()
        info = g['DocumentSummarySet']['DocumentSummary'][0]
        chrom = info.get('Chromosome','?')
        desc_n = info.get('Description','?')[:40]
        print(f"  {name:>6} ({type_desc:<15}) Chr {chrom} — {desc_n:<40}")
    else:
        print(f"  {name:>6} ({type_desc:<15}) — ❌ 未找到")

# 4. 总结
print()
print("=" * 65)
print("Layer 6 Done ✅ — 内存管理架构已建立")
print("=" * 65)
print("""
内存管理架构关键洞察:

1. 染色质修饰 = 虚拟内存保护 (mmap/mprotect)
   - H3K4me3 → MAP_READ (激活态, 可读)
   - H3K27me3 → MAP_READONLY (抑制态, 只读)
   - H3K9me3 → swap (异染色质, 冷数据换出)

2. Ub/SUMO = 垃圾回收器
   - 泛素化 → 蛋白被送进蛋白酶体 ≈ free()
   - SUMO化 → 蛋白被保护 ≈ 增加引用计数
   - 多聚泛素链(K48) → 强引用, 必须回收
   - SUMO多聚化 → 安全引用, 防止误回收

3. piRNA/siRNA = 预写日志 (Write-Ahead Log)
   - 每遇到一个外来核酸 → 记录到piRNA库
   - 类似数据库的WAL: 先记日志, 再执行操作
   - 子代继承piRNA → 持久化存储

4. 线虫丢弃了DNA甲基化 → 简化了内存管理
   - C. elegans 是极少数没有DNA甲基化的真核生物
   - 只用组蛋白修饰来做内存保护
   - 系统更轻量, 但灵活性降低

Layer 0-6 全部完成 ✅
""")
