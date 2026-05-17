"""
Layer 7: Bug追踪 — 将已知线虫突变表型映射为OS级的bug/补丁

每个已知线虫突变体 = 一个代码缺陷(bug)或补丁(patch)
"""
from Bio import Entrez
Entrez.email = 'assistant@noreply'

print("=" * 75)
print("C. elegans Operating System v1.0 — Bug Tracker")
print("映射: 已知线虫突变 → OS层缺陷报告")
print("=" * 75)

bugs = [
    # (基因, 表型, OS层, bug类型, 描述, 类比)
    ("lin-4(lf)", "异时性缺陷: L2重复重复", "Layer 4(定时器)", "逻辑错误",
     "lin-4是L2→L3的定时器。缺失时定时器不触发, L2阶段的程序反复执行",
     "while循环条件错误——计时器变量永远不会递减"),

    ("let-7(lf)", "L4→成虫过渡失败, 致死", "Layer 4(定时器)", "内核崩溃",
     "let-7是L4→成虫的定时器。缺失时无法完成最终发育转换, 系统崩溃",
     "没有let-7这个定时器信号, 最后一步永远不会执行, 系统挂死在L4阶段"),

    ("unc-86(lf)", "神经元分化失败", "Layer 3(exec)", "系统调用错误",
     "unc-86是神经前体细胞执行exec()所需的转录因子。缺失时神经元不分化",
     "exec() syscall缺失——子进程无法加载新程序, 进程是活的但不运行任何东西"),

    ("ced-3(lf)", "所有凋亡失败(131细胞不死)", "Layer 3(kill)", "内存泄漏",
     "ced-3是kill()的核心执行基因。缺失时131个应死的细胞全部存活",
     "kill() syscall返回错误但不报错——zombie进程永远不回收"),

    ("ced-3(go)", "细胞死亡提前/过多", "Layer 3(kill)", "权限错误",
     "ced-3过度活跃导致不应死的细胞也被杀",
     "kill()调用无权限检查——任何进程都能杀任何进程, 包括自己"),

    ("mec-4(d)", "机械感觉丧失", "Layer 5(read)", "驱动错误",
     "mec-4是触摸感受器的离子通道。突变时触摸信号永远读不到",
     "read()从'/dev/touch'返回空——传感器驱动坏了但系统认为一切正常"),

    ("unc-54(lf)", "瘫痪, 身体弯曲不能", "Layer 5(write)", "设备错误",
     "unc-54是肌球蛋白重链。缺失时肌肉不能收缩",
     "write()到'/dev/muscle'执行了但硬件无响应——电机坏了"),

    ("daf-2(lf)", "dauer无法退出", "Layer 4(中断)", "中断死锁",
     "daf-2是退出dauer的受体。缺失时即使食物出现也无法恢复",
     "中断信号来了但ISR(中断服务例程)无法执行——操作系统卡在suspend状态"),

    ("daf-2(go)", "从不进入dauer, 寿命缩短", "Layer 4(中断)", "临界区保护失败",
     "daf-2过度活跃导致系统从不挂起——虽然活得更短但幼体期加速",
     "ACPI S3被永久禁用——系统从不休眠, 电池耗尽更快"),

    ("lin-12(go)", "AC/VU: 两个都变VU", "Layer 2(Notch)", "分布式系统错误",
     "lin-12过度活跃导致两个细胞都激活Notch, 都选为VU, 无锚细胞",
     "分布式选主协议中两个节点都声称自己是leader——无leader"),

    ("lin-12(lf)", "AC/VU: 两个都变AC", "Layer 2(Notch)", "分布式系统错误",
     "lin-12缺失导致两个细胞都收不到Notch信号, 都选为AC",
     "分布式选主协议中两个节点都拒绝当leader——选举失败"),

    ("mig-2(lf)", "细胞迁移方向错误", "Layer 3(migrate)", "导航bug",
     "mig-2是RhoGTP酶, 指引细胞迁移方向。缺失时细胞往错误方向移动",
     "migrate()调用能执行但路径计算错误——GPS导航给了错误的目的地"),

    ("smgl-1(lf)", "SynMuv: 皮下细胞长成假性神经元", "Layer 6(mmap)", "权限提升漏洞",
     "smgl-1本是限制神经元命运的染色质抑制因子。缺失后皮下细胞获得神经命运",
     "mmap内存段的权限被放开——本应是只读的数据段变成了可执行"),

    ("prg-1(lf)", "生育力逐渐丧失(莫顿饥荒)", "Layer 6(WAL)", "日志损坏",
     "prg-1是piRNA通路核心。缺失时piRNA库逐渐丢失, 几代后转座子爆发导致不育",
     "WAL日志文件在重启后逐渐损坏——几代后数据库完全不可读"),

    ("hif-1(lf)", "低氧适应能力丧失", "Layer 4(中断)", "异常处理缺失",
     "hif-1响应低氧信号。缺失时缺氧不触发应急程序, 细胞死亡",
     "SIGTERM的信号处理函数是空的——进程收到终止信号但不做任何清理"),
]

# 按OS层分组统计
from collections import Counter
layer_counts = Counter()
type_counts = Counter()

for bug in bugs:
    layer = bug[2].replace("Layer ", "").split("(")[0]
    layer_counts[layer] += 1
    type_counts[bug[3]] += 1

print(f"\n总计bug: {len(bugs)}")
print(f"\n按OS层分布:")
for layer, count in sorted(layer_counts.items()):
    print(f"  Layer {layer}: {count}个")

print(f"\n按bug类型分布:")
for bt, count in sorted(type_counts.items(), key=lambda x: -x[1]):
    print(f"  {bt}: {count}个")

print()
print("-" * 75)
print(f"{'Bug ID':<8} {'基因':<12} {'OS层':<20} {'类型':<14} {'标题/类比':<20}")
print("-" * 75)

for i, (gene, pheno, layer, btype, desc, analogy) in enumerate(bugs, 1):
    bid = f"B-{i:02d}"
    short_pheno = pheno[:18]
    print(f"{bid:<8} {gene:<12} {layer:<20} {btype:<14} {short_pheno:<20}")
    print(f"  {'':<8} 类比: {analogy}")

print("=" * 75)

# 然后查几个典型突变体的详情
print("\n\n=== 核型bug验证 (NCBI) ===")
print("-" * 75)

verify_genes = [
    "lin-4", "let-7", "ced-3", "unc-86", "unc-54",
    "mec-4", "daf-2", "lin-12", "mig-2", "prg-1"
]

for gene in verify_genes:
    h = Entrez.esearch(db='gene', term=f'{gene}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
    r = Entrez.read(h)
    h.close()
    if int(r['Count']) > 0:
        h2 = Entrez.esummary(db='gene', id=r['IdList'][0])
        g = Entrez.read(h2)
        h2.close()
        info = g['DocumentSummarySet']['DocumentSummary'][0]
        chr_ = info.get('Chromosome','?')
        desc_ncbi = info.get('Description','?')[:50]
        summary = info.get('Summary','?')[:100]
        print(f"  {gene:>8} | Chr {chr_} | {desc_ncbi:<35}")
    else:
        print(f"  {gene:>8} | ❌ 未找到")

print()
print("=" * 75)
print("Layer 7 ✅ — Bug追踪系统已建立")
print("=" * 75)
print("""
关键发现:
1. 线虫的突变体表型几乎都可以精确映射为OS级bug
2. 每个bug的严重程度: ced-3(lf)=致命, let-7(lf)=致命, lin-4(lf)=非致命
3. 同一基因的gain/lf产生对称bug (如lin-12 go/lf = 两个AC/两个VU)
4. 最致命的bug: 系统调用核心基因 (fork/kill关键)
5. 最影响体验但不致命的: 感觉神经元驱动(mec-4)

核心洞察:
  "突变 → bug" 的映射不是比喻
  在"信息处理"层面它就是bug: 一个程序(生物体)的一个指令(基因)
  在特定条件下(环境)返回了错误结果(异常表型)
""")
