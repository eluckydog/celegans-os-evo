"""发育事件的实时"反编译"

目标: AC/VU 命运决定 (anchor cell vs ventral uterine precursor)
这是线虫中研究最透彻的Notch侧向抑制案例。
两个同等细胞, 通过lin-12/Notch信号, 一个变成AC, 另一个变成VU。

流程类似: 两个等权进程, 通过socket通信争抢master角色。
"""
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction
Entrez.email = 'assistant@noreply'

print("=" * 65)
print("AC/VU 命运决定 — 反编译一个发育事件")
print("=" * 65)

print("""
参与进程 (两个细胞):
  Z1.ppp 和 Z4.aaa  (两个生殖腺前体细胞)
  
在正常发育中:
  - Z1.ppp → AC (锚定细胞)   [master]
  - Z4.aaa → VU (子宫前体)   [slave]
  
但如果其中一个被激光杀死, 另一个自动变成AC。
意味着两个细胞都是等价的, Notch信号打破对称性。

机制:
  ┌──────────────────────────────────┐
  │两个细胞都表达:                     │
  │  lin-12 (Notch受体)               │
  │  lag-2 (DSL配体, 信号发送)         │
  │                                   │
  │竞争开始:                          │
  │  细胞A: 稍多LAG-2 → 发送信号给B    │
  │  细胞B: 接收信号 → LIN-12激活      │
  │        → 下调LAG-2表达            │
  │        → 失去竞争力               │
  │  细胞A: LAG-2积累 → 信号增强       │
  │        → 分化成AC                 │
  │  细胞B: LIN-12持续 → 分化成VU      │
  │         (Notch的胞内域NICD进入核)  │
  └──────────────────────────────────┘
""")

# 查lag-2 (DSL配体, Notch信号的"发送端")
print("【1】查 lag-2 (Notch配体, ≈ send())")
handle = Entrez.esearch(db='gene', term='lag-2[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
record = Entrez.read(handle)
handle.close()
if int(record['Count']) > 0:
    h2 = Entrez.esummary(db='gene', id=record['IdList'][0])
    g = Entrez.read(h2)
    h2.close()
    r = g['DocumentSummarySet']['DocumentSummary'][0]
    print(f"  GeneID: {record['IdList'][0]}")
    print(f"  描述: {r.get('Description','?')}")
    print(f"  染色体: {r.get('Chromosome','?')}, 位置: {r.get('MapLocation','?')}")

    # 拉lag-2 mRNA
    handle_n = Entrez.esearch(db='nucleotide', term='lag-2[Gene Name] AND Caenorhabditis elegans[Organism] AND mRNA[Filter]', retmax=2)
    nrec = Entrez.read(handle_n)
    handle_n.close()
    print(f"  mRNA记录: {nrec['Count']} 条")
    for nid in nrec['IdList'][:1]:
        h3 = Entrez.efetch(db='nucleotide', id=nid, rettype='gb', retmode='text')
        sr = SeqIO.read(h3, 'genbank')
        h3.close()
        seq_str = str(sr.seq)
        print(f"  mRNA长度: {len(seq_str)} bp")
        for f in sr.features:
            if f.type == 'CDS':
                cds = str(sr.seq[int(f.location.start):int(f.location.end)])
                from Bio.Seq import Seq
                prot = Seq(cds).translate(to_stop=True)
                print(f"  蛋白质: {len(prot)} aa")
                print(f"  N端30: {str(prot)[:30]}")
                break

print()
print("【2】查 lin-12 (Notch受体, ≈ recv()) -- 已在之前的分析完成")
print("  lin-12: 1429 aa, 114个Cys, 4555 bp mRNA")

print()
print("【3】查 lst-1 / sygl-1 (Notch通路下游靶标)")
print("  Notch激活后, 胞内域NICD进入核内,")
print("  结合转录因子LAG-1 → 激活Hes/Hey家族基因")

handle = Entrez.esearch(db='gene', term='ref-1[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
record = Entrez.read(handle)
handle.close()
if int(record['Count']) > 0:
    h2 = Entrez.esummary(db='gene', id=record['IdList'][0])
    g = Entrez.read(h2)
    h2.close()
    r = g['DocumentSummarySet']['DocumentSummary'][0]
    print(f"  ref-1 (Hes-like): {r.get('Description','?')}, Chr {r.get('Chromosome','?')}")

print()
print("=" * 65)
print("【反编译总结】AC/VU命运决定 ≈ 选主协议")
print("=" * 65)

print("""
  AC/VU命运决定 = 分布式系统中的Bully算法 (领导选举)

  ┌─────────────────────────────────────────────────────────┐
  │  算法: 领导选举 (Leader Election)                        │
  │                                                         │
  │  参与者: Z1.ppp, Z4.aaa (2个对等节点)                   │
  │                                                         │
  │  socket.send  = lag-2 (DSL配体, Notch信号发送端)         │
  │  socket.recv  = lin-12 (Notch受体)                      │
  │                                                         │
  │  while True:                                            │
  │      if my_lag2 > other_lag2:                           │
  │          send_signal(other, "我更强")                    │
  │          differentiate('AC')     # master                │
  │          break                                           │
  │      else:                                              │
  │          receive_signal()                               │
  │          differentiate('VU')     # slave                 │
  │          break                                           │
  │                                                         │
  │  初始破对称:                                             │
  │    随机涨落 + 正反馈 (哪个细胞多表达了一点lag-2)         │
  │                                                         │
  │  没有死锁: 一定是之一, 不是都AC也不是都VU                │
  │  (contrast: 两个都表达Notch, 但只有一方能赢)            │
  │                                                         │
  │  ⏱ 发育时间: 线虫受精后约7小时, 精确到分钟              │
  └─────────────────────────────────────────────────────────┘

  这还只是一个事件。线虫整个发育程序中,
  这种细胞间信号通信发生了无数次。
""")

print("=" * 65)
print("【后续建议 — 丢给工程化AI】")
print("=" * 65)
print("""
  1. 把上述"反编译"写成可执行的Python模块:
     - 拉WormBase API获取发育数据和细胞谱系
     - 可视化fork-tree (cellular pedigree)
     - 标注每个fork的基因表达切换
  
  2. 交互式"线虫OS进程浏览器":
     - 点击细胞名 → 显示其命运决定的关键基因
     - 显示Notch信号的时间线
     - 支持"敲除"特定基因后预测发育异常
  
  3. 工具化:
     - 把Seq翻译+基因分析包装成dna_tools.py (挂到C_meta)
     - 定期拉C. elegans最新注释
""")

if __name__ == "__main__":
    main()


def main():
    print("AC/VU选举算法反编译工具")
    print("运行run_evolution.py开始演化模拟")
    print("或导入本模块使用decompile_ac_vu()函数")
