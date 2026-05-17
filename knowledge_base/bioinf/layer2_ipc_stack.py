"""
Layer 2: IPC协议栈 — 信号通路映射

C. elegans 细胞间通信 ≈ IPC协议栈。
每条信号通路对应一种通信模式。
"""
from Bio import Entrez
Entrez.email = 'assistant@noreply'

print("=" * 65)
print("Layer 2: IPC协议栈 — 信号通对照表")
print("=" * 65)

print("""
┌────────────────────────────────────────────────────────────────────────┐
│                  C. elegans IPC 协议栈                                  │
├────────────┬──────────────┬──────────────┬───────────┬───────────────┤
│ IPC模式     │ 信号通路      │ 发送端(基因)  │ 接收端(基因) │ 用途          │
├────────────┼──────────────┼──────────────┼───────────┼───────────────┤
│ socket      │ Notch        │ lag-2/apx-1  │ lin-12/glp-1│ 相邻细胞选主  │
│ (1对1)      │              │              │            │               │
├────────────┼──────────────┼──────────────┼───────────┼───────────────┤
│ 广播        │ Wnt          │ mom-2/lin-44 │ pop-1/wrm-1│ 前后轴+patterning│
│ (1对多)     │              │              │            │               │
├────────────┼──────────────┼──────────────┼───────────┼───────────────┤
│ RPC         │ EGF/RTK     │ lin-3        │ let-23     │ 远程诱导     │
│ (远程1对1)  │              │              │            │ (如神经诱导)  │
├────────────┼──────────────┼──────────────┼───────────┼───────────────┤
│ 全局变量    │ TGF-β        │ daf-7        │ daf-1/daf-4│ 环境感知     │
│ (系统范围)  │              │ (神经分泌)    │            │ (dauer决定)  │
├────────────┼──────────────┼──────────────┼───────────┼───────────────┤
│ 信号量      │ Hedgehog     │ grd/wrt      │ ptc-1     │ 细胞数量控制 │
│ (配额)      │              │              │            │               │
└────────────┴──────────────┴──────────────┴───────────┴───────────────┘
""")

# 现在逐一拉取各通路的关键基因数据
pathways = {
    "Notch": {"send": "lag-2", "recv": "lin-12", "purpose": "相邻细胞, 侧向抑制(✅已分析)"},
    "Wnt": {"send": "mom-2", "recv": "pop-1", "purpose": "前后轴, 内胚层诱导"},
    "EGF": {"send": "lin-3", "recv": "let-23", "purpose": "远程信号, 神经诱导"},
    "TGF-beta": {"send": "daf-7", "recv": "daf-1", "purpose": "环境感知, dauer决定"},
}

for pathway_name, genes in pathways.items():
    print(f"\n{'='*40}")
    print(f"【{pathway_name} 通路】")
    print(f"{'='*40}")
    if genes.get('send'):
        print(f"发送端 ('send' syscall): {genes['send']}")
        h = Entrez.esearch(db='gene', term=f'{genes["send"]}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
        r = Entrez.read(h)
        h.close()
        if int(r['Count']) > 0:
            h2 = Entrez.esummary(db='gene', id=r['IdList'][0])
            g = Entrez.read(h2)
            h2.close()
            info = g['DocumentSummarySet']['DocumentSummary'][0]
            print(f"  描述: {info.get('Description','?')}")
            print(f"  染色体: {info.get('Chromosome','?')}")
            # 尝试拉蛋白长度
            try:
                hn = Entrez.esearch(db='nucleotide', term=f'{genes["send"]}[Gene Name] AND C. elegans[Organism] AND mRNA[Filter]', retmax=1)
                nr = Entrez.read(hn)
                hn.close()
                if int(nr['Count']) > 0:
                    hn2 = Entrez.efetch(db='nucleotide', id=nr['IdList'][0], rettype='gb', retmode='text')
                    from Bio import SeqIO
                    seq_r = SeqIO.read(hn2, 'genbank')
                    hn2.close()
                    for f in seq_r.features:
                        if f.type == 'CDS':
                            from Bio.Seq import Seq
                            cds = str(seq_r.seq[int(f.location.start):int(f.location.end)])
                            prot = Seq(cds).translate(to_stop=True)
                            print(f"  mRNA: {len(seq_r)}bp → 蛋白: {len(prot)}aa")
                            break
            except:
                print(f"  (蛋白长度: 无法获取)")
    
    if genes.get('recv'):
        print(f"接收端 ('recv' syscall): {genes['recv']}")
        h = Entrez.esearch(db='gene', term=f'{genes["recv"]}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
        r = Entrez.read(h)
        h.close()
        if int(r['Count']) > 0:
            h2 = Entrez.esummary(db='gene', id=r['IdList'][0])
            g = Entrez.read(h2)
            h2.close()
            info = g['DocumentSummarySet']['DocumentSummary'][0]
            print(f"  描述: {info.get('Description','?')}")
            print(f"  染色体: {info.get('Chromosome','?')}")
            try:
                hn = Entrez.esearch(db='nucleotide', term=f'{genes["recv"]}[Gene Name] AND C. elegans[Organism] AND mRNA[Filter]', retmax=1)
                nr = Entrez.read(hn)
                hn.close()
                if int(nr['Count']) > 0:
                    hn2 = Entrez.efetch(db='nucleotide', id=nr['IdList'][0], rettype='gb', retmode='text')
                    seq_r = SeqIO.read(hn2, 'genbank')
                    hn2.close()
                    for f in seq_r.features:
                        if f.type == 'CDS':
                            cds = str(seq_r.seq[int(f.location.start):int(f.location.end)])
                            prot = Seq(cds).translate(to_stop=True)
                            print(f"  mRNA: {len(seq_r)}bp → 蛋白: {len(prot)}aa")
                            break
            except:
                print(f"  (蛋白长度: 无法获取)")
    
    print(f"用途: {genes['purpose']}")
