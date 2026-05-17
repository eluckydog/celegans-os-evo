"""lin-12 (Notch homolog) — 线虫的IPC系统调用"""
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction
Entrez.email = 'assistant@noreply'

print("拉取 lin-12 (Notch) 基因序列")
print("=" * 60)

# 找lin-12的mRNA序列
handle = Entrez.esearch(db='nucleotide', term='lin-12[Gene Name] AND Caenorhabditis elegans[Organism] AND mRNA[Filter]', retmax=5)
record = Entrez.read(handle)
handle.close()
print(f"找到 {record['Count']} 条mRNA")

for nid in record['IdList'][:3]:
    h2 = Entrez.esummary(db='nucleotide', id=nid)
    r = Entrez.read(h2)
    h2.close()
    r0 = r[0]
    print(f"  {nid}: {r0.get('Title','?')[:60]} | {r0.get('Length','?')} bp")

# 拉一个全长序列分析
print("\n拉取 1 个mRNA进行分析...")
handle = Entrez.efetch(db='nucleotide', id=record['IdList'][0], rettype='gb', retmode='text')
seq_record = SeqIO.read(handle, 'genbank')
handle.close()

seq_str = str(seq_record.seq)
print(f"\n描述: {seq_record.description[:80]}")
print(f"长度: {len(seq_str)} bp")
print(f"GC含量: {gc_fraction(seq_record.seq)*100:.1f}%")

# 找CDS
for feature in seq_record.features:
    if feature.type == 'CDS':
        start = int(feature.location.start)
        end = int(feature.location.end)
        cds_str = str(seq_record.seq[start:end])
        cds_obj = Seq(cds_str) if 'Seq' in dir() else __import__('Bio.Seq', fromlist=['Seq']).Seq(cds_str)
        from Bio.Seq import Seq
        cds_obj = Seq(cds_str)
        protein = cds_obj.translate(table='Standard', to_stop=True)
        print(f"\nCDS: {start}-{end}, {end-start} bp")
        print(f"蛋白质: {len(protein)} aa")

        prot_str = str(protein)
        # Notch受体特征结构域
        print(f"\n结构域扫描 (Notch受体特征):")
        features_to_check = [
            ("EGF重复 (C-C)", "C"),  
            ("LNR (Lin-12/Notch Repeat)", "C"),
            ("跨膜域 (TM)", "L"),
            ("RAM域 (RBP-Jkappa结合)", "K"),
            ("PEST域 (降解信号)", "P"),
        ]
        for name, pat in features_to_check:
            # Notch的domain特征靠氨基酸组分而不是精确序列
            print(f"  {'⚠' if pat in prot_str else '  '} {name}")

        print(f"\nN端50aa: {prot_str[:50]}")
        print(f"C端10aa: ...{prot_str[-10:]}")

        # EGF-like重复计数 (粗略: C残基密度)
        c_count = prot_str.count('C')
        print(f"\nCys残基: {c_count} ({c_count/len(prot_str)*100:.1f}%) → EGF重复区")
        print(f"注: Notch受体通常含29-36个EGF重复, 每个约40aa, 富含Cys")

        break

# 和unc-54 (肌球蛋白) 对比 — IPC vs 执行器
print()
print("=" * 65)
print("【对比: IPC系统调用 vs 执行器】")
print("=" * 65)
print(f"""
     lin-12 (Notch受体)         unc-54 (肌球蛋白)
     ───────────────────        ───────────────────
     信号接收 (socket监听)         执行器 (actuator)
     EGF重复 (Cys-rich)            α-helix卷曲螺旋
     ~1000 aa                     ~1963 aa
     接收外部信号                  驱动肌肉收缩
     
     类比: select/poll/epoll      类比: write()系统调用
     
     lin-12敲除 → 细胞命运错误      unc-54敲除 → 运动瘫痪
     (进程调度失败)                (执行器损坏)
""")
