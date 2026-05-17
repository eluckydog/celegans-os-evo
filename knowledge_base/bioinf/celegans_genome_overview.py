"""C. elegans 全基因组概览 — 操作系统级视角"""

from Bio import Entrez
Entrez.email = 'assistant@noreply'

print("=" * 60)
print("【1】C. elegans 全基因组基本信息")
print("=" * 60)

# 查基因组组装
handle = Entrez.esearch(db='assembly', term='Caenorhabditis elegans[Organism] AND "WBcel235"[Assembly Name]', retmax=3)
record = Entrez.read(handle)
handle.close()
print(f"找到 {record['Count']} 个组装")
for aid in record['IdList'][:2]:
    h2 = Entrez.esummary(db='assembly', id=aid)
    rec = Entrez.read(h2)
    h2.close()
    r = rec['DocumentSummarySet']['DocumentSummary'][0]
    print(f"  组装名: {r.get('AssemblyName','?')}")
    print(f"  物种: {r.get('Organism','?')}")
    print(f"  染色体数: {r.get('Count','?')}")
    print(f"  总长: {r.get('TotalSequenceLength','?')} bp")

print()
print("=" * 60)
print("【2】基因统计 — '系统调用'数量")
print("=" * 60)

# 查基因库中的C. elegans基因总数
handle = Entrez.esearch(db='gene', term='Caenorhabditis elegans[Organism] AND protein coding[Gene Type]', retmax=1)
record = Entrez.read(handle)
handle.close()
print(f"预测的蛋白编码基因数: {record['Count']}")

# 也查一下全部
handle = Entrez.esearch(db='gene', term='Caenorhabditis elegans[Organism]', retmax=1)
record = Entrez.read(handle)
handle.close()
print(f"总记录(含ncRNA/假基因等): {record['Count']}")

print()
print("=" * 60)
print("【3】关键发育调控基因速查 — 操作系统的'启动脚本'")
print("=" * 60)

# 查几个经典发育基因
development_genes = [
    ('lin-4', 'microRNA pioneer'),
    ('let-7', 'microRNA'),
    ('lin-12', 'Notch homolog'),
    ('glp-1', 'Notch/Germline'),
    ('ced-3', 'caspase'),
    ('ced-4', 'apoptosis'),
    ('egl-1', 'BH3-only'),
    ('tra-1', 'sex determination'),
    ('him-8', 'meiosis'),
    ('unc-86', 'POU transcription factor'),
    ('mec-3', 'LIM homeobox'),
]

for name, desc in development_genes:
    h = Entrez.esearch(db='gene', term=f'{name.split()[0]}[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=1)
    rec = Entrez.read(h)
    h.close()
    count = rec['Count']
    status = ""
    if int(count) > 0:
        h2 = Entrez.esummary(db='gene', id=rec['IdList'][0])
        g = Entrez.read(h2)
        h2.close()
        r = g['DocumentSummarySet']['DocumentSummary'][0]
        status += f" | Chr: {r.get('Chromosome','?')}, 描述: {r.get('Description','?')[:40]}"
    print(f"  {name.split()[0]:>8} ({desc:<30}) → {count}条 {status}")
