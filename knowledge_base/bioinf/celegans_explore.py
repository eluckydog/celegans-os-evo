"""C. elegans 基因组探索"""
from Bio import Entrez
Entrez.email = 'assistant@noreply'

print("=" * 60)
print("【1】C. elegans 组装信息")
print("=" * 60)

handle = Entrez.esearch(db='assembly', term='Caenorhabditis elegans[Organism] AND latest[filter]', retmax=3)
record = Entrez.read(handle)
handle.close()
print(f"找到 {record['Count']} 个组装")
for aid in record['IdList'][:2]:
    handle2 = Entrez.esummary(db='assembly', id=aid)
    rec = Entrez.read(handle2)
    handle2.close()
    r = rec['DocumentSummarySet']['DocumentSummary'][0]
    print(f"  组装名: {r.get('AssemblyName', '?')}")
    print(f"  物种: {r.get('Organism', '?')}")
    print(f"  大小: {r.get('TotalSequenceLength', '?')} bp")
    print(f"  染色体数: {r.get('Count', '?')}")

print()
print("=" * 60)
print("【2】unc-54 基因 — 线虫版Hello World")
print("=" * 60)

handle = Entrez.esearch(db='gene', term='unc-54[Gene Name] AND Caenorhabditis elegans[Organism]', retmax=3)
record = Entrez.read(handle)
handle.close()
print(f"找到 {record['Count']} 条")
for gid in record['IdList']:
    handle2 = Entrez.esummary(db='gene', id=gid)
    rec = Entrez.read(handle2)
    handle2.close()
    r = rec['DocumentSummarySet']['DocumentSummary'][0]
    print(f"  名称: {r.get('Name', '?')} — {r.get('Description', '?')}")
    print(f"  染色体: {r.get('Chromosome', '?')}, 位置: {r.get('MapLocation', '?')}")
    print(f"  别名: {r.get('OtherAliases', '?')}")

print()
print("=" * 60)
print("【3】拉 unc-54 CDS 序列")
print("=" * 60)

handle = Entrez.esearch(db='nucleotide', term='unc-54[Gene Name] AND Caenorhabditis elegans[Organism] AND CDS[Feature Key]', retmax=5)
record = Entrez.read(handle)
handle.close()
print(f"找到 {record['Count']} 条")
for nid in record['IdList'][:3]:
    handle2 = Entrez.esummary(db='nucleotide', id=nid)
    rec = Entrez.read(handle2)
    handle2.close()
    r = rec[0]
    print(f"  {r.get('Title', '?')[:90]}")
    print(f"    长度: {r.get('Length', '?')} bp")
