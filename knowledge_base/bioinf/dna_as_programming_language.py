"""
DNA-as-Programming-Language 基础实操
把中心法则当成编译器流水线来理解
"""
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio import motifs

print("=" * 60)
print("【操作0】字符集 — 4个token, ATCG")
print("=" * 60)
# DNA序列就是一个字符串
dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print(f"序列: {dna}")
print(f"长度: {len(dna)} 碱基")
print(f"反向互补: {dna.reverse_complement()}")
print(f"转录为RNA: {dna.transcribe()}")
print()

print("=" * 60)
print("【操作1】编译 — 转录 + 翻译 (DNA → RNA → 蛋白质)")
print("=" * 60)
# 一个简单的基因
gene = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
mrna = gene.transcribe()
protein = gene.translate()
print(f"DNA:   {gene}")
print(f"mRNA:  {mrna}")
print(f"蛋白:  {protein}")
print(f"注意终止符: {gene.translate(to_stop=True)}")
print()

print("=" * 60)
print("【操作2】GC含量 — 代码的'温度'测量")
print("=" * 60)
# GC含量是DNA的一个重要指标——类似代码的复杂度度量
seq1 = Seq("ATGCGATCGATCGATCGTAGCTAGCTGATCG")
seq2 = Seq("ATATATATATATATATATATATATATATATA")
seq3 = Seq("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")
print(f"高GC序列: {gc_fraction(seq3)*100:.1f}% (稳定, 类似编译期检查多的语言)")
print(f"中等GC: {gc_fraction(seq1)*100:.1f}%")
print(f"低GC序列: {gc_fraction(seq2)*100:.1f}% (不稳定, 类似动态类型语言)")
print()

print("=" * 60)
print("【操作3】ORF — 找'函数定义'")
print("=" * 60)
# 一段长DNA里找所有可能的"函数"(ORF: open reading frame)
genomic_dna = Seq(
    "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    "TTAAGCTAGCCTAGCTAGCTAGCTAGCATGCGATCGA"
    "TAGCTAGCTAGCATGCTAGCTAGCATGCTAGCTAGCA"
)
# 6种读码框（3正向+3反向）
for frame in range(3):
    translated = genomic_dna[frame:].translate()
    print(f"读码框 {frame+1}: {translated}")
for frame in range(3):
    rev = genomic_dna.reverse_complement()
    translated = rev[frame:].translate()
    print(f"反向框 {frame+1}: {translated}")
print()

print("=" * 60)
print("【操作4】密码子表 — 指令集架构")
print("=" * 60)
# 64个密码子 → 20个氨基酸 + 1个起始 + 3个终止
# 类似CPU操作码表
from Bio.Data import CodonTable
table = CodonTable.unambiguous_dna_by_name["Standard"]
print(f"标准密码子表 (NC_000000):")
print(f"  起始密码子: {table.start_codons}")
print(f"  终止密码子: {table.stop_codons}")
# 看冗余度
codon_counts = {}
for codon, aa in table.forward_table.items():
    codon_counts[aa] = codon_counts.get(aa, 0) + 1
print(f"\n氨基酸冗余度分布:")
for aa, count in sorted(codon_counts.items(), key=lambda x: -x[1]):
    print(f"  {aa:>15}: {count}个密码子")
print()

print("=" * 60)
print("【操作5】突变 — 代码修改")
print("=" * 60)
import random
dna_list = list(str(gene))
pos = random.randint(0, len(dna_list)-1)
old_base = dna_list[pos]
new_base = random.choice([b for b in ["A","T","C","G"] if b != old_base])
dna_list[pos] = new_base
mutated = Seq("".join(dna_list))
print(f"原始: {gene}")
print(f"蛋白质: {gene.translate()}")
print(f"突变P{pos}: {old_base}→{new_base}")
print(f"突变后: {mutated}")
print(f"突变蛋白: {mutated.translate()}")
print(f"同义突变? {gene.translate() == mutated.translate()}")

print()
print("=" * 60)
print("【结论】DNA ≈ 一种4字符集的编译型语言")
print("=" * 60)
print(f"""
- token: ATCG (4个)
- 操作码: 64个密码子
- 基本指令: 翻译(3核苷酸→1氨基酸)
- 函数定义: ORF (起始→终止)
- 编译: 转录(DNA→mRNA) + 翻译(mRNA→蛋白)
- 版本控制: SNP是单字符修改, Indel是区块修改
- 内注释: 内含子(剪接掉), 外注释: 非编码区
- 动态reflection: 表观遗传(甲基化不改token但改行为)
""")
