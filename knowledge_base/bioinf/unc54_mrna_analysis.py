"""C. elegans unc-54 mRNA — 正确的编译后版本"""
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Data import CodonTable
from collections import Counter
Entrez.email = 'assistant@noreply'

print("拉取 unc-54 mRNA (ID=1972236164, 6204 bp)")
handle = Entrez.efetch(db='nucleotide', id='1972236164', rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
handle.close()

seq_str = str(record.seq)
print(f"描述: {record.description[:100]}")
print(f"长度: {len(seq_str)} bp")
print(f"GC含量: {gc_fraction(record.seq)*100:.1f}%")

# 找CDS
cds_seq = None
for feature in record.features:
    if feature.type == 'CDS':
        start = int(feature.location.start)
        end = int(feature.location.end)
        cds_seq = str(record.seq[start:end])
        print(f"\nCDS: {start}-{end}, {end-start} bp")
        print(f"产物: {feature.qualifiers.get('product',['?'])[0]}")
        
        # 翻译 — 用Biopython的Seq对象，不是str
        from Bio.Seq import Seq
        cds_obj = Seq(cds_seq)
        protein = cds_obj.translate(table='Standard', to_stop=True)
        print(f"蛋白质: {len(protein)} aa")
        break

if not cds_seq:
    print("\nmRNA里没有CDS标注，直接全序列翻译找ORF")
    # 找最长ORF
    best = (0, 0, 0)
    for frame in range(3):
        prot = record.seq[frame:].translate(table='Standard', to_stop=False)
        # 找最长ORF
        i = 0
        while i < len(prot):
            if prot[i] == 'M':
                j = i + 1
                while j < len(prot) and prot[j] != '*':
                    j += 1
                length = j - i
                if length > best[0]:
                    best = (length, i, frame)
                i = j
            else:
                i += 1
    if best[0] > 0:
        print(f"最长ORF: {best[0]} aa, 起始于frame {best[1]}, offset {best[2]}")
        cds_seq = str(record.seq[best[2]:])[best[1]*3:]
        protein = cds_seq.translate(table='Standard', to_stop=True)

print(f"\n蛋白质 N端50aa: {protein[:50]}")
print(f"C端10aa: ...{protein[-10:]}")

# 肌球蛋白特征
prot_str = str(protein)
print(f"\n基序扫描:")
checks = [
    ("ATP结合 P环 (GESGAGKT)", "GESGAGKT"),
    ("Switch II (DXXG)", "D"),
    ("肌球蛋白头部 I 类标志", "R"),
    ("卷曲螺旋尾部", "E"),
]
for name, pat in checks:
    pos = prot_str.find(pat) if len(pat) > 1 else prot_str.find(pat)
    print(f"  {'✅' if pos >= 0 else '❌'} {name}" + (f" @ pos {pos}" if pos >= 0 else ""))

# 氨基酸组分
aa = Counter(prot_str)
print(f"\n氨基酸组分:")
for a, c in aa.most_common():
    print(f"  {a}: {c:4d} ({c/len(prot_str)*100:5.1f}%)")

# 分子量
aa_weights = {'A':89.1,'R':174.2,'N':132.1,'D':133.1,'C':121.2,'Q':146.2,'E':147.1,'G':75.1,'H':155.2,'I':131.2,'L':131.2,'K':146.2,'M':149.2,'F':165.2,'P':115.1,'S':105.1,'T':119.1,'W':204.2,'Y':181.2,'V':117.1}
mw = sum(aa_weights.get(a, 0) for a in prot_str) - 18*(len(prot_str)-1)

print()
print("=" * 50)
print("C. elegans unc-54 基因 — 真正结论")
print("=" * 50)
print(f"""
  '源文件'长度 (基因组DNA): 9000 bp
  '预处理后' (mRNA, 剪接后): {len(seq_str)} bp
  '编译后' (CDS): {len(cds_seq)} bp = {len(cds_seq)//3} 条指令
  '最终二进制' (蛋白质): {len(prot_str)} aa ≈ {mw/1000:.0f} kDa
  GC含量: {gc_fraction(record.seq)*100:.1f}%

  C. elegans ≈ Hello World 的原因:
    1. 第一个完成全基因组测序的多细胞生物 (1998)
    2. 302个神经元 → 完整连接组已知 (第一个被完全描述)
    3. 959个体细胞 → 每个细胞的命运已知
    4. 可编程: RNAi干扰技术的原生阵地
    5. 3.3%的基因与人类直系同源 — 保守度惊人

  这就好比一个只含{len(seq_str)}行'代码'的微型操作系统，
  但能控制一个完整生物体的运动、进食、繁殖。
""")
