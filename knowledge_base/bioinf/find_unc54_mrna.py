"""查 unc-54 mRNA 版本"""
from Bio import Entrez
Entrez.email = 'assistant@noreply'

handle = Entrez.esearch(db='nucleotide', term='unc-54[Gene] AND Caenorhabditis elegans[Organism] AND mRNA[Filter]', retmax=5)
record = Entrez.read(handle)
handle.close()
print(f"Count: {record['Count']}")
for nid in record['IdList']:
    h2 = Entrez.esummary(db='nucleotide', id=nid)
    r = Entrez.read(h2)
    h2.close()
    r0 = r[0]
    title = r0.get('Title','?')[:60]
    length = r0.get('Length','?')
    source = r0.get('Source','?')
    print(f"  {nid}: {title} | {length} bp | {source}")
