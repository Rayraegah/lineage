from lineage import breed

matron_genes = 0x000063169218f348dc640d171b000208934b5a90189038cb3084624a50f7316c
sire_genes = 0x00005a13429085339c6521ef0300011c82438c628cc431a63298e3721f772d29
matron_cdeb = 0x000000000000000000000000000000000000000000000000000000000047ff27
blockhash = 0xf9dd4486d68b13839d2f7b345f5223f17abae39a951f2cea5b0ca0dd6dc8db83

outs = breed(matron_genes, sire_genes, matron_cdeb, blockhash)
outs2 = 0x5b174298a44b9c6521176000021c53734c9018c431a73298674a5177316c

print(hex(outs))
print(hex(outs2))