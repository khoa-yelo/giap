"""Regenerate genome.fa and proteins.faa of the tiny tblastn fixture (see PROVENANCE.md).

Run ``python make_fixture.py``: writes genome.fa and prot.faa (stored here as proteins.faa)
in the current directory.
"""

import random

CODE = {}
b = "TCAG"
aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for x in b:
    for y in b:
        for z in b:
            CODE[x + y + z] = aa[i]
            i += 1


def tr(s):
    return "".join(CODE[s[i : i + 3]] for i in range(0, len(s) - 2, 3))


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


sense = [c for c in CODE if CODE[c] != "*"]


def orf(n):
    return "ATG" + "".join(random.choice(sense) for _ in range(n - 1))


def rnd(n):
    return "".join(random.choice("ACGT") for _ in range(n))


random.seed(11)
L = 8000
g = list(rnd(L))


def put(pos, s):
    for k, c in enumerate(s):
        g[(pos + k) % L] = c


def mutate(p, frac):
    cod = [p[i : i + 3] for i in range(0, len(p), 3)]
    for k in range(1, len(cod)):
        if random.random() < frac:
            cod[k] = random.choice([c for c in sense if CODE[c] != CODE[cod[k]]])
    return "".join(cod)


p1 = orf(120)
put(1000, p1)  # plus
p2 = orf(120)
put(2000, rc(p2))  # minus [2000,2360)
p3 = orf(120)
put(3000, p3[:180] + "TAACTGA" + p3[180:])  # engineered +1 frameshift (7 nt), stops both frames
p4 = orf(100)
put(L - 150, p4)  # origin crossing plus
p6 = orf(80)
put(4000, p6)
put(5000, p6)  # two loci
p7 = orf(150)
put(5500, mutate(p7, 0.12))  # divergent ~88%
p8 = orf(120)
put(6500, p8[:180] + "A" + p8[180:])  # natural +1 frameshift
chr1 = "".join(g)
M = 3000
h = list(rnd(M))
p5 = orf(100)
r = rc(p5)
for k, c in enumerate(r):
    h[(M - 120 + k) % M] = c
with open("genome.fa", "w") as fg:
    fg.write(">chr1\n" + chr1 + "\n>pl1\n" + "".join(h) + "\n")
with open("prot.faa", "w") as f:
    for n, s in (
        ("P1", p1),
        ("P2", p2),
        ("P3", p3),
        ("P4", p4),
        ("P5", p5),
        ("P6", p6),
        ("P7", p7),
        ("P8", p8),
        ("sp|Q12345|X_ECOLI", p1[:150]),
    ):
        f.write(f">{n}\n{tr(s)}\n")
