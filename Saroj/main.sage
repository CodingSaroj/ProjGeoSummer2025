import itertools
import matplotlib.pyplot as plt

EXP = 3
ORD = 2^EXP

F = GF(ORD)
Fl = F.list()
Fx = Fl[:1]
Fg = F.gens()

zero = Fl[0]
one = Fl[-1]
g = Fg[0]

GLF2 = [[a,b,c,d] for a,b,c,d in itertools.product(Fl,Fl,Fl,Fl) if a*d - b*c != zero]

#for m in GLF2:
#    print(m[0],m[1])
#    print(m[2],m[3])
#    print()

R = GF(ORD)['x,y']
x, y = R.gens()

pols = []

base_monomials = [one,x,y,x*y,x^2,y^2]

monomials = base_monomials.copy()

for e in range(1,EXP):
    for m in base_monomials:
        monomials.append((g^e) * m)

print(monomials)

for terms in Combinations(monomials):
    p = R(sum(terms))

    if p.degree() < 2:
        continue

    f = p.factor()

    if f[0][0] == p:
        pols.append(p)

print(f'{len(pols)} irreducible polynomials found.')

def elim_symmetry(pols, txt, symm):
    sym_elim_pols = pols.copy()

    for p in pols:
        if p not in sym_elim_pols:
            continue

        eqs = symm(p)

        for eq in eqs:
            if eq != p and eq in sym_elim_pols:
                sym_elim_pols.remove(eq)
                #print(f'{p} eq. {eq} with {txt}.')

    return sym_elim_pols

pols = elim_symmetry(pols, '(x,y) => (y,x)', lambda p: [p(y,x)])
print(f'{len(pols)} irreducible polynomials found.')

pols = elim_symmetry(pols, '(x,y) => (x+i,y+j)', lambda p: [p(x+i,y+j) for i,j in itertools.product(Fl, Fl)])
print(f'{len(pols)} irreducible polynomials found.')

pols = elim_symmetry(pols, '(x,y) => M(x,y) for M in GL', lambda p: [p(a*x+b*y + f,c*x+d*y + g) for (a,b,c,d),f,g in itertools.product(GLF2, Fl, Fl)])
print(f'{len(pols)} irreducible polynomials found.')

print(pols)

if input('Print graphs(y/n): ') != 'y':
    quit()

for p in pols:
    xa = []
    ya = []

    for y in range(ORD):
        for x in range(ORD):
            if p(Fl[x],Fl[y]) == 0:
                xa.append(x)
                ya.append(y)
    
    plt.scatter(xa, ya)
    plt.grid()
    plt.title(str(p))
    plt.xlim((-0.1,ORD - 0.9))
    plt.ylim((-0.1,ORD - 0.9))
    plt.show()
