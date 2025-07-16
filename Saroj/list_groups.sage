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

zk = g

GLF2 = [[a,b,c,d] for a,b,c,d in itertools.product(Fl,Fl,Fl,Fl) if a*d - b*c != zero]

R = GF(ORD)['x,y']
x, y = R.gens()

pols = []

f = open(f'F{ORD}.txt')
lines = f.read().split('\n')[6:-1]

for line in lines:
    pols.append(sage_eval(line.replace('^','**'), preparse=False, locals={'x': x, 'y': y, f'z{EXP}': zk}))

def locus(f):
    sol_set = []

    for i,j in itertools.product(Fl,Fl):
        if f(i,j) == 0:
            sol_set.append((i,j))

    return sol_set

for conic in pols:
    pts = locus(conic)

    grp_ord = len(pts)

    Ox,Oy = pts[0]

    grp_table = [[-1 for j in pts] for i in pts]

    header = '* |'

    for i in range(grp_ord):
        header += f' {i:02d}'

    sep = '--+' + '-'*(len(header) - 3)

    print(f'{conic} has following group table:')
    print(header)
    print(sep)

    for j in range(grp_ord):
        row = f'{j:02d}|'

        for i in range(grp_ord):
            Px, Py = pts[i]
            Qx, Qy = pts[j]

            bx = Qx - Px
            by = Qy - Py

            if Px == Qx and Py == Qy:
                bx = conic.derivative(y)(Px,Py)
                by = conic.derivative(x)(Px,Py)

            l = [(Ox + t*bx, Oy + t*by) for t in Fl]

            C_int_l = [P for P in l if P in pts]

            if len(C_int_l) == 1:
                grp_table[j][i] = 0
            elif len(C_int_l) == 2:
                R = C_int_l[1] if C_int_l[0] == (Ox,Oy) else C_int_l[0]

                grp_table[j][i] = pts.index(R)

            row += f' {grp_table[j][i]:02d}'
        
        print(row)

    print()
