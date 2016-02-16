import sys

if len(sys.argv) < 4:
    print("from, to, step")
    exit()

_from = float(sys.argv[1])
_to   = float(sys.argv[2])
_step = float(sys.argv[3])

if _from >= _to:
    print("error\nfrom >= to")
    exit()

def frange(x, y, step):
    while x < y:
        yield x
        x += step

size = 200
N = 350
eta_range = frange(_from, _to+1, _step)
threshold = 1
sigma = 0

for eta in eta_range:
    with open('_eta'+str(eta), 'w') as f:
        f.write('size=' + str(size) + '\n')
        f.write('N=' + str(N) + '\n')
        f.write('eta=' + str(eta) + '\n')
        f.write('threshold=' + str(threshold) + '\n')
        f.write('sigma=' + str(sigma) + '\n')
