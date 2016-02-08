import sys

if len(sys.argv) < 4:
    print("from, step, to")
    exit()

_from = int(sys.argv[1])
_step = int(sys.argv[2])
_to   = int(sys.argv[3])
if _from >= _to:
    print("error\nfrom >= to")
    exit()

size = range(_from, _to+_step, _step)
N = 3500
eta = 1.0
threshold = 1
sigma = 0

for s in size:
    print(s)
    with open('omega_'+str(s), 'w') as f:
        f.write('size=' + str(s) + '\n')
        f.write('N=' + str(N) + '\n')
        f.write('eta=' + str(eta) + '\n')
        f.write('threshold=' + str(threshold) + '\n')
        f.write('sigma=' + str(sigma) + '\n')
