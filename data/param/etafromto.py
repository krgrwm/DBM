import sys

if len(sys.argv) < 3:
    print("from, to")
    exit()

_from = int(sys.argv[1])
_to   = int(sys.argv[2])
if _from >= _to:
    print("error\nfrom >= to")
    exit()

size = 800
N = 3500
eta = range(_from, _to+1)

for i in range(_from, _to+1):
    with open('_eta'+str(i), 'w') as f:
        f.write('size=' + str(size) + '\n')
        f.write('N=' + str(N) + '\n')
        f.write('eta=' + str(eta[i]) + '\n')
