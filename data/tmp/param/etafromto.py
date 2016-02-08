import sys

if len(sys.argv) < 3:
    print("from, to")
    exit()

_from = int(sys.argv[1])
_to   = int(sys.argv[2])
if _from >= _to:
    print("error\nfrom >= to")
    exit()

size = 2000
N = 3500
eta = range(_from, _to+1)
threshold = 1
sigma = 0

for i in range(_from, _to+1):
    with open('_eta'+str(i), 'w') as f:
        f.write('size=' + str(size) + '\n')
        f.write('N=' + str(N) + '\n')
        f.write('eta=' + str(eta[i]) + '\n')
        f.write('threshold=' + str(threshold) + '\n')
        f.write('sigma=' + str(sigma) + '\n')