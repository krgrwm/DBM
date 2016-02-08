import sys

eta_start = 0
eta_end   = 10
eta = [eta/10.0 for eta in range(eta_start, eta_end+1)]

size = 400
N = 300
threshold = 200
sigma = 0

for _eta in eta:
    with open('eta_'+str(_eta)+'sigma_'+str(sigma), 'w') as f:
        f.write('size=' + str(size) + '\n')
        f.write('N=' + str(N) + '\n')
        f.write('eta=' + str(_eta) + '\n')
        f.write('threshold=' + str(threshold) + '\n')
        f.write('sigma=' + str(sigma) + '\n')
