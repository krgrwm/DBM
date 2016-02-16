import sys

eta_start = 0
eta_end   = 10
eta_list = [eta/10.0 for eta in range(eta_start, eta_end+1)]

size = 500
N = 4000
threshold = 200

sigma_start = 0
sigma_end   = 5
sigma_list = [sigma/10.0 for sigma in range(sigma_start, sigma_end+1)]

for eta in eta_list:
    for sigma in sigma_list:
        with open('eta_'+str(eta)+'sigma_'+str(sigma), 'w') as f:
            f.write('size=' + str(size) + '\n')
            f.write('N=' + str(N) + '\n')
            f.write('eta=' + str(eta) + '\n')
            f.write('threshold=' + str(threshold) + '\n')
            f.write('sigma=' + str(sigma) + '\n')
