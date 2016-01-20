#!/usr/local/bin/python3

import sys

if len(sys.argv) < 2:
    print("data file")
    exit(1)

file = sys.argv[1] #"../data/result/test.dat.boundary"
with open(file) as f:
    header = f.readline()
    lines = f.read().splitlines()
data = list(map(lambda l: map(int, l.split()), lines))

occupied = []
c = int(len(data)/2)
R = c-2

def r2(i, j, c):
    return (i-c)**2 + (j-c)**2

for i,row in enumerate(data):
    for j,v in enumerate(row):
        if v==1 and R**2 > r2(i, j, c):
            newj = j + ((i%2)-0.5)/2.0
            occupied.append((i, newj))

print(header, end=" ")
for pv in occupied:
    print("{0} {1}".format(pv[1], pv[0]))
