#!/usr/local/bin/python3

import sys

if len(sys.argv) < 2:
    print("data file")
    exit(1)

file = sys.argv[1] #"../data/result/test.dat.boundary"
with open(file) as f:
    header = f.readline()
    lines = f.read().splitlines()
data = map(lambda l: map(int, l.split()), lines)

occupied = []
for i,row in enumerate(data):
    for j,v in enumerate(row):
        if v==1:
            newj = j + ((i%2)-0.5)/2.0
            occupied.append((i, newj))

print(header, end=" ")
for pv in occupied:
    print("{0} {1}".format(pv[1], pv[0]))
