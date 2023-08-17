#! /usr/bin/python3
import sys
import matplotlib.pyplot as plt

if len(sys.argv) == 1 or len(sys.argv) > 2:
    print("usage:", sys.argv[0], "<filename>")
    exit()

# open file and remove newlines. readlines is okay b/c
# the input file shouldn't be thaaaat big
f=open(sys.argv[1])
lines = [line.strip() for line in f.readlines()]
f.close()

# first column defines the inputs
impl_domain = [float(x.split(" ")[0]) for x in lines[1:]]
impl_range  = [float(x.split(" ")[1]) for x in lines[1:]]
plt.plot(impl_domain, impl_range, '-')
impl_range  = [float(x.split(" ")[2]) for x in lines[1:]]
plt.plot(impl_domain, impl_range, '-')
plt.show()
