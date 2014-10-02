#!/bin/sh
./testspme --model ../data/mono.mod --xyz ../data/N1000_Phi0.1.xyz \
--dim 256 --porder 6 --xi 0.8 --rmax 8.371887

./testspme --model ../data/mono.mod --xyz ../data/N1000_Phi0.1.xyz \
--dim 256 --porder 8 --xi 0.8 --rmax 8.371887

./testspme --model ../data/poly.mod --xyz ../data/N412.xyz \
--dim 256 --porder 8 --xi 0.04 --rmax 200.37188
