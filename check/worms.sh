#!/bin/sh

PROG="$1"
HEADER="# Magnetization density of HAF chain with L = 4 and T = 0.1"

H="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0"

(for h in $H; do $PROG -L 4 -T 0.1 -H $h; done) > worms.dat

echo "$HEADER" > worms-ene.dat
cat worms.dat | awk '$1=="Magnetic" { h=$4 } $1=="Energy" && $2=="Density" {print h, $4, $6}' >> worms-ene.dat

echo "$HEADER" > worms-mag.dat
cat worms.dat | awk '$1=="Magnetic" { h=$4 } $1=="Uniform" && $2=="Magnetization" {print h, $4, $6}' >> worms-mag.dat

echo "$HEADER" > worms-mag2.dat
cat worms.dat | awk '$1=="Magnetic" { h=$4 } $1=="Uniform" && $2=="Magnetization^2" {print h, $4, $6}' >> worms-mag2.dat

echo "$HEADER" > worms-smag2.dat
cat worms.dat | awk '$1=="Magnetic" { h=$4 } $1=="Staggered" && $2=="Magnetization^2" {print h, $4, $6}' >> worms-smag2.dat
