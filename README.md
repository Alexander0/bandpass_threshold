# Bandpass Threshold

This repository contains the source code for the paper:

**"Bandpass Threshold Models for Activation and Influence Propagation in Social Networks"**

## Contents

- `main.cpp` — Simulation code for the bandpass threshold model  
- `data.zip` — Example social network datasets

## Requirements

This project requires the [igraph](https://igraph.org/) C library. Please ensure it is installed before compiling the code.

## Compilation

```bash
g++ main.cpp -std=c++11 -ligraph -o main.out

## Execution

```bash
./main.out N network.txt seeds.txt

- `N` — Number of nodes in the network
- `network.txt` — File containing the network edges
- `seeds.txt` — File listing the initial seed nodes for activation
