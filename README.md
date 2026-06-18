LazyPhase
====================

<img src="logo.jpg" height=80>

### Introduction

LazyPhase is a Moltemplate- and LAMMPS-based tool for high-throughput
simulations of phase behavior of polypeptide-like polymers. A user can specify
the alphabet of side chains, their sizes and interactions.

### Dependencies

Python
Moltemplate
LAMMPS
Optional: Kokkos or other GPU interface

### Guide

For help:
```
lazyphase --help
```

Generating the sequences and preparing the high-throughput simulations:
```
lazyphase generate --alphabet SL --length 64 --num 100 --output seqs.txt
```

Run the high-throughput simulations:
```
lazyphase run --seqs seqs.txt --time 1
```
--seqs specifies the file containing the sequences, each in separate line, and
--time is the modeling time in microseconds.

### Copyright and citation

(C) Egor Vasilenko, Ivankov Lab.
If use, please cite:
