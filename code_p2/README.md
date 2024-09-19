# Project 2 - FYS4150 <!-- omit in toc -->

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)

## Introduction
This is project 2 for FYS4150, we attempt to find eigenvalues and eigenvectors using the Jacobi Rotation Algorithm.

## Features
- Scaling of equations
- Eigenvalue problems
- Unit testing

## Installation
Compiled with `c++14`.

On Unix, you can compile with `clang` or `g++` by using the following command:

```bash
g++ code_p2/task<N>.cpp -Icode_p1/include -std=c++11 -o out && ./out
```
On windows, you can compile with 
```bash
g++ task2.cpp -o task2.exe -L/ -larmadillo
```
