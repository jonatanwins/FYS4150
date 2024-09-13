# Project 1 - FYS4150 <!-- omit in toc -->

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)

## Introduction
This is project 1 for FYS4150, it is essentially a problem set exploring numerical solution techniques for the Poisson equation.

## Features
- Numerical methods for solving the Possion equation

## Installation
Compiled with `c++11`.

On Unix, you can compile with `clang` or `g++` by using the following command:

```bash
g++ code_p1/task<N>.cpp code_p1/src/utils.cpp -I code_p1/include -std=c++11 -o out  && ./out
