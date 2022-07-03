# GLS254 repository

The code repository for the paper "2D-GLS: Faster and exception-free scalar multiplication in the GLS254 binary curve".

## Project structure

* arm - Armv8 AArch64 implementation of GLS254, with runnable tests and benchmarks.
* sage - Sage implementation of the group law formualas and scalar multiplication algorithms. Reports operation counts.
* x64 - 64-bit Intel implementation of GLS254

# How to run Arm code
The project is built as follows:
```
make <command> CC = <compiler>
```
command options:
* runbench (Arm PMU cycle counters must be enabled)
* runtests

CC options:
* clang (default)
* gcc

# How to run x64 code
The project is built as follows:
```
make CC = <compiler>
```

CC options:
* clang (default)
* gcc
