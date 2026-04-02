# *LS-DivSAT*: An Efficient Local Search Solver for Diverse SAT Problem

## How to Obtain *LS-DivSAT*

*LS-DivSAT* is [publicly available on GitHub](https://github.com/LS-DivSAT/LS-DivSAT). To obtain *LS-DivSAT*, a user may use `git clone` to get a local copy of the GitHub repository:

```sh
git clone https://github.com/LS-DivSAT/LS-DivSAT.git
```

*LS-DivSAT* uses git submodules to manage most of the dependencies. Therefore, after cloning this repository from GitHub, user may use `git submodule` command to clone the dependencies.

```sh
git submodule update --init
```

## Instructions for Building *LS-DivSAT*

*LS-DivSAT* is built on top of the SAT solver [CaDiCaL](https://github.com/arminbiere/cadical.git) and [MiniSat](https://github.com/niklasso/minisat.git), as well as utilizes the libraries [nuwls](https://github.com/shaowei-cai-group/NuWLS), [dbg-macro](https://github.com/sharkdp/dbg-macro) and [clipp](https://github.com/muellan/clipp).

This project is configured by Makefile. To build this project, users may use following commands:

```sh
make
```

Note that *LS-DivSAT* should be built on a 64-bit GNU/Linux operating system.

## Instructions for Running *LS-DivSAT*

After building *LS-DivSAT*, users may run it with the following command:

```sh
./LS-DivSAT -i [CNF_PATH] -o [OUTPUT_PATH] <optional_parameters>
```

For the required parameters, we list them as follows.

| Parameter | Allowed Values | Description |
| --- | --- | --- |
| -i / -input_file_path | string | path of the input CNF file |
| -o / -output_file_path | string | path to where the solution set is saved |

For the optional parameters, we list them as follows.

| Parameter | Allowed Values | Default Value | Description | 
| --- | --- | --- | --- |
| `-seed` | integer | `1` | random seed |
| `-k` | integer | `100` | number of samples to generate |

For more infomation about optional parameters, run following command:

```sh
./LS-DivSAT -h
```

## Example Command for Running *LS-DivSAT*

```sh
./LS-DivSAT -i benchmarks/software/linux.cnf -o linux.samples -k 10
```

The command above runs *LS-DivSAT* on `benchmarks/software/linux.cnf`, generating 10 solutions. The resulting solution set will be saved in `linux.samples`.

## Implementation of *LS-DivSAT*

*LS-DivSAT* is built on top of the SAT solver [CaDiCaL](https://github.com/arminbiere/cadical.git) and [MiniSat](https://github.com/niklasso/minisat.git). The `cadical.diff` file and the `minisat_ext/minisat.patch` file contain our modifications to CaDiCaL and MiniSat. The main implementation of LS-DivSAT is in `src/DiversitySAT.cpp`.

## Benchmarks

The `benchmarks/` directory contains two subdirectories with testing benchmarks:

- Hardware Benchmarks:

    The `benchmarks/hardward/` directory contains 66 hardware benchmarks that can be used as input for testing.

- Software Benchmarks:
    
    The `benchmarks/software/` directory contains 125 software benchmarks.

## Experimental Results

The `experimental_results/` directory contains the experimental results comparing against the state-of-the-art diverse SAT solvers.

- `experimental_results/Results_of_1hour/`: contains results of LS-DivSAT and its state-of-the-art competitors on the benchmarks with 1 hour cutoff time.
- `experimental_results/Results_of_DiversekSet_time/`: contains results of LS-DivSAT at the same time when DiversekSet generates its solution sets.
- `experimental_results/Results_of_DiverSAT_time/`: contains results of LS-DivSAT at the same time when DiverSAT generates its solution sets.
- `experimental_results/Results_of_SATLIB_benchmarks/`: contains results of LS-DivSAT and its state-of-the-art competitors on the SATLIB benchmarks with 1 hour cutoff time.
