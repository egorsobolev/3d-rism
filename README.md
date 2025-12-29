# 3D RISM

**3d-rism** is a software package for solving integral equations of the theory of liquids within the **three-dimensional Reference Interaction Site Model (3D-RISM)** approximation [[1](https://doi.org/10.1021/jp971083h)].

[1] D. Beglov and B. Roux, *“An Integral Equation To Describe the Solvation of Polar Molecules in Liquid Water,”* **J. Phys. Chem. B**, 1997, 101, 7821–7826.
DOI: [10.1021/jp971083h](https://pubs.acs.org/doi/10.1021/jp971083h)

## Installation

### Dependencies

* [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS)
* [pocketfft](https://gitlab.mpcdf.mpg.de/mtr/pocketfft)
* [cpocketfft](https://github.com/egorsobolev/cpocketfft)
* [argtable2](https://argtable.sourceforge.io/)

### Download & Build

```sh
git clone --recurse-submodules https://github.com/egorsobolev/3d-rism.git
cd 3d-rism
make
make test
```

## Examples

See the test suite for usage examples.
