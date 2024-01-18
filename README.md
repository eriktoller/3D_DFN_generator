# 3D DFN generator
A generator for a 3D Discrete Fracture Network (DFN)

The repository contains two different programs:
1. `generate_DFN.py` a DFN generator
2. `load_DFN.py` a DFN viewer

The file `DFN_classes.py` contains the Python classes.

The file `demo_frac3.pkl` is a previously generated DFN with three fractures.

## DFN generator `generate_DFN.py`
A program that generates a DFN based on:
- `bbox`: the bounding box for the DFN
- `num_frac`: number of fractures
- `r_min`: minimum radius of a fracture
- `r_max`: maximum radius of a fracture

The DFN is saved to a `*.pkl`-file.

## DFN viewer `load_DFN.py`
A program that plots a previously save DFN. The user is asked to select a file when the program is started.
