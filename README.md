# 3D DFN generator
A generator for a 3D Discrete Fracture Network (DFN)

The repository contains two different programs:
1. `generate_DFN.py` a DFN generator
2. `load_DFN.py` a DFN viewer

The file `DFN_classes.py` contains the python classes.

The file `demo_frac3.pkl` is a previously generated DFN with three fractures.

## DFN generator `generate_DFN.py`
A program that generates a DFN based on:
- `bbox`: the bounding box for the DFN
- `num_frac`: numer of fractures
- `r_min`: minimum radius of a fracture
- `r_max`: maxium radius of a fracture

The DFN is save to a *.pkl'-file.

## DFN viewer `load_DFN.py`
A program the plots a privously save DFN. The user is asked to select a file when the prgoram is started.
