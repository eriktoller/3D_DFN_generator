from dataclasses import dataclass
import numpy as np
from typing import List

@dataclass
class Circle3D:
    """
    Data class for a circular fracture in 3D space.
    """
    zc3d: np.ndarray  # Center point of circle
    r: float  # Radius of circle
    n1: np.ndarray  # unit vector local x1-direction
    n2: np.ndarray  # unit vector local x2-direction
    n3: np.ndarray  # normal vector
    int_list: list  # list of intersections
    z_int: list  # list with intersection points
    ls: str = '-'  # line style
    marker: str = ''  # marker style
    color: str = 'black'  # line color
    ls_width: str = '1.5'  # line width
    alpha: float = .2  # transparency of fill
    face_color: str = 'black'  # fill color


@dataclass
class Box:
    """
    Data class for bounding box.
    """
    x1: list  # list with start and end value for x1-direction
    x2: list  # list with start and end value for x2-direction
    x3: list  # list with start and end value for x3-direction


@dataclass
class DFN:
    """
    Data class for the Discrete Fracture Network (DFN).
    """
    fractures: List[Circle3D]  # list of fractures
    name: str = 'DFN'  # name of DFN

    def num_frac(self):
        """
        :return: the numer of fractures in the DFN
        """
        return len(self.fractures)

    def num_int(self):
        """
        :return: the numer of intersections in the DFN
        """
        cint = 0
        for frac in self.fractures:
            cint = cint + len(frac.int_list)
        cint = int(cint / 2)
        return cint
