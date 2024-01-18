import pickle

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename
import matplotlib.pyplot as plt

# Input selection of data file
Tk().withdraw()  # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename()  # show an "Open" dialog box and return the path to the selected file


# Load the DFN data
pkl_file = open(filename, 'rb')
data = pickle.load(pkl_file)
DFN = data['DFN']
st = data['start_time']
et = data['end_time']
ext = data['exe_time']
fig = data['fig']

print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print('DFN loaded')
print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

# print the DFN information
print('DFN information:')
print('  Name: ' + str(DFN.name))
print('  Date and time: ' + str(st))
print('  Execution time: ' + str(ext))
print('  Number of fracture: ' + str(DFN.num_frac()))
print('  Number of intersections:' + str(DFN.num_int()))
print('  Path: ' + filename)

# Show the figure
plt.show()
