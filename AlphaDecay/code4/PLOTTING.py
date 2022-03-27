import matplotlib.pyplot as plt 
from matplotlib import rcParams
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['font.size'] =22 
rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

rcParams['figure.figsize'] = [j*1.5 for j in rcParams['figure.figsize']]
rcParams['axes.labelsize'] = 22 
rcParams['figure.subplot.bottom'] = 0.11
rcParams['figure.subplot.left'] = 0.08
rcParams['figure.subplot.right'] = 0.99
rcParams['figure.subplot.top'] = 0.95
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
