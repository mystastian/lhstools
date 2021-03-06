import matplotlib.pyplot as plt
import numpy as np

def discrete_cmap(base_cmap=None, N=10):
    """
    Author: S.Lienert, lienert@climate.unibe.ch

    Create an N-bin discrete colormap from the specified input map
    """
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

