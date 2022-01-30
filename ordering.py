"""
import numpy as np

# import scipy.cluster.hierarchy.dendrogram

from scipy.cluster.hierarchy import dendrogram, linkage

# from matplotlib import pyplot as plt

def optimalLeafOrder(allvsallmatrix):

    # TypeError: linkage() got an unexpected keyword argument 'optimal_ordering'
    Z = linkage(allvsallmatrix, "average")

#    fig = plt.figure(figsize=(5, 5))

    dn = dendrogram(Z,no_plot=True)

    newallvsallmatrix = np.zeros(allvsallmatrix.shape)

    newOrder = []

    for i in range(newallvsallmatrix.shape[0]):

        newOrder.append(dn['leaves'][i])

        for j in range(newallvsallmatrix.shape[1]):

            newallvsallmatrix[i][j] = allvsallmatrix[dn['leaves'][i]][dn['leaves'][j]]

    return newOrder
"""

