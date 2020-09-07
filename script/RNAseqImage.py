# Author: yoshi
# Date: 8/19/2020
# Updated:
# Project: covid-19
# Script: Generate mimic array images

import numpy as np
import seaborn as sns

def main():

    np.random.seed(0)
    matrix = np.random.rand(125, 125)
    sns.heatmap(matrix, yticklabels=False, xticklabels=False, cmap='Spectral')


if __name__ == '__main__':
    main()
