import numpy as np
import random
import os

def point_generate(low, high, num, sigma, width):
    L = []
    for k in range(num):
        n = random.randint(low, high)
        S = np.zeros((n,2))
        for i in range(n):
            p = random.random()
            if p <= 1.0:
                point = np.floor(np.random.uniform(0, width, size=(1, 2)))
            else:
                point = np.floor(np.random.normal(loc=width // 2, scale=sigma, size=(1,2)))
            S[i] = point
        L.append(S)
    return L
