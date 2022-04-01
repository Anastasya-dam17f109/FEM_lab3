import numpy as np
import matplotlib.pyplot as plt
import math
import itertools
import os
import os.path

# радиус кольца
R = 1
# набор узлов и получающихся КЭ
nodes = {}
mesh_elems = {}
ELEM_AMOUNT = int(input("Введите количество конечных элементов: "))
NODE_AMOUNT = 3*ELEM_AMOUNT-2

def triangulation_1D():
    point = np.array([0.0, 0.0], float)
    start_idx = 0
    # получение словаря узлов
    for i in range(NODE_AMOUNT):
        nodes.update({i: point.copy()})
        point[0] = math.cos(2*math.pi*(i+1)/NODE_AMOUNT)
        point[1] = math.sin(2 * math.pi * (i + 1) / NODE_AMOUNT)
        for j in range(2):
            if abs(point[j]) <10**(-8):
                point[j] = 0
    # составление словая конечных элементов (группировка по 3 из существующего словаря узлов)
    mesh_elems.update({0: np.array([nodes.get(start_idx), nodes.get(start_idx + 1), nodes.get(start_idx + 2)])})
    for i in range(1,ELEM_AMOUNT-1):
        start_idx = 3*i-1
        mesh_elems.update({i: np.array([nodes.get(start_idx), nodes.get(start_idx + 1), nodes.get(start_idx + 2)])})
    start_idx = 3 * (ELEM_AMOUNT-1) - 1
    mesh_elems.update({ELEM_AMOUNT-1: np.array([nodes.get(start_idx), nodes.get(start_idx + 1), nodes.get(0)])})

triangulation_1D()
print(nodes)
print(mesh_elems)