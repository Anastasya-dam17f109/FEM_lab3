import numpy as np
import matplotlib.pyplot as plt
import math
import itertools
import os
import os.path

# радиус кольца
R = 1
# кривизна конечного элемента
k = 1/R
#Q
Q = 1000
ps = 0
M = 0
# набор узлов и получающихся КЭ
nodes = {}
mesh_elems = {}
#coords =
ELEM_AMOUNT = int(input("Введите количество конечных элементов: "))
NODE_AMOUNT = 3*ELEM_AMOUNT-1
# длина конечного элемента
L = math.pi*R/ELEM_AMOUNT
# характеристики материала
E_F_ = 0.1093*10**11
E_J_ = 0.5466*10**9
E_G_ = 0.2878*10**9
# Матрица D^-1 для элемента
D_= np.array([[1.0/(E_F_),0,0],[0,1.0/(E_J_),0],[0,0,1.0/(E_G_)]])
# матрица H^-1 для элемента
H_matr = np.ones(6*6, dtype=float).reshape(6, 6)
H_matr[:3,:3] = 2*D_
H_matr[3:,3:] = 2*D_
H_matr[:3,3:] = -D_
H_matr[3:,:3] = -D_
H_matr = 2*H_matr/L
# зададим функции фи
phi1 = lambda x: 2*x**2-3*x+1
phi2 = lambda x:-4*x**2 + 4*x
phi3 = lambda x: 2*x**2 - x
# зададим производные функции фи
phi1_1 = lambda x: (4*x-3)/L
phi2_1 = lambda x:(-8*x + 4)/L
phi3_1 = lambda x: (4*x - 1)/L
# зададим функции пси
psi1 = lambda x: 1 - x
psi2 = lambda x: x

# формирование матрицы G
B1 = lambda x: np.array([[phi1_1(x), k*phi1(x),0 ],[0,0,phi1_1(x)],[-k*phi1(x),phi1_1(x),phi1(x) ]])
B2 = lambda x: np.array([[phi2_1(x), k*phi2(x),0 ],[0,0,phi2_1(x)],[-k*phi2(x),phi2_1(x),phi2(x) ]])
B3 = lambda x: np.array([[phi3_1(x), k*phi3(x),0 ],[0,0,phi3_1(x)],[-k*phi3(x),phi3_1(x),phi3(x) ]])

F1 = np.ones(6*9, dtype=float).reshape(6, 9)
F2 = np.ones(6*9, dtype=float).reshape(6, 9)
G_matr = np.ones(6*9, dtype=float).reshape(6, 9)
ksi1 = 0.5 + 1.0/(2*math.sqrt(3))
ksi2 = 0.5 - 1.0/(2*math.sqrt(3))

F1[:3,:3] = psi1(ksi1)*np.dot(D_,B1(ksi1))
F1[3:,:3] = psi2(ksi1)*np.dot(D_,B1(ksi1))
F1[:3,3:6] = psi1(ksi1)*np.dot(D_,B2(ksi1))
F1[3:,3:6] = psi2(ksi1)*np.dot(D_,B2(ksi1))
F1[:3,6:] = psi1(ksi1)*np.dot(D_,B3(ksi1))
F1[3:,6:] = psi2(ksi1)*np.dot(D_,B3(ksi1))

F2[:3,:3] = psi1(ksi2)*np.dot(D_,B1(ksi2))
F2[3:,:3] = psi2(ksi2)*np.dot(D_,B1(ksi2))
F2[:3,3:6] = psi1(ksi2)*np.dot(D_,B2(ksi2))
F2[3:,3:6] = psi2(ksi2)*np.dot(D_,B2(ksi2))
F2[:3,6:] = psi1(ksi2)*np.dot(D_,B3(ksi2))
F2[3:,6:] = psi2(ksi2)*np.dot(D_,B3(ksi2))

G_matr = (F1 + F2)*L/2.0
G_matr_t = G_matr.T
# матрица K
K = np.dot(np.dot(G_matr_t, H_matr), G_matr)

def triangulation_1D():
    point = np.array([-math.pi/2, 0.0, 0.0], float)
    start_idx = 0
    last = -math.pi/2
    step = math.pi*(1)/NODE_AMOUNT
    print(step)
    # получение словаря узлов
    for i in range(NODE_AMOUNT):
        nodes.update({i: point.copy()})
        last += step
        point[0] = last
        point[1] = math.cos(last)
        point[2] = math.sin(last)
        for j in range(2):
            if abs(point[j+1]) <10**(-8):
                point[j+1] = 0
    # составление словая конечных элементов (группировка по 3 из существующего словаря узлов)
   # print(nodes)
    #mesh_elems.update({0: np.array([nodes.get(start_idx), nodes.get(start_idx + 1), nodes.get(start_idx + 2)])})
    mesh_elems.update({0: np.array([start_idx, start_idx + 1, start_idx + 2])})
    for i in range(1,ELEM_AMOUNT-1):
        start_idx = 3*i-1
        mesh_elems.update({i: np.array([start_idx, start_idx + 1, start_idx + 2])})
    start_idx = 3 * (ELEM_AMOUNT-1) - 1
    mesh_elems.update({ELEM_AMOUNT-1: np.array([start_idx, start_idx + 1, start_idx + 2])})

# проверка, что к элементу приложена погонная сила

def isForseInElem(idx):
    node_list = mesh_elems.get(idx)
    count = 0
    #print(node_list)

    if nodes[node_list[0]][0] < 0 and nodes[node_list[2]][0]>0:
        count += 1
    if count>0:
        return True
    else:
        return False

def create_global_system():
    global gl_matr_std
    for i in range(len(mesh_elems)):
        node_list = []
        node_list = mesh_elems.get(i).copy()
        pz = 0
        if isForseInElem(i):
            pz = -Q/L

        P = np.array([ps,pz,M,4*ps, 4*pz, 4*M, ps, pz, M])*L/6
        for j in range(3):
            f_gl_std[node_list[j]*3] += P[3*j]
            f_gl_std[3*node_list[j]+1] += P[3*j+1]
            f_gl_std[3 * node_list[j] + 1] += P[3*j+2]


    for i in range( len(nodes)):
        #if abs(f_gl_std[i])<10**(-7):
           # f_gl_std[i] = 0
        if (nodes.get(i)[0] ==-math.pi/2) or (i == len(nodes)-1):
            print(i)
            f_gl_std[3*i] = 0
            f_gl_std[3*i+1] = 0
            f_gl_std[3*i+2] = 0
            gl_matr_std[3 * i+0, 3 * i+0] = 1
            #if nodes.get(i)[0] == -math.pi/2:
            gl_matr_std[3 * i + 1, 3 * i + 1] = 1
            gl_matr_std[3 * i + 2, 3 * i + 2] = 1

    for i in range(len(mesh_elems)):

        node_list = []
        for j in range(3):
            node_list.append( mesh_elems.get(i)[j])
            node_list.append(mesh_elems.get(i)[j])
            node_list.append(mesh_elems.get(i)[j])
        for j in range(9):
            i_gl = node_list[j]
            i_gl = 3*i_gl + j%3
            if (node_list[j] != 0) and (node_list[j]  != len(nodes)-1):
                for k in range(9):
                    j_gl = node_list[k]
                    j_gl = 3*j_gl + k%3
                    if (node_list[k] != 0) and (node_list[k]  != len(nodes)-1):
                        gl_matr_std[i_gl, j_gl] += K[j, k]


triangulation_1D()
f_gl_std = np.zeros(3*len(nodes))
gl_matr_std    = np.zeros(3*len(nodes)*3*len(nodes)).reshape(3*len(nodes),3*len(nodes))

create_global_system()
#print(gl_matr_std[3:,3:])
q_vec =  np.linalg.solve(gl_matr_std , f_gl_std)
u_vec = np.zeros(len(nodes))
w_vec = np.zeros(len(nodes))
teta_vec = np.zeros(len(nodes))
for i in range(len(nodes)):
    u_vec[i] = q_vec[3*i]
    w_vec[i] = q_vec[3 * i + 1]
    teta_vec[i] = q_vec[3 *i  +2]

gridsize = (3, 1)
fig2 = plt.figure(figsize=(18, 16))
ax1_1 = plt.subplot2grid(gridsize, (0, 0))
ax1_2 = plt.subplot2grid(gridsize, (1, 0))
ax1_3 = plt.subplot2grid(gridsize, (2, 0))
x = np.linspace(-math.pi/2, math.pi/2, NODE_AMOUNT)

ax1_1.plot(x, u_vec, label='график u по узлам')
ax1_2.plot(x, w_vec, label='график w по узлам')
ax1_3.plot(x, teta_vec, label='график teta по узлам')
plt.show()
#print(nodes)
#print(mesh_elems)
'''
else:
                        if k%3 ==2 and nodes.get(node_list[k])[0] == math.pi/2:
                            gl_matr_std[i_gl, j_gl] += K[j, k]
                 
            else:
                if j % 3 == 2 and nodes.get(node_list[j])[0] == math.pi / 2:
                    j_gl = node_list[k]
                    j_gl = 3 * j_gl + k % 3
                    if nodes.get(node_list[k])[0] != -math.pi / 2 or nodes.get(node_list[k])[0] != math.pi / 2:
                        gl_matr_std[i_gl, j_gl] += K[j, k]
                    else:
                        if k % 3 == 2 and nodes.get(node_list[k])[0] == math.pi / 2:
                            gl_matr_std[i_gl, j_gl] += K[j, k]
            
'''