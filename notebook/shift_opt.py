# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: qcomputing
#     language: python
#     name: qcomputing
# ---

# %matplotlib inline

# +
import math

import dimod
import numpy as np
import openjij as oj
from pyqubo import Array, Constraint, LogEncInteger, solve_qubo
# -

# # 制約条件
# $$ W_j <= \sum_i x_{ij} y_i $$

# +
W = np.array([2.1, 2.5, 1, 3.25], dtype=float)
y = np.array([0.25, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# タスクの数
M = len(W)
# 人の数
N = len(y)
# -

# # 目的変数
# $$ x \in \{0,1\}^{N \times M} $$

x = Array.create("x", shape=(N, M), vartype="BINARY")
c_ineqs = np.array(
    [LogEncInteger("y" + str(i), (0, int(np.ceil(4 * W[i])))) for i in range(M)]
)

A = np.array([1, 1, 10, 10])
B = 50
# Hy = Constraint(sum(A * (W - y @ x) ** 2), label="Hy")
Hy = Constraint(
    sum(A * (np.ceil(4 * W) + c_ineqs - np.ceil(4 * y) @ x) ** 2), label="Hy"
)
Hx = Constraint(sum(B * (1 - np.sum(x, axis=1)) ** 2), label="Hx")

H = Hx + Hy
model = H.compile()
q, offset = model.to_qubo()

sampler = oj.SQASampler(num_reads=200)
sampleset = sampler.sample_qubo(q)
decoded_sample = model.decode_sample(sampleset.first.sample, vartype="BINARY")

res_x = np.zeros(x.shape)
for i in range(N):
    for j in range(M):
        res_x[i, j] = decoded_sample.array("x", (i, j))

res_ineq = np.zeros(M)
for j in range(M):
    res_ineq[j] = sum(
        2**k * v
        for k, v in [
            (elem, decoded_sample.array("y" + str(j), elem))
            for elem in range(math.ceil(math.log2(W[j])))
        ]
    )

res_ineq

res_x.sum(axis=1)

np.ceil(4 * W)

np.ceil(4 * y) @ res_x

np.ceil(4 * W) - np.ceil(4 * y) @ res_x

res_x

W = 20
c = {0: 5, 1: 7, 2: 2, 3: 1, 4: 4, 5: 3}
w = {0: 8, 1: 10, 2: 6, 3: 4, 4: 5, 5: 3}
N = len(w)


E = 3

x = Array.create("x", shape=(N), vartype="BINARY")
y = LogEncInteger("y", (0, W))
z = LogEncInteger("z", (0, E))

# +
key1 = max(c, key=lambda k: c[k])
B = 40
A = 10 * B * c[key1]

HA1 = Constraint(A * (W - sum(w[a] * x[a] for a in range(N)) - y) ** 2, label="HA1")
HA2 = Constraint(A * (E - sum(x) - z) ** 2, label="HA2")

HB = -B * sum(c[a] * x[a] for a in range(N))

# +
H = HA1 + HA2 + HB
Q = H
model = Q.compile()
q, offset = model.to_qubo()

sampler = oj.SQASampler(num_reads=100)
sampleset = sampler.sample_qubo(q)
decoded_sample = model.decode_sample(sampleset.first.sample, vartype="BINARY")
# -

sum(
    2**k * v
    for k, v in [
        (elem, decoded_sample.array("z", elem))
        for elem in range(math.ceil(math.log2(E)))
    ]
)

# +
print()
print("[Results]")
print()
print("decoded_sample.sample:")
print(decoded_sample.sample)
print()
print("x (選ばれた宝物) :")

treasures = ["A", "B", "C", "D", "E", "F", "G"]
weight = 0
cost = 0

for k in range(N):
    if decoded_sample.array("x", k) != 0:
        print("宝物" + treasures[k])
        weight += w[k]
        cost += c[k]


sol_y = sum(
    2**k * v
    for k, v in [
        (elem, decoded_sample.array("y", elem))
        for elem in range(math.ceil(math.log2(W)))
    ]
)

print()
print("スラック変数Y = {}".format(sol_y))
print()
print("broken")
print(decoded_sample.constraints(only_broken=True))
print("合計の重さ : " + str(weight / 10) + "kg")
print("合計の価格 : $" + str(cost) + ",000")

# +
import dimod

print("[Inputs]")
print()
print("W (ナップサックの容量) : " + str(W / 10) + "kg")
print("N (宝物の数): " + str(N))
print()
print("weight list")
print(w)
print()
print("cost list")
print(c)
print()
print("A : " + str(A))
print("B : " + str(B))

H = HA + HB
Q = H
model = Q.compile()
q, offset = model.to_qubo()

sampleset = dimod.ExactSolver().sample_qubo(q)
decoded_sample = model.decode_sample(sampleset.first.sample, vartype="BINARY")
print()
print("[Results]")
print()
print("decoded_sample.sample:")
print(decoded_sample.sample)
print()
print("x (選ばれた宝物) :")

treasures = ["A", "B", "C", "D", "E", "F", "G"]
weight = 0
cost = 0

for k in range(N):
    if decoded_sample.array("x", k) != 0:
        print("宝物" + treasures[k])
        weight += w[k]
        cost += c[k]


sol_y = sum(
    2**k * v
    for k, v in [
        (elem, decoded_sample.array("y", elem))
        for elem in range(math.ceil(math.log2(W)))
    ]
)

print()
print("スラック変数Y = {}".format(sol_y))
print()
print("broken")
print(decoded_sample.constraints(only_broken=True))
print("合計の重さ : " + str(weight / 10) + "kg")
print("合計の価格 : $" + str(cost) + ",000")
# -


