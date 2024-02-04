#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018-2022 T.J.-Y. Derrien
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

## @package MachineLearning
# This package tries to make use of artificial intelligence to predict the properties of bicolor interaction for many parameters as a function of the intput parameters. We will train an AI with the several thousands theoretical data we have to make the magic happen. 
# https://scikit-learn.org/stable/tutorial/statistical_inference/settings.html

#from sklearn import datasets
#iris = datasets.load_iris()
#data = iris.data
#print(data.shape)

import numpy as np
from sklearn import linear_model

dataset_full = np.loadtxt("Database-LDA.csv")
# TODO: to improve the prediction, select the datase only with one type of similar data

dataset_in = dataset_full[:,[2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,20,21,22]] 
dataset_out = dataset_full[:,24] #15: NeMax, #16: NeFinal, #19: EnergyFinal, #20: Wpi

print("=== Shapes ===")
print((np.shape(dataset_in)))
print((np.shape(dataset_out)))

print("==== Learning ===")
alpha_val = 1E0
#reg = linear_model.Lasso(alpha=alpha_val, max_iter=1000000, random_state=True, positive=True, normalize=True, copy_X=False, fit_intercept = True, precompute=True)
reg = linear_model.Ridge(alpha=alpha_val, copy_X=False, fit_intercept=False, max_iter=None, normalize=False, random_state=True)
print((reg.fit(dataset_in,dataset_out)))

print("=== Prediction ===")
kpointgrid=16
realmeshgrid=0.43
TimeStep=0.25
lambda1=800e-9
tau1=50e-15
Efield1=5E10
polarization=0
CEP1=0.5
lambda2=400e-9
tau2=30e-15
Efield2=0
CEP2=0E0
delay=0
TheoryLevel=0
FieldInduced=0
Functional=0
FreezeHxc=0
MoveIon=0
#print(np.shape(dataset_in))
print((reg.predict([[kpointgrid, realmeshgrid, TimeStep, lambda1, tau1, Efield1, polarization, CEP1, lambda2, tau2, Efield2, CEP2, delay, TheoryLevel, FieldInduced, Functional, FreezeHxc, MoveIon]])))
#exit()
#print(reg.coef_)
##print(reg.intercept_)

