#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:14:54 2019

@author: or
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

import mne
from mne import io
from mne.datasets import sample
from mne.stats import bonferroni_correction, fdr_correction
# first crerate z transformed matrices

ket1_corr_z = np.arctan(ket1_corr)
ket2_corr_z = np.arctan(ket2_corr)
ket2vs1_real = []
for ket2,ket1 in zip(ket1_corr_z, ket2_corr_z):
    a = ket2 -ket1
    ket2vs1_real.append(a)
ket2vs1_real = np.array(ket2vs1_real)
# then average them and contrasts time 2 to time 1
ket2vs1 = np.average(ket2_corr_z, axis =0) - np.average(ket1_corr_z, axis=0)
ket2vs1.shape
#%%
# here we try to remove self correlations (diagonal) from matrix

df = pd.DataFrame(ket2vs1)
# masking on diagonal (to avoid repeated numbers)
dataCorr = df.mask(np.tril(np.ones(df.shape)).astype(np.bool))
# remove less than 0.001
dataCorr = dataCorr[abs(dataCorr) >= 0.001].stack().reset_index()
#%%
# run simple T test across matrix
corArr = np.array(dataCorr)

T, pval = stats.ttest_1samp(ket2vs1_real, 0)
alpha = 0.1

n_samples, n_tests = T.shape

# regular threshold
threshold_uncorrected = stats.t.ppf(1.0 - alpha, n_samples - 1)


# run bonferroni correction
reject_bonferroni, pval_bonferroni = bonferroni_correction(pval, alpha=alpha)
threshold_bonferroni = stats.t.ppf(1.0 - alpha / n_tests, n_samples - 1)

# run FDR correction
reject_fdr, pval_fdr = fdr_correction(pval, alpha=alpha, method='indep')
threshold_fdr = np.min(np.abs(T)[reject_fdr])

sum(reject_bonferroni)

# plot
plt.close('all')
plt.scatter(range(114), T)#, 'k', label='T-stat')
xmin, xmax = plt.xlim()
plt.hlines(threshold_uncorrected, xmin, xmax, linestyle='--', colors='k',
           label='p=0.05 (uncorrected)', linewidth=2)
plt.hlines(threshold_bonferroni, xmin, xmax, linestyle='--', colors='r',
           label='p=0.05 (Bonferroni)', linewidth=2)
plt.hlines(threshold_fdr, xmin, xmax, linestyle='--', colors='b',
           label='p=0.05 (FDR)', linewidth=2)
plt.legend()
plt.xlabel("region")
plt.ylabel("T-stat")
plt.show()

sum(reject_fdr)
len(pval)

#%%
import statsmodels

statsmodels.stats.weightstats.ttost_ind(np.average(np.array(ket1_corr,axis=0)), np.average(np.array(ket1_corr,axis=0)), low, upp)
#%%
# Ruonan suggestion:
# create contrast correlation for each subject
# run T test 
# threshold

# one option is to:
# 1. Vecotorize the matrix
# 2. contrast per subject
# 3. run T test
# 4. FDR thresholding
from nilearn.connectome import ConnectivityMeasure
correlation_measure_vectorize = ConnectivityMeasure(kind='correlation', vectorize=True,  discard_diagonal=True)

ket1_vector = []

for time_s in ket1_series:
    correlation_matrix = correlation_measure_vectorize.fit_transform([time_s])[0]
    # fisher's Z transformation
    correlation_matrix = np.arctan(correlation_matrix)
    ket1_vector.append(correlation_matrix)
ket1_vector = np.array(ket1_vector)

ket2_vector = []

for time_s in ket2_series:
    correlation_matrix = correlation_measure_vectorize.fit_transform([time_s])[0]
    # fisher's Z transformation
    correlation_matrix = np.arctan(correlation_matrix)
    ket2_vector.append(correlation_matrix)
ket2_vector = np.array(ket2_vector)

