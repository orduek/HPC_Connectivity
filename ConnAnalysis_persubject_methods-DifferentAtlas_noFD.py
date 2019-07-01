#!/usr/bin/env python
# coding: utf-8

# Here we will build a connectivity analysis pipeline for the KPE study. 
# Methods should be easily generalized for others studies. 
# 

# In[ ]:


## Use this box if you want MDSL
# using MSDL atlas - Can choose different atlas or different nodes (using ICA or something else)
# from nilearn import datasets
# atlas = datasets.fetch_atlas_msdl()

# # Loading atlas image stored in 'maps'
# atlas_filename = atlas['maps']
# # Loading atlas data stored in 'labels'
# #labels_img=yeo['thick_17']
# labels =  atlas['labels']
# coords = atlas.region_coords# grab center coordinates for atlas labels
# #coords = plotting.find_parcellation_cut_coords(labels_img=labels)
# #atlas.region_coords

# # optional set of different atlas
# #atlas_yeo_2011 = 
# #atlas_yeo = atlas_yeo_2011.thick_7


# In[ ]: Using Yeo atlas

atlas_filename = '/home/or/Downloads/1000subjects_reference_Yeo/Yeo_JNeurophysiol11_SplitLabels/MNI152/Yeo2011_17Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz'
atlas_labes = pd.read_csv('/home/or/Downloads/1000subjects_reference_Yeo/Yeo_JNeurophysiol11_SplitLabels/Yeo2011_17networks_N1000.split_components.glossary.csv')
 
# In[ ]: Using Shen Atlas


# Run either this or MSDL - not both boxes
# try to use Shen atlas
import os
import numpy as np

import pandas as pd
from nilearn import plotting
atlas_filename = '/home/or/Downloads/shenPar/shen_1mm_268_parcellation.nii.gz'
atlas_labes = pd.read_csv('/home/or/Downloads/shenPar/shen_268_parcellation_networklabels.csv')
coords = plotting.find_parcellation_cut_coords(labels_img=atlas_filename)

atlas_labes = np.array(atlas_labes)
atlas_labes.shape


# In[ ]:


# methods
def removeVars (confoundFile):
    # this method takes the csv regressors file (from fmriPrep) and chooses a few to confound. You can change those few
    import pandas as pd
    confound = pd.read_csv(confoundFile,sep="\t", na_values="n/a")
    finalConf = confound[['csf', 'white_matter',
                          'a_comp_cor_00', 'a_comp_cor_01',	'a_comp_cor_02', 'a_comp_cor_03', 'a_comp_cor_04', 
                        'a_comp_cor_05', 'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']] # can add 'global_signal' also
     # change NaN of FD to zero
    finalConf = np.array(finalConf)
    
    return finalConf


# In[ ]:


# build method for creating time series for subjects
def timeSeries(func_files, confound_files):
    total_subjects = [] # creating an empty array that will hold all subjects matrix 
    # This function needs a masker object that will be defined outside the function
    for func_file, confound_file in zip(func_files, confound_files):
        print(f"proccessing file {func_file}") # print file name
        confoundClean = removeVars(confound_file)
        confoundArray = confoundClean#confoundClean.values
        time_series = masker.fit_transform(func_file, confounds=confoundArray)
        #time_series = extractor.fit_transform(func_file, confounds=confoundArray)
        #masker.fit_transform(func_file, confoundArray)
        total_subjects.append(time_series)
    return total_subjects

# contrasting two timePoints
def contFuncs(time_series1, time_series2):
    twoMinusOneMat = []
    for scanMatrix, scanMatrix2 in zip(time_series1, time_series2):
        a = scanMatrix2 - scanMatrix
        twoMinusOneMat.append(a)
    return np.array(twoMinusOneMat)

import numpy as np
from nilearn import plotting

# create correlation matrix per subject
def createCorMat(time_series):
    # create correlation matrix for each subject
    fullMatrix = []
    for time_s in time_series:
        correlation_matrix = correlation_measure.fit_transform([time_s])[0]
        fullMatrix.append(correlation_matrix)
    return fullMatrix

# create connecotme graph per subject
def connectome_graph (fullMatrix):
    # here it is set to threshold 1%
    for matrix in fullMatrix:
        plotting.plot_connectome(matrix, coords,
                             edge_threshold="99%", colorbar=True)
        plotting.show()


# In[ ]:


# Here you set the specific methods for masking and correlation. Please see Nilearn website for more info.

from nilearn.input_data import NiftiMapsMasker
from nilearn.input_data import NiftiLabelsMasker
# in this mask we standardize the values, so mean is 0 and between -1 to 1
# masker = NiftiMapsMasker(maps_img=atlas_filename, standardize=True, smoothing_fwhm = 6,
#                          memory="/home/oad4/scratch60/shenPar_nilearn",high_pass=.01 , low_pass = .1, t_r=1, verbose=5)

# use different masker when using Yeo atlas. 
masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True,smoothing_fwhm = 6,
                        memory="/media/Data/nilearn",high_pass=.01 , low_pass = .1, t_r=1, verbose=5)
                           
from nilearn.connectome import ConnectivityMeasure
correlation_measure = ConnectivityMeasure(kind='partial correlation') # can choose partial - it might be better
#correlation_measure = ConnectivityMeasure(kind='correlation', vectorize = True, discard_diagonal = True) # can choose partial - it might be better


# In[ ]:


# now we call subjcets
# and start the real analysis
subList =  ['008','1223','1293','1307','1322','1339','1343','1387']
midSubList = ['1253','1263','1351','1364','1369','1390','1403']
totalSub = subList + midSubList
# these two functions take subject list and session number (in string) and return func file list and confound file list
def fileList(subjects, session):
    func_files = ['/media/Data/KPE_fmriPrep_preproc/kpeOutput/fmriprep/sub-%s/ses-%s/func/sub-%s_ses-%s_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz' % (sub,session,sub,session) for sub in subjects]
    return func_files

def confList(subjects, session):
    confound_files = ['/media/Data/KPE_fmriPrep_preproc/kpeOutput/fmriprep/sub-%s/ses-%s/func/sub-%s_ses-%s_task-rest_desc-confounds_regressors.tsv' % (sub,session,sub,session) for sub in subjects]
    return confound_files


# In[ ]:


# now we call for the functions for each set.
# for every time line we want to run time series and then contrast between the times
ket1_series = timeSeries(func_files=fileList(subList,'1'), confound_files=confList(subList, '1'))

ket2_series = timeSeries(func_files=fileList(subList,'2'), confound_files=confList(subList, '2'))


# In[ ]:


ket3_series = timeSeries(func_files=fileList(subList,'3'), confound_files=confList(subList, '3'))


# In[ ]:


subListKet4 =  ['008','1293','1307','1322','1339','1343','1223']
ket4_series = timeSeries(func_files=fileList(subListKet4,'4'), confound_files=confList(subListKet4, '4'))


# In[ ]:


# build correlation matrix for each time point
ket1_corr = createCorMat(ket1_series)
ket2_corr = createCorMat(ket2_series)
ket3_corr = createCorMat(ket3_series)
ket4_corr = createCorMat(ket4_series)


# In[ ]:


# start contrasting
ket2_ket1 = contFuncs(ket1_corr, ket2_corr)
ket3_ket1 = contFuncs(ket1_corr, ket3_corr)

ket3_ket1.shape


# In[ ]:
%matplotlib qt

plotting.plot_connectome(np.average(ket2_ket1, axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Ketamine 7days minux first")

plotting.plot_connectome(np.average(ket3_ket1, axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Ketamine 30days minux first")


# In[ ]:


# let do midazolam
mid1_series = timeSeries(func_files=fileList(midSubList,'1'), confound_files=confList(midSubList, '1'))

mid2_series = timeSeries(func_files=fileList(midSubList,'2'), confound_files=confList(midSubList, '2'))


# In[ ]:


mid3_no1253 = midSubList
print(mid3_no1253)
mid3_no1253.remove('1253')
mid3_series = timeSeries(fileList(mid3_no1253,'3'), confList(mid3_no1253,'3'))


# In[ ]:


# build correlation matrix for each time point
mid1_corr = createCorMat(time_series=mid1_series)
mid2_corr = createCorMat(mid2_series)
mid3_corr = createCorMat(mid3_series)


# In[ ]:


# start contrasting
mid2_mid1 = contFuncs(mid1_corr, mid2_corr)
ket3_ket1 = contFuncs(ket1_corr, ket3_corr)

mid3_mid1 = contFuncs(mid1_corr[1:], mid3_corr) # removing the first (1253) from the mid1 correlation matrix
#np.average(mid3_corr, axis=0) - np.average(mid1_corr, axis = 0) # we do so because the number do not match

mid3_mid1.shape


# In[ ]:


# plotting connectome differences
get_ipython().run_line_magic('matplotlib', 'inline')

plotting.plot_connectome(np.average(ket2_ket1, axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Ketamine 7days minus first")

plotting.plot_connectome(np.average(ket3_ket1, axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Ketamine 30days minus first")

plotting.plot_connectome(np.average(mid2_mid1, axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Midazolam 7days minus first")

#plotting.plot_connectome(mid3_mid1, coords,
 #                        edge_threshold="99%", colorbar=True, title = "midazolam 30days minux first")

plotting.plot_connectome(np.average(ket1_corr,axis=0) -np.average(mid1_corr, axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Ketamine - Midazolam 0days")

plotting.plot_connectome(np.average(ket2_corr,axis=0)-np.average(mid2_corr,axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Ketamine-Midazolan 7days", output_file='pretty_brain.png')


plotting.show()


# In[ ]:


# rearrange the array to N,N,subnumber for NBS
ket1Reshape = np.moveaxis(np.array(ket1_corr), 0,-1)
print(ket1Reshape.shape)
ket2Reshape = np.moveaxis(np.array(ket2_corr), 0,-1)
ket3Reshape = np.moveaxis(np.array(ket3_corr), 0,-1)
ket4Reshape = np.moveaxis(np.array(ket4_corr), 0,-1)
mid1Reshape = np.moveaxis(np.array(mid1_corr),0,-1)
mid2Reshape = np.moveaxis(np.array(mid2_corr),0,-1)
mid3Reshape = np.moveaxis(np.array(mid3_corr),0,-1)

#print(mid3Reshape.shape)


# In[ ]:


# now we can run NBS
# NBS is taken from: https://github.com/aestrivex/bctpy, can be installed using pip (pip install bctpy)
from bct import nbs
# we compare ket1 and ket3
pval, adj, _ = nbs.nbs_bct(ket1Reshape, ket3Reshape, thresh=3, tail='both',k=1000, paired=True, verbose = True)
# check mean p vlue
#np.mean(checkNBS[0])


# In[ ]:
# look at p values and No. of components.
print(pval.shape)
print (pval)
len(pval)
print(adj.shape)

print(adj[0:10])
ad = np.array(adj)
print(ad[:,0:10])
#bct.adjacency_plot_und(adj, coords, tube=False)


# In[ ]:


# create data frame from adjacency matrix
# here labels has two columns. Need only one
import pandas as pd
adDat = pd.DataFrame(ad, columns=atlas_labes[:,0],index=atlas_labes[:,0])
print(adDat)


# In[ ]:

%matplotlib qt
# graph adjacency matrix
import networkx as nx
#G = nx.from_numpy_matrix(np.array(ad)) 
#G = nx.DiGraph(adDat, with_labels=True)
G = nx.from_pandas_adjacency(adDat)
# drawing the adajency matrix
import matplotlib.pyplot as plt
#plt.figure(figsize=(70,50))
nx.draw(G, with_labels=True) #, font_size = 50, width = 1, alpha = 0.7, font_weight = "bold")


# Between 30 days in 1st day Ketamine has different in one components. 
# Lets graph it

# In[ ]:


# lets look at 90 days and first
pvalk2, adjk2, _ = nbs.nbs_bct(ket2Reshape, ket3Reshape, thresh=3, tail='both',k=300, paired=True, verbose = True)


# In[ ]:


# Lets check Midazolam group


# In[ ]:


pvalMid, adjMid, _ = nbs.nbs_bct(mid1Reshape, mid2Reshape, thresh=3, tail='both',k=1000, paired=True, verbose = True)


# In[ ]:


print(pvalMid.shape)
print (pvalMid)


# Midazolam has no significant change between first day and 7 or 30 days.
# Lets check differences between ket and mid.

# In[ ]:


pvalKetMid, adjKetMid, _ = nbs.nbs_bct(ket2Reshape, mid2Reshape, thresh=3, tail='both',k=1000, paired=False, verbose = True)


# In[ ]:


print(pvalKetMid.shape)
print (pvalKetMid)
# no difference between ketamine and midazolam. 


# In[ ]:


# for the sake of QA we can create histogram plots of correlation matrices

import matplotlib.pyplot as plt
import seaborn as sns




# In[ ]:


# Run analysis on all ket subjects
#color = sns.cubehelix_palette(len(ket1_corr),8)
#color = sns.palplot(sns.color_palette("RdBu_r", len(ket1_corr)))
sns.set_palette("husl") # set color pallet
correlation_vec = ConnectivityMeasure(kind='partial correlation', vectorize=True) # can choose partial - it might be better
    
# create correlation matrix for each subject
fullVec = []
for time_s in ket3_series:
    cor = correlation_vec.fit_transform([time_s])[0]
    print(cor.shape)
    plt.hist(cor, alpha = 0.3)
    fullVec.append(cor)


# In[ ]:



# loop through adj matrix
# every time find 1 index it and take actual corelation deltas from Ket3-ket1 matrix (or other of that kind)
#ketDeltaA = []
#for i in ket3_ket1:
 #   ketSubjectAdj = []#[np.zeros((39,39))]
  #  for n,k in zip(ad,i):
        
        # now row by row
   #     for rAdj, rKet in zip(n,k):
    #        a = []
            
     #       if rAdj!=0:
       #         a.insert(len(a),rKet) # inserting another element
      #      else:
        #        a.insert(len(a),0)
         #   ketSubjectAdj.insert(len(ketSubjectAdj),a)
    #ketDeltaA.append(np.array(ketSubjectAdj))#.reshape(-1)) # flatten the ,1 matrix to vector


# In[ ]:


# using adjacency matrix - take the real correaltion differences 
ketDeltaA = []
b = ad
for ket in (ket3_ket1):
    g = np.zeros([268,268])
    print(ket[b!=0])
    g[b!=0] = ket[b!=0]
    ketDeltaA.insert(len(ketDeltaA),g)
    
np.array(ketDeltaA).shape  
#print(f"First {np.array(ketDeltaA)[0,:,0:10]}")
#print(f"Sec {np.array(ketDeltaA)[1,:,0:10]}")
    
# This vectors (one per subject) will use for regression analysis, after we exclude all zeros from them (see below)


# In[ ]:


# creating same vector for midazolam group - just to help us comparing
midDeltaA = []
b = ad
for ket in (mid3_mid1):
    g = np.zeros([268,268])
    print(ket[b!=0])
    g[b!=0] = ket[b!=0]
    midDeltaA.insert(len(midDeltaA),g)
    
np.array(midDeltaA).shape  


# In[ ]:


p = ([1,2,3],[4,5,6], [7,8,9])
q = ([1,2,7],[4,5,9])
x = ([0,0,1],[0,0,1])
x= np.array(x)
p = np.array(p)
print(p)
print(p.reshape(-1))

 
    


# In[ ]:


allDeltaA = ketDeltaA + midDeltaA # combine all together
np.array(allDeltaA).shape


# In[ ]:


# We will build a simple average matrix (ket3 - ket1) according to the adjacency matrix (39,39) in order to vizualize
# results on the nodes themselves
ket3_1_average = np.average(np.array(ket3_ket1), axis = 0)

# now we replace all 1 in adj to real values and zeros will remain zeros
j = adj
j[j!=0] = ket3_1_average[j!=0] # replacing j values in real difference instead of 1
# now we can plot

plotting.plot_connectome(j, coords,
                        edge_threshold="98%", colorbar=True)#, output_file="Ket3_vsKet1.png")


plotting.plot_matrix(j, labels=atlas_labes, colorbar=True,
                     vmax=0.3, vmin=-0.3)


# In[ ]:


# this is a small snippet of code that takes the adajecnt's matrix that was created before and present it in diagonal 
# correlation matrix. 
#  change it into dataframe to add labels (using atlas labels)
import pandas as pd
d = pd.DataFrame(data=j,
                 columns=atlas_labes, index=atlas_labes)

mask = np.zeros_like(d, dtype=np.bool) # masking half
mask[np.triu_indices_from(mask)] = True


# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(15, 12))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns_plot = sns.heatmap(d, mask=mask, cmap=cmap, vmax=.1, center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5})


# In[ ]:


# now - just to make sure I understand which place in the vector (78 because its 39*2) if for which correlation
k = ad
k[k==1] = ket3_ket1[0][k==1] # replacing j values in real difference instead of 1

sub8 = pd.DataFrame(data=k,
                 columns=labels, index=labels)

mask = np.zeros_like(sub8, dtype=np.bool) # masking half
mask[np.triu_indices_from(mask)] = True


# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(15, 12))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(sub8, mask=mask, cmap=cmap, vmax=.3, center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5})


# In[ ]:


#print(ket3_ket1[1,:,0:15])

#print(k[:,38])
v = np.array(vecKetamin)
print(np.array(ketDeltaA).shape)
print(f"KDeltaArray {np.array(ketDeltaA)[0,21:39,0:10]}")
print(f"KDeltaVector {np.array(ketDeltaA)[0].reshape(-1)[0:39]}")
print(np.array(vecKetamin)[0])

#print(f"First Ket {ket1_corr[0][:,1:10]}")
#print(f"Third ket {ket3_corr[0][:,1:10]}")
#print(f"Substract {ket3_ket1[:,1:10]}")


# In[ ]:

# need to build regression model according to this vector for each subject
# change from matrix to vector
vecKetamin = [] # vector of all edges per subject
for i in ketDeltaA:
    # create a data frame 
    df = pd.DataFrame(i)
    # masking on diagonal (to avoid repeated numbers)
    dataCorr = df.mask(np.tril(np.ones(df.shape)).astype(np.bool))
    # remove less than 0.001
    dataCorr = dataCorr[abs(dataCorr) >= 0.001].stack().reset_index()
    a = np.array(dataCorr[0]).reshape(-1) # reshaping subject matrix to vector
    #print(a[a!=0]) 
    vecKetamin.append(a) # appending everything that's not zero 
np.array(vecKetamin).shape # shape of vector (subjects,length)


# In[ ]:

# # sorting correlation vector - using key=abs for absolute values
vecSort = []
for i in vecKetamin:
    print(f"Original {i}")
    a = sorted(i, key=abs, reverse=True)
    print(f"New {a}")
    vecSort.append(a)
print(np.array(vecSort).shape)
#%%
# create vector of all
vecAll = [] # vector of all edges per subject
for i in allDeltaA:
    a = np.array(i).reshape(-1) # reshaping subject matrix to vector
    #print(a[a!=0]) 
    vecAll.append(a[a!=0]) # appending everything that's not zero 
np.array(vecAll).shape # shape of


# In[ ]:


# sorting correlation vector - using key=abs for absolute values
matSort = []
vecSort = []
for i in vecKetamin:
    print(f"Original {i}")
    a = sorted(i, key=abs)
    print(f"New {a}")
    vecSort.append(a)


# In[ ]:


print(len(vecSort[0]))
print(len(set(vecSort[0])))


# In[ ]:


# lets run regression model to try and predict PCL scores with connectivity vector
# first load up pcl data
kpe_dat = pd.read_excel('/home/or/Documents/kpe_analyses/KPEIHR0009_data_all_scored.xlsx', index_col ="scr_id")


# In[ ]:


import statistics 
# subset the data frame
kpe_pcl = kpe_dat[['med_cond','pcl5_total_screen','pcl5_total_visit1', 'pcl5_total_visit7', 'pcl5_total_followup1','pcl5_total_followup2']]
# remove NaN's from FU2
for i,n in enumerate(kpe_pcl['pcl5_total_followup2']):
    if  pd.isna(n) == True:
        print('nan')
        print(kpe_pcl['pcl5_total_followup2'][i])
        kpe_pcl['pcl5_total_followup2'][i]= np.nanmean(kpe_pcl['pcl5_total_followup2'])
# remove NaN's from visit 1
for i,n in enumerate(kpe_pcl['pcl5_total_visit1']):
    if  pd.isna(n) == True:
        print('nan')
        print(kpe_pcl['pcl5_total_visit1'][i])
        kpe_pcl['pcl5_total_visit1'][i]= kpe_pcl['pcl5_total_screen'][i] # replace with screening data
kpe_pcl["PCLChange_3_1"] = kpe_pcl['pcl5_total_followup2'] - kpe_pcl["pcl5_total_visit1"]
print(kpe_pcl)
np.nanmean(kpe_pcl['pcl5_total_followup2'])


# In[ ]:


# now lets build regression model
# Y = pcl score at FU2
# X's = vector of edges from each subject
kpe_data_pcl = kpe_pcl.loc[['KPE008','KPE1293','KPE1307','KPE1322','KPE1339','KPE1343','KPE1387','KPE1223']]
print(kpe_data_pcl)


# In[ ]:


su = ['KPE008','KPE1293','KPE1307','KPE1322','KPE1339','KPE1343','KPE1387','KPE1223']

kpe_matrix_dat = pd.DataFrame(vecSort, index=su)#pd.DataFrame(vecKetamin, index=su)
print(kpe_matrix_dat.iloc[:,0:39])
# now we have flatten the matrix. 
# Remember, we now have double of everything, saw basically should split it in half. 


# In[ ]:


import statsmodels.api as sm

X = kpe_matrix_dat.iloc[:,0:3] #kpe_matrix_dat[1]
y = kpe_data_pcl["PCLChange_3_1"]

X = sm.add_constant(X)
# Note the difference in argument order
model = sm.OLS(y, X).fit()
predictions = model.predict(X) # make the predictions by the model

# Print out the statistics
print(model.summary())


# In[ ]:


fig, ax = plt.subplots(figsize=(12,8))
fig = sm.graphics.influence_plot(model, ax=ax, criterion="cooks")


# In[ ]:


allDF = kpe_data_pcl
allDF.merge(kpe_matrix_dat)
print(allDF)
#fig, ax = plt.subplots(figsize=(12,8))
#fig = sm.graphics.plot_partregress(kpe_matrix_dat[:,0], kpe_data_pcl["pcl5_total_followup2"], [kpe_matrix_dat[:,0:5],kpe_data_pcl["pcl5_total_followup2"]] , ax=ax)


# In[ ]:


# writing results as .tex file (so we can convert later to PDF or PNG)
beginningtex = """\\documentclass{report}
\\usepackage{booktabs}
\\begin{document}"""
endtex = "\end{document}"

f = open('myreg.tex', 'w')
f.write(beginningtex)
f.write(model.summary().as_latex())
f.write(endtex)
f.close()


# In[ ]:


# doing same regression model for all subjects
# first data frame
all_data_pcl = kpe_pcl.loc[['KPE008','KPE1293','KPE1307','KPE1322','KPE1339','KPE1343','KPE1387','KPE1223','KPE1263','KPE1351','KPE1364','KPE1369','KPE1390','KPE1403']]

suAll = ['KPE008','KPE1293','KPE1307','KPE1322','KPE1339','KPE1343','KPE1387','KPE1223','KPE1263','KPE1351','KPE1364','KPE1369','KPE1390','KPE1403']

all_matrix_dat = pd.DataFrame(vecAll, index=suAll)

X = all_matrix_dat.iloc[:,0:5] #kpe_matrix_dat[1]
y = all_data_pcl["pcl5_total_followup2"]

X = sm.add_constant(X)
# Note the difference in argument order
modelAll = sm.OLS(y, X).fit()
predictions = modelAll.predict(X) # make the predictions by the model

# Print out the statistics
print(modelAll.summary())


# In[ ]:


from sklearn import linear_model

lm = linear_model.LinearRegression()

X = kpe_matrix_dat #kpe_matrix_dat[1]
y = kpe_data_pcl["pcl5_total_followup2"]


#X = all_matrix_dat #kpe_matrix_dat[1]
#y = all_data_pcl["pcl5_total_followup2"]
# define the data/predictors as the pre-set feature names  
#df = pd.DataFrame(kpe_matrix_dat)


print(X.shape)
print(y.shape)
model = lm.fit(X,y)


# In[ ]:


predictions = lm.predict(X)


# In[ ]:


import bct
fig=bct.adjacency_plot_und(A = ad, coor = coords)
from mayavi import mlab 
mlab.show()


# In[ ]:


lm.score(X,y)
# This leads to perfect prediction. So lets split data for test-retest


# In[ ]:


from sklearn.model_selection import train_test_split
# create training and testing vars
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
print (X_train.shape, y_train.shape)
print (X_test.shape, y_test.shape)


# In[ ]:


# fit a model
lm = linear_model.LinearRegression()
model = lm.fit(X_train, y_train)
predictions = lm.predict(X_test)


# In[ ]:


predictions


# In[ ]:


plt.scatter(y_test, predictions)
plt.xlabel("True Values")
plt.ylabel("Predictions")


# In[ ]:


model.score(X_test, y_test)
# prediction is not great

#%%
d =  kpe_matrix_dat.iloc[:,0]
c = kpe_data_pcl['PCLChange_3_1']
ax = sns.regplot(c,d)
#ax.set_title('Correlation between highest change in connectivity and change in symptoms')
ax.set(xlabel='change in connectivity', ylabel='change in symptoms')
plt.show()

#dMid = kpe_pcl_mid["PCLChange_3_1"]
# sns.regplot(kpe_matrix_dat.iloc[:,1], d)
# plt.show()

#sns.regplot(high, dMid)
#%%
# check pearson corelation between highest change in connectivity and PCL change in symptoms
import scipy as sc
#sc.stats.pearsonr(kpe_matrix_dat.iloc[:,77],d)
sc.stats.pearsonr(d,c)
#%%
# check centrality and take the highest four
import operator

a= nx.degree_centrality(G) # using networkX to check centrality of nodes.
sorted_x = sorted(a.items(), key=operator.itemgetter(1),reverse=True)
print(sorted_x)

# should take the highest nodes in degree centrality (10?50?)
# then run analisys on them:
# 1. Seed based analysis
# 2. Regression analysis associating them with pcl change. 


#%% Graph LASSO procedure
from sklearn.covariance import GraphLassoCV

covariance_estimator = GraphLassoCV(verbose=2)

lassoTime = []
def runLasso(timeseries):
    lassoArray = []
    for i in timeseries:
        matrixLASSO = covariance_estimator.fit(i)
        lassoArray.append(matrixLASSO.covariance_)
    return lassoArray
#%%
lassoTime1 = runLasso(ket1_series)
lassoTime2 = runLasso(ket2_series)
lassoTime3 = runLasso(ket3_series)
    
# In[]:

plotting.plot_matrix(lassoTime1[1], labels=atlas_labes[:,1], colorbar=True,
                     vmax=0.8, vmin=-0.8)

plotting.plot_connectome(np.average(lassoTime1, axis=0), coords,
                         edge_threshold="99%", colorbar=True, title = "Lasso First Scan")

# In[]:
ket1ReshapeLasso = np.moveaxis(np.array(lassoTime1), 0,-1)
ket2ReshapeLasso = np.moveaxis(np.array(lassoTime2), 0,-1)
ket3ReshapeLasso = np.moveaxis(np.array(lassoTime3), 0,-1)
#print(mid3Reshape.shape)


# In[ ]:

# now we can run NBS
# NBS is taken from: https://github.com/aestrivex/bctpy, can be installed using pip (pip install bctpy)
from bct import nbs
# we compare ket1 and ket3
pvalLasso, adjLasso, _Lasso = nbs.nbs_bct(ket1ReshapeLasso, ket3ReshapeLasso, thresh=3.5, tail='both',k=1000, paired=True, verbose = True)
# check mean p vlue
#np.mean(checkNBS[0])
adLasso = adjLasso

#%% Change labels in atlas to names of netwroks
labels = pd.DataFrame(atlas_labes)
for i in labels[1]:
    print (i)
    if i==1:
        labels['network'] = 
adDatLasso = pd.DataFrame(adLasso, columns=atlas_labes[:,0],index=atlas_labes[:,0])
print(adDatLasso)


# In[ ]:

%matplotlib qt
# graph adjacency matrix
import networkx as nx
G = nx.from_numpy_matrix(np.average(ebicMat, axis=0)) 
#G = nx.DiGraph(adDat, with_labels=True)
Glasso = nx.from_pandas_adjacency(adDatLasso)
# drawing the adajency matrix
import matplotlib.pyplot as plt
#plt.figure(figsize=(70,50))
nx.draw_spring(G, with_labels=True) #, font_size = 50, width = 1, alpha = 0.7, font_weight = "bold")
#%% I graph
from igraph import *
import igraph
import cairocffi
H = igraph.Graph.Adjacency(np.average(ebicMat, axis=0).tolist())
igraph.Graph(H)

vcount = max(np.average(ebicMat, axis=0).shape)
sources, targets = np.average(ebicMat, axis=0).nonzero()
edgelist = list(zip(sources.tolist(), targets.tolist()))
g = igraph.Graph(vcount, edgelist)
%matplotlib qt
plot(g) #, layout = "fr3d")

#%% save array to txt
np.savetxt("firstLASSO.csv", np.array(np.average(lassoTime1, axis=0)), delimiter = ",")
np.savetxt("correlationMatrixFirstRun.csv", np.array(np.average(ket1_corr, axis=0)), delimiter = ",")

#%% 
# run analysis on all subjects
# build timeseries for first scan (before)
totalSub_series_1 = timeSeries(func_files=fileList(totalSub,'1'), confound_files=confList(totalSub, '1'))
totalSub_series_2 = timeSeries(func_files=fileList(totalSub,'2'), confound_files=confList(totalSub, '2'))
# In[skggm]:
# run EBIC LASSO on each subject
from inverse_covariance import QuicGraphLassoEBIC
# call estimator EBIC
estimator = QuicGraphLassoEBIC(init_method="cov", verbose=1)
# create emtpy array
ebicMat = []
ebicMat2 = []
# run loop over timeseries and create precision matrix for each subject, append to a list of matrices
for i in totalSub_series_1:
    estimator.fit(i)
    ebicMat.append(-estimator.precision_)
    

for i in totalSub_series_2:
    estimator.fit(i)
    ebicMat2.append(-estimator.precision_)

ket1_forLASSO.shape
#%%
import matplotlib.pyplot as plt
%matplotlib qt
# Display the sparse inverse covariance
plt.figure(figsize=(7.5, 7.5))
plt.imshow(
    np.triu(-estimator.precision_, 1), interpolation="nearest", cmap=plt.cm.RdBu_r
)
plt.title("Precision (Sparse Inverse Covariance) matrix")
plt.colorbar()

# And now display the corresponding graph
plotting.plot_connectome(
    -estimator.precision_,
    coords,
    title="Functional Connectivity using EBICLASSO",
    edge_threshold="99%",
    node_size=20,
)
plotting.show()

ebicMat = np.array(ebicMat)
ebicMat2 = np.array(ebicMat2)
plotting.plot_connectome(np.average(np.array(ebicMat), axis = 0),coords, title = "Average connectome for Ketamine subjects 1st scan", edge_threshold = "99%", node_size = 50, colorbar=True)

plotting.plot_connectome(np.average(np.array(ebicMat2), axis = 0),coords, title = "Average connectome for Ketamine subjects 2st scan", edge_threshold = "99%", node_size = 50, colorbar=True)
#%% Run NBS on those
allsub1_reshape= np.moveaxis(ebicMat, 0,-1)
allsub2_reshape= np.moveaxis(ebicMat2, 0,-1)

from bct import nbs

pvalALL1_2, adjAll1_2, _ = nbs.nbs_bct(allsub1_reshape, allsub2_reshape, thresh=2, tail='both',k=1000, paired=True, verbose = True)

#%%
allLASSOcv = runLasso(totalSub_series_1)
allLASSOcv = np.array(allLASSOcv)
