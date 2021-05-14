#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 17:25:43 2018

@author: kurikawa
"""
import numpy as np
import pylab as pl
import itertools as it


tinit=7000
thrs =20
ntbin=40

def gen_data_clust():
    
    coind_clust=[]
    for inet in range(5):
        for iinit in range(5):
    
            rast=np.loadtxt("Data_for_paper/rast_8080gamma28%d%d0.10.0025.dat" % (inet,iinit))
            rast_areas=np.zeros((ntbin,4))
    
            for i in range(len(rast[:,0])):
                t=int((rast[i,0]-tinit)/25.0)
                if rast[i,1]<3050:  # in L in Dec
                    rast_areas[t,0]+=1
                elif rast[i,1]<3100: # in R in Dec
                    rast_areas[t,1]+=1
                elif rast[i,1]>3200 and rast[i,1]<3440:
                    rast_areas[t,2]+=1 # in EC3
                elif rast[i,1]>3540 and rast[i,1]<3780:
                    rast_areas[t,3]+=1 # in CA1
    
            idx_clst_areas=np.zeros((ntbin,4))
            cnt=np.zeros(4)
            for t in range(5,ntbin):  # avoid initial spikes in Dec
                for iarea in range(4):
                    if rast_areas[t-1,iarea]<thrs and rast_areas[t,iarea]>=thrs:
                        cnt[iarea]+=1
                    if rast_areas[t,iarea]>=thrs:
                        idx_clst_areas[t,iarea]=cnt[iarea]
    
            t=0
            while t<ntbin-2:
                if idx_clst_areas[t,2]>0 and (idx_clst_areas[t,3]>0 or idx_clst_areas[t+1,3]>0 or idx_clst_areas[t+2,3]):
                    coind_clust.append([t,inet,iinit])
                    ttmp=t
                    spkstmp=0
                    while idx_clst_areas[ttmp,2]>0:
                        spkstmp+=rast_areas[ttmp,2]
                        ttmp+=1
                    tfwd=ttmp
                    ttmp=t-1
                    while idx_clst_areas[ttmp,2]>0:
                        spkstmp+=rast_areas[ttmp,2]
                        ttmp-=1
                    coind_clust[-1].append(spkstmp)
    
    
                    spkstmp=rast_areas[t,3]+rast_areas[t+1,3]
                    ttmp=t+2
                    while idx_clst_areas[ttmp,3]>0:
                        spkstmp+=rast_areas[ttmp,3]
                        ttmp+=1
                    tfwd=np.max([ttmp,tfwd])
                    ttmp=t-1
                    while idx_clst_areas[ttmp,3]>0:
                        spkstmp+=rast_areas[ttmp,3]
                        ttmp-=1
                    coind_clust[-1].append(spkstmp)
    
                    spkstmp=rast_areas[t,0]+rast_areas[t+1,0]+rast_areas[t+2,0]
                    if spkstmp>0:
                        ttmp=t+3
                        while idx_clst_areas[ttmp,0]>0:
                            spkstmp+=rast_areas[ttmp,0]
                            ttmp+=1
                        tfwd=np.max([ttmp,tfwd])
                        ttmp=t-1
                        while idx_clst_areas[ttmp,0]>0:
                            spkstmp+=rast_areas[ttmp,0]
                            ttmp-=1
                    t=tfwd
                    coind_clust[-1].append(spkstmp)
                else:
                    t+=1
    coind_clust=np.array(coind_clust)
    return coind_clust

def plt_bust_pre_post(coind_clust):
    pl.subplot(1,4,1)
    for i in range(len(coind_clust[:,0])):
        if coind_clust[i,0]>20:
            pl.scatter(coind_clust[i,3],coind_clust[i,5],c="b",lw=0)
        else:
            pl.scatter(coind_clust[i,3],coind_clust[i,5],c="r",lw=0)
    pl.subplot(1,4,2)
    for i in range(len(coind_clust[:,0])):
        if coind_clust[i,0]>20:
            pl.scatter(coind_clust[i,4],coind_clust[i,5],c="b",lw=0)
        else:
            pl.scatter(coind_clust[i,4],coind_clust[i,5],c="r",lw=0)
    pl.subplot(1,4,3)
    for i in range(len(coind_clust[:,0])):
        if coind_clust[i,0]>20:
            pl.scatter(coind_clust[i,3]+coind_clust[i,4],coind_clust[i,5],c="b",lw=0)
        else:
            pl.scatter(coind_clust[i,3]+coind_clust[i,4],coind_clust[i,5],c="r",lw=0)
    pl.subplot(1,4,4)
    for i in range(len(coind_clust[:,0])):
        if coind_clust[i,0]>20:
            pl.scatter(np.log(coind_clust[i,3])+np.log(coind_clust[i,4]),coind_clust[i,5],c="b",lw=0)
        else:
            pl.scatter(np.log(coind_clust[i,3])+np.log(coind_clust[i,4]),coind_clust[i,5],c="r",lw=0)

    return
            
        
def plt_prob_dec(coind_clust):
    ths=15
    prob_dec=[]
    for icond in range(3):  # condisions: ECIII, CA1 and summation of ECIII and CA1
        prob_dec.append([])
        for igamma in range(2):  # gamma is applied or not
            prob_dec[icond].append([])
            for i in range(6): # # of presynaptic spikes
                prob_dec[icond][igamma].append([])
    icond=0
    sbin=50
    max_spk=300
    ipre=3
    for i in range(len(coind_clust[:,0])):
        if coind_clust[i,0]>20:
            prob_dec[icond][0][min(int(coind_clust[i,ipre]/sbin),int(max_spk/sbin)-1)].append(coind_clust[i,5])
        else:
            prob_dec[icond][1][min(int(coind_clust[i,ipre]/sbin),int(max_spk/sbin)-1)].append(coind_clust[i,5])
    
    icond=1
    sbin=25
    max_spk=150
    ipre=4
    for i in range(len(coind_clust[:,0])):
        if coind_clust[i,0]>20:
            prob_dec[icond][0][min(int(coind_clust[i,ipre]/sbin),int(max_spk/sbin)-1)].append(coind_clust[i,5])
        else:
            prob_dec[icond][1][min(int(coind_clust[i,ipre]/sbin),int(max_spk/sbin)-1)].append(coind_clust[i,5])
    icond=2
    sbin=100
    max_spk=600
    
    for i in range(len(coind_clust[:,0])):
        if coind_clust[i,0]>20:
            prob_dec[icond][0][min(int((coind_clust[i,3]+coind_clust[i,4])/sbin),int(max_spk/sbin)-1)].append(coind_clust[i,5])
        else:
            prob_dec[icond][1][min(int((coind_clust[i,3]+coind_clust[i,4])/sbin),int(max_spk/sbin)-1)].append(coind_clust[i,5])
        
    tmpprob=np.zeros((3,2,6))    
    for icond, igamma, i in it.product(range(3),range(2),range(6)):
        if len(prob_dec[icond][igamma][i])!=0:
            tmpprob[icond,igamma,i]=len([1 for f in prob_dec[icond][igamma][i] if f>ths])/len(prob_dec[icond][igamma][i])
    sbin=[50,25,100]
    pl.figure()
    for itmp in range(3):
        pl.subplot(1,3,itmp+1)
        pl.plot(np.arange(6)*sbin[itmp],tmpprob[itmp].T)
        
    return
