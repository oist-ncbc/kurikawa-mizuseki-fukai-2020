# modify gen_fig.py to generate figures ver4 (2019.2.22)

import numpy as np
import pylab as pl
import os 
import itertools
import cmath
import math
from scipy import stats
import func as fc

"""
s-C 500-1500
R    1500-3500
delay 3500-5000


NE_EC3=240   0-240
NI_EC3=160    240-400

NE_EC5=400   405-805
NI_EC5=120   805-925

NE_CA1=480   930-1410
NI_PV =240  1410-1650
NI_OLM=240  1650-1890

NVIP=240    1895-2135
Nos=20      2135-2155
NCA3=480    2155-2635

"""
NE_EC3=240


def idx_cell_eacharea(area_name):
    if area_name=="EC3":
        (icell_s,icell_e)=(0,240)
    elif area_name=="EC5":
        (icell_s,icell_e)=(405,805)
    elif area_name=="CA1":
        (icell_s,icell_e)=(930,1410)
    elif area_name=="CA1_PV":
        (icell_s,icell_e)=(1410,1650)
    elif area_name=="CA1_OLM":
        (icell_s,icell_e)=(1650,1890)
    elif area_name=="EC3_PV":
        (icell_s,icell_e)=(240,400)
    elif area_name=="EC5_PV":
        (icell_s,icell_e)=(805,925)
    elif area_name=="CA1_VIP":
        (icell_s,icell_e)=(1895,2135)
    elif area_name=="CA3":
        (icell_s,icell_e)=(2155,2635)
    

    return (icell_s,icell_e)

def gen_often_used_data():
    rateCA1,rateCA1_PV,rateCA1_OLM,rateCA1_VIP,rateEC3,rateEC3_PV,rateEC5,rateEC5_PV,rateCA3=gf.main("rate_dyn",["normal"])
    rateCA1_dis,rateEC3_dis=gf.main("rate_dyn",["disord"])
    rateCA1_noPV,rateEC3_noPV=gf.main("rate_dyn",["noPV"])
    rateCA1_cut,rateEC3_cut,rateCA1_PV_cut,rateEC3_PV_cut=gf.main("rate_dyn",["cut"])
    
    np.save("rateCA1_all_normal",rateCA1) #etc:
    
    rate=gf.main("rate_dyn",["cut_EC3_CA1",[5000,6500]])
    np.save("rateEC3_all_cutEC3CA1_5000_6500.npy",rate)
    rate=gf.main("rate_dyn",["cut_EC3_CA1",[5000,6000]])
    np.save("rateEC3_all_cutEC3CA1_5000_6000.npy",rate)
    rate=gf.main("rate_dyn",["cut_EC3_CA1",[6500,8000]])
    np.save("rateEC3_all_cutEC3CA1_6500_8000.npy",rate)
    rate=gf.main("rate_dyn",["cut_EC3_CA1",[6000,8000]])
    np.save("rateEC3_all_cutEC3CA1_6000_8000.npy",rate)



    return

# log
"""
##### change file format #####  2019.2.28
rateCA1=np.zeros((25,480,1601))
rateEC3=np.zeros((25,240,1601))

for tmpname in ["rcl0.9","rcl1","rcl1.1","rcl1.2"]:  ## the same processes are done for "enc", "delay"
    for i in range(5):
        for j in range(5):
            fname="rateCA1_all_"+tmpname+("_%d_%d.npy" % (i,j))
            rateCA1[i*5+j]=np.load(fname)
            fname="rateEC3_all_"+tmpname+("_%d_%d.npy" % (i,j))
            rateEC3[i*5+j]=np.load(fname)
    np.save("rateCA1_all_"+tmpname+".npy",rateCA1)
    np.save("rateEC3_all_"+tmpname+".npy",rateEC3)
"""


# copy and excute for each figure
# common pre-processing process
"""
import numpy as np
import pylab as pl
%matplotlib qt
import imp
import gen_fig1 as gf
import func as fc
from scipy import stats
#common parameters
T=8000
Nnet,Ninit=5,5
sbin=10
origin=4 # 0 degree of phase in rate data 
period=int(100/sbin) 
"""

"""
# Fig1 E (just candidate)
N=len(rateCA1)
tmp_dis=np.zeros((N,2))
tmp=np.zeros((N,2))
Twin=(get_stage_duration("s-L")/sbin).astype(int)
for i in range(N):
    tmp[i,0]=np.mean(rateCA1[i,:120,Twin[0]:Twin[1]])
    tmp[i,1]=np.mean(rateCA1[i,120:240,Twin[0]:Twin[1]])
    tmp_dis[i,0]=np.mean(rateCA1_dis[i,:120,Twin[0]:Twin[1]])
    tmp_dis[i,1]=np.mean(rateCA1_dis[i,120:240,Twin[0]:Twin[1]])
    
pl.figure()
tmptmp=tmp[:,0]-tmp[:,1]
tmptmp_dis=tmp_dis[:,0]-tmp_dis[:,1]
pl.bar([0],np.mean(tmptmp))
pl.bar([1],np.mean(tmptmp_dis))
for i in range(N):
    pl.scatter(0,tmptmp[i],c="gray")
    pl.scatter(1,tmptmp_dis[i],c="gray")

    

# FigS1A
stagename="s-C"
tmp=get_stage_duration(stagename)
Twin=(tmp/sbin).astype(int)
N=len(rateCA1)  # N=Nnet*Ninit

pl.figure()
pl.subplot(8,1,1) # CA1
plt_theta_pref(np.mean(get_theta_phase_rate(rateCA1,Twin,[period,origin]),axis=0))
pl.subplot(8,1,2) # CA1 PV
plt_theta_pref(np.mean(get_theta_phase_rate(rateCA1_PV,Twin,[period,origin]),axis=0))
pl.subplot(8,1,3) # CA1 OLM
plt_theta_pref(np.mean(get_theta_phase_rate(rateCA1_OLM,Twin,[period,origin]),axis=0))

pl.subplot(8,1,4) # EC5
plt_theta_pref(np.mean(get_theta_phase_rate(rateEC5,Twin,[period,origin]),axis=0))
pl.subplot(8,1,5) # EC5 PV
plt_theta_pref(np.mean(get_theta_phase_rate(rateEC5_PV,Twin,[period,origin]),axis=0))

pl.subplot(8,1,6) # EC3
plt_theta_pref(np.mean(get_theta_phase_rate(rateEC3,Twin,[period,origin]),axis=0))
pl.subplot(8,1,7) # EC3 PV
plt_theta_pref(np.mean(get_theta_phase_rate(rateEC3_PV,Twin,[period,origin]),axis=0))

# FigS1B
stagename="s-L"
tmp=get_stage_duration(stagename)
Twin=(tmp/sbin).astype(int)
N=len(rateCA1)  # N=Nnet*Ninit
pl.figure()
cnt=1
pl.subplot(5,1,cnt) # CA3
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA3,Twin,[period,origin]),axis=0))
pl.subplot(5,1,cnt) # CA1
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA1_dis,Twin,[period,origin]),axis=0))
pl.subplot(5,1,cnt) # CA1PV
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA1_PV_dis,Twin,[period,origin]),axis=0))
pl.subplot(5,1,cnt) # EC3
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateEC3_dis,Twin,[period,origin]),axis=0))
plt_theta_pref((np.mean(get_theta_phase_rate(rateEC3,Twin,[period,origin]),axis=0))
pl.subplot(5,1,cnt) # EC3 PV
cnt+=1
plt_theta_pref(get_theta_phase_rate(rateEC3_PV_dis,Twin,[period,origin]),axis=0))

#FigS1C
stagename="s-L"
tmp=get_stage_duration(stagename)
Twin=(tmp/sbin).astype(int)
N=len(rateCA1)  # N=Nnet*Ninit
tmprate=np.zeros((N,4))

pl.figure()
for i in range(4):
    tmprate[:,i]=np.mean(rateCA1_dis[:,i*120:120*(i+1),Twin[0]:Twin[1]],axis=(1,2))
tmpmean,tmperr=np.mean(tmprate,axis=0),np.std(tmprate,axis=0)
pl.bar(range(4),tmpmean)
pl.errorbar(range(4),tmpmean,yerr=tmperr,c="r")
pl.ylim(0,7)
pl.xticks(range(4),["R","L","C","H"])
    
    
#FigS2AB
stagename="s-L"
tmp=get_stage_duration(stagename)
Twin=(tmp/sbin).astype(int)
N=len(rateCA1)  # N=Nnet*Ninit

pl.figure()
cnt=1
pl.subplot(4,1,cnt) # CA3
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA3,Twin,[period,origin]),axis=0))
pl.subplot(4,1,cnt) # CA1
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA1_cut,Twin,[period,origin]),axis=0))
pl.subplot(4,1,cnt) # CA1PV
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA1_PV_cut,Twin,[period,origin]),axis=0))
pl.subplot(4,1,cnt) # EC3
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateEC3_cut,Twin,[period,origin]),axis=0))
plt_theta_pref((np.mean(get_theta_phase_rate(rateEC3,Twin,[period,origin]),axis=0))

pl.figure()
cnt=1
pl.subplot(4,1,cnt) # CA3
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA3,Twin,[period,origin]),axis=0))
pl.subplot(4,1,cnt) # CA1
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA1_noPV,Twin,[period,origin]),axis=0))
pl.subplot(4,1,cnt) # CA1PV
cnt+=1
plt_theta_pref((np.mean(get_theta_phase_rate(rateCA1_PV_noPV,Twin,[period,origin]),axis=0))
pl.subplot(4,1,cnt) # EC3
cnt+=1
plt_theta_pref(np.mean(get_theta_phase_rate(rateEC3_noPV,Twin,[period,origin]),axis=0))
plt_theta_pref((np.mean(get_theta_phase_rate(rateEC3,Twin,[period,origin]),axis=0))


#Fig3B (+Fig3B1)
rastdata=np.loadtxt("Data_for_paper/rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.30_00.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma")
dyndata=np.loadtxt("Data_for_paper/dyn_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.30_00.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma")
pl.figure()
cnt=1
pl.subplot(5,1,cnt)
cnt+=1
pl.plot(range(0,7500), dyndata[500:8000,411:411+4*5:4],c="b")
pl.plot(range(0,7500), dyndata[500:8000,411+4*10:411+4*15:4],c="orange")
pl.xlim(0,7500)

pl.subplot(5,1,cnt)
cnt+=1
pl.plot(range(0,7500), dyndata[500:8000,412:412+4*5:4],c="b")
pl.plot(range(0,7500), dyndata[500:8000,412+4*10:412+4*15:4],c="orange")
pl.xlim(0,7500)

"""

"""
# Fig2B  ver4 
stagename="s-C"
gen_fig2Bv4(stagename,period,origin,sbin)
stagename="s-L"
gen_fig2Bv4(stagename,period,origin,sbin)

#Fig2C ver4
stagename="s-C"
gen_fig2Cv4(stagename)
stagename="s-L"
gen_fig2Cv4(stagename)
"""


def gen_FigS1C(sbin):
    rateCA1=np.load("rateCA1_all_normal.npy")
    rateCA1_dis=np.load("rateCA1_all_dis.npy")   # disorderd phase
    N=len(rateCA1)
    tmp_dis=np.zeros((N,2))
    tmp=np.zeros((N,2))
    Twin=(get_stage_duration("s-L")/sbin).astype(int)

    for i in range(N):
        tmp[i,0]=np.mean(rateCA1[i,:120,Twin[0]:Twin[1]])
        tmp[i,1]=np.mean(rateCA1[i,120:240,Twin[0]:Twin[1]])
        tmp_dis[i,0]=np.mean(rateCA1_dis[i,:120,Twin[0]:Twin[1]])
        tmp_dis[i,1]=np.mean(rateCA1_dis[i,120:240,Twin[0]:Twin[1]])
    
    pl.figure()
    diff_tmp=np.zeros((N,2))
    diff_tmp[:,0]=tmp[:,0]-tmp[:,1] # difference between L and R subgroup rates in nomral condition
    diff_tmp[:,1]=tmp_dis[:,0]-tmp_dis[:,1] # difference between L and R subgroup rates in disord condition
    pl.bar(range(2),np.mean(diff_tmp,axis=0))
    pl.plot(range(2),diff_tmp.T)
        
    print(stats.ttest_rel(diff_tmp[:,0], diff_tmp[:,1]))
    #Ttest_relResult(statistic=19.64782214394565, pvalue=2.6793871469272863e-16)
    
    return


def gen_theta_lock_behavior(args):  # args=Nnet,Ninit,T
    #Fig4A
    theta_peak=gen_theta_lock_index(args)
    tmp=np.array(theta_peak).reshape(7,Nnet*Ninit)
    pl.figure()
    pl.bar(range(7),np.mean(tmp,axis=1),alpha=0.1)
    for i in range(25):
        pl.scatter(range(7),tmp[:,i],c="gray")

def gen_rate_dyn_disconnect_condition():
    #Fig4B FigS4A
    rate=np.load("rateEC3_all_normal.npy")
    rate1=np.load("rateEC3_all_cutEC3CA1_500_1500.npy")
    rate2=np.load("rateEC3_all_cutEC3CA1_1500_3500.npy")
    rate3=np.load("rateEC3_all_cutEC3CA1_5000_6000.npy")
    rate4=np.load("rateEC3_all_cutEC3CA1_6000_8000.npy")
    tmpsbin=100
    rsbin=int(tmpsbin/sbin)
    nbin=int(800/rsbin)
    rate_g,rate1_g,rate3_g,rate2_g,rate4_g=np.zeros((2,nbin)),np.zeros((2,nbin)),np.zeros((2,nbin)),np.zeros((2,nbin)),np.zeros((2,nbin))
    for t in range(nbin):
        rate_g[0,t]=np.mean(rate[5:,:120,rsbin*t:rsbin*(t+1)])
        rate_g[1,t]=np.mean(rate[5:,120:,rsbin*t:rsbin*(t+1)])
        rate1_g[0,t]=np.mean(rate1[5:,:120,rsbin*t:rsbin*(t+1)])
        rate1_g[1,t]=np.mean(rate1[5:,120:,rsbin*t:rsbin*(t+1)])
        rate3_g[0,t]=np.mean(rate3[5:,:120,rsbin*t:rsbin*(t+1)])
        rate3_g[1,t]=np.mean(rate3[5:,120:,rsbin*t:rsbin*(t+1)])
        rate2_g[0,t]=np.mean(rate2[5:,:120,rsbin*t:rsbin*(t+1)])
        rate2_g[1,t]=np.mean(rate2[5:,120:,rsbin*t:rsbin*(t+1)])
        rate4_g[0,t]=np.mean(rate4[5:,:120,rsbin*t:rsbin*(t+1)])
        rate4_g[1,t]=np.mean(rate4[5:,120:,rsbin*t:rsbin*(t+1)])

    pl.figure()
    pl.plot(range(nbin),rate_g[0,:]-rate_g[1,:])
    pl.plot(range(nbin),rate1_g[0,:]-rate1_g[1,:])
    pl.plot(range(nbin),rate3_g[0,:]-rate3_g[1,:])
    pl.xlim(5,80)
    pl.xticks(np.arange(5,81,10),np.arange(0,75,10))

    pl.figure()
    pl.plot(range(nbin),rate_g[0,:]-rate_g[1,:])
    pl.plot(range(nbin),rate2_g[0,:]-rate2_g[1,:])
    pl.plot(range(nbin),rate4_g[0,:]-rate4_g[1,:])
    pl.xlim(5,80)
    pl.xticks(np.arange(5,81,10),np.arange(0,75,10))


def gen_suc_rate_disconnect_condition():
    #Fig4C
    # get rate of EC3
    rate1=np.load("rateEC3_all_cutEC3CA1_500_1500.npy")
    rate2=np.load("rateEC3_all_cutEC3CA1_1500_3500.npy")
    rate3=np.load("rateEC3_all_cutEC3CA1_5000_6000.npy")
    rate4=np.load("rateEC3_all_cutEC3CA1_6000_8000.npy")
    
    #rate1=main("rate_dyn",["cut_EC3_CA1",[500,1500]])
    #rate2=main("rate_dyn",["cut_EC3_CA1",[1500,3500]])
    #rate3=main("rate_dyn",["cut_EC3_CA1",[5000,6500]])
    #rate4=main("rate_dyn",["cut_EC3_CA1",[6000,8000]])
    rate=np.load("rateEC3_all_normal.npy")
    N=len(rate)
    
    suc_prob=np.zeros((5,N))
    for i in range(N):
        suc_prob[4,i]=get_suc_prob(rate[i],10)  #sbin=10
        suc_prob[0,i]=get_suc_prob(rate1[i],10)
        suc_prob[1,i]=get_suc_prob(rate2[i],10)
        suc_prob[2,i]=get_suc_prob(rate3[i],10)
        suc_prob[3,i]=get_suc_prob(rate4[i],10)

    
    pl.figure()
    pl.bar(range(5),np.mean(suc_prob,axis=1),alpha=0.2)
    for i in range(N):
        pl.scatter(range(5),suc_prob[:,i],c="gray")


    print(stats.ttest_rel(suc_prob[2], suc_prob[3]))
    #Ttest_relResult(statistic=4.697261161449587, pvalue=8.966735455349769e-05)
    print(stats.ttest_rel(suc_prob[4], suc_prob[3]))
    #Ttest_relResult(statistic=5.510180795556576, pvalue=1.149181357179599e-05)
    print(stats.ttest_rel(suc_prob[4], suc_prob[2]))
    #Ttest_relResult(statistic=1.4905963111745266, pvalue=0.14909432825561109)
    print(stats.ttest_rel(suc_prob[4], suc_prob[1]))
    #Ttest_relResult(statistic=10.22535931135964, pvalue=3.1751920498907013e-10)
    
    #pl.figure()
    #pl.plot(range(2),suc_prob[3::-1,:],c="gray")
    return


def gen_theta_pref(paras):  #paras=[Nnet,Ninit,sbin,period,origin]
    #Fig5A 5B
    [Nnet,Ninit,sbin,period,origin]=paras
    name_stages=["s-C","s-L","H2","t-C1","t-C2"]
    sum_total=gen_distr_theta("ext","suc",False,"_",name_stages,"normal",paras)
    pl.figure()
    for istage in range(len(name_stages)):
        pl.subplot(len(name_stages),1,istage+1)
        plt_theta_pref(sum_total[istage],period)
        #pl.ylim(0,5)
    sum_total=gen_distr_theta("PV","suc",False,"_",name_stages,"normal",paras)
    pl.figure()
    for istage in range(len(name_stages)):
        pl.subplot(len(name_stages),1,istage+1)
        plt_theta_pref(sum_total[istage],period)
        #pl.ylim(0,60)
    sum_total=gen_distr_theta("OLM","suc",False,"_",name_stages,"normal",paras)
    pl.figure()
    for istage in range(len(name_stages)):
        pl.subplot(len(name_stages),1,istage+1)
        plt_theta_pref(sum_total[istage],period)
        #pl.ylim(0,60)


    sum_total=gen_distr_theta("ext","suc",False,"_",name_stages,"normal",paras)
    mean_theta=gen_mean_theta_suc(sum_total,paras)
    pl.figure()
    pl.scatter(range(5),mean_theta)
    pl.scatter(range(5),mean_theta+1)

    
    return


def gen_theta_pref_compare_suc_fail(paras,trial_type):
    #Fig6C
    [Nnet,Ninit,tmpsbin,tmpperiod,tmporigin]=paras
    name_stages=["t-C2"]
    suc=gen_distr_theta("ext","suc",False,"_",name_stages,trial_type,paras)  # return value == [], there is no success trial
    fal=gen_distr_theta("ext","fal",False,"_",name_stages,trial_type,paras)
    
    pl.figure()
    pl.subplot(3,1,1)
    if len(suc)!=0: 
        tmp=suc[0]/np.sum(suc[0][:tmpperiod])    
        plt_theta_pref(tmp,tmpperiod)
    pl.subplot(3,1,2)
    if len(fal)!=0:
        tmp=fal[0]/np.sum(fal[0][:tmpperiod])    
        plt_theta_pref(tmp,tmpperiod)
    pl.subplot(3,1,3)
    if (len(suc)!=0) and (len(fal)!=0):
        tmp=suc[0]/np.sum(suc[0][:tmpperiod])-fal[0]/np.sum(fal[0][:tmpperiod])  
        plt_theta_pref(tmp,tmpperiod)

    return suc,fal
            
        

def gen_suc_rate_diff_Ach(args):
   
    flg=args[0]
    cogn_func=args[1]
    ach=args[2]
    Nnet,Ninit,sbin=args[3]
    suc_prob=np.zeros(Nnet*Ninit)
    #check success trials
    for inet in range(Nnet):
        for iinit in range(Ninit):
            #check performance
            tmp=np.load("rateEC3_all_"+cogn_func+"%g_%d_%d.npy" % (ach,inet,iinit))  # this file is generated as sbin=5
            rsbin=int(sbin/5)
            rate=np.zeros((len(tmp),int(len(tmp[0])/rsbin)+1))
            for t in range(  int(len(tmp[0])/rsbin)  ):
                rate[:,t]=np.mean(tmp[:,rsbin*t:rsbin*(t+1)],axis=1)
            suc_prob[inet*Ninit+iinit]=get_suc_prob(rate,sbin)
    return suc_prob
            
            
def gen_firing_rate_EC3_EC5_CA1():
    rateCA1=np.load("rateCA1_all_normal.npy")
    pl.figure()
    pl.subplot(4,1,1)
    pl.plot(range(801),np.mean(rateCA1[0,:120,0:],axis=0))
    pl.plot(range(801),np.mean(rateCA1[0,120:240,0:],axis=0))
    pl.plot(range(801),np.mean(rateCA1[0,240:360,0:],axis=0))
    pl.plot(range(801),np.mean(rateCA1[0,360:,0:],axis=0))
    pl.xlim(50,800)
    pl.xticks(np.arange(50,801,100),np.arange(0,751,100))

    pl.subplot(4,1,2)
    pl.plot(range(801),np.mean(rateCA1[0,:120,0:],axis=0))
    pl.plot(range(801),np.mean(rateCA1[0,120:240,0:],axis=0))
    pl.plot(range(801),np.mean(rateCA1[0,240:360,0:],axis=0))
    pl.plot(range(801),np.mean(rateCA1[0,360:,0:],axis=0))
    pl.xlim(50,250)
    pl.xticks(np.arange(50,251,50),np.arange(0,201,50))

    rateEC3=np.load("rateEC3_all_normal.npy")
    pl.subplot(4,1,3)
    pl.plot(range(801),np.mean(rateEC3[0,:120,0:],axis=0))
    pl.plot(range(801),np.mean(rateEC3[0,120:240,0:],axis=0))
    pl.xlim(50,800)
    pl.xticks(np.arange(50,801,100),np.arange(0,751,100))

    rateEC5=np.load("rateEC5_all_normal.npy")
    pl.subplot(4,1,4)
    pl.plot(range(801),np.mean(rateEC5[0,:200,0:],axis=0))
    pl.plot(range(801),np.mean(rateEC5[0,200:400,0:],axis=0))
    pl.xlim(50,800)
    pl.xticks(np.arange(50,801,100),np.arange(0,751,100))

    return




def gen_performance_diff_condition_disinhibition():
    rateEC3         =np.load("rateEC3_all_normal.npy")
    rateEC3_noPV=np.load("rateEC3_all_noPV.npy")
    rateEC3_cut   =np.load("rateEC3_all_cut.npy")
    
    N=len(rateEC3)
    suc_prob=np.zeros((N,3))
    for i in range(N):
        suc_prob[i,0]=get_suc_prob(rateEC3[i],10)  # 10 = sbin
        suc_prob[i,1]=get_suc_prob(rateEC3_noPV[i],10)
        suc_prob[i,2]=get_suc_prob(rateEC3_cut[i],10)

    pl.figure()
    pl.bar(range(2),np.mean(suc_prob[:,:2],axis=0))
    pl.plot(range(2),suc_prob[:,:2].T)
    
    pl.figure()
    pl.bar(range(2),np.mean(suc_prob[:,[0,2]],axis=0))
    pl.plot(range(2),suc_prob[:,[0,2]].T)
    
    #print(stats.ttest_rel(suc_prob[:,0], suc_prob[:,1]))
    #Ttest_relResult(statistic=5.058857051419466, pvalue=3.581134754840131e-05)

    #print(stats.ttest_rel(suc_prob[:,0], suc_prob[:,2]))
    #Ttest_relResult(statistic=9.42530943548034, pvalue=1.541603268417448e-09)
    
    
    # noPV
    rastdata=np.loadtxt("Data_for_paper/rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.30_00.015.datdec4long1.130200.0150.0180.018080gamma30_PVinact")
    pl.figure()
    tmprast=np.array([x for x in rastdata if x[1]>=0 and x[1]<120])
    pl.scatter(tmprast[:,0],tmprast[:,1],s=0.5,c="b")
    tmprast=np.array([x for x in rastdata if x[1]>=120 and x[1]<240])
    pl.scatter(tmprast[:,0],tmprast[:,1],s=0.5,c="orange")
    pl.xlim(500,8000)
    pl.ylim(0,240)
    pl.xticks(range(500,8000,2000),[0,2000,4000,6000])
    # cut
    rastdata=np.loadtxt("Data_for_paper/rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.30_00.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma")
    pl.figure()
    tmprast=np.array([x for x in rastdata if x[1]>=0 and x[1]<120])
    pl.scatter(tmprast[:,0],tmprast[:,1],s=0.5,c="b")
    tmprast=np.array([x for x in rastdata if x[1]>=120 and x[1]<240])
    pl.scatter(tmprast[:,0],tmprast[:,1],s=0.5,c="orange")
    pl.xlim(500,8000)
    pl.ylim(0,240)
    pl.xticks(range(500,8000,2000),[0,2000,4000,6000])
    return
    
def gen_Fig1EF():
    rateEC3=np.load("rateEC3_all_normal.npy")
    rateEC3_dis=np.load("rateEC3_all_dis.npy")
    N=len(rateEC3)
    tmp=np.zeros((N,2))
    tmp_dis=np.zeros((N,2))
    for i in range(N):
        tmp[i,0]=np.mean(rateEC3[i,:120,750:])
        tmp[i,1]=np.mean(rateEC3[i,120:,750:])
        tmp_dis[i,0]=np.mean(rateEC3_dis[i,:120,750:])
        tmp_dis[i,1]=np.mean(rateEC3_dis[i,120:,750:])

    pl.figure()
    pl.bar(range(2),np.mean(tmp,axis=0))
    pl.plot(range(2),tmp.T)
        
    #print(stats.ttest_rel(tmp[:,0], tmp[:,1]))
    #Ttest_relResult(statistic=6.0422073517776065, pvalue=3.071682976369521e-06)                           
    
    pl.figure()
    suc=np.zeros((N,2))
    for i in range(N):
        suc[i,0]=inf_prob(tmp[i,0],tmp[i,1])
        suc[i,1]=inf_prob(tmp_dis[i,0],tmp_dis[i,1])

    pl.bar(range(2),np.mean(suc,axis=0) )
    for i in range(N):
        pl.scatter(0,suc[i,0],c="gray")
        pl.scatter(1,suc[i,1],c="gray")
    pl.yticks([0,0.5,1.0])
    
    #print(stats.ttest_rel(suc[:,0], suc[:,1]))
    #Ttest_relResult(statistic=8.136991969689886, pvalue=2.333707331559655e-08)
    return


def main(fig,args):
    T=8000
    Nnet=5
    Ninit=5
    sbin=10
    origin=4 # 0 degree of phase in rate data 
    period=int(100/sbin)
    ADR="Data_for_paper/"
    
    if fig=="rate_dyn":  #2018.6.26
        rateCA1,rateEC3,rateEC5=[],[],[]
        rateCA1_PV,rateCA1_OLM,rateCA1_VIP,rateEC3_PV,rateEC5_PV,rateCA3=[],[],[],[],[],[]
        rateall=[]
        for i in range(len(args[1])):
            rateall.append([])
            
        for inet in range(5):
            for init in range(5):
                common_name="rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma30" % (inet,init)
                if args[0]=="normal": 
                #inet,init=4,3  # for some figure
                    fname=ADR+common_name+"_nogamma"
                elif args[0]=="Ach":
                    delay_ach,rcl_ach,enc_ach=args[1] 
                    fname1="rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_"
                    fname2="cCA1EC50.005_EC51e-05_gCAN0.02GEI_EC50.2GIE_EC50.3"
                    fname3="0.015.datdec4long30200.0150.0180.018080gamma30_nogamma"
                    fname=ADR+fname1+fname2+"%d_%d" %(inet,init)+fname3+"_examAch%g%g%g" % (delay_ach,rcl_ach,enc_ach)
                elif args[0]=="disord":
                    fname=ADR+common_name+"_disord"
                elif args[0]=="noPV":
                    fname=ADR+common_name+"_PVinact"
                elif args[0]=="cut":
                    fname="Data_for_paper/rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma" % (inet,init)
                elif args[0]=="cut_EC3_CA1":
                    Twin=args[1]
                    fname=ADR+common_name+"%d_%d" % (Twin[0],Twin[1])
                elif args[0]=="notheta":
                    fname=ADR+common_name+"_notheta%g" % args[2]
                tmp=np.loadtxt(fname)
                
                for i in range(len(args[1])):
                    celltype=args[1][i]
                    rateall[i].append(fc.get_rate(tmp,[idx_cell_eacharea(celltype),int(8000/sbin),1.0/float(sbin)]))
                    
                """ order script
                rateCA1.append(fc.get_rate(tmp,[idx_cell_eacharea("CA1"),int(8000/sbin),1.0/float(sbin)]))
                rateEC3.append(fc.get_rate(tmp,[idx_cell_eacharea("EC3"),int(8000/sbin),1.0/float(sbin)]))
                
                #rateEC3_PV.append(fc.get_rate(tmp,[idx_cell_eacharea("EC3_PV"),int(8000/sbin),1.0/float(sbin)]))    
                if args[0]=="normal":
                    rateCA1_PV.append(fc.get_rate(tmp,[idx_cell_eacharea("CA1_PV"),int(8000/sbin),1.0/float(sbin)]))
                    rateCA1_OLM.append(fc.get_rate(tmp,[idx_cell_eacharea("CA1_OLM"),int(8000/sbin),1.0/float(sbin)]))
                    rateCA1_VIP.append(fc.get_rate(tmp,[idx_cell_eacharea("CA1_VIP"),int(8000/sbin),1.0/float(sbin)]))
                    rateEC5.append(fc.get_rate(tmp,[idx_cell_eacharea("EC5"),int(8000/sbin),1.0/float(sbin)]))
                    rateEC5_PV.append(fc.get_rate(tmp,[idx_cell_eacharea("EC5_PV"),int(8000/sbin),1.0/float(sbin)]))
                    rateCA3.append(fc.get_rate(tmp,[idx_cell_eacharea("CA3"),int(8000/sbin),1.0/float(sbin)]))
                elif args[0]=="noPV" or "cut":
                    rateCA1_PV.append(fc.get_rate(tmp,[idx_cell_eacharea("CA1_PV"),int(8000/sbin),1.0/float(sbin)]))
                 """   
        
        #rateCA1=np.array(rateCA1)
        #rateEC3=np.array(rateEC3)
        return rateall
        """ 
        if args[0]=="normal":
            rateCA1_PV,rateCA1_OLM,rateCA1_VIP=np.array(rateCA1_PV),np.array(rateCA1_OLM),np.array(rateCA1_VIP)
            rateEC5,rateEC5_PV,rateEC3_PV          =np.array(rateEC5),np.array(rateEC5_PV),np.array(rateEC3_PV)
            rateCA3                                               =np.array(rateCA3)
            return rateCA1,rateCA1_PV,rateCA1_OLM,rateCA1_VIP,rateEC3,rateEC3_PV,rateEC5,rateEC5_PV,rateCA3
        elif args[0]=="noPV" or "cut" :
            rateCA1_PV=np.array(rateCA1_PV)
            return rateCA1,rateEC3,rateCA1_PV
        else:
            return rateEC3
        """
    elif fig=="FigS1C":
        gen_FigS1C(sbin)
    elif fig=="FigS2E":
        gen_performance_diff_condition_disinhibition()
    elif fig=="Fig3B1":
        gen_firing_rate_EC3_EC5_CA1()
    elif fig=="Fig4A":
        gen_theta_lock_behavior([Nnet,Ninit,T])
    elif fig=="Fig4B":
        gen_rate_dyn_disconnect_condition()
    elif fig=="Fig4C":
        gen_suc_rate_disconnect_condition()
    elif fig=="Fig5AB":
        gen_theta_pref([Nnet,Ninit,sbin,period,origin])
    elif fig=="Fig6A":
        N=Nnet*Ninit
        
        #encode"
        Achset=[0.7,0.9,1.0,1.1]
        suc_prob=np.zeros((N,len(Achset)))
        for i in range(len(Achset)):
            suc_prob[:,i]=gen_suc_rate_diff_Ach(["suc","enc",Achset[i],[Nnet,Ninit,sbin]])
        pl.figure()
        pl.bar(range(len(Achset)),np.mean(suc_prob,axis=0))
        for i in range(N):
            pl.scatter(range(len(Achset)),suc_prob[i],c="gray")
        pl.xticks(range(4),Achset)
        
        #delay"
        Achset=[0.2,0.5,0.7,0.8]
        suc_prob=np.zeros((N,len(Achset)))
        for i in range(len(Achset)):
            suc_prob[:,i]=gen_suc_rate_diff_Ach(["suc","delay",Achset[i],[Nnet,Ninit,sbin]])
        pl.figure()
        pl.bar(range(len(Achset)),np.mean(suc_prob,axis=0))
        for i in range(N):
            pl.scatter(range(len(Achset)),suc_prob[i],c="gray")
        pl.xticks(range(4),Achset)
            
        #recall"
        Achset=[0.9,1.0,1.1,1.2]
        suc_prob=np.zeros((N,len(Achset)))
        for i in range(len(Achset)):
            suc_prob[:,i]=gen_suc_rate_diff_Ach(["suc","rcl",Achset[i],[Nnet,Ninit,sbin]])
        pl.figure()
        pl.bar(range(len(Achset)),np.mean(suc_prob,axis=0))
        for i in range(N):
            pl.scatter(range(len(Achset)),suc_prob[i],c="gray")
        pl.xticks(range(4),Achset)
    elif fig=="Fig6C":
        tmpsbin=5
        tmpperiod,tmporigin=int(100/tmpsbin),int(origin*(sbin/tmpsbin))
        paras=[Nnet,Ninit,tmpsbin,tmpperiod,tmporigin]
        suc,fal           =gen_theta_pref_compare_suc_fail(paras,"enc1")
        sucenc,falenc=gen_theta_pref_compare_suc_fail(paras,"enc0.7")
        sucrcl,falrcl   =gen_theta_pref_compare_suc_fail(paras,"rcl0.9")
        
        #plot diff between normal success and abnormal fal.
        pl.figure()
        pl.subplot(2,1,1)
        tmp=suc[0]/np.sum(suc[0][:tmpperiod])-falenc[0]/np.sum(falenc[0][:tmpperiod])  
        plt_theta_pref(tmp,tmpperiod)
        pl.subplot(2,1,2)
        tmp=suc[0]/np.sum(suc[0][:tmpperiod])-falrcl[0]/np.sum(falrcl[0][:tmpperiod])  
        plt_theta_pref(tmp,tmpperiod)
        
        # plot theta preference of CA1 superimposed with differenet subgroups
        issuc=get_suc_prob_allnets(Nnet,Ninit,"rateEC3_all_enc1.npy",tmpsbin)
        Twin=(get_stage_duration("t-C2")/tmpsbin).astype(int)

        rate=np.load("rateCA1_all_enc1.npy")
        tmptmp=np.zeros((3,tmpperiod*2))
        for isubgroup in range(3):
            if isubgroup==0:  # subgroup L
                tmp=get_theta_phase_rate(rate[:,:120],Twin,[tmpperiod,tmporigin])/4.0  # factor 4 is normalization factor when comapring rate[:,:120] to rate[:,:]
            elif isubgroup==1: #subgroup C
                tmp=get_theta_phase_rate(rate[:,240:360],Twin,[tmpperiod,tmporigin])/4.0
            elif isubgroup==2: #all
                tmp=get_theta_phase_rate(rate[:,:],Twin,[tmpperiod,tmporigin])
                
            for i in [k for k in range(Nnet*Ninit) if issuc[k]==1]:
                tmptmp[isubgroup]+=tmp[i]
        tmptmp/=np.sum(issuc)
        tmptmp/=np.sum(tmptmp[2][:tmpperiod])
        
        pl.figure()
        pl.bar(range(tmpperiod*2),tmptmp[2])
        pl.bar(range(tmpperiod*2),tmptmp[0])
        pl.bar(range(tmpperiod*2),tmptmp[1],bottom=tmptmp[0])
        #trough = 0 in LFP
        pl.plot(np.arange(0,tmpperiod*2,0.1),(np.sin(2*3.14159*(np.arange(0,tmpperiod*2,0.1)/float(tmpperiod) - 0.25))+1)*0.5*np.max(tmptmp[2]) )
        pl.xticks(np.arange(0,tmpperiod*2+1,5),np.arange(0,720+1,180))
  
        
    elif fig=="theta_preference_bf_aft_Ach":
        plt_theta_pref_bf_aft_Ach()
    elif fig=="Ca_dyn": #2018.10.12
        plt_fig_Ca_dyn()
    elif fig=="theta_normal_suc": #2018.10.12
        paras=[Nnet,Ninit,sbin,period,origin]
        gen_distr_theta_normal_suc("ext",0,False,"_",paras)
        # for different types of neurons in CA1
        gen_distr_theta_normal_suc("PV",0,False,"_")
        gen_distr_theta_normal_suc("OLM",0,False,"_")
        
    elif fig=="theta_comparison": #2018.12.12
        """
        phases=np.zeros((4,5))
        # for success trials in normal condition with mixing of PV cells
        phases[0]=gen_distr_theta_normal_suc("ext",0,True,0.0)
        phases[1]=gen_distr_theta_normal_suc("ext",0,True,0.1)
        phases[2]=gen_distr_theta_normal_suc("ext",0,True,0.25)
        phases[3]=gen_distr_theta_normal_suc("ext",0,True,0.5)

        cset=["b","cyan","green","r"]
        for i in range(4):
            pl.scatter(phases[i],range(5,0,-1),c=cset[i])
            pl.scatter(phases[i]+1.0,range(5,0,-1),c=cset[i])
            pl.xlim(0,2)
        """    
        # for success and failure trials in normal condition with mixing of PV cells
        # failures in normal condition (*enc1*), in low Ach (encoding phase) condition (*enc0.7), in low Ach (recall phase) condition (rcl1)
        phases=np.zeros((4,5))
        phases[0]=gen_distr_theta_normal_suc("ext",0,False,0.0)
        phases[1]=gen_distr_theta_normal_suc("ext",1,False,0.0)

        pl.figure()
        cset=["b","cyan"]
        for i in range(2):
            pl.scatter(phases[i],range(5,0,-1),c=cset[i])
            pl.scatter(phases[i]+1.0,range(5,0,-1),c=cset[i])
            pl.xlim(0,2)
            
    elif fig=="suc_fail_norm_Ach":
        #plot distrution of success trials: args[0]="suc"  args[1]="enc"  args[2]=[Ach]
        gen_distr_theta_normal_Ach(args)
    elif fig=="suc_rate_diff_Ach":
        sucrate=np.zeros((3,4))
        for i in range(4):
            issuc=chk_if_suc(["suc","enc",[0.7,0.9,1.0,1.1][i]])
            sucrate[0,i]=np.mean(issuc)
            issuc=chk_if_suc(["suc","delay",[0.2,0.5,0.7,0.8][i]])
            sucrate[1,i]=np.mean(issuc)
            issuc=chk_if_suc(["suc","rcl",[0.9,1.0,1.1,1.2][i]])
            sucrate[2,i]=np.mean(issuc)
        pl.figure()
        pl.subplot(3,1,1)
        pl.bar(range(4),sucrate[0],align="center")
        pl.xticks(range(4), [0.7,0.9,1.0,1.1])
        pl.ylim(0,1)
        pl.subplot(3,1,2)
        pl.bar(range(4),sucrate[1],align="center")
        pl.xticks(range(4), [0.2,0.5,0.7,0.8])
        pl.ylim(0,1)
        pl.subplot(3,1,3)
        pl.bar(range(4),sucrate[2],align="center")
        pl.xticks(range(4),[0.9,1.0,1.1,1.2])
        pl.ylim(0,1)
    elif fig=="Ach_sched":
        cogn_func,ioriginal=args[0],-2
        
        fname="rateEC3_all_"+cogn_func+".npy"
        if os.path.exists(fname):
            rateEC3_all=np.load(fname)
            fname="rateCA1_all_"+cogn_func+".npy"
            rateCA1_all=np.load(fname)
            fname="theta_pref_all_"+cogn_func+".npy"
            theta_pref_all=np.load(fname)
        else:
            theta_pref_all,rateEC3_all,rateCA1_all=gen_theta_rate(cogn_func)
            np.save(fname,rateEC3_all)
            fname="rateCA1_all_"+cogn_func+".npy"
            np.save(fname,rateCA1_all)
            fname="theta_pref_all_"+cogn_func+".npy"
            np.save(fname,theta_pref_all)

        # gen scores for differnt Ach schedule   
        npara=len(theta_pref_all)
        
        tmp=[]
        for i in range(npara):
            tmp.append([])
        for i in range(npara):
            rateEC3=rateEC3_all[i]
            for inet in range(Nnet):
                for init in range(Ninit):        
                    tmp[i].append(rateEC3[6,inet,init,0]-rateEC3[6,inet,init,1])
        pl.figure()
        cset=["b","cyan","green","magenta","r"]
        for i in range(npara):
            for j in range(Nnet):
                for k in range(Ninit):
                    pl.scatter(i,tmp[i][j*Ninit+k],c=cset[j])
        pl.plot(range(npara),np.mean(tmp,axis=1))
         
        
        pl.figure()
        cset=["g","gray"]
        tmptheta=np.zeros(20)
        for j in range(2):
            if j==1:
                theta_pref=theta_pref_all[0]
            elif j==0:
                theta_pref=theta_pref_all[ioriginal]            

            for i in range(7):
                tmptheta[:10]=theta_pref[i]
                tmptheta[10:]=theta_pref[i]
                pl.subplot(7,1,i+1)
                pl.plot(np.arange(0,20),tmptheta,color=cset[j])
                if j==0:
                    pl.plot(np.arange(0,20,0.1),(1-np.sin(2*3.14159*np.arange(0,20,0.1)/10.0))*0.5*max(tmptheta)*0.25+max(tmptheta))
    return






def get_stage_duration(stagename):
    if stagename=="s-C":
        tstart,tend=500,1500
    elif stagename=="s-L":
        tstart,tend=1500,3500
    elif stagename=="s-L1":
        tstart,tend=1500,2500
    elif stagename=="s-L2":
        tstart,tend=2500,3500
    elif stagename=="H":
        tstart,tend=3500,5000
    elif stagename=="H1":
        tstart,tend=3500,4250
    elif stagename=="H2":
        tstart,tend=4250,5000
    elif stagename=="t-C":
        tstart,tend=5000,8000
    elif stagename=="t-C1":
        tstart,tend=5000,6000
    elif stagename=="t-C2":
        tstart,tend=6000,8000

    return np.array([tstart,tend])

def get_theta_phase_rate(rate,Twin,paras):
    tmpperiod,tmporigin=paras
    if len(rate.shape)!=3:  # rate[0]:Nnet*Ninit, rate[_,0]:# of neurons, rate[_,_,0]: time
        print("invalid array of rate!! \n")
        return 
    
    N=len(rate)
    phase_off_set=Twin[0]%tmpperiod
    
    tmp=np.zeros((N,tmpperiod*2))
    for itrial in range(N):
        tmptmp=np.zeros(tmpperiod*2)
        for i in range(tmpperiod):
            tmptmp[(i+phase_off_set)%tmpperiod]=np.mean(rate[itrial,:,Twin[0]+i:Twin[1]:tmpperiod])  # Hz /neuron
        tmp[itrial]=tmptmp

    tmp[:,tmpperiod:]=tmp[:,:tmpperiod] 
    thetarate=np.zeros(tmp.shape)
    thetarate[:,:-tmporigin]=tmp[:,tmporigin:]
    thetarate[:,-tmporigin:]=tmp[:,:tmporigin]

    return thetarate

def plt_theta_pref(thetadata,tmpperiod):
    #pl.bar(range(tmpperiod*2),thetadata,width=1.0)
    pl.bar(range(tmpperiod*2),thetadata)
    #trough = 0 in LFP
    pl.plot(np.arange(0,tmpperiod*2,0.1),(np.sin(2*3.14159*(np.arange(0,tmpperiod*2,0.1)/float(tmpperiod) - 0.25))+1)*0.5*np.max(thetadata) )
    pl.xticks(np.arange(0,tmpperiod*2+1,5),np.arange(0,720+1,180))
    
    return


def inf_prob(rate1,rate2):
    return np.exp(rate1)/(np.exp(rate1)+np.exp(rate2))

def gen_theta_lock_index(args):
    Nnet,Ninit,T=args
    tmpsbin=50  #[ms]
    t_rate=int(T/tmpsbin)
    inv_sbin=1.0/float(tmpsbin)
    #name_stages=["s-C","s-L1","s-L2","H1","H2","t-C1","t-C2"]
    name_stages=["t-C1","t-C2"]
      
    theta_peak=[]
    for tmpname in name_stages:
        tstage=(get_stage_duration(tmpname)/tmpsbin).astype(int)
        for inet in range(Nnet):
            for init in range(Ninit):
                fname="Data_for_paper/rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma" % (inet,init)
                tmp=np.loadtxt(fname)
                rate=fc.get_rate(tmp,[idx_cell_eacharea("EC3"),int(T*inv_sbin),inv_sbin])
        
                tmp_rate=np.mean(rate[:240,tstage[0]:tstage[1]],axis=0)  # [Hz / cell]
                cor=np.correlate(tmp_rate,tmp_rate,"full")
                theta_peak.append(np.max(cor[int(len(cor)/2)+1:int(len(cor)/2)+4])/float((tstage[1]-tstage[0])*tmpsbin))

    return theta_peak  

def get_suc_prob(rate,sbin):
    if len(rate)!=NE_EC3 or (sbin==10 and len(rate[0])!=801) or (sbin==5 and len(rate[0])!=1601):
        print("invalid rate data!! in get_suc_prob")
        return
    if sbin==10: 
        tmp1=np.mean(rate[:120,750:])
        tmp2=np.mean(rate[120:,750:])
    elif sbin==5:
        tmp1=np.mean(rate[:120,1500:])
        tmp2=np.mean(rate[120:,1500:])
    suc_prob=inf_prob(tmp1,tmp2)
    
    return suc_prob


def get_suc_prob_allnets(Nnet,Ninit,fname,tmpsbin):
    issuc=np.zeros(Nnet*Ninit)
    rate=np.load(fname)
    for i in range(Nnet*Ninit):
        if get_suc_prob(rate[i],tmpsbin)>0.6:
            issuc[i]=1
        else:
            issuc[i]=0

    return issuc

def gen_distr_theta(type_neuron,is_suc,is_mixed,mix_rate,name_stages,trial_type,paras):  
    # is_suc="suc", only success trials are used
    Nnet,Ninit,tmpsbin,tmpperiod,tmporigin=paras

    #check success trials    
    fname="rateEC3_all_"+trial_type+".npy"
    issuc=get_suc_prob_allnets(Nnet,Ninit,fname,tmpsbin)
    
    ntrial=0
    if is_suc=="suc":
        ntrial=np.sum(issuc)
    elif is_suc=="fal":
        ntrial=Nnet*Ninit-np.sum(issuc)
    elif is_suc=="all":
        ntrial=Nnet*Ninit
        
    if ntrial==0:
        print("there is no hit trial")        
        return []
    
    sum_total=np.zeros((len(name_stages),tmpperiod*2))
    
    
    for istage in range(len(name_stages)):
        tmpname=name_stages[istage]
        Twin=(get_stage_duration(tmpname)/tmpsbin).astype(int)

        if tmpsbin==10 and trial_type=="normal":
            if type_neuron=="ext":
                fname="rateCA1_all_"+trial_type+".npy"
                rate=np.load(fname)
            elif type_neuron=="PV":
                fname="rateCA1PV_all_"+trial_type+".npy"
                rate=np.load(fname)
            elif type_neuron=="OLM":
                fname="rateCA1OLM_all_"+trial_type+".npy"
                rate=np.load(fname)
        elif tmpsbin==5 and (trial_type in ["enc0.7","enc0.9","enc1","enc1.1","delay0.2","delay0.5","delay0.7","delay0.8","delay0.9","rcl0.9","rcl1","rcl1.1","rcl1.2"]):
            if type_neuron=="ext":
                fname="rateCA1_all_"+trial_type+".npy"
                rate=np.load(fname)
        else:
            print("not yet")
            return
        
        if Nnet*Ninit!=len(rate):
            print("invalid rate in gen_distr_theta!!")
            return

        if not is_mixed:
            tmp=get_theta_phase_rate(rate,Twin,[tmpperiod,tmporigin])
            
            for i in range(Nnet*Ninit):
                if (is_suc=="suc" and issuc[i]==0) or (is_suc=="fal" and issuc[i]==1) :  # when is_suc="suc", fails are skipped. when "fal", successes are skipped
                    continue
                sum_total[istage]+=tmp[i]
    
    
    return sum_total/float(ntrial)
            


    """
    else:             
        if issuc[i]==is_suc:  # if is_suc==0, fails are skipped.
            continue
        rate=np.load("rateCA1_all_enc1_%d_%d.npy" % (inet,iinit))                
        rate1=np.load("rateCA1_PV_all_enc1_%d_%d.npy" % (inet,iinit))
        tmp=np.zeros(period*2)
        for i in range(period):
            tmp[(i+phase_off_set)%period]=np.sum(rate[:,Twin[0]+i:Twin[1]:period])
        tmp[period:]=tmp[:period] 
        sum_total[istage,2]+=tmp
        tmp=np.zeros(period*2)
        for i in range(period):
            tmp[(i+phase_off_set)%period]=np.sum(rate1[:int(240*mix_rate),Twin[0]+i:Twin[1]:period])
        tmp[period:]=tmp[:period] 
        sum_total[istage,2]+=tmp
    """
        
        
    return 
    
def gen_mean_theta_suc(sum_total,paras):  # is_suc=0, only success trials are used    
    Nnet,Ninit,tmpsbin,tmpperiod,tmporigin=paras
    Nstage=len(sum_total)
    # for cal of mean phase
    tmp=np.arange(0,2*math.pi,2*math.pi/float(tmpperiod))
    mean_phase=np.zeros(5,dtype=complex)    
    wgt=np.zeros(tmpperiod,dtype=np.complex)
    for i in range(tmpperiod):
        wgt[i]=cmath.exp(1j*tmp[i])
    
    for i in range(Nstage):
        mean_phase[i]=np.mean(sum_total[i][:tmpperiod]*wgt)
        
    phases=np.zeros(Nstage)
    for k in range(Nstage):
        tmp=cmath.phase(mean_phase[k])/(2*3.14159)
        if tmp>=0: 
            phases[k]=tmp
        else:
            phases[k]=tmp+1.0
        
    return phases


def gen_fig2Bv5(stagename,period,origin,sbin):
    tmp=get_stage_duration(stagename)
    Twin=(tmp/sbin).astype(int)
    
    rateCA1=np.load("rateCA1_all_normal.npy")
    rateCA1PV=np.load("rateCA1PV_all_normal.npy")
    rateCA3=np.load("rateCA3_all_normal.npy")
    rateCA1VIP=np.load("rateCA1VIP_all_normal.npy")
    rateCA1OLM=np.load("rateCA1OLM_all_normal.npy")
    thetaCA1     =np.mean(get_theta_phase_rate(rateCA1,Twin,[period,origin]),axis=0) # averaged over networks
    thetaCA1PV =np.mean(get_theta_phase_rate(rateCA1PV,Twin,[period,origin]),axis=0)
    thetaCA3     =np.mean(get_theta_phase_rate(rateCA3,Twin,[period,origin]),axis=0)
    thetaCA1VIP=np.mean(get_theta_phase_rate(rateCA1VIP,Twin,[period,origin]),axis=0)
    thetaCA1OLM =np.mean(get_theta_phase_rate(rateCA1OLM,Twin,[period,origin]),axis=0)

    pl.figure()
    pl.subplot(5,1,1)
    pl.bar(range(period*2),thetaCA3)
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*np.max(thetaCA3) )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))

    pl.subplot(5,1,2)
    pl.bar(range(period*2),thetaCA1)
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*2.5 )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))

    pl.subplot(5,1,3)
    pl.bar(range(period*2),thetaCA1PV)
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*60 )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))
    
    pl.subplot(5,1,4)
    pl.bar(range(period*2),thetaCA1OLM)
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*60 )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))


    pl.subplot(5,1,5)
    if stagename=="s-C":
        pl.bar(range(period*2),0.2*np.array(thetaCA1VIP))  #0.2=[Ach]
    else:
        pl.bar(range(period*2),np.array(thetaCA1VIP))  #0.2=[Ach]
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*30 )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))

def gen_fig2Cv5(stagename,sbin):
    tmp=get_stage_duration(stagename)
    Twin=(tmp/sbin).astype(int)
    

    rateCA1=np.load("rateCA1_all_normal.npy")
    rateCA3=np.load("rateCA3_all_normal.npy")
    rateEC3=np.load("rateEC3_all_normal.npy")
    N=len(rateCA1)  # N=Nnet*Ninit
    
    pl.figure()
    pl.subplot(3,1,1)
    tmprate=np.zeros((N,4))
    for i in range(4):
        tmprate[:,i]=np.mean(rateCA3[:,i*120:120*(i+1),Twin[0]:Twin[1]],axis=(1,2))
    tmpmean,tmperr=np.mean(tmprate,axis=0),np.std(tmprate,axis=0)
    pl.bar(range(4),tmpmean)
    pl.errorbar(range(4),tmpmean,yerr=tmperr,c="r")
    pl.xticks(range(4),["R","L","C","H"])

    pl.subplot(3,1,2)
    tmprate=np.zeros((N,4))
    for i in range(4):
        tmprate[:,i]=np.mean(rateCA1[:,i*120:120*(i+1),Twin[0]:Twin[1]],axis=(1,2))
    tmpmean,tmperr=np.mean(tmprate,axis=0),np.std(tmprate,axis=0)
    pl.bar(range(4),tmpmean)
    pl.errorbar(range(4),tmpmean,yerr=tmperr,c="r")
    pl.ylim(0,7)
    pl.xticks(range(4),["R","L","C","H"])
    
    pl.subplot(3,1,3)
    tmprate=np.zeros((N,4))
    for i in range(4):
        tmprate[:,i]=np.mean(rateEC3[:,i*120:120*(i+1),Twin[0]:Twin[1]],axis=(1,2))
    tmpmean,tmperr=np.mean(tmprate,axis=0),np.std(tmprate,axis=0)
    pl.bar(range(4),tmpmean)
    pl.errorbar(range(4),tmpmean,yerr=tmperr,c="r")
    #pl.ylim(0,8)
    pl.xticks(range(4),["R","L","C","H"])
    
    
def gen_fig2Bv4(stagename,period,origin,sbin):
    tmp=get_stage_duration(stagename)
    Twin=(tmp/sbin).astype(int)
    
    rateCA1=np.load("rateCA1_all_normal.npy")
    rateCA1PV=np.load("rateCA1PV_all_normal.npy")
    rateCA3=np.load("rateCA3_all_normal.npy")
    rateCA1VIP=np.load("rateCA1VIP_all_normal.npy")
    thetaCA1     =np.mean(get_theta_phase_rate(rateCA1,Twin,[period,origin]),axis=0) # averaged over networks
    thetaCA1PV =np.mean(get_theta_phase_rate(rateCA1PV,Twin,[period,origin]),axis=0)
    thetaCA3     =np.mean(get_theta_phase_rate(rateCA3,Twin,[period,origin]),axis=0)
    thetaCA1VIP=np.mean(get_theta_phase_rate(rateCA1VIP,Twin,[period,origin]),axis=0)

    pl.figure()
    pl.subplot(4,1,1)
    pl.bar(range(period*2),thetaCA3)
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*np.max(thetaCA3) )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))

    pl.subplot(4,1,2)
    pl.bar(range(period*2),thetaCA1)
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*2.5 )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))

    pl.subplot(4,1,3)
    pl.bar(range(period*2),thetaCA1PV)
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*60 )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))


    pl.subplot(4,1,4)
    if stagename=="s-C":
        pl.bar(range(period*2),0.2*np.array(thetaCA1VIP))  #0.2=[Ach]
    else:
        pl.bar(range(period*2),np.array(thetaCA1VIP))  #0.2=[Ach]
    #trough = 0 in LFP
    pl.plot(np.arange(0,period*2,0.1),(np.sin(2*3.14159*(np.arange(0,period*2,0.1)/float(period) - 0.25))+1)*0.5*30 )
    pl.xticks(np.arange(0,period*2+1,5),np.arange(0,720+1,180))

    
def gen_fig2Cv4(stagename):
    tmp=get_stage_duration(stagename)
    Twin=(tmp/sbin).astype(int)
    N=len(rateCA1)  # N=Nnet*Ninit

    rateCA1=np.load("rateCA1_all_normal.npy")
    rateCA3=np.load("rateCA3_all_normal.npy")
    
    pl.figure()
    pl.subplot(2,1,1)
    tmprate=np.zeros((N,4))
    for i in range(4):
        tmprate[:,i]=np.mean(rateCA3[:,i*120:120*(i+1),Twin[0]:Twin[1]],axis=(1,2))
    tmpmean,tmperr=np.mean(tmprate,axis=0),np.std(tmprate,axis=0)
    pl.bar(range(4),tmpmean)
    pl.errorbar(range(4),tmpmean,yerr=tmperr,c="r")
    pl.xticks(range(4),["R","L","C","H"])

    pl.subplot(2,1,2)
    tmprate=np.zeros((N,4))
    for i in range(4):
        tmprate[:,i]=np.mean(rateCA1[:,i*120:120*(i+1),Twin[0]:Twin[1]],axis=(1,2))
    tmpmean,tmperr=np.mean(tmprate,axis=0),np.std(tmprate,axis=0)
    pl.bar(range(4),tmpmean)
    pl.errorbar(range(4),tmpmean,yerr=tmperr,c="r")
    pl.ylim(0,7)
    pl.xticks(range(4),["R","L","C","H"])
    
    
def plt_fig_Ca_dyn():
    
    dyndata=np.loadtxt("Data_for_paper/dyn_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.30_00.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma")
    
    pl.figure()
    pl.subplot(3,1,1)
    pl.plot(range(0,7500), dyndata[500:8000,410])
    pl.ylim(-100,-50)
    pl.subplot(3,1,2)
    pl.plot(range(0,7500), dyndata[500:8000,411])
    pl.plot(range(0,7500),np.ones(7500)*0.004)
    pl.plot(range(0,7500),np.ones(7500)*0.0003)
    pl.subplot(3,1,3)
    pl.plot(range(0,7500), dyndata[500:8000,412])
    pl.plot(range(0,7500),np.ones(7500)*2.5)
    pl.plot(range(0,7500),np.ones(7500)*0.3)
    
    return
    