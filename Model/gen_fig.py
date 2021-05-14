import numpy as np
import pylab as pl
import func as fc
import itertools

"""
        
NE_EC3=240   0-240
NE_EC=160    240-400

NE_EC5=400   405-805
NI_EC5=120   805-925

NE_CA1=480   930-1410
NE_PV =240  1410-1650
NE_OLM=240  1650-1890

NVIP=240    1895-2135
Nos=20      2135-2155
NCA1=480    2155-2635

"""
NE_EC3=240
T=8000
Nnet=5
Ninit=5

def idx_cell_eacharea(area_name):
    if area_name=="EC3":
        (icell_s,icell_e)=(0,240)
    elif area_name=="EC5":
        (icell_s,icell_e)=(405,805)
    elif area_name=="CA1":
        (icell_s,icell_e)=(930,1410)

    return (icell_s,icell_e)



def main(fig):
    if fig=="rate_dyn":  #2018.6.26
        fname="rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.34_30.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma"
        plt_rate_dyn(fname)
    elif fig=="Ach_sched":
        cogn_func,ioriginal="delay",-2
        theta_pref_all,rateEC3_all,rateCA1_all=gen_theta_rate(cogn_func)

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

#2018/7/4
def gen_theta_rate(cogn_func):
    sbin=10
    period=10  # theta period by unit
    stage_duration=[[500,1500],[1500,2500],[2500,3500],[3500,4250],[4250,5000],[5000,6500],[6500,8000]]
    theta_pref_all=[]
    rateEC3_all=[]
    rateCA1_all=[]
    
    
    if cogn_func=="enc":
        tmp_iter=itertools.product([0.7,0.9,1.0,1.1],[0.8],[1.1])
    elif cogn_func=="delay":
        tmp_iter=itertools.product([1.0],[0.2,0.5,0.7,0.8,0.9],[1.1])
    elif cogn_func=="rcl":
        tmp_iter=itertools.product([1.0],[0.8],[0.9,1.0,1.1,1.2])
    
    for enc_ach,delay_ach,rcl_ach in tmp_iter:
        if (enc_ach,delay_ach,rcl_ach,)==(1.0,0.8,1.1):
            isnormal=True
        else:
            isnormal=False
            
        ADR="Data_for_paper/"
        fname1="rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_"
        if isnormal:
            fname2="cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3"
            fname3="0.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma"
        else:
            fname2="cCA1EC50.005_EC51e-05_gCAN0.02GEI_EC50.2GIE_EC50.3"
            fname3="0.015.datdec4long30200.0150.0180.018080gamma30_nogamma"
    
        theta_pref=[]
        rateEC3=np.zeros((len(stage_duration),Nnet,Ninit,2))
        rateCA1=np.zeros((len(stage_duration),Nnet,Ninit,4))
    
        for i in range(len(stage_duration)):
            theta_pref.append([])
            for j in range(period):
                theta_pref[i].append(0)
        
        for inet,init in itertools.product(range(Nnet),range(Ninit)):
            if isnormal:
                fname=ADR+fname1+fname2+"%d_%d" %(inet,init)+fname3
            else:
                fname=ADR+fname1+fname2+"%d_%d" %(inet,init)+fname3+"_examAch%g%g%g" % (delay_ach,rcl_ach,enc_ach)
            tmp=np.loadtxt(fname)
            rate =fc.get_rate(tmp,[idx_cell_eacharea("CA1"),int(8000/sbin),1.0/float(sbin)])
            rate1=fc.get_rate(tmp,[idx_cell_eacharea("EC3"),int(8000/sbin),1.0/float(sbin)])
            for idx_stage in range(len(stage_duration)):
                tstage=stage_duration[idx_stage]
                tmpshift=int(tstage[0]/sbin)%period
                for i in range(period):
                    theta_pref[idx_stage][(i+tmpshift)%period]+=np.sum(rate[:,int(tstage[0]/sbin)+i:int(tstage[1]/sbin):period])
                for i in range(4):        
                    rateCA1[idx_stage,inet,init,i]=np.sum(rate[i*120:(i+1)*120,int(tstage[0]/sbin):int(tstage[1]/sbin)])
                for i in range(2):        
                    rateEC3[idx_stage,inet,init,i]=np.sum(rate1[i*120:(i+1)*120,int(tstage[0]/sbin):int(tstage[1]/sbin)])
        theta_pref_all.append(theta_pref)
        rateEC3_all.append(rateEC3)
        rateCA1_all.append(rateCA1)

    return theta_pref_all,rateEC3_all,rateCA1_all

def gen_theta_lock_index(fname):
    sbin=50  #[ms]
    t_rate=int(T/sbin)
    inv_sbin=1.0/float(sbin)
    stage_duration=[[10,30],[30,50],[50,70],[70,85],[85,100],[100,130],[130,160]]  # 1 unit = sbin
    
    tmp=np.loadtxt(fname)
    rate=fc.get_rate_EC(tmp,[NE_EC3,t_rate,inv_sbin])
    
    theta_peak=[]
    for tstage in stage_duration:
        tmp_rate=np.mean(rate[:240,tstage[0]:tstage[1]],axis=0)  # [Hz / cell]
        cor=np.correlate(tmp_rate,tmp_rate,"full")
        theta_peak.append(np.max(cor[int(len(cor)/2)+1:int(len(cor)/2)+4])/float((tstage[1]-tstage[0])*sbin))

    return theta_peak  


##2018.1.31
# last update  2018.6.26
def plt_rate_dyn(fname): 
    sbin=50
    inv_sbin=1.0/float(sbin)
    nbin=int(T/sbin)
    
    tmp=np.loadtxt(fname)
    rate=fc.get_rate_EC(tmp,[NE_EC3,nbin,inv_sbin])
    
    pl.plot(range(nbin+1),np.mean(rate[:int(NE_EC3/2),:],axis=0))
    pl.plot(range(nbin+1),np.mean(rate[int(NE_EC3/2):,:],axis=0))
    pl.xlim(10,160)
#    pl.savefig("rate_dyn_EC3_normal.eps")

    return

def fig_blk(): # 2018.2.1

    rate_all=[]
    for inet in range(5):
        rate_all.append([])
        for iinit in range(5):
            tmp=np.loadtxt("rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma306000_8000"
                           % (inet,iinit))
            rate=fc.get_rate_EC(tmp,[240,160,0.02])
            rate_all[inet].append(rate)

    rate_blk3=[]
    for inet in range(5):
        rate_blk3.append([])
        for inet in [0,2,3,4]:
            for iinit in range(5):
                tmp=np.loadtxt("rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma305000_6500" % (inet,iinit))
                rate=fc.get_rate_EC(tmp,[240,160,0.02])
                rate_blk3[inet].append(rate)

    rate_blk1=[]
    for inet in range(5):
        rate_blk1.append([])
        for inet in [1,2,3,4]:
            for iinit in range(5):
                tmp=np.loadtxt("rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma301500_3000" % (inet,iinit))
                rate=fc.get_rate_EC(tmp,[240,160,0.02])
                rate_blk1[inet].append(rate)
    
   
    LR_rate=np.zeros((3,Nnet,Ninit,2,80))
    for inet in [1,2,3,4]:
        tmp=np.zeros((Ninit,2,160))
        for iinit in range(Ninit):
            tmp[iinit,0]=np.mean(rate_blk1[inet][iinit][:120,:160],axis=0)
            tmp[iinit,1]=np.mean(rate_blk1[inet][iinit][120:,:160],axis=0)
            for t in range(80):
                LR_rate[0,inet,iinit,:,t]=np.sum(tmp[iinit,:,2*t:2*(t+1)],axis=1)
    for inet in [0,2,3,4]:
        tmp=np.zeros((Ninit,2,160))
        for iinit in range(Ninit):
            tmp[iinit,0]=np.mean(rate_blk3[inet][iinit][:120,:160],axis=0)
            tmp[iinit,1]=np.mean(rate_blk3[inet][iinit][120:,:160],axis=0)
            for t in range(80):
                LR_rate[1,inet,iinit,:,t]=np.sum(tmp[iinit,:,2*t:2*(t+1)],axis=1)
    for inet in [0,1,2,3,4]:
        tmp=np.zeros((Ninit,2,160))
        for iinit in range(Ninit):
            tmp[iinit,0]=np.mean(rate_all[inet][iinit][:120,:160],axis=0)
            tmp[iinit,1]=np.mean(rate_all[inet][iinit][120:,:160],axis=0)
            for t in range(80):
                LR_rate[2,inet,iinit,:,t]=np.sum(tmp[iinit,:,2*t:2*(t+1)],axis=1)
    
        for i in range(3):
            pl.plot(range(80),np.mean(np.exp(LR_rate[i,:,:,0,:])/(np.exp(LR_rate[i,:,:,0,:])+np.exp(LR_rate[i,:,:,1,:])),axis=(0,1)),label="blk%d" % i)
    pl.legend()
    pl.savefig("probL_dyn_blks.eps")
    
    
    
    
    
    rate_blk2=np.copy(rate_all)
      
    # normal condition
    rate_all=[]
    for inet in range(5):
        rate_all.append([])
        for iinit in range(5):
            tmp=np.loadtxt("rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma" % (inet,iinit))
            rate=fc.get_rate_EC(tmp,[240,160,0.02])
            rate_all[inet].append(rate)
    
        # probability of Left is infered
    
    LR_rate=np.zeros((4,Nnet,Ninit,2,80))
    for inet in [1,2,3,4]:
        tmp=np.zeros((Ninit,2,160))
        for iinit in range(Ninit):
            tmp[iinit,0]=np.mean(rate_blk1[inet][iinit][:120,:160],axis=0)
            tmp[iinit,1]=np.mean(rate_blk1[inet][iinit][120:,:160],axis=0)
            for t in range(80):
                LR_rate[0,inet,iinit,:,t]=np.sum(tmp[iinit,:,2*t:2*(t+1)],axis=1)
    for inet in [0,2,3,4]:
        tmp=np.zeros((Ninit,2,160))
        for iinit in range(Ninit):
            tmp[iinit,0]=np.mean(rate_blk3[inet][iinit][:120,:160],axis=0)
            tmp[iinit,1]=np.mean(rate_blk3[inet][iinit][120:,:160],axis=0)
            for t in range(80):
                LR_rate[1,inet,iinit,:,t]=np.sum(tmp[iinit,:,2*t:2*(t+1)],axis=1)
    for inet in [0,1,2,3,4]:
        tmp=np.zeros((Ninit,2,160))
        for iinit in range(Ninit):
            tmp[iinit,0]=np.mean(rate_blk2[inet][iinit][:120,:160],axis=0)
            tmp[iinit,1]=np.mean(rate_blk2[inet][iinit][120:,:160],axis=0)
            for t in range(80):
                LR_rate[2,inet,iinit,:,t]=np.sum(tmp[iinit,:,2*t:2*(t+1)],axis=1)
    for inet in [0,1,2,3,4]:
        tmp=np.zeros((Ninit,2,160))
        for iinit in range(Ninit):
            tmp[iinit,0]=np.mean(rate_all[inet][iinit][:120,:160],axis=0)
            tmp[iinit,1]=np.mean(rate_all[inet][iinit][120:,:160],axis=0)
            for t in range(80):
                LR_rate[3,inet,iinit,:,t]=np.sum(tmp[iinit,:,2*t:2*(t+1)],axis=1)
    
                
    tmp=np.zeros(4)
    tmpstd=np.zeros(4)
    for i in range(4):
        tmp[i]=np.mean(np.exp(3*np.mean(LR_rate[i,2:,:,0,70:],axis=2))/(np.exp(3*np.mean(LR_rate[i,2:,:,0,70:],axis=2))+np.exp(3*np.mean(LR_rate[i,2:,:,1,70:],axis=2))))
        tmpstd[i]=np.std(np.exp(3*np.mean(LR_rate[i,2:,:,0,70:],axis=2))/(np.exp(3*np.mean(LR_rate[i,2:,:,0,70:],axis=2))+np.exp(3*np.mean(LR_rate[i,2:,:,1,70:],axis=2))))
    
        
        pl.bar(range(4),tmp)
    pl.errorbar(np.arange(0.5,4.5,1.0),tmp,yerr=tmpstd,c="r")
    pl.savefig("probL_all_blk.eps")
    
    
def fig_performance():
    tmp=np.zeros((Nnet,Ninit,2))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            tmp[inet,iinit,0]=np.mean(rate_norm[inet,iinit,:120,700:800])
            tmp[inet,iinit,1]=np.mean(rate_norm[inet,iinit,120:240,700:800])
            
    for inet in range(Nnet):
        for iinit in range(Ninit):
            pl.plot(np.arange(0.4,1.5,1),tmp[inet,iinit,:],c="gray")
            pl.scatter(np.arange(0.4,1.5,1),tmp[inet,iinit,:],s=10,c="gray")
    pl.bar(range(2),np.mean(tmp,axis=(0,1)))
    pl.savefig("perf_norm.eps")
    
    # dynamics
    fname="rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma306000_8000"
    rate_blk4=fc.get_rate_all(fname,5,5)
    
    fname="rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma301500_3500"
    rate_blk1=fc.get_rate_all(fname,5,5)
    
    rate_blk1=np.array(rate_blk1)
    rate_blk4=np.array(rate_blk4)
    tmp=np.zeros((3,Nnet,Ninit,16,2))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            for it in range(16):
                tmp[0,inet,iinit,it,0] =np.mean(rate_norm[inet,iinit,:120,it*50:(it+1)*50])
                tmp[0,inet,iinit,it,1]=np.mean(rate_norm[inet,iinit,120:240,it*50:(it+1)*50])
                tmp[1,inet,iinit,it,0] =np.mean(rate_blk1[inet,iinit,:120,it*50:(it+1)*50])
                tmp[1,inet,iinit,it,1]=np.mean(rate_blk1[inet,iinit,120:240,it*50:(it+1)*50])
                tmp[2,inet,iinit,it,0] =np.mean(rate_blk4[inet,iinit,:120,it*50:(it+1)*50])
                tmp[2,inet,iinit,it,1]=np.mean(rate_blk4[inet,iinit,120:240,it*50:(it+1)*50])
                
    pl.bar(range(0,45,3),np.mean(tmp[1,:,:,1:],axis=(0,1,3)),color="b")
    pl.bar(range(1,45,3),np.mean(tmp[2,:,:,1:],axis=(0,1,3)),color="g")
    pl.savefig("theta_act_ECIII_blk1_4.eps")


def fig_disinh():
    rateall=[]
    for inet in range(Nnet):
        rateall.append([])
        for iinit in range(Ninit):
            tmp=np.loadtxt("rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma30_PVinact" % (inet,iinit))
        rateall[inet].append(fc.get_rate_EC(tmp,[2635,801,0.1]))
    rate_noPV=np.copy(rateall)

    rateall=[]
    for inet in range(Nnet):
        rateall.append([])
        for iinit in range(Ninit):
            tmp=np.loadtxt("rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma" % (inet,iinit))
            rateall[inet].append(fc.get_rate_EC(tmp,[2635,801,0.1]))
    rate_normal=np.copy(rateall)

    rateall=[]
    for inet in range(Nnet):
        rateall.append([])
        for iinit in range(Ninit):
            tmp=np.loadtxt("rast_EC5EC30.012_CA1EC50.0008_EC3_CA10.001_Gpv0_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.015.datdec4long1.130200.0150.0180.018080gamma30_nogamma" % (inet,iinit))
            rateall[inet].append(fc.get_rate_EC(tmp,[2635,801,0.1]))
    rate_disdis=np.copy(rateall)

    
    stageT=[[500,1500],[1500,3500],[3500,5000],[5000,8000]]
    idxEcell=[[[0,120],[120,240]],[[405,605],[605,805]],[[930,1050],[1050,1170],[1170,1290],[1290,1410],[1410,1650]],[[2155,2275],[2275,2395],[2395,2515],[2515,2635]]]
    phaseall=[]
    for inet in range(Nnet):
        phaseall.append([])
        for iinit in range(Ninit):
            ratetmp=rate_noPV[inet][iinit]
            phase_distr=[]
            for istage in range(len(stageT)):
                phase_distr.append([])
                listt=np.arange(stageT[istage][0],stageT[istage][1],10)*0.1
                for igroup in range(len(idxEcell)):
                    phase_distr[istage].append([])
                    for icell in range(len(idxEcell[igroup])):
                        tmpphase=np.zeros(10)
                        for t in listt:
                            iphase=int(t)%10
                            tmpphase[iphase]+=np.mean(ratetmp[idxEcell[igroup][icell][0]:idxEcell[igroup][icell][1],int(t)])
                        tmpphase=tmpphase/float(len(listt))
                        phase_distr[istage][igroup].append(tmpphase)
            phaseall[inet].append(phase_distr)
    phase_noPV=np.copy(phaseall)


    ####  perhaps, the following is different function from the just upper one
    istage=0
    phaseall=np.copy(phase_normal)
    igroup=0   # ECIII
    tmpphase=np.zeros((2,20))  
    for inet in range(Nnet):
        for iinit in range(Ninit):
            tmpphase[0,:10]+=phaseall[inet][iinit][istage][igroup][0]
            tmpphase[1,:10]+=phaseall[inet][iinit][istage][igroup][1]
    tmpphase[:,10:]=tmpphase[:,:10]
    tmpphase=tmpphase/float(Nnet*Ninit)
    phaseEC3=np.copy(tmpphase)
    
    igroup=3   # CA3
    tmpphase=np.zeros((4,20))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            tmpphase[0,:10]+=phaseall[inet][iinit][istage][igroup][0]
            tmpphase[1,:10]+=phaseall[inet][iinit][istage][igroup][1]
            tmpphase[2,:10]+=phaseall[inet][iinit][istage][igroup][2]
            tmpphase[3,:10]+=phaseall[inet][iinit][istage][igroup][3]
    tmpphase[:,10:]=tmpphase[:,:10]
    tmpphase=tmpphase/float(Nnet*Ninit)
    phaseCA3=np.copy(tmpphase)
    
    igroup=2  #CA1
    tmpphase=np.zeros((5,20))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            tmpphase[0,:10]+=phaseall[inet][iinit][istage][igroup][0]
            tmpphase[1,:10]+=phaseall[inet][iinit][istage][igroup][1]
            tmpphase[2,:10]+=phaseall[inet][iinit][istage][igroup][2]
            tmpphase[3,:10]+=phaseall[inet][iinit][istage][igroup][3]
            tmpphase[4,:10]+=phaseall[inet][iinit][istage][igroup][4]
    tmpphase[:,10:]=tmpphase[:,:10]
    tmpphase=tmpphase/float(Nnet*Ninit)
    phaseCA1=np.copy(tmpphase)
    
    
    pl.subplot(3,1,1)
    pl.bar(range(20),np.mean(phaseEC3,axis=0)*100,lw=0,color="b",label="EC3")
    pl.bar(range(20),np.mean(phaseCA3,axis=0)*100,color="r",lw=0,alpha=0.3,label="CA3")
    amp=np.max(np.mean(phaseCA1,axis=0))*80/2.0
    pl.plot(np.arange(0,20,0.1),amp*(-1*np.sin(2*3.14159*0.1*np.arange(0,20,0.1))+1))
    pl.legend()
    pl.subplot(3,1,2)
    pl.bar(range(20),np.mean(phaseCA1[:4],axis=0)*100)
    amp=np.max(np.mean(phaseCA1[:4],axis=0))*100*0.8/2.0
    pl.plot(np.arange(0,20,0.1),amp*(-1*np.sin(2*3.14159*0.1*np.arange(0,20,0.1))+1))
    pl.legend()
    pl.subplot(3,1,3)
    pl.bar(range(20),phaseCA1[4]*100,label="PV")
    amp=np.max(phaseCA1[4])*80/2.0
    pl.plot(np.arange(0,20,0.1),amp*(-1*np.sin(2*3.14159*0.1*np.arange(0,20,0.1))+1))
    pl.legend()
    pl.savefig("phase_pref_samp_center.eps")
    
    istage=0
    igroup=2
    probtmp=np.zeros((Nnet,Ninit,4))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            for icell in range(4):
                probtmp[inet,iinit,icell]=np.mean(rate_normal[inet][iinit][idxEcell[igroup][icell][0]:idxEcell[igroup][icell][1],int(stageT[istage][0]*0.1):int(stageT[istage][1]*0.1)])
    probtmp=probtmp*100  # to be [Hz]
    probtmp1=np.copy(probtmp)
    
    igroup=3
    probtmp=np.zeros((Nnet,Ninit,4))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            for icell in range(4):
                probtmp[inet,iinit,icell]=np.mean(rate_normal[inet][iinit][idxEcell[igroup][icell][0]:idxEcell[igroup][icell][1],int(stageT[istage][0]*0.1):int(stageT[istage][1]*0.1)])
    probtmp=probtmp*100  # to be [Hz]
    
    pl.subplot(2,2,1)
    pl.bar(range(4),np.mean(probtmp[:,:],axis=(0,1)))
    pl.errorbar(np.arange(4)+0.3,np.mean(probtmp[:,:],axis=(0,1)),yerr=np.std(probtmp[:,:],axis=(0,1)),color="r")
    a=np.zeros((Nnet,Ninit,4))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            a[inet,iinit]=np.exp(probtmp[inet,iinit])/np.sum(np.exp(probtmp[inet,iinit]))
    
    pl.subplot(2,2,2)
    pl.errorbar(np.arange(4)+0.3,np.mean(a[:],axis=(0,1)),yerr=np.std(a[:],axis=(0,1)),color="green")
    
    pl.subplot(2,2,3)
    pl.bar(range(4),np.mean(probtmp1[:,:],axis=(0,1)))
    pl.errorbar(np.arange(4)+0.3,np.mean(probtmp1[:,:],axis=(0,1)),yerr=np.std(probtmp1[:,:],axis=(0,1)),color="r")
    pl.ylim(0,7)
    
    pl.subplot(2,2,4)
    a=np.zeros((Nnet,Ninit,4))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            a[inet,iinit]=np.exp(probtmp1[inet,iinit])/np.sum(np.exp(probtmp1[inet,iinit]))
    pl.errorbar(np.arange(4)+0.3,np.mean(a[:],axis=(0,1)),yerr=np.std(a[:],axis=(0,1)),color="r")
    pl.savefig("info_sample_center.eps")
    
    
    
    tmp=np.zeros((3,Nnet,Ninit,2))
    for inet in range(Nnet):
        for iinit in range(Ninit):
            tmp[0,inet,iinit,0]=np.mean(rate_normal[inet][iinit][0:120,600:800])
            tmp[0,inet,iinit,1]=np.mean(rate_normal[inet][iinit][120:240,600:800])
            tmp[1,inet,iinit,0]=np.mean(rate_disdis[inet][iinit][0:120,600:800])
            tmp[1,inet,iinit,1]=np.mean(rate_disdis[inet][iinit][120:240,600:800])
            tmp[2,inet,iinit,0]=np.mean(rate_noPV[inet][iinit][0:120,600:800])
            tmp[2,inet,iinit,1]=np.mean(rate_noPV[inet][iinit][120:240,600:800])
            
    pl.subplot(1,3,1)
    pl.plot([.4,1.4],np.mean(tmp[0,:,:,:]*100,axis=1).T)
    pl.bar([0,1],np.mean(tmp[0,:,:,:]*100,axis=(0,1)))
    pl.ylim(0,15)
    pl.subplot(1,3,2)
    pl.plot([.4,1.4],np.mean(tmp[1,:,:,:]*100,axis=1).T)
    pl.bar([0,1],np.mean(tmp[1,:,:,:]*100,axis=(0,1)))
    pl.ylim(0,15)
    pl.subplot(1,3,3)
    pl.plot([.4,1.4],np.mean(tmp[2,:,:,:]*100,axis=1).T)
    pl.bar([0,1],np.mean(tmp[2,:,:,:]*100,axis=(0,1)))
    pl.ylim(0,15)
    pl.savefig("diff_EC3.eps")
