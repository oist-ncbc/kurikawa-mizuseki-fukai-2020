import numpy as np
import pylab as pl
from scipy import linalg as LA
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D

# for filtering
from scipy.signal import butter, lfilter
from scipy.signal import freqz

"""
NE_EC3 = 120 
N_EC3 = NE_EC3+60+5

NE_EC5 = 60
N_EC5  = N_EC3+ (NE_EC5+30)+5

NE_CA1 = 45 
NPV    = 45
NOLM   = 45
N_CA1  = N_EC5+ (NE_CA1+NPV+NOLM)+5

Nosc   =N_CA1+10+10
N_CA3  =Nosc + 15*3
"""

#update 2018.2.24
def get_rate_all(fname,Nnet,Ninit):

    rate_all=[]
    for inet in range(Nnet):
        rate_all.append([])
        for iinit in range(Ninit):
            tmp=np.loadtxt(fname % (inet,iinit))
            rate=get_rate_EC(tmp,[2135,801,0.1])
            rate_all[inet].append(rate)
    return rate_all


def get_rate_EC(data,args):
    ncell,T=args[0],args[1]
    inv_bin=args[2]

    rate=np.zeros((ncell,T+1))
    
    for i in range(len(data[:,0])):
        if(data[i,1]<ncell):
            rate[data[i,1],int(data[i,0]*inv_bin)]+=inv_bin*1000.0  # # of spikes per 1s, 1cell
            
    return rate




""""  following part is old ver.
def est_autocor(data):
    x=data.copy()
    n=len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    result = r/(variance*(np.arange(n, 0, -1)))
    return result



def tmp():
    rateon=np.zeros((4,300,5000))
    rateof=np.zeros((4,200,5000))
    
    for inet in range(4):
        data11=np.loadtxt("rast0.1_0.01_0.03_50 %d.dat6_off" %(inet+1))
        data10=np.loadtxt("rast-0.3_0.025_Aosc10_Aca315_%d.dat6" % (inet+11))
        rateon[inet]=get_rate_EC(data10,[300,5000,0.1])
        rateof[inet]=get_rate_EC(data11,[200,5000,0.1])

    crs_coron=np.zeros((4,3,30))
    crs_corof=np.zeros((4,3,30))

    for inet in range(4):
        crs_coron[inet][0]=est_autocor1(np.mean(rateon[inet,:80,0:200],axis=0),30)
        crs_corof[inet][0]=est_autocor1(np.mean(rateof[inet,:80,0:200],axis=0),30)
        crs_coron[inet][1]=est_autocor1(np.mean(rateon[inet,:80,200:800],axis=0),30)
        crs_corof[inet][1]=est_autocor1(np.mean(rateof[inet,:80,200:800],axis=0),30)
        crs_coron[inet][2]=est_autocor1(np.mean(rateon[inet,:80,400:],axis=0),30)
        crs_corof[inet][2]=est_autocor1(np.mean(rateof[inet,:80,400:],axis=0),30)
        
    
    return rateon,rateof,crs_coron,crs_corof


def est_autocor1(data,nbin):
    x=data.copy()
    n=len(x)
    x = x-x.mean()
    cor=np.zeros(nbin)
    var1=np.zeros((2,nbin))
    for i in range(nbin):
        cor[i]=0
        var1[0,i]=0
        var1[1,i]=0

        for j in range(n-i):
            cor[i]=cor[i]    +x[j]*x[j+i]
            var1[0,i]=var1[0,i]+x[j]*x[j]
            var1[1,i]=var1[1,i]+x[j+i]*x[j+i]

        cor[i]=cor[i]/np.sqrt(var1[0,i]*var1[1,i])
        
    return cor

def est_autocor2(data,nbin):
    x=data.copy()
    n=len(x)
    x = x-x.mean()
    stdx=np.std(x)
    x = x/stdx
    cor=np.zeros(nbin)

    for i in range(nbin):
        cor[i]=0

        for j in range(n-i):
            cor[i]=cor[i]    +x[j]*x[j+i]
        
    return cor



def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y



 # Sample rate and desired cutoff frequencies (in Hz).
#    fs = 250.0
#    lowcut = 5.0
#    highcut = 15.0

# b, a = butter_bandpass(lowcut, highcut, fs, order=order)
#        w, h = freqz(b, a, worN=2000)
#        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)







def get_rate_all(nsmp,inet):
    Nsmp=N_CA3
    cnt=0

    inv_bin_smp=0.1 #[ms-1]
    t_smp = 5000  # total time of 1 trial [ms]
    n_bin      =int(t_smp*inv_bin_smp)
    rate_tmp=np.zeros((nsmp,n_bin+1,Nsmp))

    str_tmp="rast-0.3_0.025_Aosc5_Aca315_%d.dat7" % inet
    tmp=np.loadtxt(str_tmp)

    for i in range(len(tmp[:,0])):
        if i>0 and tmp[i,0]<tmp[i-1,0]:
            cnt+=1
            if nsmp<=cnt:
                break
        if tmp[i,1]<Nsmp:
            rate_tmp[cnt,int(tmp[i,0]*inv_bin_smp),tmp[i,1]]+=1
        else:
            print tmp[i,1]

    return rate_tmp




def get_pic(rate):

    pl.subplot(3,2,1)
    pl.plot(range(101),np.mean(rate[0,:101,0:500]*10,axis=1))
    pl.plot(range(101),np.mean(rate[0,:101,500:600],axis=1))
    pl.plot(range(101),(np.sin(2*3.14159*np.arange(101)/10)+1)*0.5)
#    pl.ylim(0,2)
    pl.subplot(3,2,3)
    pl.plot(range(101),np.mean(rate[0,:101,600:1100]*10,axis=1))
    pl.plot(range(101),np.mean(rate[0,:101,1100:1200],axis=1))
    pl.plot(range(101),(np.sin(2*3.14159*np.arange(101)/10)+1)*0.5)
 #   pl.ylim(0,2)
    pl.subplot(3,2,5)
    pl.plot(range(101),np.mean(rate[0,:101,1200:1700]*10,axis=1))
    pl.plot(range(101),np.mean(rate[0,:101,1700:1800],axis=1))
    pl.plot(range(101),(np.sin(2*3.14159*np.arange(101)/10)+1)*0.5)
  #  pl.ylim(0,2)


    return


def tmp_fc(rate,trange,cellid):
    phase_lock=np.zeros(10)
    
    for i in range(trange[0],trange[1]):
        for j in range(10):
            phase_lock[0,j,0]+=np.mean(rate[0,i*10+j,cellid],axis=0)


#    for i in range(2):
#        for j in range(6):
#            phase_lock[i,:,j]=phase_lock[i,:,j]/np.max(phase_lock[i,:,j])
            
    return phase_lock

    
def get_pic1(rate,rate1,tmp,tmp1):


    pl.subplot(4,1,1)
    pl.plot(range(4000),np.mean(tmp[:4000,range(2,201,2)],axis=1),lw=3)
    pl.plot(range(4000),np.mean(tmp[:4000,range(202,401,2)],axis=1),lw=3)
    pl.plot(range(4000),np.mean(tmp1[:4000,range(2,201,2)],axis=1))
    pl.plot(range(4000),np.mean(tmp1[:4000,range(202,401,2)],axis=1))
    
    pl.subplot(4,1,2)
    pl.plot(range(120),np.mean(rate[0,:120,0:200],axis=1),lw=3)
    pl.plot(range(120),np.mean(rate[0,:120,200:400],axis=1),lw=3)
    pl.plot(range(120),np.mean(rate1[0,:120,0:200],axis=1))
    pl.plot(range(120),np.mean(rate1[0,:120,200:400],axis=1))

    pl.subplot(4,1,3)
    pl.plot(range(120),np.mean(rate[0,:120,500:580],axis=1),lw=3)
    pl.plot(range(120),np.mean(rate[0,:120,580:660],axis=1),lw=3)
    pl.plot(range(120),np.mean(rate1[0,:120,500:580],axis=1))
    pl.plot(range(120),np.mean(rate1[0,:120,580:660],axis=1))

    pl.subplot(4,1,4)
    pl.plot(range(4000),np.mean(tmp[:4000,range(1,201,2)],axis=1),lw=3)
    pl.plot(range(4000),np.mean(tmp[:4000,range(201,401,2)],axis=1),lw=3)
    pl.plot(range(4000),np.mean(tmp1[:4000,range(1,201,2)],axis=1))
    pl.plot(range(4000),np.mean(tmp1[:4000,range(201,401,2)],axis=1))

    
    pl.show()

    return

def get_pic2(tmp1,rate1):
    n_smp=len(tmp1);
    pl.subplot(4,1,1)
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(1,401,4)],axis=1),c="b")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(401,801,4)],axis=1),c="r")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(801,1121,4)],axis=1),c="g")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(1121,1201,4)],axis=1),c="magenta")
    pl.subplot(4,1,2)
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(2,401,4)],axis=1),c="b")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(402,801,4)],axis=1),c="r")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(802,1121,4)],axis=1),c="g")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(1122,1201,4)],axis=1),c="magenta")
    pl.subplot(4,1,3)
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(3,401,4)],axis=1),c="b")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(403,801,4)],axis=1),c="r")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(803,1121,4)],axis=1),c="g")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(1123,1201,4)],axis=1),c="magenta")
    pl.subplot(4,1,4)
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(4,401,4)],axis=1),c="b")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(404,801,4)],axis=1),c="r")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(804,1121,4)],axis=1),c="g")
    pl.plot(range(n_smp-1),np.mean(tmp1[:n_smp-1,range(1124,1201,4)],axis=1),c="magenta")

    #    pl.plot(range(6000),np.mean(tmp1[:6000,range(401,601,2)],axis=1),lw=3,c="g")
        


    return




def get_spike_data(filetype):
    nnet=5
    ninit=4
    nclst=4
    if filetype=="long":
        nT   =700+1
    else:
        nT   =550+1
    EC31_all=np.zeros((nnet,ninit,nclst,nT))
    EC51_all=np.zeros((nnet,ninit,nclst,nT))
    CA11_all=np.zeros((nnet,ninit,nclst,nT))
    for inet in range(nnet):
        EC31_tot=np.zeros((ninit,nclst,nT))
        EC51_tot=np.zeros((ninit,nclst,nT))
        CA11_tot=np.zeros((ninit,nclst,nT))
        for iinit in range(ninit):
            EC3=np.zeros((nclst,nT))
            EC5=np.zeros((nclst,nT))
            CA1=np.zeros((nclst,nT))
            if filetype=="normal":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.025.dat2_tmp3" % (inet,iinit),usecols=(0,1))
            elif filetype=="PVinact":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.025.dat2_PVinact" % (inet,iinit),usecols=(0,1))
            elif filetype=="disord":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.025.dat2_disord" % (inet,iinit),usecols=(0,1))
            elif filetype=="off":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.025.dat2_off" % (inet,iinit),usecols=(0,1))
            elif filetype=="noPVinh":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.025.dat2" % (inet,iinit),usecols=(0,1))
            elif filetype=="notheta1":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.0251.dat2_notheta" % (inet,iinit),usecols=(0,1))
            elif filetype=="notheta05":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.0250.5.dat2_notheta" % (inet,iinit),usecols=(0,1))
            elif filetype=="notheta02":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.0250.2.dat2_notheta" % (inet,iinit),usecols=(0,1))
            elif filetype=="notheta2":  #total activity is hold
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.0252.dat2_notheta" % (inet,iinit),usecols=(0,1))
            elif filetype=="long":
                rast=np.loadtxt("rast_EC5EC30.007_CA1EC50.001_EC3_CA10.001_Gpv0.012_cCA1EC50.005_EC51e-05_gCAN0.020.8GEI_EC50.2GIE_EC50.3%d_%d0.025.dat_long945" % (inet,iinit),usecols=(0,1))


            for i in rast:
                tmptime=int(i[0]/10.0)
                if i[1]<240:  # EC3
                    tmpidx=i[1]
                    if tmpidx<120:
                        EC3[0,tmptime]+=1
                    elif tmpidx<240:
                        EC3[1,tmptime]+=1
                elif i[1]>285+120 and i[1]<765+120:  #EC5
                    tmpidx=i[1]-(285+80)
                    if tmpidx<200:
                        EC5[0,tmptime]+=1
                    elif tmpidx<400:
                        EC5[1,tmptime]+=1
                elif i[1]>810+120 and i[1]<1050+120:
                    if tmptime>=700:
                        print iinit,tmptime
                    tmpidx=i[1]-(810+120)
                    if tmpidx<60:
                        CA1[0,tmptime]+=1
                    elif tmpidx<120:
                        CA1[1,tmptime]+=1
                    elif tmpidx<180:
                        CA1[2,tmptime]+=1
                    else:
                        CA1[3,tmptime]+=1
            EC31_tot[iinit,:,:]=EC3[:,:]
            EC51_tot[iinit,:,:]=EC5[:,:]
            CA11_tot[iinit,:,:]=CA1[:,:]
        EC31_all[inet]=EC31_tot
        EC51_all[inet]=EC51_tot
        CA11_all[inet]=CA11_tot

    return EC31_all,EC51_all,CA11_all
"""




"""
distr=np.zeros((4,12,20))
for i in range(12):
    for j in range(50,150):
        distr[0,i,j%10]+=rate[i,j]
    for j in range(200,300):
        distr[1,i,j%10]+=rate[i,j]
    for j in range(300,400):
        distr[2,i,j%10]+=rate[i,j]
    for j in range(450,550):
        distr[3,i,j%10]+=rate[i,j]
distr[:,:,10:]=distr[:,:,:10]

for i in range(50,551-10):
    for j in range(12):
        rate_ave[j,i]=np.sum(rate[j,i:i+10])
"""
