coind_clust=[]

tinit=7000
thrs =20
ntbin=40
for inet in range(5):
    for iinit in range(5):

        rast=np.loadtxt("rast_8080gamma28%d%d0.10.0025.dat" % (inet,iinit))
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
coind_clust
coind_clust=np.array(coind_clust)
coind_clust.shape
pl.scatter(coind_clust[:,3],coind_clust[:,5])
for i in range(len(coind_clust[:,0])):
    if coind_clust[i,0]>20:
        pl.scatter(coind_clust[i,3],coind_clust[i,5],c="b",lw=0)
    else
for i in range(len(coind_clust[:,0])):
    if coind_clust[i,0]>20:
        pl.scatter(coind_clust[i,3],coind_clust[i,5],c="b",lw=0)
    else:
        pl.scatter(coind_clust[i,3],coind_clust[i,5],c="r",lw=0)
pl.subplot(1,3,1)
for i in range(len(coind_clust[:,0])):
    if coind_clust[i,0]>20:
        pl.scatter(coind_clust[i,3],coind_clust[i,5],c="b",lw=0)
    else:
        pl.scatter(coind_clust[i,3],coind_clust[i,5],c="r",lw=0)
pl.subplot(1,3,2)
for i in range(len(coind_clust[:,0])):
    if coind_clust[i,0]>20:
        pl.scatter(coind_clust[i,4],coind_clust[i,5],c="b",lw=0)
    else:
        pl.scatter(coind_clust[i,4],coind_clust[i,5],c="r",lw=0)
pl.subplot(1,3,3)
for i in range(len(coind_clust[:,0])):
    if coind_clust[i,0]>20:
        pl.scatter(coind_clust[i,3]+coind_clust[i,4],coind_clust[i,5],c="b",lw=0)
    else:
        pl.scatter(coind_clust[i,3]+coind_clust[i,4],coind_clust[i,5],c="r",lw=0)
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
        pl.scatter(coind_clust[i,3]+coind_clust[i,4],coind_clust[i,5],c="b",lw=0)
    else:
        pl.scatter(np.log(coind_clust[i,3])+np.log(coind_clust[i,4]),coind_clust[i,5],c="r",lw=0)
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
