import sys,os
import math
import numpy as np

def dcoptimize(dist):
  n = len(dist)
  m = 0
  allds = []
  for i in range(n):
    for j in range(i+1,n):
      allds.append(dist[i][j])
      m += 1

  allds.sort()
  maxd = allds[m-1]

#calculate lower and upper boundaries for sigma
  counter = 0
  while allds[counter]==0:
    counter += 1

  lower =  allds[counter]
  if lower<np.percentile(np.array(allds),1):
    lower = np.percentile(np.array(allds),1)
  upper = np.percentile(np.array(allds),5)
  if upper<lower:
    upper = lower + 0.01

  hmin = 100000000000
  dc = 0
  sigma = lower

  h = 0
  while(sigma<=upper):
    
    #calculate potentials
    potentials = []
    z = 0.0
    for i in range(n):
      pt = 0.0
      for j in range(n):
        sq = (dist[i][j] / sigma) * (dist[i][j] / sigma)
        expsq = math.exp(-sq)
        if i!=j:
          pt += expsq
      potentials.append(pt)
      z += pt


    #calculate entropy
    h = 0.0
    for i in range(n):
      a = potentials[i]/z
      if a>0.0:
        b = math.log(a)
        h += a*b
    h = h*(-1)
    #store dc and hmin
    if h<hmin:
      hmin = h
      dc = sigma

    sigma += 0.005
        
        #print "dc = " + str(dc)

  return ((3/math.sqrt(2))*dc),maxd



#dist matrix and wether to include halo or not in clusters
def dclust(dist,haloflag,perc,dc):
  n = len(dist)
#Initialize structures
  rho = [0.0 for k in range(n)]
  delta = [0.0 for k in range(n)]
  neigh = [0 for k in range(n)]
  cl = [-1 for k in range(n)]
  halo = [0 for k in range(n)]
  strhoaux = [(0,0.0) for k in range(n)]

#calculate dc and maxdistance
  dc2,maxd = dcoptimize(dist)
  if dc==0:
    dc = dc2
#Calculate rho using gausian kernel
  for i in range(n):
    for j in range(i+1,n):
      sq = -(dist[i][j]/dc)*(dist[i][j]/dc)
      rho[i] += math.exp(sq)
      rho[j] += math.exp(sq)

#Sort rho with keys descentdently
  for i in range(n):
    strhoaux[i] = (i,rho[i])

  strho = sorted(strhoaux,key = lambda val:-val[1])

#Calculate delta

  for i in range(n):
    idxi,rhoidxi = strho[i]
    delta[idxi] = maxd
    for j in range(i):
      idxj,rhoidxj = strho[j]
      if delta[idxi]>dist[idxi][idxj]:
        delta[idxi] = dist[idxi][idxj]
        neigh[idxi] = idxj

  sortrho = sorted(rho)
  sortdelta = sorted(delta)

  rhothreshold = np.percentile(np.array(sortrho),perc)
  deltathreshold = np.percentile(np.array(sortdelta),75)


#Identification of cluster centers
  nclust = 0
  for i in range(n):
    if rho[i]>rhothreshold and delta[i]>deltathreshold:
      cl[i] = nclust
      #print i,nclust
      nclust+=1

#Assign clusters
  for i in range(n):
    idxi,rhoidxi = strho[i]
    if cl[idxi]==-1:
      cl[idxi] = cl[neigh[idxi]]

#Find border densities
  bord_rho = []
  for i in range(nclust):
    bord_rho.append(0)
    for j in range(n):
      if cl[j]==i:
        for k in range(n):
          if cl[k]!=cl[j] and dist[j][k]<=dc and rho[j]>bord_rho[i]:
            bord_rho[i] = rho[j]

#Generate Halo
  sum = 0
  #print n
  for i in range(n):
      #print i
    if rho[i]<bord_rho[cl[i]]:
      halo[i] = 1
    elif cl[i]>=0:
      sum += 1

  if haloflag==0:
    for i in range(n):
      if halo[i]==1:
        cl[i]=-1



  return nclust,cl









