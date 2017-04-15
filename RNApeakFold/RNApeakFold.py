import sys,os
import math
import operator
import subprocess
from numpy import *

SEQ   = '/sequence'
UBOX  = 'ubox'

fn = sys.argv[1]

d = {}
tot = 0.0

def boltzmannBasePairProb(comm,rna):
  #compute Boltzmann base pairing probability
  #cmd  = 'echo %s | RNAfold -p' % rna
  #os.system(cmd) #produces dot.ps in same directory
  auxc = comm[1:].split("|")
  auxf = ""
  for i in range(len(auxc)-2):
    auxf += auxc[i]+"_"
  auxf = auxf[:-1]+".fa"
  f = open(auxf,'w')
  f.write('>'+auxf[:-3]+'\n')
  f.write(rna+'\n')
  f.close()
  subprocess.check_call('RNAfold -p < ' + auxf + ' > /dev/null',stdin=None, stdout=None, stderr=None, shell=True)
  file = open(auxf[:-3]+'_dp.ps')
  line = file.readline()
  while line:
    if line[:len(SEQ)] == SEQ:
      rna0 = file.readline().strip()[:-1]
      break
    else:
      line = file.readline()
  #assert( rna == rna0 )
  n    = len(rna)
  line = file.readline()
  P    = {}
  #initialize P[x][y]=0
  for x in range(1,n+1):
    P[x] = {}
    for y in range(1,n+1):
      P[x][y] = 0.0
  #new get real base pairing probabilities
  while line:
    line = line.strip()
    if line=='':
      line = file.readline()
      continue
    words = line.split()
    if len(words) == 4 and words[3]==UBOX:
      x = int(words[0])
      y = int(words[1])
      p = float(words[2])**2
      P[x][y] = p; P[y][x] = p #symmetrize probabilities
    line = file.readline()
  file.close()
  #os.sys("rm dot.ps")
  return P

def getCapitals(seq):
  nseq = ""
  j = -1
  k = -1
  for i in range(len(seq)):
    if seq[i].isupper():
      if j==-1:
        j=i
      nseq += seq[i]
    elif j>-1:
      k = i
      break

  #print j,k,nseq
  return j,k,nseq

def postprocess(bps):
  #Postprocessing
  n = len(bps)
  stacked = 0
  minstack = 2
  bulge = 0
  start = -1
  started = 0
  previous = -1
  pbps = []
  for k in range(n):
    pbps.append(bps[k])

  for i in range(n):
    #print "scanning bas " + str(i) + " bprs = " + str(pbps[i])
    if(pbps[i]>i):#//hairpin
      if(started==0):
        #print "stacked starts at " + str(i)
        start = i
        started = 1
        stacked = 1
        bulge = 0
      else:
        if(bulge==0):
          previous = i-1
        else:
          previous = i-2
        if(pbps[i] == pbps[previous]-1):#stacked
          stacked+=1
          bulge = 0
          #print "Stacking = " + str(stacked)
        else:
          if((bulge==0)and(pbps[i] == pbps[i-1]-2)):
            bulge = 1
            stacked+=1
            #print "Right bulge, stacking = " + str(stacked)
          else:
            if(stacked<=minstack):
              #print "remove stem, start = " + str(start) + " stacked = " + str(stacked)
              for j in range(start,i):
                #print j,pbps[j]
                pbps[pbps[j]] = pbps[j] 
                pbps[j] = j
              #print pbps
            
            #print "stacked starts at " + str(i) 
            start = i
            stacked = 1
            started = 1
            bulge = 0
    
    if(pbps[i]==i):
      if(started>0):
        if(bulge==0):
            #print "Bulge start"
            bulge = 1
        else:
          if(stacked<=minstack):
            #print "remove stem, start = " + str(start) + " stacked = " + str(stacked)
            for j in range(start,i):
              #print i,j,pbps[j]
              pbps[pbps[j]] = pbps[j]
              pbps[j] = j
            #print pbps
          start = -1
          started = 0
          bulge = 0
          stacked = 0
    
    if(pbps[i]<i):#end hairpin
      if((started>0)and(stacked<=minstack)):
        #print "remove stem, start = " + str(start) + " stacked = " + str(stacked) 
        for j in range(start,i):
          pbps[pbps[j]] = pbps[j]
          pbps[j] = j
      bulge = 0
      start = -1
      started = 0
      stacked = 0
  return pbps

def buildDP(nP):
 L=len(nP)
 s=zeros((L,L))
 for n in xrange(1,L):
   for j in xrange(n,L):
     i=j-n
     case1=s[i+1,j-1]+ nP[i][j] #delta(seq[i],seq[j]);
     case2=s[i+1,j]
     case3=s[i,j-1]
     if i+3<=j:
       tmp=[]
       for k in xrange(i+1,j):
         tmp.append(s[i,k]+s[k+1,j])
       case4=max(tmp)
       s[i,j]=max(case1,case2,case3,case4)
     else:
       s[i,j]=max(case1,case2,case3)
 return s

def traceback(s,nP,i,j,pair):
 if i<j:
  if s[i,j]==s[i+1,j]:
   traceback(s,nP,i+1,j,pair)
  elif s[i,j]==s[i,j-1]:
   traceback(s,nP,i,j-1,pair)
  elif s[i,j]==s[i+1,j-1]+nP[i][j]: #delta(seq[i],seq[j]):
   pair.append((i,j))
   traceback(s,nP,i+1,j-1,pair)
  else:
   for k in xrange(i+1,j):
    if s[i,j]==s[i,k]+s[k+1,j]:
     traceback(s,nP,i,k,pair)
     traceback(s,nP,k+1,j,pair)
     break
 return pair

def nussinov(nP):
  n = len(nP)
  pair=traceback(buildDP(nP),nP,0,n-1,[])
  bps = []
  for k in range(n):
    bps.append(k)
  for p in pair:
    a,b = p
    bps[a] = b
    bps[b] = a
  return bps

def calculateStr(i,j,P):
  n = j-i
  nP = []
  for k in range(n):
    nP.append([])
    for l in range(n):
      nP[k].append(0.0)
  for l in range(i,j):
    for m in range(i+1,j):
      if P[l+1][m+1]>0:
        nP[l-i][m-i] = int(1000000*P[l+1][m+1])
        nP[m-i][l-i] = int(1000000*P[l+1][m+1])

  bps = nussinov(nP)  
  pbps = postprocess(bps)

  str = ""
  for k in range(n):
    if pbps[k]==k:
      str+='.'
    elif pbps[k]>k:
      str+='('
    else:
      str+=')'

  return str



f = open(fn)
line = f.readline()
while(len(line)>0):
  comm = line.strip()
  line = f.readline().strip()
  P=boltzmannBasePairProb(comm,line)
  i,j,seq = getCapitals(line)
  str = calculateStr(i,j,P)
  print comm
  print seq
  print str
  line = f.readline()
f.close()
