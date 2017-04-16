import sys,os
import networkx as nx
from eden.util import display
from sklearn import metrics
from eden.graph import Vectorizer
import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from dpcluster import dclust
import forgi.graph.bulge_graph as fgb
from random import random
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

#####fixed params#####
MINL = 10 #min length of a peak for it to be considered
MAXL = 120 #max length of a peak for it to be considered
MAXP = 800 #maximun number of peaks to be clustered in the same iteration
######################

used = []

clustalw_exe = "/Users/idotu/tools/clustalw-2.1-macosx/clustalw2"

iupac = {"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["A","C"],"B":["C","G","T"],
    "D":["A","G","T"],"H":["A","C","T"],"V":["A","C","G"],"N":["A","C","G","T"]}

letters = ['a','b','c','d','e','f','g','h','i','j','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

##########################################################
# Reads file with fasta comment \n sequence \n structure #
# Calculates graph structure of shape with sequence in   #
# unpaired regions and stems collapsed. Then it clusters #
# them using sklearn                                     #
##########################################################

def extendH(a,b,stru):
  n = len(stru)
  #print a,b,n
  ss = []
  if a==0 and b==n-1:
    return [(a,b+1)]
  aa = a
  while aa>0 and stru[aa-1]==".":
    aa -= 1
  if aa<a:
    ss.append((aa,b+1))
  bb = b
  while bb<n-1 and stru[bb+1]==".":
    bb += 1
  if bb>b:
    ss.append((a,bb+1))
  return ss


def findSubstrs(stru):
  ss = []
  hs = []
  belongs = {}
  n = len(stru)
  bg = fgb.BulgeGraph(dotbracket_str=stru)
  hairpins = []
  if len(bg.defines.keys())==1:
    return [(0,n)]
  for key in bg.defines.keys():
    if key[0]=="h":
      hairpins.append(key)
  for j in range(n):
    key = bg.get_node_from_residue_num(j+1)
    if key not in belongs:
      belongs[key] = []
    belongs[key].append(j)
          
  for h in hairpins:
    key = h
    a = -1
    b = -1
    used_keys = [h]
    last_key = h
    while key[0]!="m":
        #print key, last_key,a,b
      last_key = key
      for key1 in bg.edges[key]:
          #print key1
        if key1 not in used_keys:
          if key1[0]=="s":
            a = min(belongs[key1])
            b = max(belongs[key1])
            
            key = key1
            used_keys.append(key)
            break
          elif key1[0]=="m":
            last_key = key
            key=key1
            break
          elif key1[0]=="i":
            last_key
            key=key1
            used_keys.append(key)
            break
    #print key,last_key
      if last_key==key:
        break
    hs.append((a,b))
#print hs
  for h in hs:
    a,b = h
    ss2 = extendH(a,b,stru)
    for s in ss2:
      ss.append(s)
  return ss


def readFile(fn,op):
  f = open(fn)
  peaks = {}
  line = f.readline().strip()
  i = 0
  while(len(line)>0):
    comm = line[1:]
    line = f.readline().strip()
    #print line
    seq = ""
    while(line[0]!='.' and line[0]!='('):
      seq += line
      line = f.readline().strip()
    stru = ""
    while(len(line)>0 and line[0]!='>'):
      stru += line
      line = f.readline().strip()
    stru = stru.split()[0]
    n = len(seq)
    if op == 0:
      ss = [(0,n+1)]
    else:
      ss = findSubstrs(stru)
    m = len(ss)

    if m<len(letters):
      for j in range(m):
        a,b = ss[j]
        words = comm.split("|")
        comm2 = ""
        for k in range(len(words)):
          if k==2:
            comm2+=words[k]+letters[j]+"|"
          else:
            comm2 += words[k]+"|"
        if comm2 not in used: #and i<MAXP:
          seq2 = seq[a:b]
          stru2 = stru[a:b]
          peaks[comm2] = (seq2,stru2)
          i += 1

  f.close()
  

  nused = len(used)
  #print i,nused
  if i-nused>MAXP:
    peaks2 = {}
    f = MAXP/float(i-nused)
    for key in peaks:
      if key not in used:
        
        rd = random()
        if rd < f:
          peaks2[key] = peaks[key]
    #used.append(key)
    return peaks2
  else:
    if i < MAXP:
      return peaks
    else:
      return []

def getbps(stru):
  n = len(stru)
  bps = []
  for i in range(n):
    bps.append(i)
  pile = []
  for i in range(n):
    if stru[i]=="(":
      pile.append(i)
    elif stru[i]==")":
      k = pile.pop()
      bps[k] = i
      bps[i] = k
  #print bps
  return bps


def getGraph5(peak,stems):
  seq,stru = peak
  bps = getbps(stru)
  G = nx.Graph()
  pile = []
  open = 0
  close = 0
  single = 0
  j = 0
  opending = -1
  cpending = -1
  for i in range(len(seq)):
    #print i,pile
    if stru[i]==".":
      if single==1:
        G.add_node(j,label = seq[i])
        G.add_edge(j-1,j,label = 'bb')
        j += 1
      else:
       single = 1
       if open>0:
         G.add_node(j,label = seq[i])
         G.add_edge(j-1,j,label = 'bb')
         j += 1
         open = 0
       elif close>0:
         G.add_node(j,label = seq[i])
         G.add_edge(j-1,j,label = 'bb')
         j += 1
       else: #it is the first node
         G.add_node(j,label = seq[i])
         j += 1
    elif stru[i]=="(":
      if stems==1:
        lb = seq[i]
      else:
        lb = 'p'
      if single == 1:
        G.add_node(j,label = lb)
        G.add_edge(j-1,j,label = 'bb')
        open = 1
        pile.append((j,open))
        j += 1
        single = 0
      elif open>0:
        if bps[i-1]==bps[i]+1:
          G.add_node(j,label = lb)
          G.add_edge(j-1,j,label = 'bb')
          open += 1
          pile.append((j,open))
          j += 1
        else:#close loop
          G.add_node(j,label = lb)
          G.add_edge(j-1,j,label = 'bb')
          j += 1
          pile.append((j-1,open))
          open = 1
      elif close>0:
        G.add_node(j,label = lb)
        G.add_edge(j-1,j,label = 'bb')
        open = 1
        pile.append((j,open))
        j += 1
      else:#first nt
        open = 1
        G.add_node(j,label = lb)
        pile.append((j,open))
        j += 1       
    else: #stru[i]==")"
      if stems==1:
        lb = seq[i]
      else:
        lb = 'p'
      if single == 1:
        G.add_node(j,label = lb)
        G.add_edge(j-1,j,label = 'bb')
        k,w = pile.pop()
        G.add_edge(k,j,label = 'bp')
        j += 1
        close = 1
        single = 0
      elif close>0:
        if bps[i]==bps[i-1]-1:
          close += 1
          G.add_node(j,label = lb)
          G.add_edge(j-1,j,label = 'bb')
          k,w = pile.pop()
          G.add_edge(k,j,label = 'bp')
          j += 1
        else:
          G.add_node(j, label = lb)
          G.add_edge(j-1,j,label = 'bb')
          k,w = pile.pop()
          G.add_edge(k,j,label = 'bp')
          j += 1
          close = 1  
      else:#first nt 
        close = 1

  if DEBUG == 1:
    print seq
    print stru
    display.draw_graph(G, size=15, node_size=1500, font_size=24, node_border=True, size_x_to_y_ratio=3)
  return G


def getGraph8(peak,G1):
  seq,stru = peak
  bps = getbps(stru)
  G = G1.copy()#nx.Graph()
  pile = []
  n = len(seq)
  next_pos = 0
  last_pos = 0
  node = n
  pendingexternal = -1
  pendingstem = -1
  while(next_pos<n):
    last_pos=next_pos
    
    while(next_pos<n and stru[next_pos]=='.'):
      
      next_pos += 1
    if next_pos>last_pos:#starts with unpaired rexion
      
      G.add_node(node,label='EL')

      #Combine
      for i in range(last_pos,next_pos):
        G.add_edge(i,node,label='s')
      
      if pendingexternal!=-1:
        G.add_edge(node,pendingexternal,label='abb')
      node += 1
      pendingstem = node-1
    else:
      pendingstem = pendingexternal

    if next_pos==n:
      break
    #Start stem
    i = next_pos
    j = bps[i]
    pile.append(i)
    pile.append(j)
    pendingexternal = node
    next_pos = j + 1
    while(i<j):
      while(bps[i+1]==j-1):
        i += 1
        j -= 1
        pile.append(i)
        pile.append(j)

      G.add_node(node,label='S')

      #Combine
      for p in range(len(pile)):
        G.add_edge(pile[p],node,label='s')
      pile = []

      if pendingstem!=-1:
        G.add_edge(pendingstem,node,label = 'abb')
        pendingstem = -1

      node += 1
      if hl(i+1,j,stru):
        G.add_node(node,label='HL')
        G.add_edge(node-1,node,label='abb')
        #Combine
        for p in range(i+1,j):
          G.add_edge(p,node,label='s')
        
        
        node += 1
        i = j
      elif stru[i+1]=='.' and stru[j-1]=='.':#IL
        i += 1
        init_i = i
        init_j = j
        while(stru[i]=='.'):
          i += 1
        j = bps[i]
        G.add_node(node,label='IL')
        G.add_edge(node-1,node,label='abb')
        #Combine
        for p in range(init_i,i):
          G.add_edge(p,node,label='s')
        for p in range(j+1,init_j):
          G.add_edge(p,node,label='s')

        pendingstem = node
        node += 1

      elif stru[i+1]=='.':#LB
        i += 1
        init_i = i
        while(stru[i]=='.'):
          i += 1
        j = bps[i]
        G.add_node(node,label='LB')
        G.add_edge(node-1,node,label='abb')
        #Combine
        for p in range(init_i,i):
          G.add_edge(p,node,label='s')
        pendingstem = node
        node += 1
      else:#RB
        init_j = j
        j -= 1
        while(stru[j]=='.'):
          j -= 1
        i = bps[j]
        G.add_node(node,label='RB')
        G.add_edge(node-1,node,label='abb')
        #Combine
        for p in range(j+1,init_j):
          G.add_edge(p,node,label='s')
        pendingstem = node
        node += 1
  
  if DEBUG==1:
    print seq
    print stru
    display.draw_graph(G, size=15, node_size=1500, font_size=24, node_border=True, size_x_to_y_ratio=3)

  return G

def getGraph9(peak):
  seq,stru = peak
  bg = fgb.BulgeGraph(dotbracket_str=stru)
  
  dict = {}
  i = 0
  G = nx.Graph()
  for key2 in bg.defines.keys():
    G.add_node(i,label=key2[0])
    dict[key2] = i
    i += 1

  for key1 in bg.edges:
    for key2 in bg.edges[key1]:
      G.add_edge(dict[key1],dict[key2],label='v')

  if DEBUG==1:
    print seq
    print stru
    display.draw_graph(G, size=15, node_size=1500, font_size=24, node_border=True, size_x_to_y_ratio=3)

  return G

def getGraph11(peak,G1):
  seq,stru = peak
  bps = getbps(stru)
  bg = fgb.BulgeGraph(dotbracket_str=stru)
  n = len(seq)
  G = G1.copy()
  i = n
  dict = {}
                            
  for key2 in bg.defines.keys():
    if key2[0]=='h':
      G.add_node(i,label=key2[0])
    elif key2[0]=='s':
      G.add_node(i,label=key2[0])
    else:
      G.add_node(i,label='i')
    dict[key2] = i
    i += 1
        
  for key1 in bg.edges:
    for key2 in bg.edges[key1]:
      G.add_edge(dict[key1],dict[key2],label='v')
        
  for i in range(n):
    key = bg.get_node_from_residue_num(i+1)
    G.add_edge(i,dict[key],label='abb')

  for i in range(1,n-4):
    j = bps[i]
    if i<j:
      if bps[i+1]==i+1 and bps[j-1]==j-1:
        key = bg.get_node_from_residue_num(i+2)
        if key[0]=='i' or key[0]=='h':
          G.add_edge(i,dict[key],label='abb')
          G.add_edge(j,dict[key],label='abb')
      elif bps[i+1]==i+1:
        key = bg.get_node_from_residue_num(i+2)
        if key[0]=='i' or key[0]=='h':
          G.add_edge(i,dict[key],label='abb')
      elif bps[j-1]==j-1:
        key = bg.get_node_from_residue_num(j)
        if key[0]=='i' or key[0]=='h':
          G.add_edge(j,dict[key],label='abb')

      if j+1<n and bps[i-1]==i-1 and bps[j+1]==j+1:
        key = bg.get_node_from_residue_num(i)
        if key[0]=='i' or key[0]=='h':
          G.add_edge(i,dict[key],label='abb')
          G.add_edge(j,dict[key],label='abb')
      elif bps[i-1]==i-1:
        key = bg.get_node_from_residue_num(i)
        if key[0]=='i' or key[0]=='h':
          G.add_edge(i,dict[key],label='abb')
      elif j+1<n and bps[j+1]==j+1:
        key = bg.get_node_from_residue_num(j+2)
        if key[0]=='i' or key[0]=='h':
          G.add_edge(j,dict[key],label='abb')
    


  if DEBUG==1:
    print seq
    print stru
    display.draw_graph(G, size=15, node_size=1500, font_size=24, node_border=True, size_x_to_y_ratio=3)
        
  return G

def getGraph14(peak,all,upsc,struYes):
  seq,stru2 = peak
  if struYes==1:
    stru = stru2
  else:
    stru = ""
    for i in range(len(seq)):
      stru += "."
  bgw = 10
  bg = fgb.BulgeGraph(dotbracket_str=stru)
  n = len(seq)
  if len(bg.defines.keys())==1: #Empty structure
    G = nx.Graph()
    for i in range(n):
      G.add_node(i,label=seq[i],weight=1)
    for i in range(n-1):
      G.add_edge(i,i+1,label='bb',weight=1)
    if DEBUG==1:
      print seq
      print stru
      display.draw_graph(G, size=15, node_size=1500, font_size=24, node_border=True, size_x_to_y_ratio=3)
    return G
  
  n = len(seq)
  dict = {}
  i = 0
  G = nx.Graph()
  for key2 in bg.defines.keys():
    lb = key2[0]
    if lb!='s':
      lb = 'u'
      w = 1
    else:
      w = 1
    G.add_node(i,label=lb,weight=w)
    dict[key2] = i
    i += 1
        
  hairpins = {}
            
  for key1 in bg.edges:
    for key2 in bg.edges[key1]:
      G.add_edge(dict[key1],dict[key2],label='v',weight=bgw)

  if all==0:
    for j in range(n):
      key = bg.get_node_from_residue_num(j+1)
      if key[0]=='h' or G.number_of_nodes()==1: #or key[0]=='m' or key[0]=='i':
        if key not in hairpins:
          hairpins[key] = []
        hairpins[key].append(j)
  elif all==1:
    for j in range(n):
      key = bg.get_node_from_residue_num(j+1)
      if key[0]=='h' or key[0]=='i' or G.number_of_nodes()==1:
        if key not in hairpins:
          hairpins[key] = []
        hairpins[key].append(j)
  elif all == 2:
    for j in range(n):
      key = bg.get_node_from_residue_num(j+1)
      if key[0]=='i' or G.number_of_nodes()==1:
        if key not in hairpins:
          hairpins[key] = []
        hairpins[key].append(j)
  elif all==3:
    for j in range(n):
      key = bg.get_node_from_residue_num(j+1)
      if (key[0]=='s') or G.number_of_nodes()==1:
        if key not in hairpins:
          hairpins[key] = []
        hairpins[key].append(j)
  elif all==4:
    for j in range(n):
      key = bg.get_node_from_residue_num(j+1)
      if key not in hairpins:
        hairpins[key] = []
      hairpins[key].append(j)
  elif all==5:
    for j in range(n):
      key = bg.get_node_from_residue_num(j+1)
      if key[0]=='m' or key[0]=='f' or key[0]=='t' or G.number_of_nodes()==1:
        if key not in hairpins:
          hairpins[key] = []
        hairpins[key].append(j)
  else:
    for j in range(n):
      key = bg.get_node_from_residue_num(j+1)
      if (key[0]!='s') or G.number_of_nodes()==1:
        if key not in hairpins:
          hairpins[key] = []
        hairpins[key].append(j)

  for key in hairpins:
    
    st = min(hairpins[key])
    ed = max(hairpins[key])
    G.add_node(i,label=seq[st:st+1])
    if key[0]!='f':
      G.add_edge(dict[key],i,label='abb',weight=bgw)
    i += 1
    if key[0]!='i' and key[0]!='s':
      for j in range(st+1,ed):
        G.add_node(i,label=seq[j:j+1])
        G.add_edge(i-1,i,label='bb')
        if upsc==1:
          G.add_edge(dict[key],i,label='abb',weight=bgw)
        i += 1
      if len(hairpins[key])>1:
        G.add_node(i,label=seq[ed:ed+1])
        G.add_edge(i-1,i,label='bb')
        i += 1
      if key[0]!='t':
        G.add_edge(i-1,dict[key],label='abb',weight=bgw)
    else:
      j = st+1
      while(j in hairpins[key]):
        G.add_node(i,label=seq[j:j+1])
        G.add_edge(i-1,i,label='bb')
        i += 1
        j += 1
      if j-st>1:
        G.add_edge(i-1,dict[key],label='abb',weight=bgw)
      elif upsc==1:
        G.add_edge(dict[key],i-1,label='abb',weight=bgw)
      if st!=ed and j!=ed+1:
        G.add_node(i,label=seq[ed:ed+1])
        G.add_edge(dict[key],i,label='abb',weight=bgw)
        if upsc==1:
          G.add_edge(dict[key],i,label='abb',weight=bgw)
        i += 1
        j = ed-1
        while(j in hairpins[key]):
          G.add_node(i,label=seq[j:j+1])
          G.add_edge(i-1,i,label='bb')
          if upsc==1:
            G.add_edge(dict[key],i,label='abb',weight=bgw)
          i += 1
          j -= 1
        if ed-j>1 and upsc==0:
          G.add_edge(i-1,dict[key],label='abb',weight=bgw)

  if DEBUG==1:
    print seq
    print stru
    display.draw_graph(G, size=15, node_size=1500, font_size=24, node_border=True, size_x_to_y_ratio=3)

  return G


def peaksToGraphs(peaks,op):
  graphs = []
  dict = {}
  i = 0
  for key in peaks:
    peak = peaks[key]
    if op==0:
        G = getGraph14(peak,4,0,0)
    if op==1:
      G1 = getGraph5(peak,1)
      G = getGraph8(peak,G1)
    if op==2:
      G1 = getGraph5(peak,1)
      G = getGraph10(peak,G1)
    if op==3:
      G = getGraph9(peak)
    if op==4:
      G = getGraph14(peak,0,0,1)
    if op==5:
      G = getGraph14(peak,2,0,1)
    if op==6:
      G = getGraph14(peak,5,0,1)
    if op==7:
      G = getGraph14(peak,3,0,1)
    if op==8:
      G = getGraph14(peak,1,0,1)
    if op==9:
      G = getGraph14(peak,6,0,1)
    if op==10:
      G = getGraph14(peak,4,0,1)
    if op==11:
      G1 = getGraph5(peak,0)
      G = getGraph11(peak,G1)
    graphs.append(G)
    dict[i] = key
    i += 1
  return graphs,dict
  
def clusterGraphs(graphs,r,d,copt):
  opts = copt[1:-1]
  optl = opts.split(",")
  opt = int(optl[0])
  vectorizer = Vectorizer( r=r,d=d )
  samples = len(graphs)
  Xsp = vectorizer.transform( graphs )#sparse feature matrix
  X = Xsp.todense()#regular feature matrix
  SM=metrics.pairwise.pairwise_kernels(Xsp, metric='rbf', gamma = 1)#similarity matrix
  DM=[]#distance matrix
  for i in range(len(SM)):
    DM.append([])
    for j in range(len(SM[i])):
      val = 1.0 - SM[i][j]
      if val<0:
        DM[i].append(0.0)
      else:
        DM[i].append(val)
  if opt==0:
    nc,labels = MShift(X)
  if opt==1:
    #print DM
    nc,labels = DB_SCAN(DM,float(optl[1]),int(optl[2]))
  if opt==2:
    nc,labels = AffProp(SM)
  if opt==3:
    print DM #Matrix(X)
    return 0,[]
  if opt==4:
    nc,labels = K_Means(X)
  if opt==5:
    nc,labels = SpecClus(SM) 
  if opt==6:
    nc,labels = dclust(DM,int(optl[1]),int(optl[2]),float(optl[3]))

  return nc,labels

def SpecClus(SM):
  sc = SpectralClustering(n_clusters=8, affinity = 'precomputed')
  sc.fit_predict(SM)
  labels = sc.labels_

  n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
  return n_clusters_,labels

def K_Means(X):

  km = KMeans(n_clusters=8, precompute_distances=True, max_iter=30)
  km.fit_predict(X)
  labels = km.labels_

  n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
  return n_clusters_,labels

def Matrix(X):
  K=metrics.pairwise.pairwise_kernels(X, metric='linear')
  print K


def AffProp(SM):

  af = AffinityPropagation(preference=None, affinity='precomputed')
  af.fit_predict(SM)
  cluster_centers_indices = af.cluster_centers_indices_
  labels_ = af.labels_

  n_clusters_ = len(cluster_centers_indices)
  return n_clusters_,labels_

def MShift(X):
  # The following bandwidth can be automatically detected using
  #bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=samples/2)

  ms = MeanShift(bandwidth=None, cluster_all= False, bin_seeding=True)#bandwidth=bandwidth, bin_seeding=True)
  ms.fit_predict(X)
  labels = ms.labels_
  cluster_centers = ms.cluster_centers_

  labels_unique = np.unique(labels)
  n_clusters_ = len(labels_unique)



def DB_SCAN(X,sim,sam):
  ############# Compute DBSCAN
  db = DBSCAN(eps=sim, min_samples=sam, metric = 'precomputed').fit(X)
  core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
  core_samples_mask[db.core_sample_indices_] = True
  labels = db.labels_

  n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
  return n_clusters_,labels

def compare(finalL,L,peaks,opt,th):
  n = len(L)
  lpeaks = {}
  for key in L:
    lpeaks[key] = peaks[key]
  for key in finalL:
      lpeaks[key] = peaks[key]
  graphs,dict = peaksToGraphs(lpeaks,opt)
  
  vectorizer = Vectorizer( r=2,d=3 )
  samples = len(graphs)
  Xsp = vectorizer.transform( graphs )#sparse feature matrix
  X = Xsp.todense()#regular feature matrix
  SM=metrics.pairwise.pairwise_kernels(Xsp, metric='rbf', gamma = 1)#similarity matrix
  DM=[]#distance matrix
  for i in range(len(SM)):
    DM.append([])
    for j in range(len(SM[i])):
      val = 1.0 - SM[i][j]
      if val<0:
        DM[i].append(0.0)
      else:
        DM[i].append(val)
  avgDM = 0.0
  counts =0.0
  for i in range(len(graphs)):
    if dict[i] in L:
      for j in range(len(graphs)):
        if i!=j and dict[j] in finalL:
          avgDM += DM[i][j]
          counts += 1
  avgDM = avgDM/counts
  if avgDM>=0.0 and avgDM<= th:
    return 0
  else:
    return 1


def getClusters(nc,labels,dict,peaks,cf,clus,nclus,opt,th,usedpeaks):
  numClusters = nclus
  for i in range(nc):
    objs = []
    for j in range(len(labels)):
      if labels[j] == i:
        objs.append(dict[j])
        used.append(dict[j])
        usedpeaks[dict[j]] = peaks[dict[j]]
    newClu = 1
    for k in range(numClusters):
      newClu = compare(clus[k],objs,usedpeaks,opt,th)
      if newClu==0:
        for l in range(len(objs)):
          clus[k].append(objs[l])
        break
    if newClu==1:
      clus[numClusters] = []
      for l in range(len(objs)):
        clus[numClusters].append(objs[l])
      numClusters += 1

  return numClusters

def printResults(nc,labels,dict,peaks,cf):
  print "Number of clusters = " + str(nc)
  sum = 0
  cluster = {}
  for i in range(nc):
    print "****Cluster " + str(i) + ":"
    for j in range(len(labels)):
      if labels[j] == i:
        sum += 1
        key = dict[j]
        cluster[key] = i
        sq,st = peaks[key]
        used.append(key)
        print ">"+key
        print sq
        print st + " #S"
    print "*******************************************"
  if cf!='None':
    cpeaks = readFile(cf)
    for key in cpeaks:
      if key in peaks:
        sq,st = peaks[key]
        pos = -1
        for j in dict:
          if key ==dict[j]:
            pos = j
            break
        if key in cluster:
          print ">"+key + " - pos " + str(pos) + " --> Cluster = " + str(cluster[key])
        else:
          print ">"+key + " - pos " + str(pos) + " --> Cluster = noise"
        print sq
        print st
      else:
        sq,st = cpeaks[key]
        print ">"+key + " --> Cluster = not included"
        print sq
        print st

def equals(A,B):
  n = len(A)
  m = len(B)
  if n!=m:
    return 0
  for i in range(n):
    if A[i] not in B:
      return 1
  return 0

def getAlginedStru(stru,seq):
  n = len(seq)
  stru2 = ""
  i = 0
  for j in range(n):
    if seq[j]=="-":
      stru2 += "-"
    else:
      stru2 += stru[i]
      i += 1
  return stru2

def getLogos(seqs,strus):
  n = len(seqs)
  m = len(seqs[0])
  sql = ""
  stl = ""
  for i in range(m):
    sqd = {}
    std = {}
    for j in range(n):
      if seqs[j][i] not in sqd:
        sqd[seqs[j][i]] = 0
      sqd[seqs[j][i]] = sqd[seqs[j][i]] + 1
      if strus[j][i] not in std:
        std[strus[j][i]] = 0
      std[strus[j][i]] = std[strus[j][i]] + 1
    th = n*0.9
    sqch = "N"
    for key in sqd:
      if sqd[key]>=th:
        sqch = key
        break
    stch = ","
    for key in std:
      if std[key]>=th:
        stch = key
        break
    sql += sqch
    stl += stch
  return sql.strip("-"),stl.strip("-")


def main(fn,r,d,opt,gopt,cf,its,thClus,v):
  usedpeaks = {}
  finalClusters = {}
  numClusters = 0
  for i in range(its):
    peaks = readFile(fn,gopt)
    #print peaks
    if len(peaks)==0:
      break
    graphs,dict = peaksToGraphs(peaks,gopt)
    nc,labels = clusterGraphs(graphs,r,d,opt)
    if opt!=3:
      numClusters = getClusters(nc,labels,dict,peaks,cf,finalClusters,numClusters,gopt,thClus,usedpeaks)

  for i in range(numClusters):
    smallKeyToStru = {}
    print "Cluster " + str(i)
    ofn ="Cluster" + str(i) + ".fa"
    f = open(ofn,"w")
    for j in range(len(finalClusters[i])):
      key = finalClusters[i][j]
      sq,st = usedpeaks[key]
      words = key.split("|")
      key2 = words[0]+"|"+words[1]+"|"+words[2][:-1]
      if key2 not in smallKeyToStru:
        smallKeyToStru[key2] = st
        f.write(">"+key2+"\n")
        f.write(sq+"\n")
    f.close()
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=ofn)
    stdout, stderr = clustalw_cline()
    align = AlignIO.read(ofn[:-3]+".aln", "clustal")
    seqs = []
    strus = []
    ids = []
    for record in align:
      ids.append(record.id)
      seqs.append(record.seq)
      stru = smallKeyToStru[record.id]
      alignedstru = getAlginedStru(stru,record.seq)
      strus.append(alignedstru)
    if v==1:
      for j in range(len(seqs)):
        print seqs[j] + "\t" + ids[j]
        print strus[j]
    seqlogo,strulogo = getLogos(seqs,strus)


    print "*******************************************"
    print seqlogo
    print strulogo
    print len(finalClusters[i])
    print "*******************************************"


#######################################################
# r,d = radius,distance - computes features for every #
# pair of node neighborhhods of radius r that are at  #
# distance d in the graph.                            #
# cluster_option:                                     #
# 0 - MeanShift                                       #
# 1 - DBSCAN                                          #
# 2 - Affinity Propagation                            #
# 3 - Prints similarity matrix                        #
# 4 - KMeans                                          #
# 5 - Spectral Clustering                             #
# 6 - Density Clustering (Science)                    #
# Graph_options:                                      #
# 0 - No structure considered                         #
# 1 - graphProt minus directional                     #
# 2 - Bulge graph plus Complete Graph                 #
# 3 - Bulge graph                                     #
# 4 - Bulge graph plus sequence in HP loops           #
# 5 - Bulge graph plus sequence in IL loops           #
# 6 - Bulge graph plus sequence in external loops     #
# 7 - Bulge graph plus sequence in double stranded    #
# 8 - Bulge graph plus sequence in HP and IL          #
# 9 - Bulge graph plus sequence in single stranded    #
# 10 - Bulge graph plus sequence everywhere           #
# 11 - Same as 2 but no sequence in double stranded   #
#######################################################

if __name__ == '__main__':
  if len(sys.argv)<8:
    print "Usage: %s file r d cluster_option graph_option thClus iterations verbose debug"
    sys.exit(1)
  fn = sys.argv[1]
  r = int(sys.argv[2])
  d = int(sys.argv[3])
  opt = sys.argv[4]
  gopt = int(sys.argv[5])
  thClus = float(sys.argv[6])
  its = int(sys.argv[7])
  verbose = int(sys.argv[8])
  DEBUG = int(sys.argv[9])
  cf = "None"
  main(fn,r,d,opt,gopt,cf,its,thClus,verbose)
