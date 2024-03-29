
AAs="ACDEFGHIKLMNPQRSTVWY"
NewMat='''L 2.998 0.269 0.367 3.852 0.510 0.514 0.562
I 2.408 0.261 0.345 0.931 0.402 0.440 0.289
V 1.525 0.479 0.350 0.887 0.286 0.350 0.362
M 2.161 0.605 0.442 1.441 0.607 0.457 0.570
F 0.490 0.075 0.391 0.639 0.125 0.081 0.038
Y 1.319 0.064 0.081 1.526 0.204 0.118 0.096
G 0.084 0.215 0.432 0.111 0.153 0.367 0.125
A 1.283 1.364 1.077 2.219 0.490 1.265 0.903
K 1.233 2.194 1.817 0.611 2.095 1.686 2.027
R 1.014 1.476 1.771 0.114 1.667 2.006 1.844
H 0.590 0.646 0.584 0.842 0.307 0.611 0.396
E 0.281 3.351 2.998 0.789 4.868 2.735 3.812
D 0.068 2.103 1.646 0.182 0.664 1.581 1.401
Q 0.311 2.290 2.330 0.811 2.596 2.155 2.585
N 1.231 1.683 2.157 0.197 1.653 2.430 2.065
S 0.332 0.753 0.930 0.424 0.734 0.801 0.518
T 0.197 0.543 0.647 0.680 0.905 0.643 0.808
C 0.918 0.002 0.385 0.440 0.138 0.432 0.079
W 0.066 0.064 0.065 0.747 0.006 0.115 0.014
P 0.004 0.108 0.018 0.006 0.010 0.004 0.007'''
Parameters='''w  14 1.89 0.30 1.04 0.27 20
w  21 1.79 0.24 0.92 0.22 25
w  28 1.74 0.20 0.86 0.18 30
uw 14 1.82 0.28 0.95 0.26 20
uw 21 1.74 0.23 0.86 0.21 25
uw 28 1.69 0.18 0.80 0.18 30''' 

class Params:
     ''' define the Coils parameters '''
     def __init__(self,weighted='w',window=21):
         ''' initialize the parameters '''
         self.pow=[1.0]*7
         if weighted.lower() not in ['w','uw']:
             self.weight='uw'
         else:
             self.weight=weighted.lower()
         if self.weight=='w':
             self.pow[0]=self.pow[3]=2.5
         if window not in [14,21,28]:
             self.win=21
         else:
             self.win=window
         self.m_cc,self.sd_cc,self.m_g,self.sd_g,self.sc,self.mat=_initData(self.win,self.weight)


def _initData(window,weighted):
    ''' _initData builds the parameters starting from string '''
    from string import atof, atoi
    s=NewMat
    sl=s.split('\n')
    d={}
    for e in sl:
        v=e.split()
        for k in range(len(v[1:])):
            d[(v[0],k)]=atof(v[1:][k])
    s=Parameters
    sl=s.split('\n')
    for e in sl:
        v=e.split()
	if (v[0] == weighted) and (atoi(v[1])==window):
           vec=[]
           for n in v[2:]:
               vec.append(atof(n))        
    return vec+[d]
   


''' Original file new.mat 
% weighted
w  14 1.89 0.30 1.04 0.27 20
w  21 1.79 0.24 0.92 0.22 25
w  28 1.74 0.20 0.86 0.18 30
uw 14 1.82 0.28 0.95 0.26 20
uw 21 1.74 0.23 0.86 0.21 25
uw 28 1.69 0.18 0.80 0.18 30
%
%   a     b     c     d     e     f     g
L 2.998 0.269 0.367 3.852 0.510 0.514 0.562
I 2.408 0.261 0.345 0.931 0.402 0.440 0.289
V 1.525 0.479 0.350 0.887 0.286 0.350 0.362
M 2.161 0.605 0.442 1.441 0.607 0.457 0.570
F 0.490 0.075 0.391 0.639 0.125 0.081 0.038
Y 1.319 0.064 0.081 1.526 0.204 0.118 0.096
G 0.084 0.215 0.432 0.111 0.153 0.367 0.125
A 1.283 1.364 1.077 2.219 0.490 1.265 0.903
K 1.233 2.194 1.817 0.611 2.095 1.686 2.027
R 1.014 1.476 1.771 0.114 1.667 2.006 1.844
H 0.590 0.646 0.584 0.842 0.307 0.611 0.396
E 0.281 3.351 2.998 0.789 4.868 2.735 3.812
D 0.068 2.103 1.646 0.182 0.664 1.581 1.401
Q 0.311 2.290 2.330 0.811 2.596 2.155 2.585
N 1.231 1.683 2.157 0.197 1.653 2.430 2.065
S 0.332 0.753 0.930 0.424 0.734 0.801 0.518
T 0.197 0.543 0.647 0.680 0.905 0.643 0.808
C 0.918 0.002 0.385 0.440 0.138 0.432 0.079
W 0.066 0.064 0.065 0.747 0.006 0.115 0.014
P 0.004 0.108 0.018 0.006 0.010 0.004 0.007
'''


