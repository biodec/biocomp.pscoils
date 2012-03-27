#!/usr/bin/env python
'''
    Copyright (C) 2007 Piero Fariselli

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    psCoils.py  Copyright (C) 2007  Piero Fariselli
    contacts: Piero Fariselli
              e-mail: piero@biocomp.unibo.it, piero.fariselli@unibo.it
              Dept. of Biology
              University of Bologna
              via Irnerio 42
              40126 Bologna
              Italy
    This program comes with ABSOLUTELY NO WARRANTY; for details type psCoils.py.
    This is free software, and you are welcome to redistribute it
    under the condition of preserving this text 
'''

# Iports

import sys
import math
import string
import copy

import Params # this must be in the sys.path. Params contains the COILS paramters

from Bio.SeqFeature import FeatureLocation, SeqFeature

#
#
#---------------------------------#
def readProf(fname,aaOrder='VLIMFWYGAPSTCHRKQEND'):
    ''' readProf(fname,aaOrder) returns a dictionary profile prof[(pos,aa)]'''
    try:
        lines=open(fname,'r').readlines()
    except:
        sys.stderr.write('cannot open file '+fname+'\n')
        sys.exit()
    prof={}
    i=0
    while lines[i].find("#")>=0 and i <len(lines):
        pos=lines[i].find("ALPHABET=")
        if pos>0:
           v=lines[i][pos+len("ALPHABET="):].split()
           aaOrder=''.join(v[:20])
           
        i+=1
    if i >= len(lines):
        sys.stderr.write('Error in '+fname+'\n')
        sys.exit()
    pos=0
    while i<len(lines):
        v=lines[i].split()
        for j in range(len(aaOrder)):
            prof[(pos,aaOrder[j])]=string.atof(v[j])/100.0    
 #               prof[(pos,aaOrder[j])]=string.atoi(v[j])
        pos+=1
        i+=1
    return prof,pos
    
#---------------------------------#

#---------------------------------#
def coilProb(cScore,params):
    ''' coilProb(curScore,params) returns coiled-coils probability dependin on params'''
    t1=1/params.sd_cc
    t2=(cScore-params.m_cc)/params.sd_cc
    t4=t2*t2
    Gcc=t1*math.exp(-0.5*t4)
    t1=1/params.sd_g
    t2=(cScore-params.m_g)/params.sd_g
    t4=t2*t2
    Gg=t1*math.exp(-0.5*t4)
    return  Gg,Gcc,Gcc/(params.sc*Gg + Gcc )
#---------------------------------#

#---------------------------------#
def seqScore(mat,seq,sPos,hept_pos):
    ''' seqScore(mat,aa,pos) comoputes the sequenc based score'''
    try:
        retVal=mat[(seq[sPos],hept_pos)]
    except: 
        retVal=1
    return retVal
#---------------------------------#


#---------------------------------#
def profScore(mat,prof,sPos,hept_pos):
    ''' profScore(mat,aa,pos) comoputes the sequenc based score'''
    sum=0.0
    for aa in Params.AAs:
        prof[(sPos,aa)]
        mat[(aa,hept_pos)]
        sum+=prof[(sPos,aa)]*mat[(aa,hept_pos)]
    return sum
#---------------------------------#

    
def pred_coil(seqr, params, fScore=None):
    ''' pred_coil(seq,seqLen,params,fScore) returns the coiled coil prediction of sequence seq'''
    seqr = copy.deepcopy(seqr)
    if fScore == None:
        fScore = seqScore
    seq = seqr.seq
    seqLen = len(seqr.seq)
    hept_pos=['a','b','c','d','e','f','g']
    score=[0.0]*seqLen
    hept_seq=['x']*seqLen
    for i in range(seqLen-params.win+1):
        this_score=1.0
        actual_win=0.0
        for j in range(min(params.win,seqLen-i)):
            pos=j%7
            actual_win+=params.pow[pos]
            this_score*=math.pow( fScore(params.mat,seq,i+j,pos), params.pow[pos] )
        if actual_win > 0:
            this_score=math.pow(this_score,1/actual_win)
        else:
            this_score=0.0
        for j in range(min(params.win,seqLen-i)): 
            pos=j%7
            if this_score > score[i+j]:
                score[i+j]=this_score
                hept_seq[i+j]=hept_pos[pos]
    for i in range(seqLen):
        gg, gcc, prob = coilProb(score[i],params)
        seqf = SeqFeature(location=FeatureLocation(i,i), type="pscoils")
        seqf.qualifiers = {'gg':gg, 
            'gcc': gcc, 
            'prob': prob, 
            'score': score[i],
            'hept_seq': hept_seq[i]}
        seqr.features.append(seqf)
    return seqr  
    
#---------------------------------#
def pred_coil2(seq,prof,params,ws=0.5):
    ''' pred_coil2(seq,prof,params,ws=0.5) returns the coiled coil prediction of sequence seq'''
    hept_pos=['a','b','c','d','e','f','g']
    seqLen=len(seq)
    score=[0.0]*seqLen
    prob=[0.0]*seqLen
    gcc=[0.0]*seqLen
    gg=[0.0]*seqLen
    hept_seq=['x']*seqLen
    for i in range(seqLen-params.win+1):
        this_score=1.0
        actual_win=0.0
        for j in range(min(params.win,seqLen-i)):
            pos=j%7
            actual_win+=params.pow[pos]
            avScore=(1-ws)*profScore(params.mat,prof,i+j,pos)+ws*seqScore(params.mat,seq,i+j,pos)
            this_score*=math.pow( avScore , params.pow[pos] )
        if actual_win > 0:
            this_score=math.pow(this_score,1/actual_win)
        else:
            this_score=0.0
        for j in range(min(params.win,seqLen-i)): 
            pos=j%7
            if this_score > score[i+j]:
                score[i+j]=this_score
                hept_seq[i+j]=hept_pos[pos]
    for i in range(seqLen):
        gg[i],gcc[i],prob[i]= coilProb(score[i],params)
    return gg,gcc,prob,hept_seq,score
  
            
