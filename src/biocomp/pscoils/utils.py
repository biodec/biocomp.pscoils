from Bio import SeqIO
import sys

def readFasta(fname):
    handle = open(fname)
    dict = {}
    for seq_record in SeqIO.parse(handle, "fasta") :
        dict[seq_record.id] = seq_record.seq
    handle.close()
    return dict

def printCoil(seqr,labels='T',io=sys.stdout):
    ''' '''
    # ...
    seq = seqr.seq
    seqLen = len(seqr.seq)
    score=[0.0]*seqLen
    prob=[0.0]*seqLen
    gcc=[0.0]*seqLen
    gg=[0.0]*seqLen
    hept_seq=['x']*seqLen
    for f in seqr.features:
        for i in range(f.location.start.position,f.location.end.position+1):
            gg[i] = f.qualifiers['gg']
            gcc[i] = f.qualifiers['gcc']
            prob[i] = f.qualifiers['prob']
            hept_seq[i] = f.qualifiers['hept_seq']
            score[i] = f.qualifiers['score']
    P = prob
    Gcc = gcc
    Gg = gg
    # ...
    ccLab="C"
    loopLab="L"
    if labels !='T' :
        for i in range(seqLen):
            #print "%4d %c %c %7.3f %7.3f (%7.3f %7.3f)" % (i+1,seq[i],hept_seq[i],score[i],P[i],Gcc[i],Gg[i])
            io.write("%4d %c %c %7.3f %7.3f %7.3f %7.3f\n" % (i+1,seq[i],hept_seq[i],score[i],P[i],Gcc[i],Gg[i]))
    else:
        io.write(" Pos A Hep Score   Prob    Gcc     Gg    Pred (Loop=L Coiledcoil=C)\n")
        for i in range(seqLen):
            if P[i] > 0.5:
               clabel=ccLab
            else:
               clabel=loopLab
            io.write("%4d %c %c %7.3f %7.3f %7.3f %7.3f %s\n" % (i+1,seq[i],hept_seq[i],score[i],P[i],Gcc[i],Gg[i],clabel))
    
