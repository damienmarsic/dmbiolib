#!/usr/bin/env python
__version__='0.2.44'
last_update='2022-08-23'
author='Damien Marsic, damien.marsic@aliyun.com'

import sys,os,gzip,time,math
from glob import glob
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt

ambiguous="ryswkmbdhvn"
aa="ARNDCQEGHILKMFPSTWYV"
IUPAC=[set('ag'),set('ct'),set('gc'),set('at'),set('gt'),set('ac'),set('cgt'),set('agt'),set('act'),set('acg'),set('atgc')]
gcode={
'ggg':'G','gga':'G','ggc':'G','ggt':'G','gag':'E','gaa':'E','gac':'D','gat':'D','gcg':'A', 'gca':'A', 'gcc':'A', 'gct':'A', 'gtg':'V', 'gta':'V', 'gtc':'V', 'gtt':'V',
'agg':'R','aga':'R','agc':'S','agt':'S','aag':'K','aaa':'K','aac':'N','aat':'N','acg':'T', 'aca':'T', 'acc':'T', 'act':'T', 'atg':'M', 'ata':'I', 'atc':'I', 'att':'I',
'cgg':'R','cga':'R','cgc':'R','cgt':'R','cag':'Q','caa':'Q','cac':'H','cat':'H','ccg':'P', 'cca':'P', 'ccc':'P', 'cct':'P', 'ctg':'L', 'cta':'L', 'ctc':'L', 'ctt':'L',
'tgg':'W','tga':'*','tgc':'C','tgt':'C','tag':'*','taa':'*','tac':'Y','tat':'Y','tcg':'S', 'tca':'S', 'tcc':'S', 'tct':'S', 'ttg':'L', 'tta':'L', 'ttc':'F', 'ttt':'F'}
bpairs={'a':'t','c':'g','g':'c','t':'a'}

def aln2seq(fn,type,full,ref):
    fail=''
    seq={}
    if ref:
        ref,fail=getfasta(ref,aa,aa,False)
        if fail:
            return seq,fail
        ref=list(ref.values())[0]
    with open(fn,'r') as f:
        x=f.read().strip().split('\n')
    x=[k for k in x if k]
    for i in range(len(x)):
        if x[i].replace(' ','').isdigit():
            q=x[i].split()
            break
    else:
        fail+='\n  File format not recognized!'
    q=[int(k) for k in q]
    x=x[i+1:]
    if not x or '.' in x[0]:
        fail+='\n  No reference sequence!'
    else:
        w=x[0].split()[1:]
        wt=x[0].split()[0]
    x=x[1:]
    if len(q)!=len(w):
        fail+='\n  Number of indices is different from number of regions!'
    if not x:
        fail+='\n  No sequence found!'
    for n in x:
        l=n.split()
        y=l[0]
        z=l[1:]
        if y in seq:
            fail+='\n  Duplicate sequence name '+y+' found in file '+fn+'! All sequences must have a different name!'
        if len(z)!=len(w):
            fail+='\n  Number of regions of sequence '+y+' is different from that of reference!'
            continue
        for i in range(len(z)):
            t=''
            for j in range(len(z[i])):
                if z[i][j]=='.':
                    t+=w[i][j]
                else:
                    t+=z[i][j]
            z[i]=t
        if full:
            t=''
            t=ref[:q[0]]
            for i in range(len(q)):
                t+=z[i]
                if i==len(q)-1:
                    a=len(ref)
                else:
                    a=q[i+1]
                t+=ref[q[i]+len(z[i]):a]
            z=t
        else:
            z=','.join(z)
        if type==aa:
            z=z.upper()
        else:
            z=z.lower()
        t,req=check_seq(z.replace(',',''),type,type)
        if not t:
            fail+='\n  Sequence '+y+' in '+fn+' contains invalid characters!'
        if not req:
            fail+='\n  Sequence '+y+' in '+fn+' does not contain expected characters!'
        if y and not fail:
            seq[y]=z
    return seq,fail,wt

def check_file(filename,strict):
    try:
        f=open(filename,'r')
    except IOError:
        if strict:
            print("\n  File "+filename+" could not be found!\n")
            sys.exit()
        else:
            return False
    else:
        return True

def check_plot_format(x):
    format=''
    y=('svg','png','jpg','jpeg','pdf','ps','eps','pgf','raw','rgba','tif','tiff')
    if x not in ('Single multipage pdf',)+y:
        print("\n  File format not recognized. Options are: "+', '.join(y)+" (or single multipage pdf file with no -f argument).\n")
        sys.exit()
    if len(x)<5:
        format=x
    return format

def check_read_file(x):
    fail=''
    if check_file(x,False):
        f=open_read_file(x)
        for i in range(5):
            y=f.readline().strip()
            if y:
                break
        if y and y[0] in ('@','>'):
            y=f.readline().strip().lower()
            t,req=check_seq(y,'atgc'+ambiguous,'atgc')
            if not t:
                fail+='\n  File '+x+' contains invalid characters!'
            if not req:
                fail+='\n  File '+x+' does not contain expected characters!'
        else:
            fail+='\n  File '+x+' has incorrect format!'
        f.close()
    else:
        fail+='\n  File '+x+' not found!'
    return fail

def check_seq(seq,type,required):
    t=True
    req=False
    seq=seq.lower()
    type=type.lower()
    required=required.lower()
    for i in seq:
        if i not in type:
            t=False
            break
        if i in required:
            req=True
    return t,req

def check_sync(l1,l2):
    fail=''
    if l1[0] not in ('@','>') or l2[0] not in ('@','>'):
        fail='\n  Wrong read file format! Only fastq, gzipped fastq or fasta are accepted formats!'
    elif l1[:l1.find(' ')]!=l2[:l2.find(' ')]:
        fail='\n  R1 and R2 reads are not synchronized! Please use raw read files!'
    return fail

def complexity(seq):
    i=0
    t=[]
    while i+3<=len(seq):
        t.append(defaultdict(int))
        x=seq[i:i+3]
        i+=3
        c=[[],[],[]]
        for n in range(3):
            if x[n] in ambiguous:
                for m in IUPAC[ambiguous.index(x[n])]:
                    c[n].append(m)
            else:
                c[n].append(x[n])
        for n1 in c[0]:
            for n2 in c[1]:
                for n3 in c[2]:
                    t[-1][gcode[n1+n2+n3]]+=1
    return t

def compress(seq):
    x=''
    for n in seq:
        if x and n==x[-1]:
            continue
        x+=n
    return x

def diff(seqs):
    z=len(seqs[0])
    x=[(seqs[i],seqs[j]) for i in range(len(seqs)) for j in range(i+1, len(seqs))]
    for (a,b) in x:
        y=0
        for i in range(len(a)):
            if a[i]!=b[i]:
                y+=1
        z=min(z,y)
    return z

def entropy(matrix):
    score=0
    for n in matrix:
        H=0
        for m in n:
            if m!=0:
                H+=m*(math.log(m,2))
        H*=-1
        score+=H
    return score

def exprange(a,b,c):
    while a<b:
        yield a
        a*=c

def find_read_files():
    rfiles=glob('*.f*q.gz')
    if not rfiles:
        return {}
    x=defaultdict(int)
    for n in rfiles:
        for m in ('_R1','_R2','_1.','_2.'):
            if m in n:
                x[m]+=1
    for (a,b) in (('_R1','_R2'),('_1.','_2.')):
        if x[a]>0 and x[a]==x[b] and x[a]==len(rfiles)/2:
            break
    else:
        a,b='',''
    y=[n for n in rfiles]
    if a:
        y=[n.replace(a,'*') for n in rfiles if a in n]
    q=-1
    for i in range(1,len(y[0])):
        if len(set([k[-i:] for k in y]))>1:
            while True:
                i-=1
                if y[0][-i] in ('-','_','.'):
                    break
            q=-i
            break
    for i in range(len(y)):
        x=y[i]
        if a:
            x=y[i].replace('*',a)+' '+y[i].replace('*',b)
        y[i]=y[i][:q]+' '+x
    y.sort()
    return y

def fsize(filename):
    return os.path.getsize(filename)

def getfasta(fn,type,required,multi):
    fail=''
    seq={}
    if not check_file(fn,False):
        fail+='\n  File '+fn+' could not be found!'
        return seq,fail
    with open(fn,'r') as f:
        x=f.read().strip()
    if x[0]!='>':
        fail+='\n  '+fn+' is not a fasta file!'
        return seq,fail
    if not multi and x.count('>')>1:
        fail+='\n  File '+fn+' contains more than one sequence!'
        return seq,fail
    x=x.split('>')
    for n in x:
        if not n:
            continue
        y=n[:n.find('\n')]+' '
        y=y[:y.find(' ')]
        if not y and multi:
            fail+='\n  Sequence with no name found in file '+fn+'! All sequences must have a name!'
        z=n[n.find('\n'):]
        if y in seq:
            fail+='\n  Duplicate sequence name '+y+' found in file '+fn+'! All sequences must have a different name!'
        z=z.replace(' ','').replace('\n','')
        if type==aa:
            z=z.upper()
        else:
            z=z.lower()
        t,req=check_seq(z,type,required)
        if not t:
            fail+='\n  Sequence '+y+' in '+fn+' contains invalid characters!'
        if not req:
            fail+='\n  Sequence '+y+' in '+fn+' does not contain expected characters!'
        if y and not fail:
            seq[y]=z
    return seq,fail

def getread(f,y,counter):
    name=''
    if not y:
        l=''
        while True:
            line=f.readline().strip()
            if not line:
                line=f.readline().strip()
            if line[0]=='>':
                name=line
            if not line or (l and line[0]=='>'):
                break
            if not l and line[0]=='>':
                continue
            l+=line
    else:
        z=0
        while z<y:
            l=f.readline().strip()
            if not z:
                name=l
            if z==1:
                seq=l.lower()
            z+=1
    if not l:
        return l,f,counter,name
    counter+=1
    return seq,f,counter,name

def initreadfile(rfile):
    f=open_read_file(rfile)
    l=f.readline().strip()
    if not l or l[0] not in ('>','@'):
        f.close()
        print('  '+rfile+' does not look like a fastq or fasta file.\n')
        sys.exit()
    if l[0]=='>':
        l=f.readline().strip()
        l=f.readline().strip()
        if l and l[0]=='>':
            y=2
        else:
            y=0
    else:
        y=4
    f.seek(0)
    return f,y

def lncount(f):
    def _make_gen(reader):
        b=reader(1024*1024)
        while b:
            yield b
            b=reader(1024*1024)
    f_gen=_make_gen(f.read)
    return sum(buf.count(b'\n') for buf in f_gen)

def match(text1, text2):
    if len(text1)==len(text2):
        for l in range(len(text1)):
            if not nt_match(text1[l], text2[l]):
                return False
        return True
    else:
        return False

def nt_match(nt1, nt2):
    if nt1==nt2:
        return True
    elif (nt2 in ambiguous and nt1 in IUPAC[ambiguous.index(nt2)]) or (nt1 in ambiguous and nt2 in IUPAC[ambiguous.index(nt1)]):
        return True
    else:
        return False

def open_read_file(x):
    if x[-2:]=='gz':
        f=gzip.open(x,'rt')
    else:
        f=open(x,'r')
    return f

def plot_end(fig,name,format,mppdf):
    plt.legend()
    fig.subplots_adjust(bottom=0.15)
    fig.tight_layout()
    if format:
        g=name+'.'+format
        plt.savefig(g,dpi=600)
        print('  Figure was saved into file: '+g+'\n')
    else:
        mppdf.savefig()
    plt.close()

def plot_start(x,y,z):
    colors=plt.get_cmap(x,y)
    fig=plt.figure(figsize=(12,6.75))
    plt.title(z,size=15,weight='roman')
    return colors,fig

def pr2(f,t):
    print(t)
    f.write(t+'\n')

def progress_check(c,show,t):
    if c in show:
        k=show[c]
        print('\r  '+t+' '*(8-len(k))+k+'%',end='')

def progress_end():
    print('\b\b\b\b\b\b100.0%\n')

def progress_start(nr,t):
    y=np.arange(0,nr,nr/1000)
    x=[round(n) for n in y]
    z=np.arange(0,100,0.1)
    y=[str(round(n,1)) for n in z]
    show=dict(zip(x,y))
    print('  '+t+'     0.0%',end='')
    return show

def readcount(R,fail):
    if R[-3:]=='.gz':
        f=gzip.open(R,'r')
    elif R[-1]=='q':
        f=open(R,'rb')
    else:
        f=open(R,'r')
    if R[-3:]=='.gz' or R[-1]=='q':
        nr=lncount(f)//4
    else:
        nr=f.read().count('>')
    f.close()
    if nr<2:
        fail+='\n  File '+R+' is not a read file!'
    return(nr,fail)

def rename(name):
    if glob(name):
        t=str(time.time())
        n=name[:name.rfind('.')]+'-'+t[:t.find('.')]+name[name.rfind('.'):]
        os.rename(name,n)
        print('\n  Existing '+name+' file was renamed as '+n+'\n  Creating new '+name+' file...\n')

def revcomp(seq):
    rs=(seq[::-1]).lower()
    x=''.join([bpairs.get(rs[i], 'X') for i in range(len(seq))])
    return x

def shortest_probe(seqs,lim,host,t):
    if lim<1:
        lim=1
    fail=''
    q=-1
    x=min([len(k) for k in seqs])
    y=set([k[-x:] for k in seqs])
    if len(y)!=len(seqs):
        fail='\n  Duplicate '+t+' found! '+t[0].upper()+t[1:]+'s must all be different when trimmed to their maximal common size!'
    if host and len([k for k in y if k in host+host[:x-1]]):
        fail+='\n  '+t[0].upper()+t[1:]+' found in the host genome!'
    if not fail:
        q=lim
        while True:
            y=set([k[-q:] for k in seqs])
            if len(y)==len(seqs) and max([k.count(p) for k in seqs for p in y])==1 and not len([k for k in y if k in host+host[:q-1]]):
                break
            q+=1
    return q,fail

def sortfiles(x):
    x.sort()
    y=[k for k in x if not any(i.isdigit() for i in k[:k.rfind('--')])]
    x=[k for k in x if not k in y]
    z=[(int(''.join([n for n in k if n.isdigit()])),k) for k in x]
    for n in sorted(z):
        y.append(n[1])
    return y
    
def transl(seq):
    seq=seq.lower()
    x=''.join([gcode.get(seq[3*i:3*i+3],'X') for i in range(len(seq)//3)])
    return x













