#!/usr/bin/env python
__version__='0.2.6'
last_update='2022-06-14'
author='Damien Marsic, damien.marsic@aliyun.com'

import sys,glob,os,gzip,time,math,argparse
import numpy as np

ambiguous="ryswkmbdhvn"
aa="ARNDCQEGHILKMFPSTWYV"
IUPAC=[set('ag'),set('ct'),set('gc'),set('at'),set('gt'),set('ac'),set('cgt'),set('agt'),set('act'),set('acg'),set('atgc')]

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
        y=n[:n.find('\n')]
        y=y[:max(len(y),y.find(' '))]
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
            z+=1
    if not l:
        return l,f,counter,name
    counter+=1
    return l.lower(),f,counter,name

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
    counter=0
    return f,y,counter

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

def override(func):
    class OverrideAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            func()
            parser.exit()
    return OverrideAction

def plot_end(fig,counter,x,format,mppdf):
    plt.legend()
    fig.subplots_adjust(bottom=0.15)
    fig.tight_layout()
    if format:
        counter+=1
        g=x+str(counter)+'.'+format
        plt.savefig(g,dpi=600)
        print('  Figure was saved into file: '+g+'\n')
    else:
        mppdf.savefig()
    plt.close()
    return counter

def plot_start(x,y):
    colors=plt.get_cmap(x)
    fig=plt.figure(figsize=(12,6.75))
    plt.title(y,size=15,weight='roman')
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
    if glob.glob(name):
        t=str(time.time())
        n=name[:name.rfind('.')]+'-'+t[:t.find('.')]+name[name.rfind('.'):]
        os.rename(name,n)
        print('\n  Existing '+name+' file was renamed as '+n+'\n  Creating new '+name+' file...\n')

def shortest_probe(seqs,lim,t):
    if lim<1:
        lim=1
    fail=''
    q=-1
    x=min([len(k) for k in seqs])
    y=set([k[-x:] for k in seqs])
    if len(y)!=len(seqs):
        fail='\n  Duplicate '+t+' found! '+t[-1:].upper()+t[1:]+'s must all be different when trimmed to their maximal common size!'
    if not fail:
        q=lim
        while len(y)!=len(seqs) or len(set([k.count(p) for k in seqs for p in y]))>1:
            q+=1
            y=set([k[-q:] for k in seqs])
    return q,fail

def version():
    print('\n  Project: '+sys.argv[0][max(sys.argv[0].rfind('\\'),sys.argv[0].rfind('/'))+1:-3]+'\n  version: '+__version__+'\n  Latest update: '+last_update+'\n  Author: '+author+'\n  License: GNU General Public v3 (GPLv3)\n')
