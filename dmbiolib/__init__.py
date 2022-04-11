#!/usr/bin/env python
__version__='0.1.0'
last_update='2022-04-11'
author='Damien Marsic, damien.marsic@aliyun.com'

import sys,glob,os,gzip,time,math

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

def check_format(x):
    format=''
    if x not in ('Single multipage pdf','svg','png','jpg','jpeg','pdf','ps','eps','pgf','raw','rgba','tif','tiff'):
        print("\n  File format not recognized. Options are svg, png, jpg, pdf, ps, eps, pgf, raw, rgba, tif, tiff (or single multipage pdf file with no -f argument).\n")
        sys.exit()
    if len(x)<5:
        format=x
    return format

def check_seq(filename,type,required,strict,limit):
    if filename[-3:]=='.gz':
        f=gzip.open(filename,'rt')
    else:
        f=open(filename,'r')
    t=True
    req=False
    count=0
    new=''
    for line in f:
        if limit and count==limit+1:
            break
        if type==aa:
            l=line.upper().strip()
        else:
            l=line.lower().strip()
        if l and l[0] in ('>','@'):
            count+=1
            new=l[0]
            continue
        if not new:
            continue
        if new=='@':
            new=''
        for i in l:
            if i not in type and not i.isdigit() and i!=' ':
                t=False
                break
            if i in required:
                req=True
        if t==False:
            if strict:
                break
            f.close()
            return False
    else:
        f.close()
        if req==True:
            return True
        if not strict and count>0:
            return False
        if count==0:
            print("\n  File "+filename+" does not contain any fasta or fastq sequence!\n")
            sys.exit()
    f.close()
    if limit and count==limit+1:
        return True
    if t==False:
        print("\n  File "+filename+" contains invalid characters!\n")
    elif req==False:
        print("\n  File "+filename+" does not contain any '"+required+"' nucleotide!\n")
    sys.exit()

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

def endplot(fig,counter,x,format,mppdf):
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

def findreadfiles(pattern,seqtype,exclude):
    if pattern=='auto':
        x='*.f*'
    elif check_file(pattern,False):
        x=pattern
    else:
        x=pattern+'*'
    rfiles=glob.glob(x)
    if len(rfiles)>1:
        rfiles=[x for x in rfiles if (len(x)>5 and x[-5:] in ('fasta','fastq','fa.gz','fq.gz')) or (len(x)>3 and x[-3:] in ('.fa','.fq')) or (len(x)>8 and x[-8:] in ('fasta.gz','fastq.gz'))]
    rfiles=[x for x in rfiles if check_seq(x,seqtype[0],seqtype[1],False,5)]
    if not rfiles:
        print("\n  No read file found. Check pattern, or move to correct directory, or include path in pattern.\n")
        sys.exit()
    for n in exclude:
        if n and n in rfiles:
            rfiles.remove(n)
    rfiles=dict.fromkeys(rfiles,0)
    for rfile in rfiles:
        if rfile[-3:]=='.gz':
            f=gzip.open(rfile,'r')
        elif rfile[-1]=='q':
            f=open(rfile,'rb')
        else:
            f=open(rfile,'r')
        if rfile[-3:]=='.gz' or rfile[-1]=='q':
            nr=lncount(f)//4
        else:
            nr=f.read().count('>')
        f.close()
        rfiles[rfile]=nr
    print('  Read file(s) (number of reads):')
    for n in rfiles:
        print('  '+n+' ('+str(rfiles[n])+')')
    print()
    return rfiles

def getread(f,y,z,counter):
    if not y:
        l=''
        while True:
            line=f.readline().strip()
            if not line:
                line=f.readline().strip()
            if not line or (l and line[0]=='>'):
                break
            if not l and line[0]=='>':
                continue
            l+=line
    else:
        while z<y:
            l=f.readline().strip()
            z+=1
    if not l:
        return l,f,z,counter
    l=l.lower()
    z=0
    counter+=1
    return l,f,z,counter

def initreadfile(rfile):
    if rfile[-3:]=='.gz':
        f=gzip.open(rfile,'rt')
    else:
        f=open(rfile,'r')
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
    z=y-2
    counter=0
    return f,y,z,counter

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

def override(func):
    class OverrideAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            func()
            parser.exit()
    return OverrideAction

def pr2(f,t):
    print(t)
    f.write(t+'\n')

def readfasta(filename,type):
    f=open(filename,'r')
    seq=1
    for line in f:
        if line and line[0]=='>':
            if seq!=1:
                print('\n  File '+filename+' contains more than one sequence !\n')
                sys.exit()
            seq=''
            continue
        if seq!=1:
            if type==aa:
                l=line.upper().strip()
            else:
                l=line.lower().strip()
            l=''.join(x for x in l if not x.isdigit() and x!=' ')
            seq+=l
    f.close()
    if not seq or seq==1:
        if type==aa:
            x='Protein'
        else:
            x='DNA'
        print('\n  '+x+' sequence could not be found in '+filename+' !\n')
        sys.exit()
    return seq

def readmultifasta(filename,type):
    f=open(filename,'r')
    temp=f.read().strip()
    f.close()
    seq={}
    if not '>' in temp:
        print('\n  '+filename+' is not a fasta file!\n\n')
        sys.exit()
    while temp[0]!='>':
        temp=temp[1:]
    temp=temp.split('>')
    for n in temp:
        x=n[:n.find('\n')]
        x=x[:max(len(x),x.find(' '))]
        y=n[n.find('\n'):]
        if x in seq:
            print('\n  Duplicate name '+x+' found in parental sequences! All sequences must have a different name!\n\n')
            sys.exit()
        y=y.replace(' ','').replace('\n','')
        if type==aa:
            y=y.upper()
        else:
            y=y.lower()
        if x:
            seq[x]=y
    return seq

def rename(name):
    if glob.glob(name):
        t=str(time.time())
        n=name[:name.rfind('.')]+'-'+t[:t.find('.')]+name[name.rfind('.'):]
        os.rename(name,n)
        print('\n  Existing '+name+' file was renamed as '+n+'\n  Creating new '+name+' file...\n')

def startplot(x,y):
    colors=plt.get_cmap(x)
    fig=plt.figure(figsize=(12,6.75))
    plt.title(y,size=15,weight='roman')
    return colors,fig

def version():
    print('\n  Project: '+sys.argv[0][max(sys.argv[0].rfind('\\'),sys.argv[0].rfind('/'))+1:-3]+'\n  version: '+__version__+'\n  Latest update: '+last_update+'\n  Author: '+author+'\n  License: GNU General Public v3 (GPLv3)\n')
