import math,itertools

def filter_seq(tol,seqs,method="log",method2="log",min1=None,max1=None,min2=None,max2=None):
    if len(seqs)==2:
        seq1,seq2=seqs
    elif len(seqs)==3:
        seq1,seq2,seq3=seqs
    else:
        raise Exception()
    import itertools
    last=None
    last2=None
    l1,l2,l3=[],[],[]
    if method=="log":
        seq1=filter(lambda x:x>0,seq1)
    if method2=="log":
        seq2=filter(lambda x:x>0,seq2)
    min1=min1 if min1!=None else min(seq1)
    max1=max1 if max1!=None else max(seq1)
    min2=min2 if min2!=None else min(seq2)
    max2=max2 if max2!=None else max(seq2)
    tot1=float(max1-min1) if method=="lin" else float(math.log(max1)-math.log(min1))
    tot2=float(max2-min2) if method2=="lin" else float(math.log(max2)-math.log(min2))
    for es in itertools.izip(*seqs):
        if len(es)==2:
            e1,e2=es
        elif len(es)==3:
            e1,e2,e3=es

        d1,d2=-1,-1
        if method=="log":
            #oktoprint1=last is None or last == 0 or abs(math.log(e1)-math.log(last))/tot1>tol
            if e2!=0:
                d1=1. if last is None or last == 0 else abs(math.log(e1)-math.log(last))/tot1
            else:
                d1=None
        elif method=="lin":
            #oktoprint1=last is None or last == 0 or abs(e1-last)/tot1>tol
            d1=1. if last is None or last == 0 else abs(e1-last)/tot1
        if method2=="log":
            #oktoprint2=last is None or last == 0 or abs(math.log(e2)-math.log(last2))/tot2>tol
            if e2!=0:
                d2=1. if last is None or last == 0 else abs(math.log(e2)-math.log(last2))/tot2
            else:
                d2=None
        elif method2=="lin":
            #oktoprint2=last is None or last == 0 or abs(e2-last2)/tot2>tol
            d2=1. if last is None or last == 0 else abs(e2-last2)/tot2>tol
        #if oktoprint1 or oktoprint2:
        if d1!=None and d2!=None and math.sqrt(d1*d1+d2*d2)>tol:
            l1.append(e1)
            l2.append(e2)
            if len(es)==3:
                l3.append(e3)
            last=e1
            last2=e2
    if last!=e1:
        l1.append(e1)
        l2.append(e2)
        if len(es)==3:
            l3.append(e3)
    if len(seqs)==2:
        return l1,l2
    elif len(seqs)==3:
        assert len(l1)==len(l2)==len(l3)
        return l1,l2,l3

def remove_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def plot_vector_field(xs,ys,f,alen=None):
    from matplotlib import pyplot as plt
    import numpy
    fx,fy=numpy.zeros((len(xs),len(ys))),numpy.zeros((len(xs),len(ys)))
    for i,x in enumerate(xs):
        for j,y in enumerate(ys):
            ff=f(x,y)
            if alen!=None:
                lf=math.sqrt(ff[0]**2+ff[1]**2)/alen
                ff=(ff[0]/lf,ff[1]/lf)

            fx[j,i]=ff[0]
            fy[j,i]=ff[1]

    plt.quiver(xs,ys,fx,fy)

colors_soft=["#7ac36a","#5a9bd4","#faa75b","#9e67ab","#737373","#f15a60","#ce7058","#d77fb4"]
colors_bright=["#008c48","#185aa9","#f47d23","#662c91","#010202","#ee2e2f","#a21d21","#b43894"]
