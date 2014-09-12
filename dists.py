import collections

class Dist(collections.MutableMapping):
    def __init__(self):
        self.d={}
        self.s=0
    def __getitem__(self,key):
        return self.d.get(key,0)
    def __iter__(self):
        return self.d.__iter__()
    def __len__(self):
        return len(self.d)

    def __setitem__(self,key,val):
        assert val>=0
        self.s=self.s - self[key] + val
        self.d[key]=val
    def __delitem__(self,key):
        del self.d[key]

    def moment(self,n):
        if self.s==0:
            return 0
        m=0.
        for key,val in self.iteritems():
            m+=val*key**n
        return m/float(self.s)

    def get_cumdist(self,start=0,normalized=True):
        assert start not in self.d
        if normalized:
            p=[1.]
        else:
            p=[self.s]
        vals=[start]+sorted(self.d.keys())
        if normalized:
            for val in vals[1:]:
                p.append(p[-1]-self.d[val]/float(self.s))
        else:
            for val in vals[1:]:
                p.append(p[-1]-self.d[val])
        return CumDist(vals,p)

class CumDist(object):
    def __init__(self,vals,p):
        self.vals=vals
        self.p=p



def get_moment(vals,probs,T,m):
    s=0
    for i in range(1,len(vals)):
        p=probs[i-1]-probs[i]
        v=vals[i]
        s+=p*v**m
    s+=probs[i]*T**m
    return s
