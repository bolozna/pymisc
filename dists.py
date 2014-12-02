import collections,itertools,random

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

    def get_moment(self,n):
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
        p[-1]=0
        return CumDist(vals,p)

class CumDist(object):
    """Cumulative probability distribution.

    P_{<}(t)
    """
    def __init__(self,vals,ps,maxval=None):
        self.vals=vals
        self.ps=ps
        self.maxval=maxval

    def rescale(self,newavg):
        """Returns a rescaled distribution

        Parameters
        ----------
        newavg (float) : The average of the new distribution  
        """
        new_vals=[]
        c=newavg/float(self.get_moment(1))
        for val in self.vals:
            new_vals.append(c*val)
        return CumDist(new_vals,copy(self.ps))

    def get_moment(self,moment):
        s=0
        for i in range(1,len(self.ps)):
           s+=(self.ps[i-1]-self.ps[i])*(self.vals[i]**moment)
        if self.ps[-1]!=0:
            s+=self.ps[-1]*(self.maxval**moment)
        return s

    def sample(self):
        u=random.random()
        i=bisect(dist.values,u)
        return keys[i] #+1 -1 ???

    def get_residual_waiting_dist(self):
        """Residual waiting time distribution.

        Method returning the residual waiting time distribution for a renewall
        process that is produced assuming that this distribution is the inter-
        event time distribution.

        The cumulative distribution of the residual waiting time distribution
        can be calculated with the following formula:
        G(\tau_R) = 1-\int_0^{\tau_R} frac{F(\tau -)}{\tau_F} d\tau,
        where F is the inter-event time distribution and \tau_F is the average of F.        

        Reference: Feller 1971 Chapter XI
        """
        pass


def get_moment(vals,probs,T,m):
    s=0
    for i in range(1,len(vals)):
        p=probs[i-1]-probs[i]
        v=vals[i]
        s+=p*v**m
    s+=probs[i]*T**m
    return s


if __name__=="__main__":
    d=Dist()
    d[10]=1
    d[100]=1
    dc=d.get_cumdist()
    assert dc.get_moment(1)==d.get_moment(1)
    assert dc.get_moment(2)==d.get_moment(2)
