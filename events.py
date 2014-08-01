import scipy.stats
import random
import math

class KaplanMeierEstimator(object):
    def __init__(self):
        self.censor_times={}
        self.event_times={}

    def add_event(self, time, count=1):        
        self.event_times[time]=self.event_times.get(time,0)+count

    def add_events(self,times):
        for time in times:
            self.add_event(time)

    def add_censored(self, time, count=1):
        self.censor_times[time]=self.censor_times.get(time,0)+count

    def add_censoreds(self,times):
        for time in times:
            self.add_censored(time)

    def get_cumulative(self):
        t_set=set(self.event_times.keys()).union(set(self.censor_times.keys()))
        t_list=sorted(t_set)
        del t_set # just to save some memory

        events_cum,events_sum=[],0
        censored_cum,censored_sum=[],0
        for t in t_list:
            events_sum+=self.event_times.get(t,0)
            events_cum.append(events_sum)
            censored_sum+=self.censor_times.get(t,0)
            censored_cum.append(censored_sum)
        return t_list,events_cum,censored_cum
        

    def get_estimator(self):
        """Returns the Kaplan-Meier estimator as tuple (t,s).

        t : an ordered array of times
        s : the estimator value at each time
        """
        t_list,events_cum,censored_cum=self.get_cumulative()

        ndts=events_cum[-1]+censored_cum[-1]
        Si=[1]
        Ti=[0]
        for i,t in enumerate(t_list):
            if i==0:
                ni=ndts
            else:
                ni=ndts-events_cum[i-1]-censored_cum[i-1]
            di=self.event_times.get(t,0)
            if ni>0:
                if di>0:
                    Ti.append(t)
                    si=(ni-di)/float(ni)
                    Si_prev= Si[-1] if len(Si)>0 else 1
                    Si.append(Si_prev*si)                

        return Ti,Si

    def _get_naive_estimator(self):
        """For testing."""
        t_list=sorted(self.event_times.keys())

        events_cum,events_sum=[],0
        for t in t_list:
            events_sum+=self.event_times[t]#*(self.endTime/float(self.endTime-t))
            events_cum.append(events_sum)
            
        ts,ns=[0],[1]
        for t,n in zip(t_list,events_cum):
            ts.append(t)
            ns.append(1-n/float(events_sum))
        return ts,ns


    def get_naive_estimator(self):
        """Returns an estimator where the sensored events are simply discarded.
        """
        t_list,events_cum,censored_cum=self.get_cumulative()

        ndts=events_cum[-1]
        Si=[1]
        Ti=[0]
        for i,t in enumerate(t_list):
            if i==0:
                ni=ndts
            else:
                ni=ndts-events_cum[i-1]
            di=self.event_times.get(t,0)
            if ni>0:
                if di>0:
                    Ti.append(t)
                    si=(ni-di)/float(ni)
                    Si_prev= Si[-1] if len(Si)>0 else 1
                    Si.append(Si_prev*si)                

        return Ti,Si

class IntereventTimeEstimator(KaplanMeierEstimator):
    def __init__(self,endTime,mode='censorlast'):
        """ Constructor.
        
        Parameters
        ----------
        endTime : float
           The last time point in the observation period
        mode : string
           'censorlast' : the last iet is censored when the observation period ends
           'censorall' : in addition to the last iet, the first one is censored
           'periodic' : periodic boundary conditions
        
        """
        super(IntereventTimeEstimator, self).__init__()
        self.endTime=endTime
        assert mode in ["censorlast","censorall","periodic"]
        self.mode=mode

    def add_time_seq(self,seq):
        if len(seq)!=0:
            for i,time in enumerate(seq):
                if i!=0:
                    dt=time-last
                    assert dt>=0
                    self.add_event(dt)
                    if self.mode=='censorall':
                        self.add_event(dt)
                elif self.mode=='censorall':
                    self.add_censored(time)                                  
                elif self.mode=='periodic':
                    firstTime=time
                last=time
            if self.mode=='periodic':
                self.add_event(firstTime+self.endTime-last)
            else:
                self.add_censored(self.endTime-last)


    def add_time_duration_seq(self,seq):
        if len(seq)!=0:
            for i,(time,duration) in enumerate(seq):
                if i!=0:
                    dt=time-last
                    assert dt>=0
                    self.add_event(dt)
                    if self.mode=='censorall':
                        self.add_event(dt)
                elif self.mode=='censorall':
                    self.add_censored(time)                                  
                elif self.mode=='periodic':
                    firstTime=time
                last=time+duration
            if self.mode=='periodic':
                self.add_event(firstTime+self.endTime-last)
            else:
                self.add_censored(self.endTime-last)
                if self.mode=='censorall':
                    self.add_censored(self.endTime-last)


    def get_length_bias_estimator(self):
        """Returns an estimator correcting for the length bias of windows cencored data. 

        This estimator only takes into account the non-censored data. For the non-censored
        data the data points are distributed as p'(\tau) \sim (T-t)p(\tau), where p'(\tau)
        is the observed distribution of the inter-event times (i.e., the uncensored 
        distribution), p(\tau) is the original unobserved distribution which produced the
        underlying event sequence, and T is the time window length.
        """
        t_list=sorted(self.event_times.keys())

        events_cum,events_sum=[],0
        for t in t_list:
            events_sum+=self.event_times[t]*(self.endTime/float(self.endTime-t))
            events_cum.append(events_sum)

        #t_list,events_cum,censored_cum=self.get_cumulative()

        ts,ns=[0],[1]
        for t,n in zip(t_list,events_cum):
            ts.append(t)
            ns.append(1-n/float(events_sum))
        return ts,ns


def edgestotimeseqs(edges, issorted=True):
    """Generator transforming edges to time sequences.
    """
    if not issorted:
        edges=sorted(edges)
    current=None
    for event in edges:
        if (event[0],event[1]) != current:
            if current!=None:
                l.sort()
                yield l
            l=[]
        l.append(event[2])
        current=(event[0],event[1])
    l.sort()
    yield l


def iets(events):
    """Generator for inter-event times.
    """
    for i,event in enumerate(events):
        if i!=0:
            yield event-lastevent
        lastevent=event

def normalize(events,form="timeseqs"):
    """Normalizes times in an event list in place.

    Normalization is done such that the event times are between 0 and 1. For 
    edge the network is in addition made undirected in a way that the smaller
    node in the edge is always given before the larger one. E.g., (2,1,t) is
    transformed into (1,2,t).
    """
    mint,maxt=None,None

    if form=="edges":
        for fr,to,t in events:
            if maxt==None or t>maxt:
                maxt=t
            if mint==None or t<mint:
                mint=t
        for i in range(len(events)):
            events[i][2]=events[i][2]-mint
            if events[i][0]>events[i][1]:
                temp=events[i][0]
                events[i][0]=events[i][1]
                events[i][1]=temp
    if form=="timeseqs":
        for timeseq in events:
            for t in timeseq:
                if maxt==None or t>maxt:
                    maxt=t
                if mint==None or t<mint:
                    mint=t
        for timeseq in events:
            for i,t in enumerate(timeseq):
                timeseq[i]=timeseq[i]-mint

    return maxt-mint


def random_timeseq(tdist,trdist,endtime,starttime=0):
    """ Generate a random time sequence.

    Parameters
    ----------
    tdist : function
      A function returning realizations from the inter-event time distribution
    trdist : function
      A function returning realizations from the residual waiting time distribution
    endtime : float
      The end time of the observation period
    starttime : float
      The start time of the observation period

    Returns
    -------
    A sequence of times of the events.
    """
    l=[]
    t1=trdist()
    if t1+starttime>endtime:
        return l
    else:
        l.append(t1+starttime)
    while True:
        ti=tdist()
        if l[-1]+ti>endtime:
            break
        l.append(l[-1]+ti)

    return l

def random_timeseq_burnin(tdist,endtime,starttime=0,burninfactor=10):
    t=0
    dt=endtime-starttime
    bt=dt*burninfactor
    l=[]
    while True:
        ti=tdist()
        t=t+ti
        if t>bt:
            if starttime+t-bt > endtime:
                return l
            l.append(starttime+t-bt)
            break
    while True:
        ti=tdist()
        if l[-1]+ti>endtime:
            break
        l.append(l[-1]+ti)

    return l

def random_timeseq_empirical(tdist_cum,trdist_cum,T):
    """Returns a random inter-even time distribution given cumulative distributions of 
    inter-event times and residual inter-event times. If the cumulative distribution of 
    inter-event times doesn't go to one all the rest are assumed to be larger than T.
    """
    pass


def exprv(rate):
    p=random.random()
    return -math.log(p)/float(rate)

def random_timeseq_exp(rate,starttime,endtime):
    #return random_timeseq(lambda :scipy.stats.expon.rvs(rate),lambda :scipy.stats.expon.rvs(rate),endtime,starttime)
    return random_timeseq(lambda :exprv(rate),lambda :exprv(rate),endtime,starttime)

def plaw(exp,mint=1.):
    """Generate a value from power-law distribution.

    Probability density is:
    p(\tau) = 0, when \tau < \tau_m
    p(\tau) = (\alpha - 1) \tau_m^{\alpha - 1} \tau^{-\alpha}, \tau > \tau_m    
    """
    p=random.random()
    return mint * (1. - p)**(1./float(1.-exp))

def plaw_residual(exp,mint=1.):
    """Generate a value from the distribution of the residual of power-law iet.

    Probability density is:
    p(\tau_R) = \frac{\alpha - 2}{\alpha - 1} \tau_m^{-1}, when \tau_R < \tau_m
    and
    p(\tau_R) = \frac{\alpha - 2}{\alpha - 1} \tau_m^{\alpha - 2} * \tau_R^{1-\alpha}, when \tau_R > \tau_m

    Cumulative probability distribution is:
    P(\tau_R) = \frac{\alpha - 2}{\alpha -1} \tau_m^{-1} \tau_R, when \tau_R < \tau_m
    and
    P(\tau_R) = 1 - \frac{\tau_m^{\alpha - 2}}{\alpha - 1} \tau_R^{2 - \alpha}, when \tau_R > \tau_m
    
    Inverse cumulative probability distribution:
    P^{-1}(p) = \frac{\alpha - 1}{\alpha - 2} \tau_m p, when p < \frac{\alpha - 2}{\alpha - 1}
    and
    P^{-1}(p) = \tau_m ((\alpha - 1)(1-p))^{frac{1}{2 - \alpha}}, when p < \frac{\alpha - 2}{\alpha - 1}
    """
    p=random.random()
    if p < float(exp - 2.)/float(exp - 1.):
        return float(exp-1.)/float(exp-2.)*mint*p
    else:
        return mint*((exp -1.)*(1.-p))**(1./float(2.-exp))

def random_timeseq_plaw(exp,mint,endtime,starttime=0,burnin=None):
    #lambda :scipy.stats.pareto.rvs(exp)
    if burnin==None:
        return random_timeseq(lambda :plaw(exp,mint),lambda :plaw_residual(exp,mint),endtime,starttime)
    else:
        return random_timeseq_burnin(lambda :plaw(exp,mint),endtime,starttime,burninfactor=burnin)



if __name__=="__main__":
    km=KaplanMeierEstimator()
    km.add_event(3)
    km.add_event(11,2)
    km.add_censored(9)
    km.add_censored(12,6)
    t,s=km.get_estimator()
    print t
    print s


    tn,sn=km.get_naive_estimator()
    print tn
    print sn

    print km._get_naive_estimator()
