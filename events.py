import scipy.stats
import scipy.optimize
from scipy.special import binom
import random
import math
import itertools
import dists

class KaplanMeierEstimator(object):
    def __init__(self,censored_times=None,event_times=None):
        if censored_times==None:
            self.censor_times={}
        else:
            self.censor_times=censored_times
        if event_times==None:
            self.event_times={}
        else:
            self.event_times=event_times

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
        

    def get_estimator(self,variance=False,correctionFactor=1.0):
        """Returns the Kaplan-Meier estimator as tuple (t,s).

        t : an ordered array of times
        s : the estimator value at each time
        v : Variance of the estimator (Greenwood's formula)
        """
        t_list,events_cum,censored_cum=self.get_cumulative()

        ndts=events_cum[-1]+censored_cum[-1]
        Si=[1]
        Ti=[0]

        Vmod=0
        Vi=[1*Vmod]

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
                    if (ni-di)!=0:
                        Vmod+=di/float(ni*(ni-di))
                        Vi.append(correctionFactor*Si[-1]**2 * Vmod)
                    else:
                        Vi.append(None)

        if variance:
            return Ti,Si,Vi
        else:
            return Ti,Si

    @staticmethod
    def get_confidence_intervals(Si,Vi,confLevel,logcorrect=False,correctionFactor=1.0):
        assert len(Si)==len(Vi)
        Si_up=[1]
        Si_down=[1]
        z=scipy.stats.norm.ppf(1-(1.-confLevel)/2.)
        if logcorrect:
            def loginv(w):
                return math.exp(w)/float(1+math.exp(w))
            for i in range(1,len(Si)):
                w=math.log(Si[i]/(1.0-Si[i]))
                wvar=(1./float(Si[i]*(1-Si[i])))**2*correctionFactor*Vi[i]
                Si_up.append(loginv(w+z*math.sqrt(wvar)))
                Si_down.append(loginv(w-z*math.sqrt(wvar)))
        else:
            for i in range(1,len(Si)):
                Si_up.append(Si[i]+z*math.sqrt(correctionFactor*Vi[i]))
                Si_down.append(Si[i]-z*math.sqrt(correctionFactor*Vi[i]))
        return Si_up,Si_down

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

class IntereventTimeEstimator(object):
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
        #super(IntereventTimeEstimator, self).__init__()
        self.endTime=endTime
        assert mode in ["censorlast","censorall","periodic"]
        self.mode=mode

        self.observed_iets={}
        self.forward_censored_iets={}
        self.backward_censored_iets={}
        self.empty_seqs=0
        self.nseqs=0

    def add_time_seq(self,seq):
        self.nseqs+=1
        if len(seq)!=0:
            for i,time in enumerate(seq):
                if i!=0:
                    dt=time-last
                    assert dt>=0
                    self.observed_iets[dt]=self.observed_iets.get(dt,0)+1
                elif self.mode=='censorall':
                    self.backward_censored_iets[time]=self.backward_censored_iets.get(time,0)+1
                elif self.mode=='periodic':
                    firstTime=time
                last=time
            if self.mode=='periodic':
                dt=firstTime+self.endTime-last
                self.observed_iets[dt]=self.observed_iets.get(dt,0)+1
            else:
                dt=self.endTime-last
                self.forward_censored_iets[dt]=self.forward_censored_iets.get(dt,0)+1
        else:
            self.empty_seqs+=1

    def read_seqs(self,filename):
        """Reads a file containing event sequences.
        """
        f=open(filename,'r')
        for line in f:
            seq=map(float,line.split())
            self.add_time_seq(seq)

    """                
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
    """

    def get_estimator(self,variance=False):
        if self.mode=='periodic':
            censored_times=None
            event_times=self.observed_iets
        elif self.mode=='censorlast':
            censored_times=self.forward_censored_iets
            event_times=self.observed_iets
        elif self.mode=='censorall':
            censored_times=self.backward_censored_iets.copy()
            for key,val in self.forward_censored_iets.iteritems():
                censored_times[key]=censored_times.get(key,0)+val
            event_times=self.observed_iets.copy()
            for key,val in self.observed_iets.iteritems():
                event_times[key]=event_times[key]*2

        km_estimator=KaplanMeierEstimator(censored_times=censored_times,event_times=event_times)
        return km_estimator.get_estimator(variance=variance)

    def get_naive_estimator(self):
        km_estimator=KaplanMeierEstimator(event_times=self.observed_iets)
        return km_estimator.get_naive_estimator()


    def get_length_bias_estimator(self):
        """Returns an estimator correcting for the length bias of windows cencored data. 

        This estimator only takes into account the non-censored data. For the non-censored
        data the data points are distributed as p'(\tau) \sim (T-t)p(\tau), where p'(\tau)
        is the observed distribution of the inter-event times (i.e., the uncensored 
        distribution), p(\tau) is the original unobserved distribution which produced the
        underlying event sequence, and T is the time window length.
        """
        t_list=sorted(self.observed_iets.keys())

        events_cum,events_sum=[],0
        for t in t_list:
            events_sum+=self.observed_iets[t]*(self.endTime/float(self.endTime-t))
            events_cum.append(events_sum)

        #t_list,events_cum,censored_cum=self.get_cumulative()

        ts,ns=[0],[1]
        for t,n in zip(t_list,events_cum):
            ts.append(t)
            ns.append(1-n/float(events_sum))
        return ts,ns

    """
    def get_vardi_estimator(self):
        assert self.mode=="censorall"
        #setup
        nx=sum(self.observed_iets.itervalues())
        ny=sum(self.backward_censored_iets.itervalues())
        nz=sum(self.forward_censored_iets.itervalues())
        nw=self.empty_seqs
        ts=list(set(itertools.chain(self.observed_iets.iterkeys(),self.forward_censored_iets.iterkeys(),self.backward_censored_iets.iterkeys())))
        ts.append(self.endTime)
        ts.append(self.endTime+1)
        ts.sort()
        
        #step a
        p_old=[1./float(len(ts)) for key in ts]
        #sum_p_old=1.

        while True:
            #step b
            r=[]
            sumi=0
            for k in range(len(ts)):
                t=ts[k]
                sum_p_old=sum(map(lambda j:p_old[j],range(k,len(ts))))
                sumi=sumi+(self.forward_censored_iets.get(t,0)+self.backward_censored_iets.get(t,0))/float(sum_p_old)
                r.append(self.observed_iets.get(t,0)+p_old[k]*sumi/float(sum_p_old))
            r.append(p_old[-1]*(sumi+nw/float(p_old[-1]  )))

            f=lambda mu:sum( (r[k]*ts[k]/float((nx+nz)*mu+(ny+nw)*ts[k]) for k in range(len(ts))) )-1
            #mu_new=scipy.optimize.newton(f,ts[-1],tol=10**-10)
            #mu_new=scipy.optimize.ridder(f,0,ts[-1])
            mu_new=scipy.optimize.bisect(f,ts[0],ts[-1])
            
            #step c
            p_new=[]
            for k,t in enumerate(ts):
                p_new.append(r[k]*mu_new/float((nx+nz)*mu_new+(ny+nw)*t ))

            #step d
            if sum(map(lambda x,y:abs(x-y),p_old,p_new))>10**-6:
                p_old=p_new
                #sum_p_old=sum(p_old)
            else:
                cump=[1]
                if ts[0]!=0:
                    ts.insert(0,0)
                for pval in p_new:
                    cump.append(cump[-1]-pval)
                print mu_new
                return ts,cump
    """


    def get_npmle_estimator(self,return_mu=False):
        """Modified RT algorithm by Soon et al. (1996)
        """
        assert self.mode=="censorall"
        #setup
        nx=sum(self.observed_iets.itervalues())
        ny=sum(self.backward_censored_iets.itervalues())
        nz=sum(self.forward_censored_iets.itervalues())
        nw=self.empty_seqs
        ts=list(set(itertools.chain(self.observed_iets.iterkeys(),self.forward_censored_iets.iterkeys(),self.backward_censored_iets.iterkeys())))
        if nw!=0:
            ts.append(self.endTime)
        ts.sort()
        
        #step a
        p_old=[1./float(len(ts)) for key in ts]
        #sum_p_old=1.
        v_old=nw/float(ny+nw)

        while True:
            #step b
            r=[]
            sumi=0
            sum_p_old=sum(map(lambda j:p_old[j],range(len(ts))))
            for k in range(len(ts)-1):
                t=ts[k]
                sum_p_old-=p_old[k-1] if k>0 else 0.
                sumi=sumi+(self.forward_censored_iets.get(t,0)+self.backward_censored_iets.get(t,0))/float(sum_p_old)

                #sumi=0.
                #for i in range(k):
                #    sum_p_old=sum(map(lambda j:p_old[j],range(i,len(ts))))
                #    sumi+=(self.forward_censored_iets.get(ts[i],0)+self.backward_censored_iets.get(ts[i],0))/float(sum_p_old)

                r.append(self.observed_iets.get(t,0)+p_old[k]*sumi)
                
            r.append(p_old[-1]*sumi)

            if v_old!=0:
                rh1=v_old*nw/float(v_old )
            else:
                rh1=0

            f=lambda mu:sum( (r[k]*ts[k]/float((nx+nz)*mu+(ny+nw)*ts[k]) for k in range(len(ts))) )-1+rh1/float(ny+nw)
            #mu_new=scipy.optimize.newton(f,ts[-1],tol=10**-10)
            #mu_new=scipy.optimize.ridder(f,0,ts[-1])
            #mu_new=scipy.optimize.bisect(f,ts[0],100*ts[-1])
            mu_new=scipy.optimize.bisect(f,0,100*ts[-1])
            
            #step c
            p_new=[]
            for k,t in enumerate(ts):
                p_new.append(r[k]*mu_new/float((nx+nz)*mu_new+(ny+nw)*t ))
            v_new=rh1*mu_new/float(ny+nw)

            #step d
            if sum(map(lambda x,y:abs(x-y),p_old,p_new))>10**-4:
                #print sum(map(lambda x,y:abs(x-y),p_old,p_new)),sum(r),rh1,nx+nz,ny+nw,sum_p_old
                p_old=p_new
                v_old=v_new
            else:
                #print mu_new,v_new
                #print ts,p_new
                cump=[1]
                if ts[0]!=0:
                    ts.insert(0,0)
                for pval in p_new:
                    cump.append(cump[-1]-pval)
                if return_mu:
                    return ts,cump,mu_new
                else:
                    return ts,cump


    def estimate_moment(self,moment,method="naive",central=False):
        """Returns an estimate for a moment of the inter-event time distribution.

        Choose one of the following methods.
        'naive' : Moment of the observed inter-event time distribution
        'lowerbound' : Moment of a distribution where one has the observed inter-event 
        times, the censored inter-event times and time window widths equal to the number
        of empty time sequences.
        """ 
        s,n=0,0
        if central:
            u=self.estimate_moment(1,method=method,central=False)
        else:
            u=0.0
        if method=="naive":
            for iet,num in self.observed_iets.iteritems():
                s+=num*((iet-u)**moment)
                n+=num
            if n!=0:
                if central:
                    if n>1:
                        return s/float(n-1)
                    else:
                        return None
                else:
                    return s/float(n)
            else:
                return None
        elif method=="lowerbound":
            assert self.mode in ["censorall","censorlast"]
            mult=2 if self.mode == "censorall" else 1
            for iet,num in self.observed_iets.iteritems():
                s+=mult*num*((iet-u)**moment)
                n+=mult*num
            for iet,num in itertools.chain(self.forward_censored_iets.iteritems(),self.backward_censored_iets.iteritems()):
                s+=num*((iet-u)**moment)
                n+=num
            s+=self.empty_seqs*((self.endTime-u)**moment)
            n+=self.empty_seqs
            if central:
                if n>1:
                    return s/float(n-1)
                else:
                    return None
            else:
                return s/float(n)
        elif method=="lbias":
            for iet,num in self.observed_iets.iteritems():
                s+=num/(1-iet/float(self.endTime))*((iet-u)**moment)
                n+=num/(1-iet/float(self.endTime))
            if n!=0:
                if central:
                    if n>1:
                        return s/float(n-1)
                    else:
                        return None
                else:
                    return s/float(n)
            else:
                return None
        elif method=="lbias_lbound":
            if central:
                raise Exception("Not implemented.")
            s2,n2=0,0
            assert self.mode in ["censorall","censorlast"]
            mult=2 if self.mode == "censorall" else 1
            for iet,num in self.observed_iets.iteritems():
                s+=2*num/(1-iet/float(self.endTime))*((iet-u)**moment)
                n+=2*num/(1-iet/float(self.endTime))
            for iet,num in itertools.chain(self.forward_censored_iets.iteritems(),self.backward_censored_iets.iteritems()):
                s2+=num*((iet-u)**moment)
                n2+=num
            s2+=self.empty_seqs*((self.endTime-u)**moment)
            n2+=self.empty_seqs
            if n!=0:
                lbias= s/float(n)
                p_lbias=n/float(n+n2)
            else:
                lbias= 0
                p_lbias=0
            lbound=s2/float(n2)
            p_lbound=n2/float(n+n2)
            return lbias*p_lbias+lbound*p_lbound
        elif method=="poisson":
            if central:
                raise Exception("Not implemented.")
            if len(self.observed_iets)!=0:
                rate=len(self.observed_iets)/float(self.endTime*self.nseqs)
            else:
                rate=1./float(self.endTime*self.nseqs)
            return math.factorial(moment)/float(rate**moment)
        elif method=="npmle":
            if central:
                raise Exception("Not implemented.")
            if len(self.observed_iets)!=0:
                ts,ps=self.get_npmle_estimator()
                cd=dists.CumDist(ts,ps,maxval=self.endTime)
                return cd.get_moment(moment)
            else:
                return self.endTime
        elif method=="npmle_mod":
            if central:
                raise Exception("Not implemented.")
            if len(self.observed_iets)>1:
                ts,ps=self.get_npmle_estimator()
                cd=dists.CumDist(ts,ps,maxval=ts[-1])
                #cd.ps[-1]=0
                return cd.get_moment(moment)
            elif len(self.observed_iets)==0:
                return None #self.endTime**moment
            elif len(self.observed_iets)==1:
                return None #(self.endTime/2.)**moment
        elif method=="km":
            if central:
                raise Exception("Not implemented.")

            if len(self.observed_iets)!=0:
                ts,ps=self.get_estimator()
                cd=dists.CumDist(ts,ps,maxval=ts[-1])
                return cd.get_moment(moment)
            else:
                return self.endTime
        elif method=="km_mod":
            if central:
                raise Exception("Not implemented.")

            niets=sum(self.observed_iets.itervalues())
            if niets>1 :
                ts,ps=self.get_estimator()
                cd=dists.CumDist(ts,ps,maxval=ts[-1])
                #cd.ps[-1]=0
                return cd.get_moment(moment)
            elif niets==0:
                return self.endTime**moment
            elif niets==1:
                return (self.endTime/2.)**moment

        else:
            raise Exception("Invalid parameter value for 'method': "+method)

    def estimate_moment_number_of_data_points(self,method):
        """Returns the number of data points used by the estimator specified in the parameter 'method'.

        This function is useful when constructing new estimators based on moments. For example, the
        sample variance 's' can be calculated by the following code:

        >>> m1=est.estimate_moment(1,method="naive")
        >>> m2=est.estimate_moment(1,method="naive")
        >>> n=est.estimate_moment_number_of_data_points("naive")
        >>> s=n/float(n-1)*(m2-m1**2)
        """
        if method=="naive":
            n=0
            for iet,num in self.observed_iets.iteritems():
                    n+=num
            return n
        elif method=="lowerbound":
            n=0
            for iet,num in self.observed_iets.iteritems():
                n+=num
            for iet,num in itertools.chain(self.forward_censored_iets.iteritems(),self.backward_censored_iets.iteritems()):
                n+=num
            n+=self.empty_seqs
            return n
        else:
            raise Exception("Invalid parameter value for 'method': "+method)

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

def create_scaled_sequence(uldist,uldist_residual,avgdist,n,T):
    """Create time sequences that follow an underlying distribution.
    """
    for i in range(n):
        navg=avgdist()
        ndist=lambda :navg*uldist()
        ndist_residual=lambda :navg*uldist_residual()
        yield random_timeseq(ndist,ndist_residual,T,0)

def get_burstiness(m1,m2):
    if m1==0:
        return None
    if abs(m2-m1*m1)<max(m2,m1*m1)/10**6:
        s=0
    else:
        s=math.sqrt(m2-m1*m1)
    return (s-m1)/float(s+m1)

def get_central_moment(n,moments=[]):
    un=0
    for j in range(0,n+1):
        un+=binom(n,j)*(-1)**(n-j)*moments[j-1]*moments[0]**(n-j)
    return un

def get_burstiness_less_biased(m1,m2,n,normalcorrection=True,u2=None,u3=None,u4=None,u5=None,u6=None):
    """If all central moments are given solve equation (2) in:
    'Signigicance tests for coefficients of variation and variability profile', Sokal & Braumann (1980)
    """
    if m1==0:
        return None
    if abs(m2-m1*m1)<max(m2,m1*m1)/10**6:
        s=0
        return -1
    else:
        if u2==None:
            s=math.sqrt(float(n)/float(n-1)*(m2-m1*m1))
        else:
            s=math.sqrt(u2)

    if normalcorrection and u3==None:
        cv=(1+1./float(4*n))*s/float(m1)
        #elif m3!=None and m4==None:
        #thetahat=s**2/float(m1**2)
        #gamma=(m3-3*m2*m1+2*m1**3)/float(s**3)
        #thetadilde=thetahat-(thetahat**(3./2.))/float(n)*(3*thetahat**(0.5)-2*gamma)
        #cv=math.sqrt(thetadilde)
    elif u2!=None and u3!=None and u4!=None and u5!=None and u6!=None:
        vp=s/float(m1)
        gamma1=u3/float(s**3)
        gamma2=u4/float(s**4)

        correction1=(1-1./float(4*(n-1)) +vp**2/n - vp*gamma1/float(2*n) - gamma2/float(8*n)) 
        correction2=1./float(2*(n-1)**2) - (3*n-5)/float(4*n**2*(n-1))*vp*gamma1 - 3*(n-5)/float(16*n**2*(n-1))*gamma2 + 1./float(8*n**2)*vp*u5/(u2**2.5) + 1./float(16*n**2)*(u6/float(u2**3)-15) - 1./float(n**2) * vp**3 *gamma1 - (3*n**2 - 6*n +5)/float(8*n**2*(n-1)**2)*gamma1**2 + 1./float(2*n**2)*vp**2*gamma2
        cv=vp/float(correction1+correction2)

        target=s/float(m1) 
        ev=lambda vp:vp*(1-1./float(4*(n-1)) +vp**2/n - vp*gamma1/float(2*n) - gamma2/float(8*n)+1./float(2*(n-1)**2) - (3*n-5)/float(4*n**2*(n-1))*vp*gamma1 - 3*(n-5)/float(16*n**2*(n-1))*gamma2 + 1./float(8*n**2)*vp*u5/(u2**2.5) + 1./float(16*n**2)*(u6/float(u2**3)-15) - 1./float(n**2) * vp**3 *gamma1 - (3*n**2 - 6*n +5)/float(8*n**2*(n-1)**2)*gamma1**2 + 1./float(2*n**2)*vp**2*gamma2)-target

        evp=lambda vp:1-1./float(4*(n-1)) +3*vp**2/n - 2*vp*gamma1/float(2*n) - gamma2/float(8*n)+1./float(2*(n-1)**2) - 2*(3*n-5)/float(4*n**2*(n-1))*vp*gamma1 - 3*(n-5)/float(16*n**2*(n-1))*gamma2 + 2*1./float(8*n**2)*vp*u5/(u2**2.5) + 1./float(16*n**2)*(u6/float(u2**3)-15) - 4*1./float(n**2) * vp**3 *gamma1 - (3*n**2 - 6*n +5)/float(8*n**2*(n-1)**2)*gamma1**2 + 3*1./float(2*n**2)*vp**2*gamma2
        evpp=lambda vp: 6*vp/n - 2*gamma1/float(2*n) - 2*(3*n-5)/float(4*n**2*(n-1))*gamma1 + 2*1./float(8*n**2)*u5/(u2**2.5)  - 12*1./float(n**2) * vp**2 *gamma1  + 6*1./float(2*n**2)*vp*gamma2

        try:
            cv=scipy.optimize.newton(ev,target)
        except RuntimeError:
            cv=scipy.optimize.newton(evp,target)
            assert evpp(cv)<0
            #from matplotlib import pyplot as plt
            #t=lambda x:ev(x/1000.)
            #plt.plot(map(lambda x:x/1000.,range(10000)),map(t,range(10000)))
            #raise Exception()
    else:
        cv=s/float(m1)        

    return float(cv-1)/float(cv+1)


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


    for mode in ['censorlast','censorall','periodic']:
        iet_est=IntereventTimeEstimator(12,mode=mode)
        iet_est.add_time_seq([3,4,5,9])
        print mode,
        print ": ",
        print iet_est.get_naive_estimator(),
        print ", ",
        print iet_est.get_estimator(variance=True)

    iet_est=IntereventTimeEstimator(12,mode='censorall')
    iet_est.add_time_seq([3,4,5,9])
    print iet_est.get_npmle_estimator()
    print iet_est.get_length_bias_estimator()
