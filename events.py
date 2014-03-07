
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
    def __init__(self,endTime):
        super(IntereventTimeEstimator, self).__init__()
        self.endTime=endTime

    def add_time_seq(self,seq):
        for i,time in enumerate(seq):
            if i!=0:
                dt=time-last
                assert dt>=0
                self.add_event(dt)
            last=time
        self.add_censored(self.endTime-last)

    def add_time_duration_seq(self,seq):
        for i,(time,duration) in enumerate(seq):
            if i!=0:
                dt=time-last
                assert dt>=0
                self.add_event(dt)
            last=time+duration
        self.add_censored(self.endTime-last)


def edgestotimeseqs(edges):
    """Generator transforming edges to time sequences.
    """
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




if __name__=="__main__":
    km=KaplanMeierEstimator()
    km.add_event(3)
    km.add_event(11,2)
    km.add_censored(9)
    km.add_censored(12,6)
    t,s=km.get_estimator()
    print t
    print s
