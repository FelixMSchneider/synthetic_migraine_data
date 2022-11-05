import pylab as plt
import numpy as np
import datetime

def lowpass(data, freq, df, corners=4, zerophase=True):
    """
    Butterworth-Lowpass Filter.

    Filter data removing data over certain frequency ``freq`` using ``corners``
    corners.
    The filter uses :func:`scipy.signal.iirfilter` (for design)
    and :func:`scipy.signal.sosfilt` (for applying the filter).

    :type data: numpy.ndarray
    :param data: Data to filter.
    :param freq: Filter corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / order.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    from  scipy.signal import iirfilter,sosfilt, zpk2sos 
    import warnings
    
    fe = 0.5 * df
    f = freq / fe
    # raise for some bad scenarios
    if f > 1:
        f = 1.0
        msg = "Selected corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    z, p, k = iirfilter(corners, f, btype='lowpass', ftype='butter',
                        output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        firstpass = sosfilt(sos, data)
        return sosfilt(sos, firstpass[::-1])[::-1]
    else:
        return sosfilt(sos, data)
 

class Syndata(): 
    def __init__(self, modT=False, randseed=43 ):
        self.npts=365*2 
        self.starttime=datetime.datetime(2009, 1, 1, 0, 0)
        self.dt=1 # 1 day
        self.dtime=datetime.timedelta(self.dt) # 1 day
        self.T0=30
        self.threshold=1.1 
        self.modT=modT
        self.randseed=randseed

    def get_w_T(self,  regphase=365, startmod=365*2,DYfac=2 ):
        """
        calulates T(t) and w(t)
        if self.modT = False T and w are constant
        if self.modT = True T(t) is derived by get_lin_mod_T
        """
        if self.modT:
            self.T=self.get_lin_mod_T(regphase, startmod,DYfac)
        else:
            self.T = np.full_like(self.times, self.T0 )
        self.w=2*np.pi/self.T

    def run(self, regphase=365, startmod=365*2,DYfac=2 ):
        self.prng = np.random.RandomState(self.randseed)
        self.get_timeaxis()
        self.get_w_T(regphase, startmod,DYfac)
        self.get_syndata()
        


    def get_timeaxis(self):
       timeaxis=[]
       iaxis=[]
       for i in range(self.npts):
            timeaxis.append(self.starttime + i*self.dtime)
            iaxis.append(i*self.dt)
       self.timeaxis=timeaxis
       self.times=iaxis
         
    @staticmethod
    def lp_filter(data,dt,cornerT):
        """ 
        Butterworth lowpass filter       
        """
        dt_sec=dt.total_seconds()
        cornerT=cornerT*24*60*60
        _df=1/dt_sec
        data=lowpass(data, 1/cornerT,_df, zerophase=True)
        return data
     
    def get_lin_mod_T(self, regphase=365, startmod=365*2,DYfac=2 ):
        """
        returns T(t) linearly modulated  from T0 to DYfac*T0
        regphase: length of time interval in which T is linear modulated
        startmod: t1 (start of modulation)
        """
        T0=self.T0
        taxis=self.times
        dt=taxis[1]-taxis[0]
        Tmod=[]
        Tn=T0
    
        #m=dy/dx
        DY=T0*DYfac
        m=DY/regphase
    
        stopmod=startmod+regphase
        for t in taxis:
            if t> startmod and t <= stopmod:
                Tn=Tn+m*dt
            Tmod.append(Tn)
        return np.array(Tmod)

    def get_biorythm(self):
        """
        returns biorythm which has angular frequency w(t)
        """
        di=self.times[1]-self.times[0]
        return np.sin(np.cumsum(self.w*di))

    def get_noise(self, lp=6):
        noise=self.prng.randn(len(self.times))
        if lp:
            noise=self.lp_filter(noise,self.dtime,lp)
        return noise

    def get_syndata(self, fac=1):
        threshold=self.threshold
        self.biorythm=self.get_biorythm()
        self.noise_lp=self.get_noise()
        self.trigger=self.biorythm + fac *self.noise_lp
        noise_lp=self.noise_lp
        biorythm=self.biorythm

        self.syndata = np.heaviside(np.maximum(self.trigger, threshold)-threshold, 0)
        
    def plot(self):
        fig=plt.figure(figsize=(8,6))

        ax1=fig.add_subplot(211)
        ax2=fig.add_subplot(212, sharex=ax1)

        threshold=self.threshold
        timeaxis=self.timeaxis
        ax1.plot([timeaxis[0], timeaxis[-1]], [threshold, threshold], "k--")
        ax1.plot(self.timeaxis,self.biorythm )
        ax1.plot(self.timeaxis,self.trigger )
        ax2.plot(self.timeaxis, self.syndata ,"o")

        plt.savefig("syndata.png")  

if __name__ == "__main__":
    syn=Syndata()
    syn.modT=True
    syn.run(regphase=100, startmod=300,DYfac=2)
    syn.plot()
