import pylab as plt
import numpy as np
from matplotlib import cm

from get_syndata import Syndata


import math
def _nearest_pow_2(x):
    """
    Find power of two nearest to x

    >>> _nearest_pow_2(3)
    2.0
    >>> _nearest_pow_2(15)
    16.0

    :type x: float
    :param x: Number
    :rtype: int
    :return: Nearest power of 2 to x
    """
    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b

def plot_spectrogram(syn, logscale=True, local_normalize=True):

    threshold=syn.threshold
    biorythm=syn.biorythm
    trigger=syn.trigger
    timeaxis=syn.timeaxis
    syndata=syn.syndata 
    T=syn.T
    
    
    
    fig = plt.figure(figsize=(10,7))
    ax0 = fig.add_subplot(311)
    ax1 = fig.add_subplot(312,sharex=ax0)
    ax2 = fig.add_subplot(313,sharex=ax0)
    ax=ax2
    
    ax1.plot([timeaxis[0], timeaxis[-1]], [threshold, threshold], "k--")
    
    ax1.plot(timeaxis, biorythm)
    ax1.plot(timeaxis, trigger)
    ax0.plot(timeaxis, syndata)
    
    
    samp_rate = syn.dt
    wlen=200
    data=syndata
    
    npts = len(data)
    nfft = int(_nearest_pow_2(wlen * samp_rate))
    nlap = int(nfft * float(0.9))
    data = data - data.mean()
    end = npts / samp_rate
    
    
    from matplotlib import mlab
    mult=2
    mult = int(_nearest_pow_2(mult))
    mult = mult * nfft
    specgram, freq, time = mlab.specgram(data, Fs=samp_rate, NFFT=nfft,
                                             pad_to=mult, noverlap=nlap)
    clip=[0., 1.0]
    dbscale=False
    if dbscale:
        specgram = 10 * np.log10(specgram[1:, :])
    else:
        specgram = np.sqrt(specgram[1:, :])
    
    freq = freq[1:]
    
    vmin, vmax = clip
    
    if vmin < 0 or vmax > 1 or vmin >= vmax:
        msg = "Invalid parameters for clip option."
        raise ValueError(msg)
    
    from matplotlib.colors import Normalize
    
    # calculate half bin width
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0
   
    zorder=None 
    # argument None is not allowed for kwargs on matplotlib python 3.3
    kwargs = {k: v for k, v in (('cmap', cm.viridis), ('zorder', zorder))
              if v is not None}
    
    # pcolor expects one bin more at the right end
    freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
    time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
    # center bin
    time -= halfbin_time
    freq -= halfbin_freq

    # Log scaling for frequency values (y-axis)
    if logscale: ax.set_yscale('log')

    A=specgram
    
    if local_normalize:
        A=(A/A.max(axis=0))
    
    _range = float(A.max() - A.min())
    vmin = A.min() + vmin * _range
    vmax = A.min() + vmax * _range
    norm = Normalize(vmin, vmax, clip=True)
    
    ax.pcolormesh(np.array(timeaxis)[np.array(time).astype(int)], 1/freq, A, norm=norm, **kwargs)
    ax.set_ylim(2,0.5*5*365/2)
    ax.plot(timeaxis, T, "k--", alpha=0.2, linewidth=2)
    ax.set_xlim(timeaxis[0], timeaxis[-1] )
    
    
    
    ax.set_ylabel("T [Days]")
    
    plt.savefig("spectrogram_from_migattack.png")




syn=Syndata()
syn.npts=365*5
syn.modT=True
syn.randseed=10
syn.run()

plot_spectrogram(syn)
 
