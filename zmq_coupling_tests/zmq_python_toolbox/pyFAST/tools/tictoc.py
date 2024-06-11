import numpy as np
import time

def pretty_time(t):
    # fPrettyTime: returns a 6-characters string corresponding to the input time in seconds.
    #   fPrettyTime(612)=='10m12s'
    # AUTHOR: E. Branlard
    if(t<0):
        s='------';
    elif (t<1) :
        c=np.floor(t*100);
        s='{:2d}.{:02d}s'.format(0,int(c))
    elif(t<60) :
        s=np.floor(t);
        c=np.floor((t-s)*100);
        s='{:2d}.{:02d}s'.format(int(s),int(c))
    elif(t<3600) :
        m=np.floor(t/60);
        s=np.mod( np.floor(t), 60);
        s='{:2d}m{:02d}s'.format(int(m),int(s))
    elif(t<86400) :
        h=np.floor(t/3600);
        m=np.floor(( np.mod( np.floor(t) , 3600))/60);
        s='{:2d}h{:02d}m'.format(int(h),int(m))
    elif(t<8553600) : #below 3month
        d=np.floor(t/86400);
        h=np.floor( np.mod(np.floor(t), 86400)/3600);
        s='{:2d}d{:02d}h'.format(int(d),int(h))
    elif(t<31536000):
        m=t/(3600*24*30.5);
        s='{:4.1f}mo'.format(m)
        #s='+3mon.';
    else:
        y=t/(3600*24*365.25);
        s='{:.1f}y'.format(y)
    return s


class Timer(object):
    """ Time a set of commands, as a context manager
    usage:

        with Timer('A name'):
            cmd1
            cmd2
    """
    def __init__(self, name=None, writeBefore=False, silent=False, nChar=40):
        self.name        = name
        self.writeBefore = writeBefore
        self.silent=silent
        self.nChar=nChar
        self.sFmt='{:'+str(nChar+1)+'s}'

    def ref_str(self):
        s='[TIME] '
        if self.name:
            s+=self.sFmt.format(self.name[:self.nChar])
        return s

    def __enter__(self):
        if self.silent:
            return
        self.tstart = time.time()
        if self.writeBefore:
            s=self.ref_str()
            print(s,end='')

    def __exit__(self, type, value, traceback):
        if self.silent:
            return
        if self.writeBefore:
            print('Elapsed: {:6s}'.format(pretty_time(time.time() - self.tstart)))
        else:
            s=self.ref_str()
            print(s+'Elapsed: {:6s}'.format(pretty_time(time.time() - self.tstart)))
