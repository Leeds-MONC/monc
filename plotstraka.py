#!/usr/bin/env python2.7

import netCDF4
import pylab as plt
import numpy as np
import sys, getopt, os, glob, re

def natsort(list_):
    MAXINTLENGTH=6
    # decorate
    tmp=[]
    for i in list_:
        iall=re.findall('\d+', i)
        iall.reverse()
        n=0
        for j in range(len(iall)):
            n+=(10**MAXINTLENGTH)**j*int(iall[j])
        tmp.append((n,i))
    tmp.sort()
    # undecorate
    return [ i[1] for i in tmp ]

def usage():
    print 'straka_plots.py -[hi] -f < filenames >'

def main(argv, case='Straka'):

    l_interactive=False
    try:
        opts, args = getopt.getopt(argv[1:], "hf:i",
                                   ["help", "filenames=", "interactive"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    print opts, args
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        if opt in ("-i", "--interactive"):
            l_interactive=True
        elif opt in ("-f", "--filenames"):
            filenames = natsort(glob.glob(r'%s' % arg))

    nfiles=len(filenames)
    plt.figure(figsize=(8,2*nfiles))
    axDict={}
    for ifile,filename in enumerate(filenames):
        short_name=os.path.basename(filename)
        print short_name

        ncid = netCDF4.Dataset(filename)

        x =ncid.variables['x'][:]
        y =ncid.variables['y'][:]
        z =ncid.variables['z'][:]

        # Center y
        y = y - 0.5*y[-1]

        dx = np.zeros_like(x)
        dx[:-1] = np.diff(x)
        dx[-1]=dx[-2]

        dy = np.zeros_like(y)
        dy[:-1] = np.diff(y)
        dy[-1]=dy[-2]

        dz = np.zeros_like(z)
        dz[:-1] = np.diff(z)
        dz[-1]=dz[-2]

        yoff = y + 0.5*dy
        xoff = x + 0.5*dx
        zoff = z - 0.5*dz

        u =ncid.variables['u'][:]
        v =ncid.variables['v'][:]
        w =ncid.variables['w'][:]

        theta_pert = ncid.variables['th'][:]

        ncid.close()

        # Generate meshes X-Z
        XXZ,ZXZ=np.meshgrid(xoff,zoff)

        XXoffZ,ZXoffZ=np.meshgrid(xoff,z)
        XXZoff,ZXZoff=np.meshgrid(x,zoff)
        # Generate meshes Y-Z
        YYZ,ZYZ=np.meshgrid(yoff,zoff)

        YYoffZ,ZYoffZ=np.meshgrid(yoff,z)
        YYZoff,ZYZoff=np.meshgrid(y,zoff)
                                                                                                      49,1          40%

        u_th = np.zeros_like(u)
        u_th[:,:,:-1] = 0.5*(u[:,:,:-1] + u[:,:,1:])

        v_th = np.zeros_like(v)
        v_th[:,:-1,:] = 0.5*(v[:,:-1,:] + v[:,1:,:])

        w_th = np.zeros_like(w)
        w_th[:-1,:,:] = 0.5*(w[:-1,:,:] + w[1:,:,:])


        if case=='Straka': # Some configurable plotting options
            XDIM=0
            JSAMP=12
            KSAMP=4
            umin=-100.
            umax=100.
            wmin=-30
            wmax=30
            arrow_scale=1000
            Xoff=YYZ
            Zoff=ZYZ
            X=YYZ[::KSAMP,::JSAMP]
            Z=ZYZ[::KSAMP,::JSAMP]
            uwind=v_th[XDIM,:,:]
            wwind=w_th[XDIM,:,:]
            v1=theta_pert[XDIM,:,:]
            v1label='theta perturbation (K)'
            v1_min=-16.
            v1_max=0.

        titles={}

        if ifile==0:
            axDict[ifile]=plt.subplot(nfiles,1,ifile+1)
        else:
            axDict[ifile]=plt.subplot(nfiles,1,ifile+1, sharex=axDict[0], sharey=axDict[0])
        ax=axDict[ifile]
        titles[ax]='%s' % v1label

        plt.contour(Xoff, Zoff, v1.T,
                    levels=[-1,-2,-3,-4,-5,
                             -6,-7,-8,-9,-10,
                             -11,-12,-13,-14,-15])
        plt.xlabel('y (km)')
        plt.ylabel('z (km)')
        plt.axis('tight')
        plt.axis(ymin=0.0, xmin=0.0, xmax=16000.)
        plt.text(.99,.9, short_name, horizontalalignment='right',transform=ax.transAxes,size=8, color='gray')


    if l_interactive:
        plt.show()
    else:
        plt.savefig('straka.png')


if __name__=='__main__':
    main(sys.argv, case='Straka')
                                                                                                      158,1         Bot

                                                                                                      91,1          80%



