# hkbin
# Copyright (C) 2020-2021  Sebastiaan L. Zoutendijk
# Based on pyGravSphere
# Copyright (C) 2019-2020  Anna Genina, Justin Read
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from GSpro import profiles
from scipy.integrate.quadrature import simps as integrator
from GSpro import fitting_funcs as fits
import numpy as np
import matplotlib.pyplot as plt

# Changes relative to pyGravSphere are enclosed by `#SLZ` and `#/SLZ`.

myfontsize = 18
mylinewidth = 1
ranalmin = 0.001
ranalmax = 500.
ranalpnts = 10000
ranal = np.logspace(np.log10(ranalmin),np.log10(ranalmax),ranalpnts) #fine bins for integration

figx = 5
figy = 5

y_sigLOSmax = 45
ymin_Sigstar = 1e1
ymax_Sigstar = 1e5
yMlow = 1e6
yMhigh = 1e10



def get_surfden_bins(R,ms,Nbin,maxdatrad,maxdatfitrad,p0in,p0in_min,p0in_max, Mstar_rlim, outdir, gal_num): #photometric positions, photometric (number)masses as input, number of stars per bin, Returns bins.
    cnt = 0 #bin id
    jsum = 0.0 #stars in bin
    norm = np.zeros(len(R))
    rbin_phot_t = np.zeros(len(R))
    surfden_t = np.zeros(len(R)) #temporary (number)mass keeper
    index = np.argsort(R)  # sort by position
    for i in range(len(R)): # for each star
        #SLZ: add any remaining stars to the last bin instead of discarding them
        surfden_t[cnt] = surfden_t[cnt] + ms[index[i]] # add masses to bin
        jsum = jsum + ms[index[i]] # add masses to get one number
        rbin_phot_t[cnt] = R[index[i]] # update max position
        if ((jsum >= Nbin) and (np.sum(ms)-np.sum(norm)-jsum) >= Nbin) or (i == len(R)-1): # gone above Nbin and enough stars for new bin
        #/SLZ
            norm[cnt] = jsum # total number in bin
            if (cnt == 0):
                area = np.pi*rbin_phot_t[cnt]**2.0 # if this is bin zero
            else:
                area = np.pi*(rbin_phot_t[cnt]**2.0-rbin_phot_t[cnt-1]**2.0)         #if not bin zero

            surfden_t[cnt] = surfden_t[cnt]/area #divide number by area - get surface density
            jsum = 0.0
            cnt = cnt + 1 # move on to next bin

    surfdenerr_t = surfden_t / np.sqrt(norm) # Poisson error in each bin
    rbin_phot_t = rbin_phot_t[:cnt] # there are only cnt bins
    surfden_t = surfden_t[:cnt]
    surfdenerr_t = surfdenerr_t[:cnt]
    rbin_phot = rbin_phot_t[rbin_phot_t < maxdatrad] # bins only go as far as (50 kpc)
    surfden = surfden_t[rbin_phot_t < maxdatrad]
    surfdenerr = surfdenerr_t[rbin_phot_t < maxdatrad]
    rbin_photfit = rbin_phot_t[rbin_phot_t < maxdatfitrad]
    surfdenfit = surfden_t[rbin_phot_t < maxdatfitrad]
    surfdenerrfit = surfdenerr_t[rbin_phot_t < maxdatfitrad]  # fitting purposes, exclude some really outer radius

    pfits = fits.tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_photfit,surfdenfit,surfdenerrfit) #get 3 plummer fit. p0 - priors and initial values

    Mstar_rad = rbin_phot
    norm = np.max(profiles.threeplummass(np.linspace(0,Mstar_rlim,100),\
                                pfits[0],pfits[1],\
                                pfits[2],\
                                pfits[3],pfits[4],pfits[5]))   # what is the maximum mass (at Mstar_rlim)?
    Mstar_prof = profiles.threeplummass(Mstar_rad,pfits[0],pfits[1], \
                               pfits[2],\
                               pfits[3],pfits[4],pfits[5])  / norm  # mass at each photometric bin, normalized
    Mstar_surf = profiles.threeplumsurf(ranal,pfits[0],pfits[1],\
                               pfits[2],\
                               pfits[3],pfits[4],pfits[5]) / norm   # fitted surface density at each photometric bin

    Mcum_surf = 0.0
    i = 1

    print 'Norm ::', norm,pfits[0]+pfits[1]+pfits[2]

    while (Mcum_surf < (pfits[0]+pfits[1]+pfits[2])/2.0/norm): #  cum_mass/tot_mass = 0.5  compute the half-mass radius

        Mcum_surf =   Mcum_surf + \
                    2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1]) #  Sum 2 pi r M dr
        #Alternatively : use Simpson's rule  Mcum_surf =  Mcum_surf + integrator(Mstar_surf[:i] * 2 * np.pi * ranal[:i],ranal[:i])
        i = i + 1
    Rhalf = ranal[i-1]
    print 'Rhalf calculated: ', Rhalf

    return rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits


def calc_virial_moments(Rhalf,nmonte,R,vz,vzerr,ms,pfits,maxdatrad,Nbin, outdir, gal_num):
    #Improve estimator using fitted surfden:
    rint = np.logspace(np.log10(Rhalf/100.0),\
                       np.log10(Rhalf*1000.0),10000)
    index = np.argsort(R) # sort all kinematics positions by distance

    #First calculate and subtract the mean vz:  ms ~ number contribution

    vzmean = np.sum(vz*ms)/np.sum(ms)   #com velocity
    vzmeanerr = 0.

    #And now the 2nd and 4th moments:
    cnt = 0
    jsum = 0.0
    norm = np.zeros(len(R))
    vlos4med = np.zeros(len(R))
    vlos2med = np.zeros(len(R))
    #SLZ: calculate moments of verr
    verr4med = np.zeros(len(R))
    verr2med = np.zeros(len(R))
    #/SLZ
    rbin_tmp = np.zeros(len(R))
    for i in range(len(R)):       #fill bins with specified stars per bin
        #SLZ: add remaining stars to last bin, use mean instead of max radius, calculate moments of verr
        vlos4med[cnt] = vlos4med[cnt] + \
            (vz[index[i]]-vzmean)**4.*ms[index[i]]    # number weighted 4th moment
        vlos2med[cnt] = vlos2med[cnt] + \
            (vz[index[i]]-vzmean)**2.*ms[index[i]]    # number weighted 2nd moment
        verr4med[cnt] = verr4med[cnt] + \
            vzerr[index[i]]**4.*ms[index[i]]    # number weighted 4th moment of verr
        verr2med[cnt] = verr2med[cnt] + \
            vzerr[index[i]]**2.*ms[index[i]]    # number weighted 2nd moment of verr
        rbin_tmp[cnt] = rbin_tmp[cnt] + \
            R[index[i]]*ms[index[i]]            # number weighted mean radius
        jsum = jsum + ms[index[i]]
        if ((jsum >= Nbin) and (np.sum(ms)-np.sum(norm)-jsum) >= Nbin) or (i == len(R)-1): # gone above Nbin and enough stars for new bin
        #/SLZ
            norm[cnt] = jsum
            jsum = 0.0
            cnt = cnt + 1
    vlos4med = vlos4med[:cnt]
    vlos2med = vlos2med[:cnt]
    #SLZ: calculate moments of verr
    verr4med = verr4med[:cnt]
    verr2med = verr2med[:cnt]
    #/SLZ
    norm = norm[:cnt]
    vlos4med = vlos4med / norm     # normalize with tot stars
    vlos2med = vlos2med / norm
    #SLZ: calculate moments of verr, use mean instead of max radius
    verr4med = verr4med / norm
    verr2med = verr2med / norm
    rbin_tmp = rbin_tmp[:cnt] / norm
    #/SLZ

    #SLZ: estimate unbiased intrinsic moments using h and k estimators
    m2 = vlos2med    # sample moments
    m4 = vlos4med
    m2err = verr2med
    m4err = 3*verr4med
    h2 = norm/(norm-1) * m2    # h statistics: estimators of population moments
    h4 = (3*(3-2*norm)*norm**2*m2**2 + norm**2*(norm**2-2*norm+3)*m4) / ((norm-3)*(norm-2)*(norm-1)*norm)
    k2 = norm/(norm-1) * m2    # k statistics: estimators of population cumulants
    k4 = (norm**2 * ((norm+1)*m4 - 3*(norm-1)*m2**2)) / ((norm-1)*(norm-2)*(norm-3))
    k2err = m2err
    k4err = m4err - 3*m2err**2
    k2int = k2 - k2err    # intrinsic k statistics: estimators of intrinsic population cumulants
    k4int = k4 - k4err
    m2int = (norm-1)/norm * k2int    # estimators of intrinsic sample moments
    m4int = ((norm-1)*(norm-2)*(norm-3)*k4int/norm**2 + 3*(norm-1)*m2int**2) / (norm+1)
    h2int = norm/(norm-1) * m2int    # intrinsic h statistics: estimators of intrinsic population moments
    h4int = (3*(3-2*norm)*norm**2*m2int**2 + norm**2*(norm**2-2*norm+3)*m4int) / ((norm-3)*(norm-2)*(norm-1)*norm)

    vlos2var = ((norm-1)*h4 - (norm-3)*h2**2) / (norm*(norm-1)) # slightly biased

    vlos4med = h4int
    vlos2med = h2int
    vlos4err = np.full(vlos4med.shape, np.nan) # not enough data to estimate
    vlos2err = np.sqrt(vlos2var) # slightly biased
    #/SLZ

    #Demand positive:

    print 'Dispersion profile:', np.sqrt(vlos2med)


    #SLZ: relax requirement on positive vlos4med
    #/SLZ

    #SLZ: disable vlos4med fit
    #/SLZ


    #Cut back to maxdatrad:
    rbin_tmp_full = rbin_tmp
    vlos4err_full = vlos4err
    vlos2err_full = vlos2err
    vlos4med_full = vlos4med
    vlos2med_full = vlos2med
    vzmean_full = vzmean
    vzmeanerr_full = vzmeanerr
    vlos4err = vlos4err[rbin_tmp < maxdatrad]
    vlos2err = vlos2err[rbin_tmp < maxdatrad]
    vlos4med = vlos4med[rbin_tmp < maxdatrad]
    vlos2med = vlos2med[rbin_tmp < maxdatrad]
    rbin_tmp = rbin_tmp[rbin_tmp < maxdatrad]

    #SLZ: disable vlos4med fit
    #/SLZ


    #plt.figure()
    #plt.loglog()
    #plt.errorbar(rbin_tmp_full,vlos4med_full,vlos4err_full,color='k')
    #plt.errorbar(rbin_tmp,vlos4med,vlos4err,color='b')
    #tvl4 = fits.tvl4func(rint,rbin_tmp,vlos4med,\
    #                pfits_powline[0],pfits_powline[1],\
#                    pfits_powline[2],router,gamout)
#    plt.plot(rint,tvl4,'r')
#    plt.show()


    #Plot the 1st, 2nd and 4th moments:
    if (len(rbin_tmp_full) > 1):
        #Calculate sigLOS(Rhalf) and output:
        if (np.max(rbin_tmp_full) > Rhalf):
            j=0L
            while (rbin_tmp_full[j] < Rhalf):
                j=j+1
            print 'sigLOS(Rhalf) [km/s]:', np.sqrt(vlos2med_full[j])

        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)

        plt.loglog()
        plt.errorbar(rbin_tmp_full,vlos4med_full,vlos4err_full,color='black')
        plt.errorbar(rbin_tmp,vlos4med,vlos4err,color='blue')
        #SLZ: disable vlos4med fit
        #/SLZ
        plt.xlim([Rhalf/10.0,Rhalf*100.0])
        plt.ylim([1.0,y_sigLOSmax**4.0*100.0])
        plt.xlabel(r'$R\,[{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\langle v_{\rm los}^4\rangle\,({\rm km}^4\,{\rm s}^{-4})$',\
                       fontsize=myfontsize)
        plt.savefig(outdir+'Galaxy_%s_output_vlos4.pdf' % gal_num,bbox_inches='tight')
        plt.close()
        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)

        plt.loglog()
        plt.errorbar(rbin_tmp_full,np.sqrt(vlos2med_full),\
                         vlos2err_full/2.0/np.sqrt(vlos2med_full),color='k')
        plt.errorbar(rbin_tmp,np.sqrt(vlos2med),\
                         vlos2err/2.0/np.sqrt(vlos2med),color='b')
        plt.xlim([Rhalf/10.0,Rhalf*100.0])
        plt.ylim([1.0,y_sigLOSmax])

        plt.xlabel(r'$R\,[{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\sigma_{\rm LOS}\,({\rm km}\,{\rm s}^{-1})$',\
                       fontsize=myfontsize)
        plt.savefig(outdir+'Galaxy_%s_output_vlos2.pdf' % gal_num,bbox_inches='tight')
        plt.close()
    #SLZ: do not calculate VSPs, not enough data for ultra-faint dwarf galaxies
    vs1bin = np.nan
    vs2bin = np.nan
    vs1err = np.nan
    vs2err = np.nan
    #/SLZ

    #Output also 2nd moment (pruning any would-be NaN values):
    rbin_kin = rbin_tmp[vlos2med > 0]
    sigpmean = np.sqrt(vlos2med[vlos2med > 0])
    sigperr = vlos2err[vlos2med > 0]/2.0/\
        np.sqrt(vlos2med[vlos2med > 0])

    #SLZ: do not calculate Richardson & Fairbairn estimators, not enough data
    #/SLZ

    mean_disp = np.sum(sigpmean)/np.float(len(sigpmean))
    print 'Mean dispersion:', mean_disp
    print 'Mean dispersion error:', np.sqrt(np.sum((sigpmean-mean_disp)**2.0)/\
        np.float(len(sigpmean)-1))/np.sqrt(len(sigpmean))

    #plt.errorbar(rbin_kin,sigpmean, sigperr)
    #plt.show()

    #SLZ: no VSPs to plot
    #/SLZ

    return rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err
