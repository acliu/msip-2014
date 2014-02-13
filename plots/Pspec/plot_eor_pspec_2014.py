#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os, glob

#BINS = [(0.05,0.1),(0.1,0.2),(0.2,0.4),(0.4,0.8),(0.8,1.6)]
BINS = None

def plot(k3pk, style, edge=False, **kwargs):
    ks = k3pk.keys(); ks.sort()
    k3pk = n.array([k3pk[k] for k in ks])
    ks = n.array(ks)
    pk = k3pk / (ks**3/(2*n.pi**2))
    if not BINS is None and edge:
        k3pk_bin = {}
        for b in BINS:
            dsum,dwgt = 0.,0.
            ksum = 0.
            #for k,d in zip(ks, pk):
            for k,d in zip(ks, k3pk):
                #print k, b, k < b[0] or k >= b[1]
                if k < b[0] or k >= b[1]: continue
                w = 1./d**2
                dsum += d*w; ksum += k*w; dwgt += w
            #k,d = ksum/dwgt, dsum/dwgt
            k = (b[0] + b[1])/2
            d = dsum/dwgt
            #k3pk_bin[k] = k**3/(2*n.pi**2)*d
            k3pk_bin[k] = d
        ks = k3pk_bin.keys(); ks.sort()
        k3pk = n.array([k3pk_bin[k] for k in ks])
        ks = n.array(ks)
    print ks
    print k3pk
    #k_edges = n.array([ks/n.sqrt(1.9),ks*n.sqrt(1.9)]).transpose().flatten()
    dk = ks[1] - ks[0]
    k_edges = n.array([n.concatenate([[ks[0]-dk/2],ks[:-1]+dk/2]),ks+dk/2]).transpose().flatten()
    k3pk_edges = n.array([k3pk,k3pk]).transpose().flatten()
    k_edges = n.concatenate([[k_edges[0]], k_edges])
    k3pk_edges = n.concatenate([[1e10], k3pk_edges])
    if edge:
    #if False:
        #p.plot(k_edges, n.sqrt(k3pk_edges), style, **kwargs)
        p.plot(k_edges, k3pk_edges, style, **kwargs)
        #if samp_var:
        #    p.plot(k_edges[:7], k3pk_edges[:7], style, linewidth=8)
    else:
        p.plot(ks, k3pk, style, **kwargs)
    #if not upperlimit: p.errorbar(kpl, pk.real, yerr=err, fmt=color+'.', capsize=0)
    #else: p.plot(kpl, pk.real + err, color+'-')
    ##p.errorbar(kpl, pk.imag, yerr=err, fmt=color+'x')
    #p.subplot(122)
    #k0 = n.abs(kpl).argmin()
    #pkfold = pk[k0:].copy()
    #errfold = err[k0:].copy()
    #pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
    #pkneg,errneg = pk[k0-1:0:-1].copy(), err[k0-1:0:-1].copy()
    #pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
    #errfold[1:] = n.sqrt(1./(1./errpos**2 + 1./errneg**2))
    ##p.errorbar(k, k3*pk, yerr=k3*err, fmt=color+'.', capsize=0)
    #if not upperlimit: p.errorbar(k[k0:], k3[k0:]*pkfold, yerr=k3[k0:]*errfold, fmt=color+'.', capsize=0)
    #else: p.plot(k[k0:], k3[k0:]*(pkfold + errfold), color+'-')
    #kpl = kpl
    #k3pk = k3*pk
    #k3err = k3*err
    ##p.plot(kpl, k3*pk+k3*err, color+'.-')
    ##p.plot(kpl, k3pk+k3err, color+'.-')
    #for _k,_k3pk,_k3err in zip(kpl,k3pk,k3err):
    #    print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    #print '-'*20
    #for _k,_k3pk,_k3err in zip(k[k0:],k3[k0:]*pkfold,k3[k0:]*errfold):
    #    print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)

o = optparse.OptionParser()
PAPER_32 = {
    0.016: 1263175125.40380,
    0.031: 1433167925.58217,
    0.056: 559467861.45039,
    0.082: 14012.80011,
    0.108: 2743.92294,
    0.135: 3631.23328,
    0.162: 3951.62416,
    0.188: 4853.14656,
    0.215: 5489.85212,
    0.242: 5121.77880,
    0.269: 6853.34330,
    0.295: 12067.94883,
    0.322: 11281.72628,
    0.349: 13315.13109,
    0.376: 23695.60404,
    0.403: 29216.86240,
    0.429: 26172.44106,
    0.456: 27728.75816,
    0.483: 56819.39294,
    0.510: 97287.16351,
}
GMRT2013 = {
    0.1 : 2e5,
    0.13: 4e5,
    0.17: 1e5,
    0.19: 2e5,
    0.27: 2e5,
    0.31: 4e5,
    0.40: 5e5,
    0.50: 248**2,
    0.60: 3e5,
}
MWA2013 = {
    0.3: (5e3)**2,
    0.4: (9e3)**2,
    0.9: (15e3)**2,
}

HERA37 = {}
HERA127 = {}
HERA331 = {}
HERA547 = {}

# Kinda dumb, but I'm going to read in these files "by hand" (ACL)



midHERA127 = n.load('../pspec_errs/hera127_mid_z9.0_nf0.37_B0.014.npz')
#sillyInf = midHERA127['err'][0]
HERA127['mid'] = {}
for k,d in zip(midHERA127['ks'],midHERA127['err']):
    if d == n.inf:
        HERA127['mid'][k] = 10**20
    else:
        HERA127['mid'][k] = d * 2. # 2-sigma errors.

midHERA331 = n.load('../pspec_errs/hera331_mid_z9.0_nf0.37_B0.014.npz')
HERA331['mid'] = {}
for k,d in zip(midHERA331['ks'],midHERA331['err']):
    if d == n.inf:
        HERA331['mid'][k] = 10**20
    else:
        HERA331['mid'][k] = d * 2. # 2-sigma errors.

optHERA331 = n.load('../pspec_errs/hera331_opt_z9.0_nf0.37_B0.014.npz')
HERA331['opt'] = {}
for k,d in zip(optHERA331['ks'],optHERA331['err']):
    if d == n.inf:
        HERA331['opt'][k] = 10**20
    else:
        HERA331['opt'][k] = d * 2. # 2-sigma errors.

cosmoSig = {}
for k,d in zip(optHERA331['ks'],optHERA331['delta2']):
    cosmoSig[k] = d # Cosmological signal

#print optHERA331['ks']
#print midHERA331['ks']
#print optHERA331['err']
#print HERA127['mid']
#assert False

#for f in glob.glob('../dillon/*csv'):
#    print os.path.basename(f)[:-len('.csv')]
#    try: H = eval(os.path.basename(f)[:-len('.csv')])
#    except(NameError): continue
#    d1 = n.loadtxt(f,skiprows=7,delimiter=',',usecols=(0,1))
#    d2 = n.loadtxt(f,skiprows=5,delimiter=',',usecols=(3,4))
#    d3 = n.loadtxt(f,skiprows=5,delimiter=',',usecols=(6,7))
#    for cut,dat in zip(['1','2','3'],[d1,d2,d3]):
#        ks,ds = dat[:,0], dat[:,1]
#        H[cut] = {}
#        for k,d in zip(ks,ds):
#            H[cut][k] = d / 1.4 # fudge factor rescaling 8.5 to 7.3

# XXX hack to approximate HERA-61 from HERA-127, noting that Jonnie has:
#  37 2.18 2.76 10.87
#  61 4.02 5.03 18.30
#  91 6.31 7.78 26.71
# 127 8.91 10.88 35.65
# which shows that in total, significance of 61 is approx 61/127 * significance of 127
# we'll use this to rationalize what is about to happen
HERA61 = HERA127.copy()
for cut in HERA61:
    HERA61[cut] = HERA127[cut].copy()
    for k in HERA61[cut]:
        HERA61[cut][k] = HERA127[cut][k] * 127./61

#BINS = None
##d = [RS_VS_KPL[k] for k in kpl]
#colors = 'kbcm' * 10
#for sep in RS_VS_KPL:
#    if not sep == 'total': continue
#    dsum, dwgt = {}, {}
#    ks = RS_VS_KPL[sep].keys(); ks.sort()
#    for k in ks:
#        _d,_n = RS_VS_KPL[sep][k]
#        #print k, _d, _n, '->',
#        dsum[k] = dsum.get(k,0) + _d / _n**2
#        dwgt[k] = dwgt.get(k,0) + 1 / _n**2
#        #print dsum[k] / dwgt[k], 1./n.sqrt(dwgt[k])
#        
#    kpl = dsum.keys(); kpl.sort()
#    d = [dsum[k]/dwgt[k] for k in kpl]
#    nos = [1./n.sqrt(dwgt[k]) for k in kpl]
#    d,kpl,nos = n.array(d, dtype=n.complex), n.array(kpl), n.array(nos)
#    dual_plot(kpl, d, 2*nos, color=colors[0]) # 2-sigma error bars
#    colors = colors[1:] + colors[0]

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

import glob
#re_z = re.compile(r'power_21cm.*_z(\d+\.\d+).*\.dat')

#for filename in glob.glob('lidz_mcquinn_k3pk/*7.3*dat'):

#for filename in glob.glob('lidz_mcquinn_k3pk/*hmass*dat'):
#for cnt,filename in enumerate(glob.glob('lidz_mcquinn_k3pk/*hmass*.[36]*dat')):
#for cnt,filename in enumerate(glob.glob('../lidz_mcquinn_k3pk/power_21cm_z7.3*dat')):
#    print 'Reading', filename
#    d = n.array([map(float, L.split()) for L in open(filename).readlines()])
#    ks, pk = d[:,0], d[:,1]
    #ks = n.concatenate([[.01], ks])
    #pk = n.concatenate([[pk[0]], pk])
    #z_file = float(re_z.match(os.path.basename(filename)).groups()[0])
    #z = C.pspec.f2z(.164)
#    z = 7.3
#    print 'Setting redshift at', z
#    k3pk = ks**3 / (2*n.pi**2) * pk
#    label = 'Lidz et al. (2008)'
#    if cnt == 0: plot(dict(zip(ks,k3pk*mean_temp(z)**2)), 'k:', label=label, linewidth=1.5)
#    else: plot(dict(zip(ks,k3pk*mean_temp(z)**2)), 'k:', linewidth=1.5)
    
#tau_h = 100. + 15 # ns
#k_h = C.pspec.dk_deta(C.pspec.f2z(.164)) * tau_h

#p.subplot(121)
##p.plot([k_h, k_h], [1e-10, 1e20], 'k--')
##p.plot([-k_h, -k_h], [1e-10, 1e20], 'k--')
#p.gca().set_yscale('log', nonposy='clip')
#p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
#p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$')
#p.ylim(1e5,3e16)
#p.grid()


##p.subplot(122)
##p.plot([k_h, k_h], [1e-10, 1e20], 'k--')
#plot(PAPER_32, 'ko', label='PAPER-32 (Parsons et al. 2013)')
#plot(GMRT2013, 'yv', label='GMRT (Paciga et al. 2013)')
#plot(MWA2013, 'co', label='MWA (Dillon et al. 2013)')
#CUT = '2'
#for color,H in zip(('m','orange','darkviolet','orangered'),['HERA37', 'HERA127', 'HERA331', 'HERA547']):
#for cnt,(color,H) in enumerate(zip(('dodgerblue','red','black'),['HERA127', 'HERA331', 'HERA547'])):
#for cnt,(color,H) in enumerate(zip(('dodgerblue','red'),['HERA127', 'HERA331'])):
#    CUT = str(cnt+1)
#    Hdict = eval(H)
#    plot(Hdict[CUT], color, linestyle='-', edge=True, linewidth=2, label=H.replace('HERA','HERA-').replace('547','568'))

plot(HERA127['mid'],'dodgerblue', linestyle='-', edge=True, linewidth=2, label='HERA-127, foreground avoidance')
plot(HERA331['mid'],'red', linestyle='-', edge=True, linewidth=2, label='HERA-331, foreground avoidance')
plot(HERA331['opt'],'black', linestyle='-', edge=True, linewidth=2, label='HERA-331, foreground modeling')
plot(cosmoSig,'black', linestyle='--', label='Mesinger et al. 2011', linewidth=1.5)

p.setp(p.gca().get_xticklabels(), fontsize=20)
p.setp(p.gca().get_yticklabels(), fontsize=20)
p.gca().set_yscale('log', nonposy='clip')
p.gca().set_xscale('log', nonposy='clip')
p.tick_params(axis='both', which='major', length=8)
p.tick_params(axis='both', which='minor', length=4)
p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize=20)
#p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$')
p.ylabel(r'$\Delta^2(k)\ [{\rm mK}^2]$', fontsize=20)
p.ylim(3e-1**2,1e3)
p.xlim(.01, 2.0)
p.legend(loc='best', fontsize=14)
p.grid()
#p.savefig('eor_pspec_2014.png', bbox='tight')
p.show()
