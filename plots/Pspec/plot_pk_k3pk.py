#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os, glob

ONLY_POS_K = True

def dual_plot(kpl, pk, err, umag=16., f0=.164, color='', upperlimit=False):
    z = C.pspec.f2z(f0)
    kpr = C.pspec.dk_du(z) * umag
    k = n.sqrt(kpl**2 + kpr**2)
    k3 = n.abs(k**3 / (2*n.pi**2))
    print 'k [h Mpc^-1], P(k) [K^2], err (2sigma)'
    for _k,_pk,_err in zip(kpl,pk,err):
        print '%6.3f, %9.5f, %9.5f' % (_k, _pk.real/1e6, _err/1e6)
    print '-'*20
    #p.subplot(121)
    #pk = n.abs(pk)
    #p.errorbar(kpl, pk, yerr=err, fmt=color+'.-')
    #p.errorbar(kpl, n.abs(pk), yerr=err, fmt='mx')
    #if not upperlimit: p.errorbar(kpl, pk.real, yerr=err, fmt=color+'.', capsize=0)
    #else: p.plot(kpl, pk.real + err, color+'-')
    #p.errorbar(kpl, pk.imag, yerr=err, fmt=color+'x')
    #p.subplot(122)
    k0 = n.abs(kpl).argmin()
    pkfold = pk[k0:].copy()
    errfold = err[k0:].copy()
    pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
    pkneg,errneg = pk[k0-1:0:-1].copy(), err[k0-1:0:-1].copy()
    pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
    errfold[1:] = n.sqrt(1./(1./errpos**2 + 1./errneg**2))
    #p.errorbar(k, k3*pk, yerr=k3*err, fmt=color+'.', capsize=0)
    if not upperlimit: p.errorbar(k[k0:], k3[k0:]*pkfold, yerr=k3[k0:]*errfold, fmt=color+'.', capsize=0)
    else: p.plot(k[k0:], k3[k0:]*(pkfold + errfold), color+'-')
    kpl = kpl
    k3pk = k3*pk
    k3err = k3*err
    #p.plot(kpl, k3*pk+k3*err, color+'.-')
    #p.plot(kpl, k3pk+k3err, color+'.-')
    for _k,_k3pk,_k3err in zip(kpl,k3pk,k3err):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    print '-'*20
    for _k,_k3pk,_k3err in zip(k[k0:],k3[k0:]*pkfold,k3[k0:]*errfold):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)

o = optparse.OptionParser()
#opts,args = o.parse_args(sys.argv[1:])
#args = glob.glob('data/pspec_t1_c110-149_sep*,1.npz')
args = glob.glob('data/pspec_v009_c110-149_sep*,1.npz')
#dspec = 'data/pspec_t1_c110-149_sep0,1_nocov.npz' # XXX need to update this
#args = sys.argv[1:]

FG_VS_KPL_NOS = 168.74e6
FG_VS_KPL = { # K^2 # XXX need to update this as well
    '-0.054':   5.37262e+13,
    '-0.027':   7.15304e+14, 
    ' 0.000':   3.50958e+15, 
    ' 0.027':   4.12396e+14, 
    ' 0.054':   2.60795e+13,
    #-0.0536455587089:   5.37262e+13,
    #-0.0268227793545:   7.15304e+14, 
    #0.0:                3.50958e+15, 
    #0.0268227793545:    4.12396e+14, 
    #0.0536455587089:    2.60795e+13,
}

RS_VS_KPL = {} # K^2
dsum, dwgt = {}, {}
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    RS_VS_KPL[filename] = {}
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    for _kpl, _pk, _err in zip(kpl, pk, err):
        RS_VS_KPL[filename][_kpl] = (_pk, _err)
        dsum[_kpl] = dsum.get(_kpl, 0) + _pk / _err**2
        dwgt[_kpl] = dwgt.get(_kpl, 0) + 1 / _err**2
#RS_VS_KPL = {}
if True:
    RS_VS_KPL['total'] = {}
    for _kpl in dsum:
        RS_VS_KPL['total'][_kpl] = (dsum[_kpl] / dwgt[_kpl], 1./n.sqrt(dwgt[_kpl]))

if False: # put in raw delay spec
    f = n.load(dspec)
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    if True: #if 'I' in sep: # Add foregrounds
        for cnt,k in enumerate(kpl):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            pk[cnt] += FG_VS_KPL[k]
            err[cnt] = 2*n.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
    pk *= .76 # psa747 calibration of Pic A = 370.6 Jy @ 160 MHz (which includes resolution effects)
    err *= .76 
    pk *= 2.35 # Use power**2 beam, which is a 1.69/0.72=2.35 penalty factor
    err *= 2.35
    if False: # For aggressive fringe-rate filtering, change beam area
        pk *= 1.90 # ratio of power**2 beams for filtered * unfiltered beams: 0.306 / 0.162
        err *= 1.90
    dual_plot(kpl, pk, 2*err, color='c', upperlimit=True) # 2-sigma error bars

BINS = None
#d = [RS_VS_KPL[k] for k in kpl]
colors = 'kbcm' * 10
for sep in RS_VS_KPL:
    if not sep == 'total': continue
    dsum, dwgt = {}, {}
    ks = RS_VS_KPL[sep].keys(); ks.sort()
    for k in ks:
        _d,_n = RS_VS_KPL[sep][k]
        #print k, _d, _n, '->',
        dsum[k] = dsum.get(k,0) + _d / _n**2
        dwgt[k] = dwgt.get(k,0) + 1 / _n**2
        #print dsum[k] / dwgt[k], 1./n.sqrt(dwgt[k])
        
    kpl = dsum.keys(); kpl.sort()
    d = [dsum[k]/dwgt[k] for k in kpl]
    nos = [1./n.sqrt(dwgt[k]) for k in kpl]
    #d = [RS_VS_KPL[k][0] for k in kpl]
    #nos = [RS_VS_KPL[k][1] for k in kpl]
    if True: #if 'I' in sep: # Add foregrounds
        for cnt,k in enumerate(kpl):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            d[cnt] += FG_VS_KPL[k]
            nos[cnt] = 2*n.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
    d,kpl,nos = n.array(d, dtype=n.complex), n.array(kpl), n.array(nos)
    # Currently calibrated to Pictor A @ 160 MHz = 424 Jy
    # To recalibrate to new Pic A, must multiply by square of ratio of fluxes
    #d *= 0.774 # Recalibrate to Pic A from Perley et al. 1997
    #d *= 1.125 # Recalibrate to Pic A from Slee 1995
    d *= .76 # psa747 calibration of Pic A = 370.6 Jy @ 160 MHz (which includes resolution effects)
    nos *= .76 
    d *= 2.35 # Use power**2 beam, which is a 1.69/0.72=2.35 penalty factor
    nos *= 2.35
    if False: # For aggressive fringe-rate filtering, change beam area
        d *= 1.90 # ratio of power**2 beams for filtered * unfiltered beams: 0.306 / 0.162
        nos *= 1.90
    if True: # compensate for noise attenuation in lstbin outlier cut
        nos *= 1.586
    #for _kpl,_pk,_nos in zip(kpl,d,nos): print _kpl, _pk, _nos
    print sep, colors[0]
    dual_plot(kpl, d, 2*nos, color=colors[0]) # 2-sigma error bars
    colors = colors[1:] + colors[0]

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

if False:
    import glob
    re_z = re.compile(r'power_21cm.*_z(\d+\.\d+).*\.dat')

    for filename in glob.glob('lidz_mcquinn_k3pk/*7.3*dat'):
    #for filename in glob.glob('lidz_mcquinn_k3pk/*7.0*dat'):
    #for filename in glob.glob('lidz_mcquinn_k3pk/hmass/power*hmass*dat'):
        print 'Reading', filename
        d = n.array([map(float, L.split()) for L in open(filename).readlines()])
        ks, pk = d[:,0], d[:,1]
        z_file = float(re_z.match(os.path.basename(filename)).groups()[0])
        z = C.pspec.f2z(.164)
        k3pk = ks**3 / (2*n.pi**2) * pk
        #p.subplot(122)
        p.plot(ks, k3pk * mean_temp(z)**2, 'm-')
    
tau_h = 100. + 15 # ns
k_h = C.pspec.dk_deta(C.pspec.f2z(.164)) * tau_h
#p.subplot(121)
#p.plot([k_h, k_h], [1e-10, 1e20], 'k--')
#p.plot([-k_h, -k_h], [1e-10, 1e20], 'k--')
#p.gca().set_yscale('log', nonposy='clip')
#p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
#p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$')
#p.ylim(1e5,3e16)
#p.grid()


#p.subplot(122)
p.plot([k_h, k_h], [1e-10, 1e20], 'k--')
#p.plot([.1,.13,.17,.19,.27,.31,.4,.5,.6][:-1], [2e5,4e5,1e5,2e5,2e5,4e5,5e5,248**2,3e5][:-1], 'yv', label='GMRT2013')
#p.plot(ks, (ks/0.15)**3 * 484, 'c--', label='noise_sim')
#p.plot(ks, (ks/0.15)**3 * 484 * 55./75, 'c--', label='noise_sim')


p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$')
p.ylim(3e-1,3e10)
p.xlim(0, 0.6)
p.grid()
p.show()
