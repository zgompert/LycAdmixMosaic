import moments
import scipy
import numpy
import matplotlib
import matplotlib.pyplot as pyplot
import sys
import re
import pylab

## vcf
vv = "/uufs/chpc.utah.edu/common/home/gompert-group5/projects/LycAdmix/WingMap/Genetics/vars/split/GNP.vcf";

## downsample to 50
n = 50

ns = numpy.array([int(n)])

## read data from filtered vcf and create downsampled/projected 1d sfs 
dd = moments.Misc.make_data_dict_vcf(vv,"GNP_ids.txt")
pop_id = ['GNP']
fs = moments.Spectrum.from_data_dict(dd, pop_id, projections=ns, polarized=False)
fs_proj = fs 

of = "dd_"+pop_id[0]+"_2Ne.txt"
pngf = "dd_"+pop_id[0]+"_2Ne.png"

ofile = open(of,"w")

ns =  fs_proj.sample_sizes

## define single population two values for Ne, sudden shift at time T
def Ne2(params, ns):
    nu, T = params
    # Create equilibrium neutral model
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    # Instantaneous change
    fs.integrate([nu], T)
    return fs

func = Ne2

# Parameters are: (nu and T)
upper_bound = [100, 10]
lower_bound = [0.1, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [5, 0.4]
np = len(p0)
p0 = moments.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
# powell method for apporximate solution
popt = moments.Inference.optimize_log_powell(p0, fs_proj, func, lower_bound=lower_bound, upper_bound=upper_bound,verbose=5, maxiter=10)

# optimization
ofile = open(of,"w")
print(ofile)

## initial defaults
ll_opt = -99999999999
theta = 1.0
mpopt = p0
popt = p0
for x in range(20):
    theta = 1.0
    popt = moments.Misc.perturb_params(popt, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
    # powell method for apporximate solution
    popt = moments.Inference.optimize_log_powell(popt, fs_proj, func, lower_bound=lower_bound, upper_bound=upper_bound,verbose=5, maxiter=30)
    # BFGS for final fit
    popt = moments.Inference.optimize_log(popt, fs_proj, func, lower_bound=lower_bound, upper_bound=upper_bound,verbose=5, maxiter=30)
    model = func(popt, ns)
    ll_model = moments.Inference.ll_multinom(model,fs_proj)
    theta = moments.Inference.optimal_sfs_scaling(model,fs_proj)
    ## write current
    ofile.write("{0}".format(ll_model))
    ofile.write(" {0}".format(theta))
    for a in range(np):
        ofile.write(" {0}".format(popt[a]))
        
    ofile.write("\n")

    if(ll_model > ll_opt):
        ll_opt = ll_model
        mpopt = popt		
        theta_opt = moments.Inference.optimal_sfs_scaling(model,fs_proj)

## write max
ofile.write("{0}".format(ll_opt))
ofile.write(" {0}".format(theta_opt))
for a in range(np):
    ofile.write(" {0}".format(mpopt[a]))

ofile.write("\n")
ofile.close()


pylab.figure(figsize=(8,6))
moments.Plotting.plot_1d_comp_multinom(model,fs_proj,show=False)
pylab.savefig(pngf, dpi=400)

