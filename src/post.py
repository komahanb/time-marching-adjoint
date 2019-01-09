import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import axes3d
import numplot.core as npl
import numpy as np

# These are the "Colors 20" colors as RGB.    
colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(colors)):    
    r, g, b = colors[i]    
    colors[i] = (r / 255., g / 255., b / 255.)
    
def exact_solution(t):
    return 50 - 2*np.exp(-0.196*t)

def get_error(sets, xkey, ykey, exact):
    """
    """
    error_data = {}
    keys = sets.keys()
    for key in keys:
        data = sets[key]
        u = exact(data[xkey])

        # Create a data map for error
        error = {}
        error[xkey] = data[xkey]
        error['error'] = abs(u - data[ykey])

        # Add this data to map
        error_data[key] = error
        
    return error_data

def get_rmse(sets):
    """
    """
    rmse = {}
    keys = sets.keys()
    for key in keys:
        error = sets[key]['error']
        val = np.sqrt(np.sum(error**2)/(len(error)))
        rmse[key] = val
    return rmse

def plot_solution(sets, xkey, ykey, name):
    '''
    Plot x vs u for different data sets
    '''    
    # Plot
    plt.figure()    
    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')   
    #ax.annotate('192,000 dof', xy=(4, 500))
    #ax.annotate('2 million dof', xy=(25, 10000))
    #plt.axis([5, 45, -.1, 1.0])
    #plt.xticks(np.arange(5, 50, step=5))
        
    # plot solutions on same graph
    keys = sets.keys()
    cidx = 0
    for key in keys:
        cidx += 2
        data = sets[key]
        plt.semilogy(data[xkey], data[ykey], '-' , lw=3, mec='black', label=key, color=colors[cidx])
        
    # Axis formatting
    plt.legend(loc='upper right')
    plt.xlabel(xkey)
    plt.ylabel(ykey)

    plt.savefig(name, bbox_inches='tight', pad_inches=0.05)
        
    return

def plot_refinement(timer, name):
    '''
    Plot computational time for jacobi, seidel and sor
    '''
    plt.figure()
    
    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')   
    plt.axis([100, 10000, 1.0e-5, 1.0e-1])

    #plt.loglog(timer['npts'], timer['imp-euler'] , '-s' , lw =2, mec='black', label='Implicit Euler', color=colors[8])
    #plt.loglog(timer['npts'], timer['cni']       , '-o' , lw =2, mec='black', label='Crank Nicolson', color=colors[10])               
    plt.loglog(timer['npts'], timer['dirk2'] , '-v' , lw =2, mec='black', label='DIRK2', color=colors[2])
    plt.loglog(timer['npts'], timer['dirk3'] , '-D' , lw =2, mec='black', label='DIRK3', color=colors[14])
    plt.loglog(timer['npts'], timer['dirk4'] , '-*' , lw =2, mec='black', label='DIRK4', color=colors[16])    

    plt.ylabel('RMSE') 
    plt.xlabel('1/h')
    plt.legend(loc='lower left',framealpha=0.0)
    
    plt.savefig(name, bbox_inches='tight', pad_inches=0.05)
    
    return

def find_rmse(files):    
    sets = {}
    fnames = files.keys()
    for fname in fnames:
        sets[fname] = npl.Map(files[fname])
    
    # plot solution vs time
    #plot_solution(sets, 'time', 'u', 'solution.pdf')
    
    # plot error vs time
    error_set = get_error(sets, 'time', 'u', exact_solution)
    #plot_solution(error_set, 'time', 'error', 'error.pdf')
    
    # Get RMSE
    rmse = get_rmse(error_set)
    
    return rmse

def get_files(n):
    files = {}
#    files['bdf1']  = 'fode-bdf1-'  + n + '.dat'
#    files['bdf2']  = 'fode-bdf2-'  + n + '.dat'
#    files['bdf3']  = 'fode-bdf3-'  + n + '.dat'
    files['dirk2'] = 'fode-dirk2-' + n + '.dat'
    files['dirk3'] = 'fode-dirk3-' + n + '.dat'
    files['dirk4'] = 'fode-dirk4-' + n + '.dat'
    #files['abm1']  = 'fode-abm1-'  + n + '.dat'
    #files['abm2']  = 'fode-abm2-'  + n + '.dat'
    #files['abm3']  = 'fode-abm3-'  + n + '.dat'
    return files

def get_numerical_orders(data, methodkey):
    """
    """
    orders = {}
    # unpack data based on size
    npts = data.keys()
    for pt in npts:
        # unpack data based on method
        methods = data[pt].keys()
        for method in methods:
            print pt, method, data[pt][method]            
    return orders

#####################################################################
# Plot solution
#####################################################################

sizes = ['10000', '5000', '2500', '1250']

RMSE = {}
for size in sizes:
    RMSE[size] = find_rmse(get_files(size))    
print RMSE

# Get order of convergence plots
orders = get_numerical_orders(RMSE)
print orders
stop

######################################################################

#timer = npl.Map("timing.dat")
#plot_time(timer, 'timing-study.pdf')

refine = npl.Map('spatial-error-case1.dat')
keys = refine.keys()
xkey = refine['npts'][-3:]
for key in keys:
    if key != 'npts' and key != 'h':
        val = refine[key][-3:]
        print key, np.log(val[2]/val[0])/np.log(xkey[0]/xkey[2])
plot_refinement(refine, 'spatial-error.pdf')

stop

refine = npl.Map('spatial-error-case2.dat')
keys = refine.keys()
xkey = refine['npts'][-3:]
for key in keys:
    if key != 'npts' and key != 'h':
        val = refine[key][-3:]
        print key, np.log(val[2]/val[0])/np.log(xkey[0]/xkey[2])
plot_refinement(refine, 'spatial-error-case2.pdf')

stop

######################################################################
nx     = 502
ny     = 3001
tindex = 0

exp_euler = npl.Map("case2-transport-explicit-euler.dat")
exp_euler_new = {}
exp_euler_new['time']  = exp_euler['time'].reshape((nx,ny))
exp_euler_new['x']     = exp_euler['x'].reshape((nx,ny))
exp_euler_new['state'] = exp_euler['state'].reshape((nx,ny))

imp_euler = npl.Map("case2-transport-implicit-euler.dat")
imp_euler_new = {}
imp_euler_new['time']  = imp_euler['time'].reshape((nx,ny))
imp_euler_new['x']     = imp_euler['x'].reshape((nx,ny))
imp_euler_new['state'] = imp_euler['state'].reshape((nx,ny))

cni = npl.Map("case2-transport-cni.dat")
cni_new = {}
cni_new['time']  = cni['time'].reshape((nx,ny))
cni_new['x']     = cni['x'].reshape((nx,ny))
cni_new['state'] = cni['state'].reshape((nx,ny))

dirk2 = npl.Map("case2-transport-implicit-dirk2.dat")
dirk2_new = {}
dirk2_new['time']  = dirk2['time'].reshape((nx,ny))
dirk2_new['x']     = dirk2['x'].reshape((nx,ny))
dirk2_new['state'] = dirk2['state'].reshape((nx,ny))

dirk3 = npl.Map("case2-transport-implicit-dirk3.dat")
dirk3_new = {}
dirk3_new['time']  = dirk3['time'].reshape((nx,ny))
dirk3_new['x']     = dirk3['x'].reshape((nx,ny))
dirk3_new['state'] = dirk3['state'].reshape((nx,ny))

dirk4 = npl.Map("case2-transport-implicit-dirk4.dat")
dirk4_new = {}
dirk4_new['time']  = dirk4['time'].reshape((nx,ny))
dirk4_new['x']     = dirk4['x'].reshape((nx,ny))
dirk4_new['state'] = dirk4['state'].reshape((nx,ny))

solutions = {}
solutions['exp_euler'] = exp_euler_new
solutions['imp_euler'] = imp_euler_new
solutions['cni']       = cni_new
solutions['dirk2']     = dirk2_new
solutions['dirk3']     = dirk3_new
solutions['dirk4']     = dirk4_new

# Case 2
plot_time_snaps(exp_euler_new, [1000,2000,3000], exact2, 'case2-timesnaps-expeuler.pdf')
plot_time_snaps(imp_euler_new, [1000,2000,3000], exact2, 'case2-timesnaps-impeuler.pdf')
plot_time_snaps(cni_new      , [1000,2000,3000], exact2, 'case2-timesnaps-cni.pdf')
plot_time_snaps(dirk2_new    , [1000,2000,3000], exact2, 'case2-timesnaps-dirk2.pdf')
plot_time_snaps(dirk3_new    , [1000,2000,3000], exact2, 'case2-timesnaps-dirk3.pdf')
plot_time_snaps(dirk4_new    , [1000,2000,3000], exact2, 'case2-timesnaps-dirk4.pdf')

plot_space_snaps(exp_euler_new, [125,250], [15.0, 25.0], exact2, 'case2-spacesnaps-expeuler.pdf')
plot_space_snaps(imp_euler_new, [125,250], [15.0, 25.0], exact2, 'case2-spacesnaps-impeuler.pdf')
plot_space_snaps(cni_new      , [125,250], [15.0, 25.0], exact2, 'case2-spacesnaps-cni.pdf')
plot_space_snaps(dirk2_new    , [125,250], [15.0, 25.0], exact2, 'case2-spacesnaps-dirk2.pdf')
plot_space_snaps(dirk3_new    , [125,250], [15.0, 25.0], exact2, 'case2-spacesnaps-dirk3.pdf')
plot_space_snaps(dirk4_new    , [125,250], [15.0, 25.0], exact2, 'case2-spacesnaps-dirk4.pdf')

stop

time1 = dirk4_new['time'][250,:]
plt.plt(time1, dirk4_new['state'][250,:])
plt.show()

#plot_space_snaps(imp_euler_new, [125,250], [15.0, 25.0], exact1, 'case2-spacesnap-impeuler.pdf')
#plot_space_snaps(cni_new      , [125,250], [15.0, 25.0], exact1, 'case2-spacesnap-cni.pdf')

stop

# Case 2
exp_euler = npl.Map("case2-transport-explicit-euler.dat")
exp_euler_new = {}
exp_euler_new['time'] = exp_euler['time'].reshape((nx,ny))
exp_euler_new['x'] = exp_euler['x'].reshape((nx,ny))
exp_euler_new['state'] = exp_euler['state'].reshape((nx,ny))

imp_euler = npl.Map("case2-transport-implicit-euler.dat")
imp_euler_new = {}
imp_euler_new['time'] = imp_euler['time'].reshape((nx,ny))
imp_euler_new['x'] = imp_euler['x'].reshape((nx,ny))
imp_euler_new['state'] = imp_euler['state'].reshape((nx,ny))

cni = npl.Map("case2-transport-cni.dat")
cni_new = {}
cni_new['time'] = cni['time'].reshape((nx,ny))
cni_new['x'] = cni['x'].reshape((nx,ny))
cni_new['state'] = cni['state'].reshape((nx,ny))

solutions = {}
solutions['exp_euler'] = exp_euler_new
solutions['imp_euler'] = imp_euler_new
solutions['cni'] = cni_new

stop

Z = jacobi['T'].reshape((nx,ny))

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

## plt.contour(X,Y,Z)
## plt.xlabel("x")
## plt.ylabel("y")
## plt.show()
## stop

ax.plot_surface(X, Y, Z, cmap="autumn_r", lw=0.5, rstride=1, cstride=1)
#ax.contour(X, Y, Z, 10, lw=3, cmap="autumn_r", linestyles="solid", offset=-1)
#ax.contour(X, Y, Z, 10, lw=3, colors="k", linestyles="solid")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("T")
plt.show()
stop

sol20 = npl.Map("linear-solution-npts-20.dat")
sol30 = npl.Map("linear-solution-npts-30.dat")
sol40 = npl.Map("linear-solution-npts-40.dat")
sol50 = npl.Map("linear-solution-npts-50.dat")

solution = {}
solution['sol10'] = sol10
solution['sol20'] = sol20
solution['sol30'] = sol30
solution['sol40'] = sol40
solution['sol50'] = sol50
plot_solution(solution, 'linear_solution.pdf')

linear_summary = npl.Map("linear_summary.dat")
nonlinear_summary1 = npl.Map("nonlinear_summary-case-1.dat")
nonlinear_summary2 = npl.Map("nonlinear_summary-case-2.dat")
nonlinear_summary3 = npl.Map("nonlinear_summary-case-3.dat")

# plot nonlinear solutions vs x
nonsol1 =  npl.Map("nonlinear-solution-npts-20-case-1.dat")
nonsol2 =  npl.Map("nonlinear-solution-npts-20-case-2.dat")
nonsol3 =  npl.Map("nonlinear-solution-npts-20-case-3.dat")
plot_nonlinear_solution(nonsol1, 'nonlinear_solution-case1.pdf')
plot_nonlinear_solution(nonsol2, 'nonlinear_solution-case2.pdf')
plot_nonlinear_solution(nonsol3, 'nonlinear_solution-case3.pdf')

stop

#spatial_rmse = npl.Map("imp-euler-spatial-error.dat")
#plot_rmse(spatial_rmse, "imp-euler-spatial-error.pdf")

stop

# t = np.arange(10, 40, step=.1)
# plt.plot(t, exact1(15,t))
# plt.plot(t, exact1(25,t))
# plt.show()
# stop



#plot_space_snaps(exp_euler_new, [125,250], [15.0, 25.0], exact1, 'case1-spacesnap-expeuler.pdf')

#stop

## avals = exp_euler_new['time'][250,:]
## xvals = exp_euler_new['x'][250,:]

## #plt.plot(avals, exp_euler_new['state'][::500,:])
## #plt.show()
## print avals
## print xvals
## stop
## for aa in avals:
##     print aa

