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

def plot_refinement(sets, xkey, ykey, name):
    '''
    Plot computational time for jacobi, seidel and sor
    '''
    plt.figure()
    
    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')   
    #plt.axis([100, 10000, 1.0e-5, 1.0e-1])

    # plot solutions on same graph
    keys = sets.keys()
    cidx = 0
    for key in keys:
        cidx += 2
        data = sets[key]
        plt.loglog(data[xkey], data[ykey], '-' ,
                   lw=3, mec='black',
                   label=key,
                   color=colors[cidx])
    # Axis formatting
    plt.legend(loc='lower left')
    plt.xlabel(xkey)
    plt.ylabel(ykey)

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

def invertmap(sets):
    # plot solutions on same graph
    keys = sets.keys()
    for key in keys:
        methods = sets[key].keys()
    xvals = []
    for key in keys:
        xvals.append(int(key))
    rmse = {}
    for key in methods:
        points = sets.keys()
        data = {}
        data['spacing'] = xvals
        yvals = []
        for point in points:
            yvals.append(sets[point][key])
        data['error'] = yvals
        rmse[key] = data
    return rmse

def get_files(n):
    files = {}
    files['bdf1']  = 'fode-bdf1-'  + n + '.dat'
    files['bdf2']  = 'fode-bdf2-'  + n + '.dat'
    files['bdf3']  = 'fode-bdf3-'  + n + '.dat'
    files['dirk2'] = 'fode-dirk2-' + n + '.dat'
    files['dirk3'] = 'fode-dirk3-' + n + '.dat'
    files['dirk4'] = 'fode-dirk4-' + n + '.dat'
    files['abm1']  = 'fode-abm1-'  + n + '.dat'
    files['abm2']  = 'fode-abm2-'  + n + '.dat'
    files['abm3']  = 'fode-abm3-'  + n + '.dat'
    return files

def get_numerical_orders(data):
    """
    """
    orders = {}
    # unpack data based on size
    methods = data.keys()
    for method in methods:
        spacings = data[method]['spacing']
        rmse = data[method]['error']
        p = np.log(rmse[0]/rmse[1])/np.log(1.0*spacings[0]/spacings[1])
        orders[method] = p
    return orders

#####################################################################
# Plot solution
#####################################################################

sizes = ['1250', '2500', '5000', '10000']

RMSE = {}
for size in sizes:
    RMSE[size] = find_rmse(get_files(size))    
print RMSE

rmse = invertmap(RMSE)
plot_refinement(rmse, 'spacing', 'error', 'convergence.pdf')

# Get order of convergence plots
orders = get_numerical_orders(rmse)
print orders
