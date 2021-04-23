# MOLUSC Survivors Plotter
# Mackenna Wood, UNC Chapel Hill
import numpy as np
import math as math
import matplotlib.pyplot as plt
from astropy.table import Table
import scipy.optimize
from scipy.ndimage.filters import gaussian_filter
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Georgia']
rcParams['font.weight'] = 'bold'
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.linewidth'] = 1

# Input files
survivors_file = 'output_kept.csv'
all_file = 'output_all.csv'
# Output files
out_file = 'output_corner.pdf'  # writeout file for the corner plots
out_file2 = 'output_dtct_lims.pdf' # writeout file for the detection limit plots
out_file3 = 'output_srv.pdf' # writeout file for the survivor plots

# Convenience functions
def jup_mass_to_sol(jupiter_mass):
    return 9.457e-4*jupiter_mass

def sol_mass_to_jup(solar_mass):
    return solar_mass/9.457e-4

def period_to_a(per):
    # Returns the semi major axis (in AU) for a given period (in days) assuming an equal mass binary with solar mass
    # stars
    G = 39.478  # Gravitational constant in AU^3/years^2*M_solar
    return (np.divide(per, 365)**2 * G * 2/(4*np.pi**2))**(1 / 3)

def a_to_period(a):
    # Returns period (days) of an equal mass, solar-mass binary with a given semi-major axis
    G = 39.478 # Gravitational constant in AU^3/years^2*M_solar
    return np.sqrt((4*np.pi**2 * a**3) / (2*G)) * 365
  
  def corner(file_in, file_out=None, given_params='auto', n_gen=5000000, smoothing=False, color='blue'):
    # Creates a corner plot showing period, mass ratio, eccentricity and inclination
    # Plot is a 4x4 array of subplots, with the bottom left triangle showing contour plots of 2D parameter spaces
    # The diagonal subplots show histograms of each of the parameters. The top right triangle lists the surviving
    # fraction (i.e. len(file_in)/n_gen
    # Inputs:
    #       file_in - this is the "kept" csv output by MOLUSC
    #       file_out - optional, name of the file to write the plot to
    #       given_params - optional, can be "all", "auto" or a list of parameters to plot e.g. ['P','e']. "all" is
    #           equivalent to ['P','e','q','cos_i']. "auto" chooses parameters based on their ranges, keeping those that
    #           have not been limited to single values
    #       n_gen - number of companions generated to produce the output, used to calculate survivor fraction
    #       smoothing - boolean, decides whether to apply gaussian filter smoothing
    #       color - decides what color scheme to make the plot, can be "blue", "purple", "green", or "gray"
    # Outputs:
    #       corner plot, saved to file_out if applicable
    fs = 18  # fontsize
    rcParams['xtick.labelsize'] = 'large'
    rcParams['ytick.labelsize'] = 'large'

    # Arrange Data
    t = Table.read(file_in, format='ascii.csv')

    P = t['period(days)']
    mass_ratio = t['mass ratio']
    e = t['eccentricity']
    cos_i = abs(t['cos_i'])
    a  = t['semi-major axis(AU)']
    # Get data ranges
    P_range = [min(P), max(P)]
    e_range = [min(e), max(e)]
    a_range =  [min(a), max(a)]
    i_range = [min(cos_i), max(cos_i)]
    m_range = [min(mass_ratio), max(mass_ratio)]

    # Choose color scheme
    if color == 'blue':  # default option
        c1 = '#2a567b'  # border color
        c2 = '#9dc9ee'  # face color
        colors2 = ['white', '#9dc9ee', '#2a567b']  # contour colors for 2-level
        colors3 = ['white', '#9dc9ee', '#4d9de0', '#2a567b']  # contour colors
    elif color == 'purple':
        c1 = '#483b57'
        c2 = '#c6bcd3'
        colors2 = ['white', '#c6bcd3', '#483b57']
        colors3 = ['white', '#c6bcd3', '#836c9f', '#483b57']
    elif color == 'green':
        c1 = '#26724a'
        c2 = '#94d5b2'
        colors2 = ['white', '#94d5b2', '#26724a']
        colors3 = ['white', '#94d5b2', '#3bb273', '#26724a']
    elif color == 'gray':
        c1 = '#3b3b3b'
        c2 = '#b9b9b9'
        colors2 = ['white', '#b9b9b9', '#3b3b3b']
        colors3 = ['white', '#b9b9b9', 'gray', '#3b3b3b']
    else:
        print('Unrecognized color. Using blue')
        c1 = '#2a567b'  # border color
        c2 = '#9dc9ee'  # face color
        colors2 = ['white', '#9dc9ee', '#2a567b']  # contour colors for 2-level
        colors3 = ['white', '#9dc9ee', '#4d9de0', '#2a567b']  # contour colors

    # Determine which parameters to plot
    params = []
    data_list = []

    if given_params == 'all':
        params = ['P (days)', 'e', 'cos(i)', 'q']
        data_list = [P, e, cos_i, mass_ratio]
    elif given_params == 'auto':
        if P_range[1]-P_range[0] > 1:
            params.append('P (days)')
            data_list.append(P)
        if e_range[1]-e_range[0] > 0.05:
            params.append('e')
            data_list.append(e)
        if i_range[1]-i_range[0] > 0.1:
            params.append('sin(i)')
            data_list.append(cos_i)
        if m_range[1] - m_range[0] > 0.05:
            params.append('q')
            data_list.append(mass_ratio)
    else:
        list_params = given_params
        for x in list_params:
            if x == 'P':
                params.append('P (days)')
                data_list.append(P)
            elif x == 'e':
                params.append('e')
                data_list.append(e)
            elif x == 'cos(i)' or x == 'cos_i' or x == 'cos i':
                params.append('cos(i)')
                data_list.append(cos_i)
            elif x == 'mass ratio' or x == 'q':
                params.append('q')
                data_list.append(mass_ratio)

    N = len(params)
    if N < 2:
        print('More than one unfixed parameter is required to make a corner plot.')
        return -1

    # Create data structure
    data = np.vstack([np.transpose(x) for x in data_list]).transpose()

    # Create nxn diagonal array of subplots
    fig, axes = plt.subplots(N, N, figsize=(9, 7))

    # Changes to apply to all axes
    plt.subplots_adjust(wspace=.05, hspace=.05)  # Change spacing between axes
    for i in range(N):
        for j in range(N):
            # Change axes line width
            for axis in ['top', 'bottom', 'left', 'right']:
                axes[i, j].spines[axis].set_linewidth(2)
            # Change tick params
            axes[i, j].tick_params(axis='both', direction='in', top=True, right=True, labelbottom=False, labelleft=False, width=1.5, length=4)
            axes[i, j].tick_params(which='minor', length=0)
            # Add Axis Labels
            if i == (N-1):  # x axis labels
                axes[i, j].tick_params(axis='x', labelbottom=True, labeltop=False, labelrotation=45)
            if j == 0:  # y axis labels
                axes[i, j].tick_params(axis='y', labelleft=True, labelright=False, labelrotation=45)

    # Diagonal: Plot Survivor Histograms
    n_bins = 15     # number bins for each histogram
    for i in range(N):
        # For each parameter I need to choose log or linear, set the bins and x scale appropriately and set x bounds
        # For all parameters I plot the histogram, and the 16th, 50th and 84th percentile lines, and set the y bounds
        bins = n_bins  # default
        if params[i] == 'P (days)':
            # Period, log normal
            bins = np.logspace(np.log10(P_range[0]), np.log10(P_range[-1]), n_bins)
            axes[i, i].set_xscale('log')
            axes[i, i].set_xbound(P_range[0], P_range[1])
            if math.ceil(np.log10(a_range[1])) - math.floor(np.log10(a_range[0])) >= 15:
                ticks = np.arange(math.floor(np.log10(P_range[0])), math.ceil(np.log10(P_range[1])) + 1, 4)[1:-1]
            elif math.ceil(np.log10(a_range[1])) - math.floor(np.log10(a_range[0])) >= 10:
                ticks = np.arange(math.floor(np.log10(P_range[0])), math.ceil(np.log10(P_range[1])) + 1, 3)[1:-1]
            elif math.ceil(np.log10(a_range[1])) - math.floor(np.log10(a_range[0])) > 5:
                ticks = np.arange(math.floor(np.log10(P_range[0])), math.ceil(np.log10(P_range[1])) + 1, 2)[1:-1]
            else:
                ticks = np.arange(math.floor(np.log10(P_range[0])), math.ceil(np.log10(P_range[1]))+1, 1)[1:-1]
            ticks = [10.**x for x in ticks]
            character = 'P'
        elif params[i] == 'e':
            # Eccentricity, uniform
            bins = np.linspace(e_range[0], e_range[1], n_bins)
            axes[i, i].set_xbound(e_range[0], e_range[1])
            ticks = np.linspace(math.floor(e_range[0]), math.ceil(e_range[1]), 5)[1:-1]
            character = 'e'
        elif params[i] == 'a (AU)':
            # Separation, log normal
            bins = np.logspace(np.log10(a_range[0]), np.log10(a_range[-1]), n_bins)
            axes[i, i].set_xscale('log')
            axes[i, i].set_xbound(a_range[0], a_range[1])
            if math.ceil(np.log10(a_range[1])) - math.floor(np.log10(a_range[0])) > 15:
                ticks = np.arange(math.floor(np.log10(a_range[0])), math.ceil(np.log10(a_range[1])) + 1, 4)[1:-1]
            elif math.ceil(np.log10(a_range[1])) - math.floor(np.log10(a_range[0])) > 11:
                ticks = np.arange(math.floor(np.log10(a_range[0])), math.ceil(np.log10(a_range[1])) + 1, 3)[1:-1]
            elif math.ceil(np.log10(a_range[1])) - math.floor(np.log10(a_range[0])) > 6:
                ticks = np.arange(math.floor(np.log10(a_range[0])), math.ceil(np.log10(a_range[1])) + 1, 2)[1:-1]
            else:
                ticks = np.arange(math.floor(np.log10(a_range[0])), math.ceil(np.log10(a_range[1]))+1, 1)[1:-1]
            ticks = [10.**x for x in ticks]
            character = 'a'
        elif params[i] == 'cos(i)':
            # sin(i), uniform
            bins = np.linspace(i_range[0], i_range[1], n_bins)
            axes[i, i].set_xbound(i_range[0], i_range[1])
            ticks = np.linspace(math.floor(i_range[0]), math.ceil(i_range[1]), 5)[1:-1]
            character = 'cos(i)'
        elif params[i] == 'q':
            # mass ratio, uniform
            bins = np.linspace(m_range[0], m_range[1], n_bins)
            axes[i, i].set_xbound(round(m_range[0]), round(m_range[1]))
            ticks = np.linspace(math.floor(m_range[0]), math.ceil(m_range[1]), 5)[1:-1]
            character = 'M'

        # Plot for all
        n, _ = np.histogram(data[:, i], bins=bins)
        widths = [bins[i + 1] - bins[i] for i in range(len(bins) - 1)]
        axes[i, i].bar(bins[:-1], n/sum(n), widths, align='edge', edgecolor=c1, linewidth=3, color=c2)
        axes[i, i].bar(bins[:-1], n/sum(n), widths, align='edge', linewidth=0, color=c2)
        lims = axes[i, i].get_ylim()
        # plot median and 16th and 84th percentiles
        axes[i, i].vlines([np.median(data[:, i]), np.percentile(data[:, i], 16.0), np.percentile(data[:, i], 84.0)],
                          lims[0], lims[1], linestyle='--', color='k')
        axes[i, i].set_ybound(lims[0], lims[1])

        # Set x ticks
        axes[i, i].set_xticks(ticks)

        # Add title printing the median and 16th and 84th percentile
        # upper = np.percentile(data[:, i], 84.0) - np.median(data[:, i])
        # lower = np.median(data[:, i]) - np.percentile(data[:, i], 16.0)
        # title = character + r'= %.3f $\pm _{%.2f}^{%.2f}$' %(np.median(data[:, i]), lower, upper)
        # axes[i, i].set_title(title)

    # Add the x-axis label on the lowest diagonal, and y-axis label on highest diagonal
    axes[N-1, N-1].set_xlabel(params[N-1], fontsize=fs)
    axes[0, 0].set_ylabel('Frequency', fontsize=fs)

    # Lower Triangle: Plot Parameters
    if len(P) <= 5000:  # For plots with less than 5000 survivors use a scatterplot
        for i in range(N):
            axes[i, 0].set_ylabel(params[i], fontsize=fs)
            for j in range(i):
                # Plot Data
                # 1. Determine if either axis uses a log scale. Period and semi-major axis need log scales
                if params[i] == 'P (days)' or params[i] == 'a (AU)':
                    axes[i, j].set_yscale('log')
                if params[j] == 'P (days)' or params[j] == 'a (AU)':
                    axes[i, j].set_xscale('log')
                # 2. Plot the scatterplot
                axes[i, j].scatter(data[:,j],data[:,i], color=colors3[1], marker='.', s=0.3, alpha=.9)
                # 3. Add x-axis labels
                axes[(N-1), j].set_xlabel(params[j], fontsize=fs)
    else:  # For more than 5000 survivors use a 2d histogram
        for i in range(N):
            if i > 0:  # Set y-axis labels
                axes[i, 0].set_ylabel(params[i], fontsize=fs)
            for j in range(0, i):
                # 1. Determine if either axis uses a log scale. Period and semi-major axis need log scales
                if params[i] == 'P (days)' or params[i] == 'a (AU)':
                    ys = np.log10(data[:, i])
                else:
                    ys = data[:, i]
                if params[j] == 'P (days)' or params[j] == 'a (AU)':
                    xs = np.log10(data[:, j])
                else:
                    xs = data[:, j]
                if not smoothing:
                    # 2. Get the histogram counts and bins
                    b = 30  # number of bins in histogram
                    counts, x_bins, y_bins = np.histogram2d(xs, ys, bins=b)
                else:
                    # 2. Get the histogram counts and bins
                    b = 50  # number of bins in histogram
                    smoothing_factor = 1.1
                    counts, x_bins, y_bins = np.histogram2d(xs, ys, bins=b)
                    # 2b. Put the counts through a gaussian filter
                    counts = gaussian_filter(counts, smoothing_factor)
                # 3. Check the number of unique contours and set color scale accordingly
                ls = np.unique([np.percentile(counts, 39.3), np.percentile(counts, 86), np.percentile(counts, 98)])
                if len(ls) == 2:  c = colors2
                else:  c = colors3
                # 4. Plot the filled contour and contour lines
                if params[i] == 'P (days)' or params[i] == 'a (AU)':
                    ys = np.power(10, y_bins[:-1])
                    axes[i, j].set_yscale('log')
                else:
                    ys = y_bins[:-1]
                if params[j] == 'P (days)' or params[j] == 'a (AU)':
                    xs = np.power(10, x_bins[:-1])
                    axes[i, j].set_xscale('log')
                else:
                    xs = x_bins[:-1]
                axes[i, j].contourf(xs, ys, counts.transpose(), levels=ls, extend='both', colors=c)
                axes[i, j].contour(xs, ys, counts.transpose(), levels=ls, colors='k', alpha=0.8)
                # 5. Add x-axis labels
                axes[(N-1), j].set_xlabel(params[j], fontsize=fs)
                # 6. Adjust x-axis ticks
                axes[i, j].set_xticks(axes[j, j].get_xticks())

    # Upper Triangle:  Write survivor fraction, make invisible
    axes[0, N-2].text(0, 1, ('Surviving Fraction: %.3f' %(len(t)/n_gen)), fontsize=20)
    for i in range(N):
        for j in range(0, i):
            axes[j, i].axis('off')

    if file_out is not None:
        plt.savefig(file_out, bbox_inches='tight', pad_inches=0.25)
    plt.show()

    return 0


def detection_limits(file_in, star_mass, file_out=None, mark_P=None):
    # Plot 1,2 and 3 sigma mass detection limits
    # Inputs:
    #       file_in - this is the "kept" csv output by MOLUSC
    #       star_mass - primary star mass in solar masses
    #       file_out - optional, file name to write image to, does not write out image if not provided
    #       mark_P - optional, if included will mark vertical lines at the period provided
    # Output:
    #       plot with detection limits
    # Some extra image parameters
    rcParams['xtick.labelsize'] = 'small'
    rcParams['ytick.labelsize'] = 'small'
    rcParams['axes.labelsize'] = 10

    fig, ax = plt.subplots(figsize=(3.352242, 2.514181))  # this is sized to fit in one column of AAS journal format
    # Read in survivor data
    survivors = Table.read(file_in, format='ascii.csv')
    # Need to split into logarithmic period bins and take the 95th percentile of the masses in each bin
    survivors['mass (M_Jup)'] = survivors['mass ratio']*star_mass*1047.35

    # Bin
    bins = np.arange(0., 10., .5)
    n, _ = np.histogram(np.log10(survivors['period(days)']), bins=bins)
    # Create a subchart for each bin. Find the 95 percentile mass in that bin, and the median period
    period = []
    ninety_five = []
    ninety = []
    sixty = []
    n = []
    for j in range(1, len(bins)):
        survivor_bin = survivors[[bins[j-1] <= np.log10(x) < bins[j] for x in survivors['period(days)']]]
        if len(survivor_bin) > 0:
            ninety_five.append(np.percentile(survivor_bin['mass (M_Jup)'], 95))
            ninety.append(np.percentile(survivor_bin['mass (M_Jup)'], 90))
            sixty.append(np.percentile(survivor_bin['mass (M_Jup)'], 60))
            period.append(np.median(survivor_bin['period(days)']))
            n.append(len(survivor_bin))
    # Plot
    ax.plot(period, ninety_five, marker='o', ms=4, color='#4d9de0', label='95th percentile')
    ax.plot(period, ninety, marker='o', ms=4,  color='#e15554', label='90th')
    ax.plot(period, sixty, marker='o', ms=4, color='#836c9f', label='60th')
    ax.set_ylim(0.6, 1.2*ax.get_ylim()[1])
    if mark_P:
        ax.vlines([mark_P], ax.get_ylim()[0], ax.get_ylim()[1], ls=':', colors='k', lw=2)
    ax.vlines([6.4e5], ax.get_ylim()[0], ax.get_ylim()[1], ls='--', colors='k', lw=2)  # DS Tuc B

    # Add a secondary axis, showing mass in solar masses
    secax = ax.secondary_yaxis('right', functions=(jup_mass_to_sol, sol_mass_to_jup))
    secax.set_ylabel(r'Mass ($M_\odot$)')
    # Add a secondary x axis, showing semi-major axis for a solar mass companion
    triax = ax.secondary_xaxis('top', functions=(period_to_a, a_to_period))
    triax.set_xlabel('a (AU)')
    triax.tick_params(axis='x', which='minor', bottom=False, top=False)
    # Correct axis scales, labels, and ticks
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Period (days)')
    ax.set_ylabel(r'Mass $(M_{Jup})$')
    plt.xticks()
    plt.yticks()
    # Add legend
    plt.legend(fontsize=6)

    # Write out
    if file_out:
        plt.savefig(file_out, bbox_inches='tight', pad_inches=0.2)

    plt.show()

    return

def survivor(survivors_file, all_file, param, file_out=None):
  # Creates three panel plot with histograms showing generated, and surviving parameter distribution and bar plot showing survivor fraction
  # Inputs:
  #     survivors_file - filename for kept companions. This is the _kept.csv generated by  MOLUSC
  #     all_file - filename for all companions. This is the _all.csv generated by MOLUSC
  #     param  - orbital parameter to make the plot for. Can be either period (P)  or mass ratio (m)
  #     file_out - optional, filename to write image to when complete
  # Load kept and survivors
  t = Table.read(survivors_file, format='ascii.csv')
  t_all = Table.read(all_file, format='ascii.csv')
  # Figure Options
  rcParams['xtick.labelsize'] = 'small'
  rcParams['ytick.labelsize'] = 'small'
  rcParams['axes.labelsize'] = 10
  rcParams['hatch.linewidth'] = 1.5

  fig, axes = plt.subplots(3, 1, sharex='col', figsize=(3.352242, 3.352242)) #  sized for one column of AAS journal format

  if param == 'P' or param == 'period':
      # Plot two histograms, top one a histogram of period (with log bins), bottom one a plot of survivorship fraction
      bins = np.logspace(0, 10, 20)
      n, _ = np.histogram(t['period(days)'], bins=bins)
      N, _ = np.histogram(t_all['period(days)'], bins=bins)
      scale = 'log'
      axes[2].set_xlabel('Period(days)')
  elif param == 'm' or param == 'mass ratio' or param == 'mass_ratio':
      # Plot two histograms, top one a histogram of mass, bottom one a plot of survivorship fraction
      bins = np.linspace(0, 1, 20)
      n, _ = np.histogram(t['mass ratio'], bins=bins)
      N, _ = np.histogram(t_all['mass ratio'], bins=bins)
      scale = 'linear'
      axes[2].set_xlabel('Mass Ratio')

  widths = [bins[i + 1] - bins[i] for i in range(len(bins) - 1)]

  # Top plot, generated histogram
  steps = np.concatenate((N/sum(N), [(N/sum(N))[-1]]))  # steps will stop at the start of the last bar unless repeated
  axes[0].step(bins, steps, where='post', c='k', linewidth=2)
  axes[0].bar(bins[:-1], N / sum(N), widths, align='edge', linewidth=1, color='#4d9de0', edgecolor='#4d9de0')
  # middle plot, surviving histogram
  steps = np.concatenate((n/sum(N), [(n/sum(N))[-1]]))
  axes[1].step(bins, steps, where='post', c='k', linewidth=2)
  axes[1].bar(bins[:-1], n / sum(N), widths, align='edge', linewidth=1, color='#4d9de0', edgecolor='#4d9de0')
  # Bottom plot, surviving fraction
  steps = np.concatenate((n/N, [(n/N)[-1]]))
  axes[2].step(bins, steps, where='post', c='k', linewidth=2)
  axes[2].bar(bins[:-1], n / N, widths, align='edge', linewidth=1, color='#4d9de0', edgecolor='#4d9de0')

   # Set axis labels and limits
  axes[0].set_ylabel('Generated\nFrequency')
  axes[1].set_ylabel('Survivor\nFrequency')
  axes[2].set_ylabel('Surviving\nFraction')
  axes[0].set_xlim(bins[0], bins[-1])
  axes[1].set_xlim(bins[0], bins[-1])
  axes[2].set_xlim(bins[0], bins[-1])

  # Add a secondary x axis, showing semi-major axis for a solar mass companion
  secax = axes[0].secondary_xaxis('top', functions=(period_to_a, a_to_period))
  secax.set_xlabel('a (AU)')
  secax.tick_params(axis='x', which='minor', bottom=False, top=False)

  # Set to log scale if plotting period, or linear for mass ratio
  for a in axes:
      a.set_xscale(scale)

  if file_out:
      plt.savefig(file_out, bbox_inches='tight', pad_inches=0.2)
  plt.show()
  return
