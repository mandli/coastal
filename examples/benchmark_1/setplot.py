
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy

# Plot customization
import matplotlib

# Markers and line widths
matplotlib.rcParams['lines.linewidth'] = 2.0
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['lines.markersize'] = 8

# Font Sizes
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.labelsize'] = 15
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12

# DPI of output images
# matplotlib.rcParams['savefig.dpi'] = 300 # Publication quality
matplotlib.rcParams['savefig.dpi'] = 100

import matplotlib.pyplot as plt
import clawpack.pyclaw.solution as solution


# Color and linestyles
rgb_converter = lambda triple: [float(rgb) / 255.0 for rgb in triple]
water_color = rgb_converter((67,183,219))
bathy_linestyle = '-'
surface_linestyle = '-'


def add_legend(axes,label,location=0,color='r',linestyle='-'):
    r""""""
    
    # Create new line for legend
    line = plt.Line2D((0,0),(0,1),color=color,linestyle=linestyle)
    
    # Append extra legend entry to list of handles and labels
    handles,labels = axes.get_legend_handles_labels()
    handles.append(line)
    labels.append(label)
    
    # Add legend to axes
    axes.legend(handles,labels,loc=location)


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    # Read in plotdata
    with open("plot.data", 'r') as plot_data:
        line = plot_data.readline().split()
        xlimits = [float(value) for value in line]
        file_format = plot_data.readline().strip()
        dry_tolerance = float(plot_data.readline())

    # Load bathymetry
    b = solution.Solution(0, path=plotdata.outdir, read_aux=True).state.aux[0,:]

    def bathy(cd):
        return b
        

    def eta(cd):
        return cd.q[0,:] + bathy(cd)
        

    def u(cd):
        index = numpy.nonzero(cd.q[0, :] > dry_tolerance)
        u = numpy.zeros(cd.q[0, :].shape)
        u[index] = cd.q[1, index] / cd.q[0, index]
        return u


    plotdata.clearfigures()
    plotdata.format = file_format
    plotdata.output_controller.file_format = file_format

    # ========================================================================
    #  Generic Helper Functions
    # ========================================================================
    def pcolor_afteraxes(cd):
        pass
        
    def contour_afteraxes(cd):
        pass
        
    def bathy_ref_lines(cd):
        pass

    # ========================================================================
    #  Plot Limits
    # ========================================================================
    # These are now passed into the function besides this simple one
    depth_limits = [-1.1, 1.0]
    speed_limits = [-1.0, 1.0]
    momentum_limits = [-1.0, 1.0]

    # ========================================================================
    #  Surface Elevation
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name="Depth and Speed",figno=0)
    plotfigure.show = True
    
    # Custom axes for this to look nice
    def twin_axes(cd):
        fig = plt.gcf()
        fig.clf()
        
        # Get x coordinate values
        x = cd.patch.dimensions[0].centers
        
        # Create axes for each plot, sharing x axis
        depth_axes = fig.add_subplot(211)
        vel_axes = fig.add_subplot(212) #,sharex=depth_axes)     # the velocity scale
        
        # Water
        depth_axes.fill_between(x, bathy(cd), eta(cd),color=water_color)
        # Plot bathy
        depth_axes.plot(x, bathy(cd), 'k', linestyle=bathy_linestyle)
        # Plot surface
        depth_axes.plot(x, eta(cd), 'k', linestyle=surface_linestyle)
        
        # Remove ticks from top plot
        num_ticks = len(depth_axes.xaxis.get_ticklocs())
        depth_axes.xaxis.set_ticklabels(["" for n in xrange(num_ticks)])
        
        # ax1.set_title('')
        depth_axes.set_title('Depth and Speed at t = %3.2f' % (cd.t))
        depth_axes.set_xlim(xlimits)
        depth_axes.set_ylim(depth_limits)
        depth_axes.set_ylabel('Depth (m)')
        
        # Velocity
        vel_axes.plot(x, u(cd), 'b', linestyle=surface_linestyle, label="Velocity")

        # Add legend
        vel_axes.legend(loc=4)
        vel_axes.set_title('')
        vel_axes.set_ylabel('Velocity (m/s)')
        vel_axes.set_xlabel('x (m)')
        vel_axes.set_xlim(xlimits)
        vel_axes.set_ylim(speed_limits)
        
        # Add axis labels (not sure why this needs to be done)
        locs = vel_axes.xaxis.get_ticklocs()
        vel_axes.xaxis.set_ticklabels([str(loc) for loc in locs])
        
        # This does not work on all versions of matplotlib
        # try:
        #     plt.subplots_adjust(hspace=0.1)
        # except:
        #     pass
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.afteraxes = twin_axes

    # ========================================================================
    #  Momentum
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name="momentum",figno=134)
    plotfigure.show = True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Momentum"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = momentum_limits
    
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = 'b-'
    plotitem.show = True

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = [1,2]  # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
