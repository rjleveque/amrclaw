
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
"""

import os

outdir = os.path.abspath('_output')

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.clawutil.data import ClawData
    import os
    import glob

    plotdata.clearfigures() # clear any old figures,axes,items data

    print "plotdata.outdir = ",plotdata.outdir

    slice_dirs = glob.glob(outdir + '/slice_*')
    print "Found slice directories: ", slice_dirs

    slice_data = ClawData()
    slice_data.read(outdir + '/slices.data', force=True)


    # make a figure for each slice found, using the function defined below:

    for k,v in slice_data.iteritems():
        if k[0]=='z':
            make_plot_figure('xy',k,v,plotdata)
        elif k[0]=='y':
            make_plot_figure('xz',k,v,plotdata)
        elif k[0]=='x':
            make_plot_figure('yz',k,v,plotdata)
        


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                    type='each_gauge')
    #plotfigure.show = False
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q'

    # Plot q as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata


def make_plot_figure(plane, k, v, plotdata):

    import os
    from clawpack.visclaw import colormaps

    title = '%s slice (%s = %s)' % (plane,k[0],v)
    outdir = os.path.join(plotdata.outdir, 'slice_%s%s' % (plane,k[1]))

    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name=title)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.title = title
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.outdir = outdir
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.1
    plotitem.pcolor_cmax = 1.
    plotitem.add_colorbar = True

    plotitem.amr_celledges_show = [0]  
    plotitem.amr_patchedges_show = [1]

