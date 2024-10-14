# ----------------------------------------------------------------
# Plot function
# 
# ----------------------------------------------------------------
# Author: ravi_ram
# ----------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

# combined plot for weighted mean evaluation graphs and porkchop plot
def plot(df, title, start_date, arrival_date, xlist, ylist,
                  xy_contour_data_1, xy_contour_data_2, clevels,
                  xy_tof_data, tlevels):

    # subplots [3X3] grid.
    # 1-col width X 1-row height  [3 optimal graphs]
    # 2-col width X 3-rows height [1 porkchop plot]
    fig = plt.figure(figsize=(10, 6), constrained_layout=True)
    fig.suptitle('Optimal Travel', fontsize=10)
    gs = fig.add_gridspec(3,3)
    
    # setup for grids, ticks, legends
    def setup_optimal(ax):   
        # set grid lines
        ax.grid(True, which='major', color='k', linestyle='-', lw=0.2, alpha=0.5)
        ax.grid(True, which='minor', color='k', linestyle='--', lw=0.2, alpha=0.5)
        # set ticks
        ax.minorticks_on()
        ax.xaxis.set_tick_params(labelsize=7, rotation=90)
        ax.yaxis.set_tick_params(labelsize=7)    
        # legends
        #ax.legend(loc='best', prop={'size': 7})    
        ax.legend(fontsize=7)
        
        # vertical line       
        ax.axvline(pd.to_datetime(arrival_date), color='m',
                   linestyle='--', lw=0.9)
        ax.text(pd.to_datetime(arrival_date), 0.75, arrival_date, color='r',
                ha='right', va='top', rotation=90, size=7,
                transform=ax.get_xaxis_transform())       
        return ax
    # end function 

    # optimal - c3
    ax1 = fig.add_subplot(gs[0, 0])
    df.plot(x='Arrival Date', y=['C3 Departure', 'C3 Arrival'],
            rot=90, grid=True, x_compat=True,
            lw=0.5, fontsize=8, marker='.', ms=3, ax=ax1)    
    ax1.set_title('C3', fontsize=7)
    ax1.set_ylabel(r'C3 $(km^{2}/s^{2})$', fontsize=7)
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    # optimal - delta-v
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    df.plot(x='Arrival Date', y=['Departure delV', 'Arrival delV', 'Total delV'],
            rot=90, grid=True, x_compat=True,
            lw=0.5, fontsize=8, marker='.', ms=3, ax=ax2)    
    ax2.set_title('dV', fontsize=7)
    ax2.set_ylabel(r'dV (km/s)', fontsize=7)
    plt.setp(ax2.get_xticklabels(), visible=False)

    # optimal - phase angle
    ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)
    df.plot(x='Arrival Date', y=['Departure Phase', 'Arrival Phase'],
            rot=90, grid=True, x_compat=True,
            lw=0.5, fontsize=8, marker='.', ms=3, ax=ax3)    
    ax3.set_title('Phase Angle', fontsize=7)    
    ax3.set_ylabel(r'γ $^\circ$', fontsize=7)
        
    # set ticks and grid lines for optimal graphs
    for ax in [ax1, ax2, ax3]:
        setup_optimal(ax)       


    # 2nd [2-cols width X 3-rows height] porkchop
    ax4 = fig.add_subplot(gs[:, 1:3]) 
    ax4.set_title(title, fontsize=7)

    # pork
    def plot_pork(ax):
        # countour text format ( include 'days')
        def tp_fmt(x):
            s = f"{x:.1f}"
            return rf"{s} days"
        # end function
        
        # find the minimum value with the corresponding dep, arr dates
        # draw ∆V_dep contours cmap=rainbow, hsv, 
        cp1 = ax.contour(xlist, ylist, xy_contour_data_1, clevels, cmap="spring")
        ax.clabel(cp1, inline=True, fontsize=7)      

        # find the minimum value with the corresponding dep, arr dates
        # draw ∆V_arr contours
        cp2 = ax.contour(xlist, ylist, xy_contour_data_2, clevels, cmap="winter")
        ax.clabel(cp2, inline=True, fontsize=7)
        # get legend handles
        h1,_ = cp1.legend_elements()
        h2,_ = cp2.legend_elements()
        ax.legend([h1[0], h2[0]], ['Departure', 'Arrival'])  
        
        # draw time-of-flight contours
        tp = ax.contour(xlist, ylist, xy_tof_data, tlevels, colors='k',linestyles=':')
        ax.clabel(tp, inline=True, fmt=tp_fmt, fontsize=7)
        # optimal arrival horizontal line
        ad = pd.to_datetime(arrival_date, dayfirst=True, exact=True, format='%d-%m-%Y')
        ax.text(0.75, ad, arrival_date, color='k',
                ha='right', va='bottom', rotation=0, size=7,
                transform=ax.get_yaxis_transform())        
        ax.axhline(ad, color='r', linestyle='-.', lw=1.0)
        # optimal arrival vertical line
        sd = pd.to_datetime(start_date, dayfirst=True, exact=True, format='%d-%m-%Y')
        ax.axvline(sd, color='r', linestyle='-.', lw=1.0)
        ax.text(sd, 0.25, start_date, color='k',
                ha='right', va='bottom', rotation=90, size=7,
                transform=ax.get_xaxis_transform())
        
        # intersection
        ax.scatter(sd, ad, s=20, c='k', marker='o', alpha=0.7, linewidths=1.0)
        
        # major grid
        x_tick_spacing, y_tick_spacing = 5, 3
        ax.xaxis.set_major_locator(ticker.MultipleLocator(x_tick_spacing))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(y_tick_spacing))
        # fontsize of major, minor ticks label
        ax.xaxis.set_tick_params(labelsize=7, rotation=90)
        ax.yaxis.set_tick_params(labelsize=7)
        ax.set_xlabel("Dep Date (yyyy-mm-dd)", fontsize=8)        
        ax.set_ylabel("Arr Date (yyyy-mm-dd)", fontsize=8)

        # Customize the major grid
        ax.grid(which='major', linestyle='dashdot', linewidth='0.5', color='gray')
        # Customize the minor grid
        ax.grid(which='minor', linestyle='dotted', linewidth='0.5', color='gray')       
        # Turn on the minor ticks (minor grid)
        ax.minorticks_on()        
        # Turn off the display of all ticks.
        ax.tick_params(which='both',
                        top='off',
                        left='off',
                        right='off',
                        bottom='off')
        
        # set aspect ratio and layout
        ax.set_aspect(1)         
        return
    # end function    

    # set grids and ticks
    plot_pork(ax4)
 
    # show
    plt.savefig("plot.png", format="png")
    plt.show()
    # return
    return

