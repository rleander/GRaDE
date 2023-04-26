import numpy as np

# Voor gumbelplotjes
def add_return_period(tretlist,ax,ymax,yrange):
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    lbl = ["%d"%int(tr) for tr in tretlist]
    pos = [-np.log(-np.log(1.-1./tr)) for tr in tretlist]
    for tr in tretlist:
        xpos = -np.log(-np.log(1.-1./tr))
        ax2.plot([xpos,xpos],[0,ymax],color='green',linestyle='--', alpha=0.6)
    ax2.set_xticks(pos)
    ax2.set_yticks(np.arange(0,ymax,yrange))
#   ax2.yaxis.set_minor_locator(MultipleLocator(1000))
    ax2.set_xticklabels(lbl, rotation='vertical')
    return ax2

def add_gum_axis(tretlist,ax,ymax,yrange):
    ax2 = add_return_period(tretlist,ax,ymax,yrange)
    grid2 = ax2.grid()
    grid1 = ax.grid(which='major')
    grid1 = ax.grid(which='minor', color='grey', linewidth=0.1)
    ylbl1 = ax.set_ylabel('Discharge [m3 s-1]')
    xlbl1 = ax.set_xlabel('Standardised Gumbel Variate [-]')
    xlbl2 = ax2.set_xlabel('Return period [year]')

def add_return_period_vert(tretlist,ax,xmax,xrange):
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    lbl = ["%d"%int(tr) for tr in tretlist]
    pos = [-np.log(-np.log(1.-1./tr)) for tr in tretlist]
    for tr in tretlist:
        ypos = -np.log(-np.log(1.-1./tr))
        ax2.plot([0,xmax],[ypos,ypos],color='green',linestyle='--', alpha=0.6)
    ax2.set_yticks(pos)
    ax2.set_xticks(np.arange(0,xmax,xrange))
    ax2.set_yticklabels(lbl, rotation='horizontal')
    return ax2

def add_gum_axis_vert(tretlist,ax,xmax,xrange):
    ax2 = add_return_period_vert(tretlist,ax,xmax,xrange)
    grid2 = ax2.grid()
    grid1 = ax.grid(which='major')
    grid1 = ax.grid(which='minor', color='grey', linewidth=0.1)
    ylbl1 = ax.set_xlabel('Discharge [m3 s-1]')
    xlbl1 = ax.set_ylabel('Standardised Gumbel Variate [-]')
    xlbl2 = ax2.set_ylabel('Return period [year]')



def weissman(xi,**kwargs):
    if 'tretlist' in kwargs:
        fx = 1./kwargs['tretlist']
    elif 'flist' in kwargs:
        fx = 1.-kwargs['flist']
    elif 'fxlist' in kwargs:
        fx = kwargs['fxlist']
    else:
        return
    nx = len(xi)
    if 'k' in kwargs: 
        k = kwargs['k']
    elif 'fxth' in kwargs: 
        k = int(0.3 + (nx + 0.4)*(1.-kwargs['fxth']))
    sigma_hat = np.mean(xi[-k:])-xi[-k]
    retvals = xi[-k] + sigma_hat * np.log(1./fx*(k/nx))
    return retvals 

          
