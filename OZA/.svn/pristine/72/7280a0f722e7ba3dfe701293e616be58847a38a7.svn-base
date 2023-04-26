import numpy as np

# Voor gumbelplotjes
def add_return_period(tretlist,ax):
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    lbl = ["%d"%int(tr) for tr in tretlist]
    pos = [-np.log(-np.log(1.-1./tr)) for tr in tretlist]
    ax2.set_xticks(pos)
    ax2.set_xticklabels(lbl, rotation='vertical')
    return ax2

def add_gum_axis(tretlist,ax):
    ax2 = add_return_period(tretlist,ax)
    grid2 = ax2.grid()
    grid1 = ax.grid()
    ylbl1 = ax.set_ylabel('Afvoer [m3 s-1]')
    xlbl1 = ax.set_xlabel('Standardised Gumbel Variate [-]')
    xlbl2 = ax2.set_xlabel('Herhalingstijd [jaar]')
    return ax2
