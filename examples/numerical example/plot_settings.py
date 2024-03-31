#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "sans serif"
plt.rcParams['font.serif'] = "cm"

matplotlib.rcParams['xtick.major.pad']='10'
matplotlib.rcParams['ytick.major.pad']='7'
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({'legend.fontsize': 15})


framealpha = 0.95


# plotting styles

# horizontal line 
horizontal = {
            'value':1e-2,
            'label':r'$E \cdot L = {0:3.2f}$',
            'color':'black',
            'lw':1,
            }

def add_horizontal_line(ax,
                        value=None,
                        label=None):
    #
    if value is None:
        value = horizontal['value']
    if label is None:
        label = horizontal['label'].format(value)
    #
    ax.axhline(value,
               color=horizontal['color'],
              label=label,
              lw=horizontal['lw'],
              )


# vertical line that indicates the breakdown
breakdown = {
            'label':r'$\Delta t_b/T \approx {0:3.2f}$',
            'color':'gray',
            'lw':1,
            }

def add_vertical_breakdown_line(ax,
                                t_breakdown):
    #
    ax.axvline(t_breakdown,
        color=breakdown['color'],
        label=breakdown['label'].format(t_breakdown),
        lw=breakdown['lw'],
        )


# analytical solution
analytical = {
        'label':'Exact',
        'color':'black',
        'lw':3,
        }

# perturbative solutions
# 0 -> Gaussian
# 1 -> K = 2
# 2 -> K = 8
'''
perturbative_dashes = {
                0:[4,4], 
                1:[1.5,1.5],
                2:[1,1,1,3,4,3],
                     }

perturbative_color = {
        0:'crimson',
        1:'dodgerblue',
        2:'limegreen',
            }

perturbative_label = {
        0:r'GP',
        1:r'NPP, $K = 2$',
        2:r'NPP, $K=8$',
            }

perturbative_lw={
            0:3,
            1:3,
            2:3,
            }
''';

perturbative = {'dashes':[ [4,4], 
                            [1.5,1.5],
                            [1,1,1,3,4,3],
                        ],
                'color': [  'crimson',
                            'dodgerblue',
                            'limegreen',
                        ],
                'label': [ r'GP',
                        r'NPP, $K = 2$',
                        r'NPP, $K=8$',
                        ],
                'lw': [3,3,3],
                     }


x0_vertical = {'color':'gray',
            'ls':'--'}


#
dt_plot = {'times':[5e-2, # time for middle column
                    2e-1, # time for right column
                    ],
            'color':['gray','gray'],
            'label': [r'$\Delta t/T =0.05$',
                    r'$\Delta t/T =0.2$'],
            'dashes': [[5,5],
                        [2,2],
                    ],
            }
