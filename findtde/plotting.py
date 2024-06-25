#!/usr/bin/env python3
"""Python module of plotting functions for findTDE calculations."""
import os
from pathlib import Path
import random as rand
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
# import seaborn as sns

from findtde.utils import idw, idw_heatmap

base_path = Path.cwd()
bin_path, inp_path, perfect_path = base_path / 'bin', base_path / 'inp', base_path / 'perfect'

# interpolation and plotting functions
txt_positions_reg = ['top center', 'bottom center']
txt_positions_extra = ['top center', 'bottom center', 'top right', 'top left', 'bottom right', 'bottom left']


def improve_text_position(x, txt_positions=txt_positions_reg):
    """
    Function to improve the text position for Plotly scatter plots.
    More efficient if the x values are sorted.
    """
    return [txt_positions[(i % len(txt_positions)-(rand.randint(1,2)))] for i in range(len(x))]


# Plotly hover text formatting for heatmaps
hover_se = '<br><b>Polar</b>: %{x:.2f}'+\
            '<br><b>Azimuthal</b>: %{y:.2f}'+\
            '<br><b>Ef</b>: %{z:.2f} eV'

hover_tde = '<br><b>Polar</b>: %{x:.2f}'+\
            '<br><b>Azimuthal</b>: %{y:.2f}'+\
            '<br><b>TDE</b>: %{color:.2f} eV'


# Change Plotly default template to simple white and modify for 
pl_paper_theme = pio.templates['simple_white']
pl_paper_theme.layout.xaxis.ticks = 'inside'
pl_paper_theme.layout.yaxis.ticks = 'inside'
pl_paper_theme.layout.xaxis.mirror = 'ticks'  # True | "ticks" | False | "all" | "allticks"
pl_paper_theme.layout.yaxis.mirror = 'ticks'  # True | "ticks" | False | "all" | "allticks"
pl_paper_theme.layout.font.size = 32
# pl_paper_theme.layout.xaxis.title.standoff = 20
pl_paper_theme.layout.xaxis.title.font.size = 40
# pl_paper_theme.layout.yaxis.title.standoff = 20
pl_paper_theme.layout.yaxis.title.font.size = 40
#pl_paper_theme.layout.coloraxis.colorbar.title.standoff = 20
pio.templates.default = pl_paper_theme


def generate_tde_line_plot(tde_data_df, im_write=False, im_name='tde_lineplot.png'):
    """
    Given a DataFrame of TDE data, produces line plot for each calculated direction with KE on x-axis
    and difference in final energy on y-axis.
    """
    fig = px.line(tde_data_df, #log_y=True,
                  markers=True, color_discrete_sequence=px.colors.qualitative.Dark24) # Alphabet, Light24, & Dark24 have most

    fig.update_traces(connectgaps=True, marker_size=9, line = dict(width=4, dash='dash'))
    fig.update_layout(autosize=False, width=1300, height=600,
                      xaxis_title=r'$KE_{i} \text{ [eV]}$',
                      yaxis_title=r'$\Delta E \text{ [eV]}$',
                      legend_title='Pseudo',
                      font_size=20
                     )
    fig.show()
    
    if im_write == True:
        fig.write_image(os.path.join(base_path, im_name))
    elif im_write == False:
        pass
    
    return


def generate_tde_scatter_plot(tde_sph_arr, tde_directions, txt_show=False, im_write=False, im_name='tde_scatter_plot.png'):
    """
    Given an Nx3 array of TDE data, produces line plot for each calculated direction with polar angle on x-axis,
    azimuthal angle on y-axis, TDE value as color, and direction pseudos as text.
    """
    fig = go.Figure(data=go.Scatter(x=tde_sph_arr[:, 0],
                                    y=tde_sph_arr[:, 1],
                                    text=list(tde_directions) if txt_show == True else None,
                                    mode='markers+text',
                                    marker=dict(color=tde_sph_arr[:, 2],
                                                colorscale='thermal_r',
                                                size=40, colorbar=dict(thickness=20) #, title=r'$E_{d} \; [eV]$', title_side='right')
                                               )
                                   ))
    
    if txt_show == True:
        fig.update_traces(textposition=improve_text_position(list(tde_directions), txt_positions=txt_positions_reg))
    elif txt_show == False:
        pass

    fig.add_vline(x=90, line=dict(color='black', width=3, dash='dash'))
    fig.add_vline(x=210, line=dict(color='black', width=3, dash='dash'))

    fig.update_layout(# title=('Interpolated SEs from Sph. Pert. of Ga #34 with Latt. Dir. TDEs'),
                      # xaxis_title=r'$\phi \; [^{\circ}]$',
                      # yaxis_title=r'$\theta \; [^{\circ}]$',
                      xaxis=dict(tickmode='linear', dtick=30),
                      yaxis=dict(tickmode='linear', dtick=30),
                      xaxis_range=[80, 160],
                      # xaxis_range=[80, 220],
                      yaxis_range=[0, 180],
                      yaxis_autorange='reversed',
                      # font_size=32,  # 16
                      autosize=False, width = 1200, height = 900  # width = 1600, height = 900
                     )

    fig.show()
    
    if im_write == True:
        fig.write_image(os.path.join(base_path, im_name), scale=2)
    elif im_write == False:
        pass
    
    return


def generate_tde_heatmap_plot(polars, azimuthals, energies, im_write=False, im_name='tde_heatmap.png'):
    """
    Given an Nx3 array of TDE data, produces line plot for each calculated direction with polar angle on x-axis,
    azimuthal angle on y-axis, TDE value as color, and direction pseudos as text.
    """
    fig = go.Figure(data=go.Heatmap(x=polars,
                                    y=azimuthals,
                                    z=energies,
                                    hovertemplate=hover_tde,
                                    colorbar_title=r'$E_{d} \; [eV]$',
                                    colorbar_title_side='right',
                                    colorscale='electric_r'
                                   ))

    fig.add_vline(x=90, line=dict(color='black', width=3, dash='dash'))
    fig.add_vline(x=210, line=dict(color='black', width=3, dash='dash'))

    fig.update_layout(# title=('Interpolated Ga #34 with Latt. Dir. TDEs'),
                      # xaxis_title=r'$\phi \; [^{\circ}]$',
                      # yaxis_title=r'$\theta \; [^{\circ}]$',
                      xaxis=dict(tickmode='linear', dtick=30),
                      yaxis=dict(tickmode='linear', dtick=30),
                      xaxis_range=[90, 210],
                      # xaxis_range=[90, 150],
                      yaxis_range=[0, 180],
                      yaxis_autorange='reversed',
                      font_size=16,  # 24
                      autosize=False, width = 1400, height = 900  # width = 1600, height = 900
                     )
    
    fig.show()
    
    if im_write == True:
        fig.write_image(os.path.join(base_path, im_name))
    elif im_write == False:
        pass
    
    return
