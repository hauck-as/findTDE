#!/usr/bin/env python3
"""Python module of plotting functions for findTDE calculations."""
import os
from pathlib import Path
from typing import Any, Union
from pymatgen.util.typing import PathLike
import random as rand
import numpy as np
from numpy.typing import ArrayLike
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


def improve_text_position(
    x: list[Any],
    txt_positions: list[str] = txt_positions_reg
) -> list[str]:
    """
    Function to improve the text position for Plotly scatter plots. More efficient if the x values are sorted.

    Args
    ---------
        x (list[Any]):
            List of any type of data to be used for annotating on a plot.
        txt_positions (list(str)):
            List of strings corresponding to possible text position options. Defaults to txt_positions_reg,
            defined prior to function definition as ['top center', 'bottom center']. Another option,
            txt_positions_extra, is included: ['top center', 'bottom center', 'top right', 'top left',
            'bottom right', 'bottom left']. Another option may be specified.

    Returns
    ---------
        List of text positions associated for each x value.
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


def generate_tde_line_plot(
    tde_data_df: pd.DataFrame,
    fig_name: PathLike | None = 'tde_lineplot.png'
) -> go.Figure():
    """
    Given a DataFrame of TDE data, produces line plot for each calculated direction with KE on x-axis and difference
    in final energy on y-axis.

    Args
    ---------
        tde_data_df (DataFrame):
            DataFrame of findTDE data created by find_tde_analysis. Column names are direction pseudonyms and indices
            are kinetic energy values. Values are the energy differences from displacement calculations referenced to
            a perfect supercell single-point energy calculation.
        fig_name (PathLike or None):
            Path/filename+extension to be used to save figure, or None to not save figure. Defaults to
            'tde_lineplot.png'.

    Returns
    ---------
        Figure showing which KE values for each direction create nonzero energies (suggesting a defect was generated).
    """
    fig = px.line(
        tde_data_df,
        #log_y=True,
        markers=True,
        color_discrete_sequence=px.colors.qualitative.Dark24  # Alphabet, Light24, & Dark24 have most
    )

    fig.update_traces(
        connectgaps=True,
        marker_size=9,
        line = dict(
            width=4,
            dash='dash'
        )
    )
    
    fig.update_layout(
        autosize=False,
        width=1300,
        height=600,
        xaxis_title=r'$KE_{i} \text{ [eV]}$',
        yaxis_title=r'$\Delta E \text{ [eV]}$',
        legend_title='Pseudo',
        font_size=20
    )
    
    if fig_name is not None:
        fig.write_image(fig_name, scale=2)
    
    return fig


def generate_tde_scatter_plot(
    tde_sph_arr: ArrayLike,
    tde_directions: ArrayLike,
    txt_show: bool = False,
    symmetry_lines: dict = {'x': None, 'y': None},
    xaxis_title: str = r'$\Large{\phi \; \left[ ^{\circ} \right]}$',
    xaxis_range: list = [0, 360],
    yaxis_title: str = r'$\Large{\theta \; \left[ ^{\circ} \right]}$',
    yaxis_range: list = [0, 180],
    colorbar_title: str = r'$\Large{E_{d} \; \left[ \text{eV} \right]}$',
    fig_name: PathLike | None = 'tde_scatter.png'
) -> go.Figure():
    """
    Given an Nx3 array of TDE data, produces scatter plot for each calculated direction with polar angle on x-axis,
    azimuthal angle on y-axis, TDE value as color, and direction pseudos as text.

    Args
    ---------
        tde_sph_arr (ArrayLike):
            2D array of spherical directions and corresponding TDE values. Each row has a length of 3 with values for
            the polar angle (phi, in the xy-plane around the z-axis, in units of degrees), the azimuthal angle (theta,
            from the z-axis, in units of degrees), and the TDE value (in units of eV).
        tde_directions (ArrayLike):
            Collection of data (typically of type dict_keys, must be able to be converted to a list) corresponding to
            the direction pseudonyms for each TDE calculation.
        txt_show (bool):
            Boolean determining if directions should be displayed on the scatter plot.
        symmetry_lines (dict):
            Dictionary of symmetries to be displayed on the plot. Must have 'x' and 'y' as the keys and either None or
            lists of integers/floats as the values. Defaults to {'x': None, 'y': None} (no symmetries displayed).
        xaxis_range (list):
            List of 2 numbers indicating the minimum and maximum x-values. Oriented such that 0 degrees is along the a
            lattice vector of the crystal. Defaults to [0, 360] (360 degrees in the xy-plane around the z-axis).
        yaxis_range (list):
            List of 2 numbers indicating the minimum and maximum y-values. Oriented such that 0 degrees (top of plot)
            is along the positive c-axis of the crystal and 180 degrees (bottom of plot) is along the negative c-axis
            of the crystal. Defaults to [0, 180].
        fig_name (PathLike or None):
            Path/filename+extension to be used to save figure, or None to not save figure. Defaults to
            'tde_scatter.png'.

    Returns
    ---------
        Figure showing TDE values on a colorscale on scatter plot points indicating the displacement directions.
    """
    fig = go.Figure(
        data=go.Scatter(
            x=tde_sph_arr[:, 0],
            y=tde_sph_arr[:, 1],
            text=list(tde_directions) if txt_show == True else None,
            mode='markers+text',
            marker=dict(
                color=tde_sph_arr[:, 2],
                colorscale='thermal_r',
                size=40,
                colorbar=dict(
                    thickness=20,
                    title=colorbar_title,
                    title_side='right'
                )
            )
        )
    )
    
    if txt_show == True:
        fig.update_traces(textposition=improve_text_position(list(tde_directions), txt_positions=txt_positions_reg))

    if symmetry_lines['x'] is not None:
        for i in symmetry_lines['x']:
            fig.add_vline(x=i, line=dict(color='black', width=3, dash='dash'))
    if symmetry_lines['y'] is not None:
        for i in symmetry_lines['y']:
            fig.add_hline(y=i, line=dict(color='black', width=3, dash='dash'))

    fig.update_layout(
        autosize=False,
        width=1600,
        height=900,
        xaxis=dict(
            title=xaxis_title,
            range=xaxis_range,
            tickmode='linear',
            dtick=30
        ),
        yaxis=dict(
            title=yaxis_title,
            range=yaxis_range,
            autorange='reversed',
            tickmode='linear',
            dtick=30
        )
    )
    
    if fig_name is not None:
        fig.write_image(fig_name, scale=2)
    
    return fig


def generate_tde_heatmap_plot(
    polars: ArrayLike,
    azimuthals: ArrayLike,
    energies: ArrayLike,
    symmetry_lines: dict = {'x': None, 'y': None},
    xaxis_title: str = r'$\Large{\phi \; \left[ ^{\circ} \right]}$',
    xaxis_range: list = [0, 360],
    yaxis_title: str = r'$\Large{\theta \; \left[ ^{\circ} \right]}$',
    yaxis_range: list = [0, 180],
    colorbar_title: str = r'$\Large{E_{d} \; \left[ \text{eV} \right]}$',
    fig_name: PathLike | None = 'tde_heatmap.png'
) -> go.Figure():
    """
    Given an Nx3 array of TDE data, produces an interpolated scatter plot heatmap for each calculated direction with
    polar angle on x-axis, azimuthal angle on y-axis, TDE value as color, and direction pseudos as text.

    Args
    ---------
        polars (ArrayLike):
            Polar angles used for interpolations from idw_heatmap, in units of degrees.
        azimuthals (ArrayLike):
            Azimuthal angles used for interpolations from idw_heatmap, in units of degrees.
        energies (ArrayLike):
            Interpolated energies calculated from idw_heatmap, in units of eV.
        symmetry_lines (dict):
            Dictionary of symmetries to be displayed on the plot. Must have 'x' and 'y' as the keys and either None or
            lists of integers/floats as the values. Defaults to {'x': None, 'y': None} (no symmetries displayed).
        xaxis_range (list):
            List of 2 numbers indicating the minimum and maximum x-values. Oriented such that 0 degrees is along the a
            lattice vector of the crystal. Defaults to [0, 360] (360 degrees in the xy-plane around the z-axis).
        yaxis_range (list):
            List of 2 numbers indicating the minimum and maximum y-values. Oriented such that 0 degrees (top of plot)
            is along the positive c-axis of the crystal and 180 degrees (bottom of plot) is along the negative c-axis
            of the crystal. Defaults to [0, 180].
        fig_name (PathLike or None):
            Path/filename+extension to be used to save figure, or None to not save figure. Defaults to
            'tde_heatmap.png'.

    Returns
    ---------
        Figure showing TDE values on a colorscale on scatter plot points indicating the displacement directions.
    """
    fig = go.Figure(
        data=go.Heatmap(
            x=polars,
            y=azimuthals,
            z=energies,
            hovertemplate=hover_tde,
            colorbar=dict(
                thickness=20,
                title=colorbar_title,
                title_side='right'
            ),
            colorscale='electric_r'
        )
    )

    if symmetry_lines['x'] is not None:
        for i in symmetry_lines['x']:
            fig.add_vline(x=i, line=dict(color='black', width=3, dash='dash'))
    if symmetry_lines['y'] is not None:
        for i in symmetry_lines['y']:
            fig.add_hline(y=i, line=dict(color='black', width=3, dash='dash'))
    
    fig.update_layout(
        autosize=False,
        width=1600,
        height=900,
        xaxis=dict(
            title=xaxis_title,
            range=xaxis_range,
            tickmode='linear',
            dtick=30
        ),
        yaxis=dict(
            title=yaxis_title,
            range=yaxis_range,
            autorange='reversed',
            tickmode='linear',
            dtick=30
        )
    )
    
    if fig_name is not None:
        fig.write_image(fig_name, scale=2)
    
    return fig
