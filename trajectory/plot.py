"""Routines for plotting segments"""

import pandas as pd
import matplotlib.pyplot as plt

def sel_axis(df, axis):
    t = df[df.axis == axis].copy()

    t = pd.concat([t.iloc[0:1], t]).reset_index()
    t.at[0, 't'] = 0
    t.at[0, 'v_f'] = 0
    t = t.set_index('t')
    t = t[['v_f', 'ss', 'del_t']]
    return t


def plot_axis(df, axis, ax=None):

    t = sel_axis(df, axis)

    ax = t[['v_f']].plot(ax=ax)

    return ax


def plot_segment_list(df, ax=None):

    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(20,3))

    for axn in df.axis.unique():
        ax = plot_axis(df, axn, ax=ax)

    # Map sub-segment names ( a,c,d) to matplotlib colors.
    # The color here is actually the one that ends the segment, but it is
    # interpreted as the one that starts it.
    cm = {'a': ('b',3),
          'c': ('r',3),
          'd': ('g',3)}

    for idx, row in sel_axis(df, 0).iterrows():
        #rectangle = plt.Rectangle((idx,0), row.del_t, row.v_f, fc=cm[row.ss], alpha=0.05)
        #plt.gca().add_patch(rectangle)

        ax.axvline(x=idx, color='k', alpha=.5, lw=1)