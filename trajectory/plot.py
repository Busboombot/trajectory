"""Routines for plotting segments"""

import matplotlib.pyplot as plt
import pandas as pd
from warnings import warn


def sel_axis(df, axis):
    t = df[df.axis == axis].copy()

    t = pd.concat([t.iloc[0:1], t]).reset_index()

    t.at[0, 't'] = 0
    t.at[0, 'v_f'] = 0
    t = t.set_index('t')
    t = t[['v_f', 'ss', 'del_t']]

    return t


def plot_axis(df, axis, ax=None):
    df_ = df[df.axis == axis].copy()

    discontinuities = [0]
    def _f():
        t = 0

        for i, (idx, r) in enumerate(df_.iterrows()):
            if i == 0:
                yield {'t': t, 'v': r.v_i}
            elif abs(last_row.v_f-r.v_i)>1: # Discontinuity limit
                # If there is no discontinuity, we don't need to yield this part,
                # because both v_0 and v_1 go through the same point, so the line
                # segment would be zero length
                a = f"{last_row.seg}/{last_row.axis}"
                b = f"{r.seg}/{r.axis}"
                # warn(f"Discontinuty {a}@{last_row.v_f} -> {b}@{r.v_i}")
                discontinuities[0] =discontinuities[0] +1

                yield {'t': t, 'v': r.v_i}

            t = t + r.del_t
            yield {'t': t, 'v': r.v_f}
            last_row = r

    t = pd.DataFrame(list(_f())).set_index('t')

    ax = t[['v']].plot(ax=ax, marker='+')

    for t, row in t.iterrows():
        # rectangle = plt.Rectangle((idx,0), row.del_t, row.v_f, fc=cm[row.ss], alpha=0.05)
        # plt.gca().add_patch(rectangle)
        ax.axvline(x=t, color='k', alpha=.5, lw=.5, linestyle='dotted')

    if discontinuities[0]:
        warn(f"Found {discontinuities[0]} discontinuities in axis {axis}")

    return ax

def plot_trajectory(df, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(18, 3))

    for axn in df.axis.unique():
        ax = plot_axis(df, axn, ax=ax)

    return ax

def plot_params_df(*args):
    from trajectory.gsolver import ACDBlock

    import numpy as np

    cols = ['t', 'seg', 'axis', 'x', 'v_i', 'v_f', 'del_t']

    def p_to_s(p, seg, axis):
        a = pd.Series(index=cols, dtype=np.float64)

        a.seg = seg
        a.axis = axis
        a.x = p.x
        a.del_t = p.t
        a.v_i = p.v_0
        a.v_f = p.v_1

        return a

    seg_0 = []
    segments = []

    for i,a in enumerate(args):

        if isinstance(a, ACDBlock):
            # Individual Params
            seg_0.append(a)
        else: # Assume tuple of params -- a segment
            segments.append(a)

    if len(seg_0):
        segments = [seg_0] + segments


    rows = []
    for i, seg in enumerate(segments):
        for j, p in enumerate(seg):
            rows.append(p_to_s(p, i, j))

    df = pd.DataFrame(rows)
    df['seg'] = df.seg.astype(int)
    df['axis'] = df.axis.astype(int)
    df['t'] = df.groupby('axis').del_t.cumsum()
    df = df.fillna(0)

    return df


def plot_params(*args, ax=None):

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(18, 2))

    df = plot_params_df(*args)
    plot_trajectory(df, ax=ax)

    return ax
