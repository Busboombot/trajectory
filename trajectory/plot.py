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

    def _f():
        t = 0
        for i, (idx, r) in enumerate(df_.iterrows()):
            if i == 0:
                yield {'t': t, 'v': r.v_i, 'ss': r.ss}
            elif last_row.v_f != r.v_i:
                a = f"{last_row.seg}/{last_row.axis}{last_row.ss}"
                b = f"{r.seg}/{r.axis}{r.ss}"
                warn(f"Discontinuty {a}@{last_row.v_f} -> {b}@{r.v_i}")
                yield {'t': t, 'v': r.v_i, 'ss': r.ss}

            t = t + r.del_t
            yield {'t': t, 'v': r.v_f, 'ss': r.ss}
            last_row = r

    t = pd.DataFrame(list(_f())).set_index('t')

    ax = t[['v']].plot(ax=ax)

    for t, row in t.iterrows():
        # rectangle = plt.Rectangle((idx,0), row.del_t, row.v_f, fc=cm[row.ss], alpha=0.05)
        # plt.gca().add_patch(rectangle)
        ax.axvline(x=t, color='k', alpha=.5, lw=.5, linestyle='dotted')

    return ax

def plot_segment_list(df, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(20, 3))

    for axn in df.axis.unique():
        ax = plot_axis(df, axn, ax=ax)


def plot_params_df(*args):
    from trajectory.planner import Block

    import numpy as np

    cols = ['t', 'seg', 'axis', 'x', 'v_i', 'v_f', 'ss', 'del_t', 'v0m', 'v1m', 'calc_x', 'err']

    def p_to_s(p, seg, axis):
        a = pd.Series(index=cols, dtype=np.float64)
        c = pd.Series(index=cols, dtype=np.float64)
        d = pd.Series(index=cols, dtype=np.float64)

        a.seg = c.seg = d.seg = seg
        a.axis = c.axis = d.axis = axis

        a.x = p.x_a
        c.x = p.x_c
        d.x = p.x_d

        a.ss = 'a'
        c.ss = 'c'
        d.ss = 'd'

        a.del_t = p.t_a
        c.del_t = p.t_c
        d.del_t = p.t_d

        a.v_i = p.v_0
        a.v_f = p.v_c
        c.v_i = p.v_c
        c.v_f = p.v_c
        d.v_i = p.v_c
        d.v_f = p.v_1

        yield a
        yield c
        yield d

    seg_0 = []
    segments = []

    for i,a in enumerate(args):

        if isinstance(a, Block):
            # Individual Params
            seg_0.append(a)
        else: # Assume tuple of params -- a segment
            segments.append(a)

    if len(seg_0):
        segments = [seg_0] + segments


    rows = []
    for i, seg in enumerate(segments):
        for j, p in enumerate(seg):
            rows.extend(p_to_s(p, i, j))

    df = pd.DataFrame(rows)
    df['seg'] = df.seg.astype(int)
    df['axis'] = df.axis.astype(int)
    df['t'] = df.groupby('axis').del_t.cumsum()
    df = df.fillna(0)

    return df


def plot_params(*args, ax=None):

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(20, 3))

    df = plot_params_df(*args)
    plot_segment_list(df, ax=ax)

    return ax
