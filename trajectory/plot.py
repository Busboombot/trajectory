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
