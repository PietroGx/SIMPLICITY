import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'long_paper_figures'))
from long_shedders_plots import PIE_COLORS


def _color_for(k):
    return PIE_COLORS.get(k.traits.get('cluster'), 'grey')


def plot_fig4_tree(ax, bt_tree, title):
    """One scenario's phylogenetic tree, branches/tips colored by the same
    long/standard/mixed/founder labels as Figure 3's pies/PCA."""
    if bt_tree is None:
        ax.axis('off')
        ax.text(0.5, 0.5, f"{title}\n(Data Missing)", ha='center', va='center',
                transform=ax.transAxes, fontsize=8)
        return

    bt_tree.drawTree()
    bt_tree.plotTree(ax, colour=_color_for, width=1.2)
    bt_tree.plotPoints(ax, colour=_color_for, size=20, zorder=3)

    ax.set_title(title, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    ax.set_xlabel("Time", fontsize=7)
