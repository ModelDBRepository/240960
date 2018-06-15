__author__ = 'milsteina'
import matplotlib as mpl
import matplotlib.lines as mlines
import scipy.stats as stats
import matplotlib.gridspec as gridspec
from matplotlib import cm
from function_lib import *

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.size'] = 12.
mpl.rcParams['font.sans-serif'] = 'Myriad Pro'
mpl.rcParams['text.usetex'] = False


def plot_synaptic_param_distribution(cell, syn_type, param_name, scale_factor=1., param_label=None,
                                 ylabel='Peak conductance', yunits='uS', svg_title=None):
    """
    Takes a cell as input rather than a file. No simulation is required, this method just takes a fully specified cell
    and plots the relationship between distance and the specified synaptic parameter for all spines. Used while
    debugging specification of synaptic parameters.
    :param cell: :class:'HocCell'
    :param syn_type: str
    :param param_name: str
    :param scale_factor: float
    :param param_label: str
    :param ylabel: str
    :param yunits: str
    :param svg_title: str
    """
    colors = ['k', 'r', 'c', 'y', 'm', 'g', 'b']
    dend_types = ['basal', 'trunk', 'apical', 'tuft']

    if svg_title is not None:
        remember_font_size = mpl.rcParams['font.size']
        mpl.rcParams['font.size'] = 20
    fig, axes = plt.subplots(1)
    maxval, minval = None, None
    for i, sec_type in enumerate(dend_types):
        syn_list = []
        distances = []
        param_vals = []

        for branch in cell.get_nodes_of_subtype(sec_type):
            for spine in branch.spines:
                syn_list.extend(spine.synapses)
            syn_list.extend(branch.synapses)
        for syn in [syn for syn in syn_list if syn_type in syn._syn]:
            if syn.node.type == 'spine_head':
                this_distance = cell.get_distance_to_node(cell.tree.root, syn.node.parent.parent, syn.loc)
            else:
                this_distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
            distances.append(this_distance)
            if sec_type == 'basal':
                    distances[-1] *= -1
            param_vals.append(getattr(syn.target(syn_type), param_name) * scale_factor)
        if param_vals:
            axes.scatter(distances, param_vals, color=colors[i], label=sec_type)
            if maxval is None:
                maxval = max(param_vals)
            else:
                maxval = max(maxval, max(param_vals))
            if minval is None:
                minval = min(param_vals)
            else:
                minval = min(minval, min(param_vals))

    axes.set_ylabel(ylabel + ' (' + yunits + ')')
    if (maxval is not None) and (minval is not None):
        buffer = 0.1 * (maxval - minval)
        axes.set_ylim(minval - buffer, maxval + buffer)
    axes.set_xlabel('Distance to soma (um)')
    axes.set_xlim(-200., 525.)
    axes.set_xticks([-150., 0., 150., 300., 450.])
    plt.legend(loc='best', scatterpoints=1, frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    if param_label is not None:
        plt.title(param_label, fontsize=mpl.rcParams['font.size'])
    clean_axes(axes)
    axes.tick_params(direction='out')
    if not svg_title is None:
        if param_label is not None:
            svg_title = svg_title+' - '+param_label+'.svg'
        else:
            svg_title = svg_title+' - '+syn_type+'_'+param_name+' distribution.svg'
        fig.set_size_inches(5.27, 4.37)
        fig.savefig(data_dir+svg_title, format='svg', transparent=True)
    plt.show()
    plt.close()
    if svg_title is not None:
        mpl.rcParams['font.size'] = remember_font_size


def plot_mech_param_distribution(cell, mech_name, param_name, scale_factor=10000., param_label=None,
                                 ylabel='Conductance density', yunits='pS/um2', svg_title=None):
    """
    Takes a cell as input rather than a file. No simulation is required, this method just takes a fully specified cell
    and plots the relationship between distance and the specified mechanism parameter for all dendritic segments. Used
    while debugging specification of mechanism parameters.
    :param cell: :class:'HocCell'
    :param mech_name: str
    :param param_name: str
    :param scale_factor: float
    :param ylabel: str
    :param yunits: str
    :param svg_title: str
    """
    colors = ['k', 'r', 'c', 'y', 'm', 'g', 'b']
    dend_types = ['basal', 'trunk', 'apical', 'tuft']

    if svg_title is not None:
        remember_font_size = mpl.rcParams['font.size']
        mpl.rcParams['font.size'] = 20
    fig, axes = plt.subplots(1)
    maxval, minval = None, None
    for i, sec_type in enumerate(dend_types):
        distances = []
        param_vals = []
        for branch in cell.get_nodes_of_subtype(sec_type):
            for seg in [seg for seg in branch.sec if hasattr(seg, mech_name)]:
                distances.append(cell.get_distance_to_node(cell.tree.root, branch, seg.x))
                if sec_type == 'basal':
                    distances[-1] *= -1
                param_vals.append(getattr(getattr(seg, mech_name), param_name) * scale_factor)
        if param_vals:
            axes.scatter(distances, param_vals, color=colors[i], label=sec_type)
            if maxval is None:
                maxval = max(param_vals)
            else:
                maxval = max(maxval, max(param_vals))
            if minval is None:
                minval = min(param_vals)
            else:
                minval = min(minval, min(param_vals))
    axes.set_xlabel('Distance to soma (um)')
    axes.set_xlim(-200., 525.)
    axes.set_xticks([-150., 0., 150., 300., 450.])
    axes.set_ylabel(ylabel+' ('+yunits+')')
    if (maxval is not None) and (minval is not None):
        buffer = 0.1 * (maxval - minval)
        axes.set_ylim(minval-buffer, maxval+buffer)
    if param_label is not None:
        plt.title(param_label, fontsize=mpl.rcParams['font.size'])
    plt.legend(loc='best', scatterpoints=1, frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    clean_axes(axes)
    axes.tick_params(direction='out')
    if not svg_title is None:
        if param_label is not None:
            svg_title = svg_title+' - '+param_label+'.svg'
        else:
            svg_title = svg_title+' - '+mech_name+'_'+param_name+' distribution.svg'
        fig.set_size_inches(5.27, 4.37)
        fig.savefig(data_dir + svg_title, format='svg', transparent=True)
    plt.show()
    plt.close()
    if svg_title is not None:
        mpl.rcParams['font.size'] = remember_font_size


def plot_sum_mech_param_distribution(cell, mech_param_list, scale_factor=10000., param_label=None,
                                 ylabel='Conductance density', yunits='pS/um2', svg_title=None):
    """
    Takes a cell as input rather than a file. No simulation is required, this method just takes a fully specified cell
    and plots the relationship between distance and the specified mechanism parameter for all dendritic segments. Used
    while debugging specification of mechanism parameters.
    :param cell: :class:'HocCell'
    :param mech_param_list: list of tuple of str
    :param scale_factor: float
    :param ylabel: str
    :param yunits: str
    :param svg_title: str
    """
    colors = ['k', 'r', 'c', 'y', 'm', 'g', 'b']
    dend_types = ['basal', 'trunk', 'apical', 'tuft']

    if svg_title is not None:
        remember_font_size = mpl.rcParams['font.size']
        mpl.rcParams['font.size'] = 20
    fig, axes = plt.subplots(1)
    maxval, minval = None, None
    for i, sec_type in enumerate(dend_types):
        distances = []
        param_vals = []
        for branch in cell.get_nodes_of_subtype(sec_type):
            for seg in branch.sec:
                this_param_val = 0.
                this_distance = None
                for mech_name, param_name in mech_param_list:
                    if hasattr(seg, mech_name):
                        if this_distance is None:
                            this_distance = cell.get_distance_to_node(cell.tree.root, branch, seg.x)
                            if sec_type == 'basal':
                                this_distance *= -1
                        this_param_val += getattr(getattr(seg, mech_name), param_name) * scale_factor
                if this_distance is not None:
                    distances.append(this_distance)
                    param_vals.append(this_param_val)
        if param_vals:
            axes.scatter(distances, param_vals, color=colors[i], label=sec_type)
            if maxval is None:
                maxval = max(param_vals)
            else:
                maxval = max(maxval, max(param_vals))
            if minval is None:
                minval = min(param_vals)
            else:
                minval = min(minval, min(param_vals))
    axes.set_xlabel('Distance to soma (um)')
    axes.set_xlim(-200., 525.)
    axes.set_xticks([-150., 0., 150., 300., 450.])
    axes.set_ylabel(ylabel+' ('+yunits+')')
    buffer = 0.1 * (maxval - minval)
    axes.set_ylim(minval-buffer, maxval+buffer)
    if param_label is not None:
        plt.title(param_label, fontsize=mpl.rcParams['font.size'])
    plt.legend(loc='best', scatterpoints=1, frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    clean_axes(axes)
    axes.tick_params(direction='out')
    if not svg_title is None:
        if param_label is not None:
            svg_title = svg_title+' - '+param_label+'.svg'
        else:
            mech_name, param_name = mech_param_list[0]
            svg_title = svg_title+' - '+mech_name+'_'+param_name+' distribution.svg'
        fig.set_size_inches(5.27, 4.37)
        fig.savefig(data_dir + svg_title, format='svg', transparent=True)
    plt.show()
    plt.close()
    if svg_title is not None:
        mpl.rcParams['font.size'] = remember_font_size
