import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Times New Roman'
plt.rcParams['font.size'] = '18'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['mathtext.rm'] = 'sans'
plt.rcParams['mathtext.it'] = 'sans:italic'
plt.rcParams['mathtext.default'] = 'it'
ax = fig.add_subplot(111)
ax.plot([1,2,3], '--rs', markersize=11, markerfacecolor='none', linewidth=2, label='EPORSS')
ax.plot([1,2,3], '--bx', markersize=11, linewidth=2, label='Greedy')
ax.plot([1,2,3], '--g^', markersize=11, markerfacecolor='none', linewidth=2,
         label='Modified greedy')
ax.plot([1,2,3], '--co', markersize=11, markerfacecolor='none', linewidth=2, label='SATURATE')

def export_legend(ax, filename="legend.pdf"):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    ax2.axis('off')
    legend = ax2.legend(*ax.get_legend_handles_labels(), frameon=True, loc='lower center', ncol=10, framealpha=1)
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent()
    # print(bbox.extents)
    # bbox = bbox.from_extents(*(bbox.extents + np.array([-100,-100,100,100])))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(ax)
plt.show()