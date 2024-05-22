from matplotlib.colors import LinearSegmentedColormap, Normalize

SPEI_COLORS = ['#8B1A1A', '#DE2929', '#F3641D', '#FDC404',
          '#9AFA94', '#03F2FD', '#12ADF3', '#1771DE', '#00008B']
SPEI_CMAP_N = 1000
SPEI_CMAP = LinearSegmentedColormap.from_list('spei_cmap', SPEI_COLORS, N=SPEI_CMAP_N)
SPEI_V_MIN = -2.0
SPEI_V_MAX = 2.0
SPEI_NORM = Normalize(vmin=SPEI_V_MIN, vmax=SPEI_V_MAX)
