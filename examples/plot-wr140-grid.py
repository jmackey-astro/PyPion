
from pypion.example_plots_yt import yt_example_plots

ytp = yt_example_plots("/mnt/massive-stars/data/wr140/mide-n128/", "wr140_mhd_compton_mide_d3l7n128", "img", "wr140_mhd_compton_mide_d3l7n128")
for i in range(0,ytp.get_nframes(),5):
  if i==0: continue
  #ytp.plot_grid_plot(i)
  #ytp.plot_Bfield_XY(i, 1e-19, 1e-12)
  #ytp.plot_Bfield_Bmag(i, 1e-6, 1.0, 8)
  #ytp.plot_Bfield_both(i, 1e-19, 1e-12)
  ytp.plot_Xray_intensity_XY(i, 1e-4, 1e1, 16)

quit()

