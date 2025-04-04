
from pypion.example_plots_yt import yt_example_plots
ytp = yt_example_plots("/mnt/massive-stars/data/ostar-xray/mhd_d2n0128l3_s7/", "Ostar_mhd_d2n0128l3_s7", "img", "Ostar_mhd_d2n0128l3_s7", quantities=["density"])
for i in range(0,ytp.get_nframes(),100):
  if i==0: continue
  ytp.plot_density_2d(i, 1e-27, 1e-22, 5.0e3, 2.0e8)
quit()

#ytp = yt_example_plots("/mnt/massive-stars/data/ostar-xray/mhd_d3n0128l3_s7/", "Ostar_mhd_d3n0128l3_s7", "img", "Ostar_mhd_d3n0128l3_s7", quantities=["density"])
#for i in range(0,ytp.get_nframes(),100):
#  if i==0: continue
#  ytp.plot_Bfield_both(i, 1e-27, 1e-22)
#quit()

