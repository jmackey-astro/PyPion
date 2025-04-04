
from pypion.example_plots_yt import yt_example_plots

ytp = yt_example_plots("/mnt/massive-stars/data/ostar-xray/mhd_d3n0128l3_s7/", "Ostar_mhd_d3n0128l3_s7", "img", "Ostar_mhd_d3n0128l3_s7", quantities=["density"])
for i in range(0,ytp.get_nframes(),100):
  if i==0: continue
  #ytp.plot_grid_plot(i)
  #ytp.plot_Bfield_XY(i, 1e-27, 1e-22)
  #ytp.plot_Bfield_Bmag(i, 1.0e-9, 3.0e-5, 1)
  #ytp.plot_Bfield_both(i, 1e-27, 1e-22)
  #ytp.plot_density_slice_both(i, 1e-27, 1e-22)
  #ytp.plot_Xray_intensity(i, 1e-18, 1e-15)
  #ytp.plot_projected_quantity(i, "density", ["density"], [0,0,1], [1,0,0])
  #ytp.plot_projected_quantity(i, "density", ["density"], [1,0,0], [0,1,1])
  ytp.plot_projected_quantity(i, "xray_0.3", ["xray_emission"], [1,0,0], [0,1,1])
quit()

