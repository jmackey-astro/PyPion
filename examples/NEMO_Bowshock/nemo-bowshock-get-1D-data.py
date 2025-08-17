from pypion.example_plots_yt import yt_example_plots



ytp = yt_example_plots("/home/tony/Desktop/multi-ion-bowshock/sim-output/silo/", "Ostar_mhd-nemo_d2n0128l3", "img", "Ostar_mhd-nemo_d2n0128l3", quantities=["density"])




for i in range(0, ytp.get_nframes(), 1):


  ytp.get_1ddata_from_2ddata(i, 1.0, 60, 10)