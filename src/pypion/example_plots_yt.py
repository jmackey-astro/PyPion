from silo_to_yt import *
from read_torus_fits import *
import sys
import pickle
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import pandas as pd
import sigfig
from yt.visualization.volume_rendering.api import PointSource
yt.set_log_level("ERROR")
#print(plt.style.available)
plt.style.use('classic')

from matplotlib import ticker
from matplotlib import cm
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colorbar import Colorbar
import sys
import pathlib
home = str(pathlib.Path.home())
#from misc.make_movies import make_movies
from datetime import datetime


class YTPlotFunction():

    # Constructor
    def __init__(self, data_path, file_base, img_dir, sim_name, quantities=["density"]):

        self.evolution = make_snapshots(data_path, file_base) # list of all snapshots
        self.data_path = data_path # path to data
        self.file_base = file_base # base of filename
        self.img_dir = img_dir # path to save images
        self.sim_name = sim_name
        if not os.path.exists(self.img_dir):
            os.makedirs(self.img_dir)
            print("Created images directory:{} ".format(self.img_dir))
        else:
            print("Images directory already exists")
        self.quantities = quantities # list of quantities to include in dataset

    def get_nframes(self):
        return self.evolution.shape[0]

    ###########################################################################
    # Calculate yt dataset of synchrotron emissivity
    def add_synchrotron_emission(self, ds):
        from unyt import K, g, cm, s
        def _Isync(field, data):
            Isync = data["magnetic_field_magnitude"]**(3/2) * data["pressure"] * data["NG_Mask"]
            return Isync * (K*g**(7/4)/(cm**(3/4)*s**(7/2)))**-1 * cm

        # add synchrotron emission field to dataset
        ds.add_field(("gas", "Isync"), function=_Isync, units="auto", sampling_type="cell", force_override=True) 


    ###########################################################################
    def plot_synchrotron_emission(self, i, north_vec, norm, **kwargs):

        """
        Function that creates synchrotron emission maps using the
        yt.OffAxisProjectionPlot function. Converts the yt figure
        into a mpl figure and returns the mpl fig object.

        Parameters
        ----------
        i : int

            index of snapshot to plot

        north_vec : list

            list of three floats that define the north vector of the
            projection plot
        
        norm : list

            list of three floats that define the normal vector to the
            projection plot

        **kwargs : dict
        
            dictionary of keyword arguments to pass to the yt.OffAxisProjectionPlot
        
        Returns
        -------
        fig : matplotlib.figure.Figure

        """

        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["magnetic_field", "pressure", "NG_Mask"])
        time = (ds.current_time.to("Myr")).value
        print(f"Successfully Loaded dataset: {str(ds)} at time {time:12.3e} Myr")
        self.add_synchrotron_emission(ds)
        print(f"Added Synchrotron emission field")
        print(f"Plotting Isync for snapshot {i}...")
        prj = yt.OffAxisProjectionPlot(ds, normal=norm, north_vector=north_vec, fields=("gas", "Isync")) # create projection plot
        prj.set_cmap(("gas", "Isync"), "gist_heat")
        prj.hide_colorbar("Isync")
        prj.set_figure_size(5)
        prj.zoom(kwargs.get("zoom", 1))
        prj.set_zlim(("gas", "Isync"), 1e10, 5e15)
        # prj.set_font({"family": "sans-serif", "size": 8})

        fig = prj.export_to_mpl_figure((1,1), cbar_mode=kwargs.get("cbarmode", "single")) # export to matplotlib figure
        ax = fig.axes[0] 
        ax.set_xlabel("x (pc)") 
        ax.set_ylabel("y (pc)")
        st = r"$t=$ " + f"{time:.2f} Myr" 
        ax.text(0.8, 0.9, st, color="black", fontsize=12, transform=ax.transAxes, 
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
        return fig 


    ###########################################################################
    def plot_projected_quantity(self, i, plot_quantity, ds_quantities, north_vec, norm, **kwargs):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], ds_quantities, start_time=self.start_time)
        time = ds.current_time.value
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting {plot_quantity} for snapshot {i}...")
        prj = yt.OffAxisProjectionPlot(ds, normal=norm, north_vector=north_vec, fields=("gas", plot_quantity))
        prj.set_cmap(("gas", plot_quantity), kwargs.get("cmap", "viridis"))
        prj.set_figure_size(kwargs.get("figsize", 5))
        # prj.set_zlim(("gas", plot_quantity), kwargs.get("zlim", None))
        print(f"Saving image {plot_quantity}_{i}.png...")
        phase = ((self.period + (time*u.s).to(u.yr))/self.period).value

        fig = prj.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("x (AU)")
        ax.set_ylabel("y (AU)")
        st = r"$\phi$ = " + f"{phase:.2f}"
        ax.text(0.1, 0.9, st, color="black", fontsize=8,
                transform=ax.transAxes,
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
        
        return fig
    
    ###########################################################################
    def temp_quantity_plotter(self, plot_quantity, ds_quantities, north_vec, norm, **kwargs):
        self.plot_indices = self.three_slice_indices() # list of indices for pre-periastron, periastron, and post-periastron snapshots
        num = 0
        for index in self.plot_indices:
            fig = self.plot_projected_quantity(index, plot_quantity, ds_quantities, north_vec, norm, **kwargs)
            fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_{plot_quantity}_{index}.png"), dpi=300, bbox_inches="tight")
            num += 1


    ###########################################################################
    def plot_smr_plot(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["windtracer"], start_time=self.start_time)
        time = ds.current_time.value
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting SMR for snapshot {i}...")
        slc = yt.SlicePlot(ds, "z", "windtracer")
       
        slc.set_cmap("windtracer", "gray")
        slc.set_figure_size(5)
        slc.zoom(2)
        slc.set_log("windtracer", False)
        slc.set_zlim("windtracer", -10, 0)
        slc.annotate_grids(linewidth=1, edgecolors="black")
        slc.annotate_cell_edges(line_width=0.00001, alpha=0.5, color="black")
        # slc.set_font({"family": "times new roman"})
        slc.set_font({"family": "mpl-default", "size": 11, "weight": "bold"})
        fig = slc.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("$\mathrm{x}$ (pc)")
        ax.set_ylabel("y (pc)")

        ax.legend(loc="upper right", frameon=True, framealpha=1, facecolor="white", fontsize=14)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_smr_{i}.png"), dpi=300, bbox_inches="tight")

    ###########################################################################
    def plot_Bfield_XY(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["density", "magnetic_field"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        # Coordinates are in (z,y,x) ordering.
        slc = yt.SlicePlot(ds, "z", "density", width=(4.5, "pc"),origin="native")
        slc.set_cmap("density", "viridis")
        slc.set_figure_size(5)
        slc.set_log("density", True)
        slc.set_zlim(("gas", "density"), zmin=(10e-28, "g/cm**3"), zmax=(1e-22, "g/cm**3"))
        slc.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), color="white", factor=1, density=1.5)
        #slc.annotate_grids()
        #slc.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        fig = slc.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("x (pc)")
        ax.set_ylabel("y (pc)")
        num = "{:05.4f}".format(time.value) + " Myr"
        ax.text(0.725, 0.95, num, transform=ax.transAxes, fontsize=16)

        num = str(i).zfill(5)

        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bfield_contours_XY_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)


    ###########################################################################
    def plot_Bfield_ZX(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["density", "magnetic_field"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        # Coordinates are in (z,y,x) ordering.
# plot density
        slc = yt.SlicePlot(ds, "y", "density", width=(4.5, "pc"),origin="native") #, center=[0.5,0.5,0.7857142857])
        slc.swap_axes()
        slc.set_cmap("density", "viridis")
        slc.set_figure_size(5)
        slc.set_log("density", True)
        slc.set_zlim(("gas", "density"), zmin=(10e-28, "g/cm**3"), zmax=(1e-22, "g/cm**3"))
        #slc.set_log("density", False)
        #slc.set_zlim(("gas", "density"), zmin=-27.5, zmax=-22.0)
# plot |B|
#        slc = yt.SlicePlot(ds, "y", "magnetic_field_magnitude", width=(4.5, "pc"),origin="native") #, center=[0.5,0.5,0.7857142857])
#        slc.set_cmap("magnetic_field_magnitude", "magma")
#        slc.set_figure_size(5)
#        slc.set_log("magnetic_field_magnitude", True)
#        slc.set_zlim(("gas", "magnetic_field_magnitude"), zmin=(1e-8, "G"), zmax=(1e-4, "G"))
# add streamlines
        slc.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=1, density=1.5)
        #slc.hide_colorbar("density")
        #slc.annotate_grids()
        #slc.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        fig = slc.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("x (pc)")
        ax.set_ylabel("z (pc)")
        num = "{:05.4f}".format(time.value) + " Myr"
        ax.text(0.725, 0.95, num, transform=ax.transAxes, fontsize=16)
  
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bfield_contours_ZX_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)


    ###########################################################################
    def plot_Bfield_Bmag(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["magnetic_field"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 2),
                axes_pad = 0.075,
                label_mode = "all",
                share_all = True,
                cbar_mode="single",
                cbar_location="right",
                cbar_size="7%",
                cbar_pad="28%")


        # first plot: XZ plane, on the left
        slc1 = yt.SlicePlot(ds, "y", "magnetic_field_magnitude", width=(4.5, "pc"),origin="native")
        #slc1.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc1.swap_axes()
        slc1.set_cmap("magnetic_field_magnitude", "magma")
        slc1.set_zlim(("gas", "magnetic_field_magnitude"), zmin=(1e-9, "G"), zmax=(3e-5, "G"))
        slc1.set_log("magnetic_field_magnitude", True)
        slc1.set_figure_size(7)
        #slc1.hide_colorbar()
        st = r"$t=$ " + f"{time:.3f}" 
        slc1.annotate_text([5.5e18,0,-4e18], st, text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14
        # add streamlines
        slc1.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=1, density=1.3)
        #slc1.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        plot1 = slc1.plots['magnetic_field_magnitude']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        slc1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="in", width=1, length=2)
        plot1.axes.yaxis.set_label_coords(-0.1,0.5)
        plot1.axes.xaxis.set_label_coords(0.5,-0.1)
      
        # second plot: XY plane, on the right
        slc2 = yt.SlicePlot(ds, "z", "magnetic_field_magnitude", width=(4.5, "pc"),origin="native")
        #slc2.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc2.set_cmap("magnetic_field_magnitude", "magma")
        slc2.set_figure_size(7)
        slc2.set_log("magnetic_field_magnitude", True)
        #slc2.zoom(16)
        slc2.set_zlim(("gas", "magnetic_field_magnitude"), zmin=(1e-9, "G"), zmax=(3e-5, "G"))
        slc2.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), color="white", factor=1, density=1.3)

        #slc2.hide_colorbar()
        plot2 = slc2.plots['magnetic_field_magnitude']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        plot2.axes.yaxis.tick_right()
        slc2._setup_plots()
        plot2.axes.yaxis.set_label_position("right")
        plot2.axes.yaxis.set_label_coords(1.1,0.5)
        plot2.axes.xaxis.set_label_coords(0.5,-0.1)
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="in", width=1, length=2)
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bmag_contours_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_Bfield_both(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["density", "magnetic_field"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 2),
                axes_pad = 0.075,
                label_mode = "all",
                share_all = True,
                cbar_mode="single",
                cbar_location="right",
                cbar_size="7%",
                cbar_pad="28%")


        # first plot: XZ plane, on the left
        slc1 = yt.SlicePlot(ds, "y", "density", width=(4.5, "pc"),origin="native")
        #slc1.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc1.swap_axes()
        slc1.set_cmap("density", "viridis")
        slc1.set_zlim(("gas", "density"), zmin=(3e-28, "g/cm**3"), zmax=(1e-22, "g/cm**3"))
        slc1.set_log("density", True)
        slc1.set_figure_size(7)
        #slc1.hide_colorbar()
        st = r"$t=$ " + f"{time:.3f}" 
        slc1.annotate_text([5.5e18,0,-4e18], st, text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14
        # add streamlines
        slc1.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=1, density=1.3)
        #slc1.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        #fig1 = slc1.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        #ax1 = fig1.axes[0] 
        #mm = ax1.images
        #cb = mm[-1].colorbar
        #cb.remove()
        #ax1.set_xlabel("x (pc)") 
        #ax1.set_ylabel("y (pc)")
        #ax1.text(0.8, 0.9, st, color="black", fontsize=12, transform=ax1.transAxes, 
        #            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))

        plot1 = slc1.plots['density']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]


        #new_ax = fig.add_axes((0.075,0.75,0.85,0.03))
        #cb1 = ax1.images.colorbar
        #cb1.remove()
        #cb2 = plot1.figure.colorbar(ax1.images, cax=ax1, orientation='horizontal', location='top') #, format='%.1E'
        #new_ax.tick_params(which="major", direction='inout', width=1, length=8)
        #new_ax.tick_params(which="minor", direction='in', width=1, length=2)
        #new_ax.tick_params(labelsize=14)
        
        #cb1.set_label("$\\rho$ (g cm$^{-3}$)", loc='right', fontsize=14)
        slc1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="in", width=1, length=2)
        plot1.axes.yaxis.set_label_coords(-0.1,0.5)
        plot1.axes.xaxis.set_label_coords(0.5,-0.1)
      
        # second plot: XY plane, on the right
        slc2 = yt.SlicePlot(ds, "z", "density", width=(4.5, "pc"),origin="native") #, center=[0.5,0.5,0.7857142857])
        #slc2.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc2.set_cmap("density", "viridis")
        slc2.set_figure_size(7)
        slc2.set_log("density", True)
        #slc2.zoom(16)
        slc2.set_zlim(("gas", "density"), zmin=(3e-28, "g/cm**3"), zmax=(1e-22, "g/cm**3"))
        slc2.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), color="white", factor=1, density=1.3)

        #slc2.annotate_grids()
        #slc2.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        #slc2.hide_colorbar()
        plot2 = slc2.plots['density']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        plot2.axes.yaxis.tick_right()
        #plot2.axes.yaxis.set_label_position("right")
        #new_ax = plot2.figure.add_axes((0.025,0.94,0.75,0.025))
        #cb2 = plot2.figure.colorbar(plot2.image, new_ax, extend='both', orientation='horizontal', format='%.0E')
        slc2._setup_plots()
        plot2.axes.yaxis.set_label_position("right")
        #plot2.axes.yaxis.set_label("right")
        plot2.axes.yaxis.set_label_coords(1.1,0.5)
        plot2.axes.xaxis.set_label_coords(0.5,-0.1)
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="in", width=1, length=2)
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bfield_contours_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_density_2d(self, i):
        # Plots density and temperature on a log scale, with T reflected onto
        # the lower half-plane.
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["density", "temperature"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho-2d for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2,1),
                axes_pad = 0.0,
                label_mode = "all",
                share_all = False,
                cbar_mode="each",
                cbar_location="right",
                cbar_size="5%",
                cbar_pad="1%")


        # first plot: upper half-plane, density
        slc1 = yt.SlicePlot(ds, "theta", "density")
        slc1.swap_axes()
        slc1.set_cmap("density", "viridis")
        slc1.set_zlim(("gas", "density"), zmin=(3e-28, "g/cm**3"), zmax=(2e-22, "g/cm**3"))
        slc1.set_log("density", True)
        slc1.set_figure_size(7)
        st = r"$t=$ " + f"{time:.3f}" 
        slc1.annotate_text([5.75e18,-10.4e18,0], st, text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14

        plot1 = slc1.plots['density']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        slc1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="inout", width=1, length=3)
        plot1.axes.set_yticks([0,1,2])
      
        # second plot: lower half-plane, temperature
        slc2 = yt.SlicePlot(ds, "theta", "temperature")
        slc2.swap_axes()
        slc2.flip_vertical() # reflect the domain
        slc2.set_cmap("temperature", "plasma")
        slc2.set_figure_size(7)
        slc2.set_log("temperature", True)
        slc2.set_zlim(("gas", "temperature"), zmin=(5000.0, "K"), zmax=(2e8, "K"))
        plot2 = slc2.plots['temperature']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        slc2._setup_plots()
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="inout", width=1, length=3)
        plot2.axes.set_yticks([1,2])
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_rho2d_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_density_slice_both(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["density"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho slices for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 2),
                axes_pad = 0.075,
                label_mode = "all",
                share_all = True,
                cbar_mode="single",
                cbar_location="right",
                cbar_size="7%",
                cbar_pad="28%")


        # first plot: XZ plane, on the left
        slc1 = yt.SlicePlot(ds, "y", "density", width=(4.5, "pc"),origin="native")
        #slc1.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc1.swap_axes()
        slc1.set_cmap("density", "viridis")
        slc1.set_zlim(("gas", "density"), zmin=(3e-28, "g/cm**3"), zmax=(1e-22, "g/cm**3"))
        slc1.set_log("density", True)
        slc1.set_figure_size(7)
        #slc1.hide_colorbar()
        st = r"$t=$ " + f"{time:.3f}" 
        slc1.annotate_text([5.5e18,0,-4e18], st, text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14

        plot1 = slc1.plots['density']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        slc1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="in", width=1, length=2)
        plot1.axes.yaxis.set_label_coords(-0.1,0.5)
        plot1.axes.xaxis.set_label_coords(0.5,-0.1)
      
        # second plot: XY plane, on the right
        slc2 = yt.SlicePlot(ds, "z", "density", width=(4.5, "pc"),origin="native") #, center=[0.5,0.5,0.7857142857])
        #slc2.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc2.set_cmap("density", "viridis")
        slc2.set_figure_size(7)
        slc2.set_log("density", True)
        #slc2.zoom(16)
        slc2.set_zlim(("gas", "density"), zmin=(3e-28, "g/cm**3"), zmax=(1e-22, "g/cm**3"))
        #slc2.hide_colorbar()
        plot2 = slc2.plots['density']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        plot2.axes.yaxis.tick_right()
        #plot2.axes.yaxis.set_label_position("right")
        #new_ax = plot2.figure.add_axes((0.025,0.94,0.75,0.025))
        #cb2 = plot2.figure.colorbar(plot2.image, new_ax, extend='both', orientation='horizontal', format='%.0E')
        slc2._setup_plots()
        plot2.axes.yaxis.set_label_position("right")
        #plot2.axes.yaxis.set_label("right")
        plot2.axes.yaxis.set_label_coords(1.1,0.5)
        plot2.axes.xaxis.set_label_coords(0.5,-0.1)
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="in", width=1, length=2)
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_rho_slices_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_Xray_intensity(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["xray_emission"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 2),
                axes_pad = 0.075,
                label_mode = "all",
                share_all = True,
                cbar_mode="single",
                cbar_location="right",
                cbar_size="7%",
                cbar_pad="28%")


        # first plot: XZ plane, on the left
        prj1 = yt.ProjectionPlot( ds, "y", ("gas", "xray_0.3"),  width=(4.5, "pc"),origin="native", method="integrate", buff_size=(1024, 1024))
        prj1.swap_axes()
        prj1.set_cmap("xray_0.3", "magma")
        prj1.set_zlim(("gas", "xray_0.3"), zmin=(3e-18, "erg/cm**2/s/arcmin**2"), zmax=(3e-15, "erg/cm**2/s/arcmin**2"))
        prj1.set_log("xray_0.3", True)
        prj1.set_figure_size(7)
        st = r"$t=$ " + f"{time:.3f}" 
        prj1.annotate_text([5.5e18,0,-4e18], st, text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14
        # add streamlines
        #prj1.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=1, density=1.3)
        #slc1.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        plot1 = prj1.plots['xray_0.3']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        prj1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="in", width=1, length=2)
        plot1.axes.yaxis.set_label_coords(-0.1,0.5)
        plot1.axes.xaxis.set_label_coords(0.5,-0.1)
      
        # second plot: XY plane, on the right
        prj2 = yt.ProjectionPlot( ds, "z", ("gas", "xray_0.3"),  width=(4.5, "pc"),origin="native", method="integrate", buff_size=(1024, 1024))
        prj2.set_cmap("xray_0.3", "magma")
        prj2.set_zlim(("gas", "xray_0.3"), zmin=(3e-18, "erg/cm**2/s/arcmin**2"), zmax=(3e-15, "erg/cm**2/s/arcmin**2"))
        prj2.set_log("xray_0.3", True)
        prj2.set_figure_size(7)

        plot2 = prj2.plots['xray_0.3']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        plot2.axes.yaxis.tick_right()
        prj2._setup_plots()
        plot2.axes.yaxis.set_label_position("right")
        plot2.axes.yaxis.set_label_coords(1.1,0.5)
        plot2.axes.xaxis.set_label_coords(0.5,-0.1)
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="in", width=1, length=2)
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Proj_Xray003_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_Xray_IR(self, i, irpath):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["xray_emission"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting X-ray + IR for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 2),
                axes_pad = 0.075,
                label_mode = "all",
                share_all = True,
                cbar_mode="single",
                cbar_location="right",
                cbar_size="7%",
                cbar_pad="28%")


        # first plot: XZ plane, on the left
        prj1 = yt.ProjectionPlot( ds, "y", ("gas", "xray_0.3"),  width=(4.5, "pc"), origin="native", method="integrate", buff_size=(1024, 1024))
        prj1.swap_axes()
        prj1.set_cmap("xray_0.3", "magma")
        prj1.set_zlim(("gas", "xray_0.3"), zmin=(3e-18, "erg/cm**2/s/arcmin**2"), zmax=(3e-15, "erg/cm**2/s/arcmin**2"))
        prj1.set_log("xray_0.3", True)
        prj1.set_figure_size(7)
        st = r"$t=$ " + f"{time:.3f}" 
        prj1.annotate_text([4.1e18,0,-2.75e18], st, text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14
        # add streamlines
        #prj1.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=1, density=1.3)
        #slc1.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        plot1 = prj1.plots['xray_0.3']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        prj1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="in", width=1, length=2)
        plot1.axes.yaxis.set_label_coords(-0.14,0.5)
        plot1.axes.xaxis.set_label_coords(0.5,-0.1)
        plot1.axes.set_xlim(-1,0.8)
        plot1.axes.set_ylim(-1.5,1.5)
        plot1.axes.set_xticks([-1,-0.5,0,0.5])
      
        # second plot: XY plane, on the right
        prj2 = yt.ProjectionPlot( ds, "z", ("gas", "xray_0.3"),  width=(4.5, "pc"), origin="native", method="integrate", buff_size=(1024, 1024))
        prj2.set_cmap("xray_0.3", "magma")
        prj2.set_zlim(("gas", "xray_0.3"), zmin=(3e-18, "erg/cm**2/s/arcmin**2"), zmax=(3e-15, "erg/cm**2/s/arcmin**2"))
        prj2.set_log("xray_0.3", True)
        prj2.set_figure_size(7)

        plot2 = prj2.plots['xray_0.3']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        plot2.axes.yaxis.tick_right()
        prj2._setup_plots()
        plot2.axes.yaxis.set_label_position("right")
        plot2.axes.yaxis.set_label_coords(1.15,0.5)
        plot2.axes.xaxis.set_label_coords(0.5,-0.1)
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="in", width=1, length=2)
        plot2.axes.set_xlim(-1,0.8)
        plot2.axes.set_ylim(-1.5,1.5)
        plot2.axes.set_xticks([-1,-0.5,0,0.5])

        #prj = ds.proj(("gas", ""xray_0.3""), 0)
        #frb = prj.to_frb((4.5, "pc"), 1024)
        ax = grid[0].axes
        ff = "../torus/"+irpath+"/Ostar_"+irpath+"_24um_XZ_t90.fits"
        ir = PlottingFITS(ff)
        #ir = PlottingFITS("../torus/hires/Ostar_mhdBY_d3n0384l3_s7_01453056_24um_t90.fits")
        data = ir.data()
        x_axis = ir.extents()[0] / 206264.80749673 / 2
        y_axis = ir.extents()[1] / 206264.80749673 / 2
        ext = [-y_axis, y_axis, -x_axis, x_axis]
        im1 = ax.contour(np.flip(data,axis=0),extent=ext,levels=[50,100,200,300,400,500,600],colors="white")

        ax = grid[1].axes
        ff = "../torus/"+irpath+"/Ostar_"+irpath+"_24um_XY_t90.fits"
        ir = PlottingFITS(ff)
        #ir = PlottingFITS("../torus/hires/Ostar_mhdBY_d3n0384l3_s7_01453056_24um_t90.fits")
        data = ir.data()
        x_axis = ir.extents()[0] / 206264.80749673 / 2
        y_axis = ir.extents()[1] / 206264.80749673 / 2
        ext = [-y_axis, y_axis, -x_axis, x_axis]
        im1 = ax.contour(np.flip(data,axis=0),extent=ext,levels=[50,100,200,300,400,500,600],colors="white")

        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Proj_Xray003_IR_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)


    ###########################################################################
    def volume_rendering(self, **kwargs):
        from yt.units import cm
        # Creating volume renderings at each time step in ts
        i=0
        for file in self.evolution[40:]:
            ds = get_ds(file)
            sc = yt.create_scene(ds)
            # Print the time of the current scene
            print(f"Time: {ds.current_time.to('s')}")
            time = np.float64(ds.current_time.to('Myr').value)
            print(type(time))
            # if kwargs['trajectory_file'] is not None:
                # star1_x, star1_y, star2_x, star2_y, star1_z, star2_z = self.get_star_position(kwargs['trajectory_file'], [time])
            # identifying the source
            source = sc[0]
            source.set_field('density')
            source.set_log(True)

            # building transfer function
            bounds = (1e-28, 1e-21)
            tf = yt.ColorTransferFunction(x_bounds=np.log10(bounds), nbins=500)
            # Automatically add a number of layers
            tf.add_layers(8, w=0.002, colormap='viridis')
            
            source.tfh.tf = tf
            source.tfh.bounds = bounds
            
            cam = sc.camera
            cam
            cam.zoom(2.2)
            cam.resolution = (4096, 4096)
            cam.switch_orientation(normal_vector=[-1,0,0], north_vector=[0.2,1,0])

            # colors = np.random.random([1, 4])
            # colors[:, 3] = 1.0
        

            # points = PointSource(np.array([star1_x*cm, star1_y*cm, star1_z*cm]), colors=colors, radii=5e12)
            # sc.add_source(points)

            sc.save(os.path.join(self.img_dir, f"sim-vol-img_{i+1}"), sigma_clip=6.0)
            print(f"Saving image {os.path.join(self.img_dir, f'sim-vol-img_{i+1}.png')}")
            del sc, cam
            i+=1

        #make_movies(self.img_dir, self.img_dir, "sim-vol-img.mp4")



###########################################################################
###########################################################################

