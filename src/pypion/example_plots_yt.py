import sys
import os
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from .pion2yt import pion2yt
import yt
from yt.visualization.volume_rendering.api import PointSource
yt.set_log_level("ERROR")
#print(plt.style.available)
#plt.style.use('classic')

from matplotlib import ticker
from matplotlib import cm
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colorbar import Colorbar
import sys
import pathlib
home = str(pathlib.Path.home())
#from misc.make_movies import make_movies
from datetime import datetime


class yt_example_plots(pion2yt):

    # Constructor
    def __init__(self, data_path, file_base, img_dir, sim_name, quantities=["density"]):

        #self.p2y = pion2yt

        self.evolution = self.make_snapshots(data_path, file_base) # list of all snapshots
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
    def plot_projected_quantity(self, i, plot_quantity, ds_quantities, north_vec, norm, **kwargs):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], ds_quantities)
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting {plot_quantity} for snapshot {i}...")
        prj = yt.OffAxisProjectionPlot(ds, normal=norm, north_vector=north_vec, fields=("gas", plot_quantity))
        prj.set_cmap(("gas", plot_quantity), kwargs.get("cmap", "viridis"))
        prj.set_figure_size(kwargs.get("figsize", 5))
        # prj.set_zlim(("gas", plot_quantity), kwargs.get("zlim", None))

        fig = prj.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        #ax.set_xlabel("x (AU)")
        #ax.set_ylabel("y (AU)")
        st = r"$t=$ = " + f"{time:.5f}"
        ax.text(0.1, 0.9, st, color="black", fontsize=8,
                transform=ax.transAxes,
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
        num = str(i).zfill(5)
        print(f"Saving image {self.sim_name}_{plot_quantity}_{num}.png...")
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_{plot_quantity}_{num}.png"), dpi=300, bbox_inches="tight")
        #return fig


    ###########################################################################
    def temp_quantity_plotter(self, plot_quantity, ds_quantities, north_vec, norm, **kwargs):
        self.plot_indices = self.three_slice_indices() # list of indices for pre-periastron, periastron, and post-periastron snapshots
        num = 0
        for index in self.plot_indices:
            fig = self.plot_projected_quantity(index, plot_quantity, ds_quantities, north_vec, norm, **kwargs)
            fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_{plot_quantity}_{index}.png"), dpi=300, bbox_inches="tight")
            num += 1


    ###########################################################################
    def plot_grid_plot(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["density"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting SMR grid for snapshot {i}...")
        slc = yt.SlicePlot(ds, "z", "density")
       
        slc.set_cmap("density", "GnBu")
        slc.set_figure_size(5)
        slc.zoom(1)
        slc.set_log("density", True)
        #slc.set_zlim("density", )
        slc.annotate_grids(linewidth=1, edgecolors="black")
        slc.annotate_cell_edges(line_width=0.00001, alpha=0.6, color="black")
        # slc.set_font({"family": "times new roman"})
        slc.set_font({"family": "Sans", "size": 14, "weight": "normal"})
        fig = slc.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        t = ax.get_xlabel()
        ax.set_xlabel(t, weight='bold')  # not working because it is mathmode

        #ax.legend(loc="upper right", frameon=True, framealpha=1, facecolor="white", fontsize=14)
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_grid_{num}.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_grid_{num}.pdf"), dpi=300, bbox_inches="tight")

    ###########################################################################
    def plot_Bfield_XY(self, i, dmin, dmax):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["density", "magnetic_field"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        # Coordinates are in (z,y,x) ordering.
        #slc = yt.SlicePlot(ds, "z", "density", width=(4.5, "pc"),origin="native")
        slc = yt.SlicePlot(ds, "z", "density", origin="native")
        slc.set_cmap("density", "viridis")
        slc.set_figure_size(5)
        slc.set_log("density", True)
        slc.set_zlim(("gas", "density"), zmin=(dmin, "g/cm**3"), zmax=(dmax, "g/cm**3"))
        slc.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), color="white", factor=1, density=1.5)
        #slc.annotate_grids()
        #slc.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        fig = slc.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        #ax.set_xlabel("x (pc)")
        #ax.set_ylabel("y (pc)")
        #num = "{:05.4f}".format(time.value) + " Myr"
        num = "{:05.4f}".format(time)
        ax.text(0.725, 0.95, num, transform=ax.transAxes, fontsize=16)

        num = str(i).zfill(5)

        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bfield_contours_XY_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)


    ###########################################################################
    def plot_Bfield_ZX(self, i, dmin, dmax):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["density", "magnetic_field"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        # Coordinates are in (z,y,x) ordering.
        # plot density
        #slc = yt.SlicePlot(ds, "y", "density", width=(4.5, "pc"),origin="native")
        slc = yt.SlicePlot(ds, "y", "density", origin="native")
        slc.swap_axes()
        slc.set_cmap("density", "viridis")
        slc.set_figure_size(5)
        slc.set_log("density", True)
        slc.set_zlim(("gas", "density"), zmin=(dmin, "g/cm**3"), zmax=(dmax, "g/cm**3"))
        # add streamlines
        slc.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=1, density=1.5)
        #slc.hide_colorbar("density")
        #slc.annotate_grids()
        #slc.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        fig = slc.export_to_mpl_figure((1,1), cbar_mode="single") # export to matplotlib figure
        ax = fig.axes[0]
        #ax.set_xlabel("x (pc)")
        #ax.set_ylabel("z (pc)")
        num = "{:05.4f}".format(time.value) + " Myr"
        ax.text(0.725, 0.95, num, transform=ax.transAxes, fontsize=16)
  
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bfield_contours_ZX_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_Bfield_Bmag(self, i, bmin, bmax, zoom):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["magnetic_field"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho+B for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075, 0.075, 0.85, 0.85),
                nrows_ncols = (1, 2),
                axes_pad = 0.075,
                label_mode = "all",
                share_all = True,
                cbar_mode="single",
                cbar_location="right",
                cbar_size="7%",
                cbar_pad="28%")


        # first plot: XZ plane, on the left
        #slc1 = yt.SlicePlot(ds, "y", "magnetic_field_magnitude", width=(4.5, "pc"),origin="native")
        slc1 = yt.SlicePlot(ds, "z", "magnetic_field_magnitude", origin="native")
        #slc1.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc1.swap_axes()
        slc1.set_cmap("magnetic_field_magnitude", "magma")
        slc1.set_zlim(("gas", "magnetic_field_magnitude"), zmin=(bmin, "G"), zmax=(bmax, "G"))
        slc1.set_log("magnetic_field_magnitude", True)
        slc1.zoom(zoom)
        slc1.set_figure_size(7)
        #slc1.hide_colorbar()
        st = r"$t=$ " + f"{time:.3f}" 
        slc1.annotate_text([0.325, 0.25], st, coord_system="figure", text_args={"color": "black", "fontsize": 12}, inset_box_args={'boxstyle': 'square', 'pad': 0.3, 'facecolor': 'white', 'linewidth': 1, 'edgecolor': 'black', 'alpha': 0.5}) # , fontsize=14
        # add streamlines
        slc1.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=zoom, density=1.3)
        #slc1.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        plot1 = slc1.plots['magnetic_field_magnitude']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        slc1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="in", width=1, length=2)
        #plot1.axes.yaxis.set_label_coords(-0.1,0.5)
        #plot1.axes.xaxis.set_label_coords(0.5,-0.1)
      
        # second plot: XY plane, on the right
        #slc2 = yt.SlicePlot(ds, "z", "magnetic_field_magnitude", width=(4.5, "pc"),origin="native")
        slc2 = yt.SlicePlot(ds, "z", "magnetic_field_magnitude", origin="native")
        #slc2.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc2.set_cmap("magnetic_field_magnitude", "magma")
        slc2.set_figure_size(7)
        slc2.set_log("magnetic_field_magnitude", True)
        slc2.zoom(zoom)
        slc2.set_zlim(("gas", "magnetic_field_magnitude"), zmin=(bmin, "G"), zmax=(bmax, "G"))
        slc2.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), color="white", factor=zoom, density=1.3)

        #slc2.hide_colorbar()
        plot2 = slc2.plots['magnetic_field_magnitude']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        plot2.axes.yaxis.tick_right()
        slc2._setup_plots()
        plot2.axes.yaxis.set_label_position("right")
        #plot2.axes.yaxis.set_label_coords(1.1,0.5)
        #plot2.axes.xaxis.set_label_coords(0.5,-0.1)
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="in", width=1, length=2)
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bmag_contours_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_Bfield_both(self, i, dmin, dmax):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["density", "magnetic_field"])
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
        #slc1 = yt.SlicePlot(ds, "y", "density", width=(4.5, "pc"),origin="native")
        slc1 = yt.SlicePlot(ds, "y", "density",origin="native")
        #slc1.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc1.swap_axes()
        slc1.set_cmap("density", "viridis")
        slc1.set_zlim(("gas", "density"), zmin=(dmin, "g/cm**3"), zmax=(dmax, "g/cm**3"))
        slc1.set_log("density", True)
        slc1.set_figure_size(7)
        #slc1.hide_colorbar()
        st = r"$t=$ " + f"{time:.3f}" 
        slc1.annotate_text([0.325, 0.25], st, coord_system="figure", text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14
        # add streamlines
        slc1.annotate_streamlines(("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), color="white", factor=1, density=1.3)
        #slc1.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"),lim=(0.5,0.65))

        plot1 = slc1.plots['density']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        #cb1.set_label("$\\rho$ (g cm$^{-3}$)", loc='right', fontsize=14)
        slc1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="in", width=1, length=2)
        #plot1.axes.yaxis.set_label_coords(-0.1,0.5)
        #plot1.axes.xaxis.set_label_coords(0.5,-0.1)
      
        # second plot: XY plane, on the right
        #slc2 = yt.SlicePlot(ds, "z", "density", width=(4.5, "pc"),origin="native")
        slc2 = yt.SlicePlot(ds, "z", "density", origin="native")
        #slc2.set_font({'family': 'serif', 'style': 'normal', 'weight': 'normal', 'size': 12})
        slc2.set_cmap("density", "viridis")
        slc2.set_figure_size(7)
        slc2.set_log("density", True)
        #slc2.zoom(16)
        slc2.set_zlim(("gas", "density"), zmin=(dmin, "g/cm**3"), zmax=(dmax, "g/cm**3"))
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
        #plot2.axes.yaxis.set_label_coords(1.1,0.5)
        #plot2.axes.xaxis.set_label_coords(0.5,-0.1)
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="in", width=1, length=2)
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Bfield_contours_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def plot_density_2d(self, i, dmin, dmax, Tmin, Tmax):
        # Plots density and temperature on a log scale, with T reflected onto
        # the lower half-plane.
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["density", "temperature"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting rho-2d for snapshot {i}...")
        
        fig = plt.figure()
        grid = ImageGrid(fig, (0.075, 0.075, 0.85, 0.85),
                nrows_ncols = (2, 1),
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
        slc1.set_zlim(("gas", "density"), zmin=(dmin, "g/cm**3"), zmax=(dmax, "g/cm**3"))
        slc1.set_log("density", True)
        slc1.set_figure_size(7)
        st = r"$t=$ " + f"{time:.3f}" 
        slc1.annotate_text([0.85, 0.525], st, coord_system="figure", text_args={"color":"black", "fontsize":12}, inset_box_args={'boxstyle':'square', 'pad':0.3, 'facecolor':'white', 'linewidth':1, 'edgecolor':'black', 'alpha':0.5}) # , fontsize=14

        plot1 = slc1.plots['density']
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]

        slc1._setup_plots()
        plot1.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot1.axes.tick_params(which="minor", direction="inout", width=1, length=3)
        #plot1.axes.set_yticks([0,1,2])
      
        # second plot: lower half-plane, temperature
        slc2 = yt.SlicePlot(ds, "theta", "temperature")
        slc2.swap_axes()
        slc2.flip_vertical() # reflect the domain
        slc2.set_cmap("temperature", "plasma")
        slc2.set_figure_size(7)
        slc2.set_log("temperature", True)
        slc2.set_zlim(("gas", "temperature"), zmin=(Tmin, "K"), zmax=(Tmax, "K"))
        plot2 = slc2.plots['temperature']
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        slc2._setup_plots()
        plot2.axes.tick_params(which="major", direction="inout", width=1, length=8)
        plot2.axes.tick_params(which="minor", direction="inout", width=1, length=3)
        #plot2.axes.set_yticks([1,2])
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_rho2d_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    ###########################################################################
    def plot_Vfield_Bfield_2d(self, i, Vmin, Vmax, Bmin, Bmax):
        '''
        Plots a 2D slice of the velocity and magnetic field magnitudes from a simulation snapshot.

        The upper panel shows the velocity magnitude with velocity vectors overlaid.
        The lower panel shows the magnetic field magnitude with magnetic field vectors overlaid.

        Parameters:
        -----------
        i : int
            Index of the simulation snapshot to load and plot.
        Vmin : float
            Minimum value for the velocity magnitude color scale.
        Vmax : float
            Maximum value for the velocity magnitude color scale.
        Bmin : float
            Minimum value for the magnetic field magnitude color scale.
        Bmax : float
            Maximum value for the magnetic field magnitude color scale.

        Returns:
        --------
        None
        '''

        # Define a derived field: magnitude of velocity vector
        def _V_magnitude(field, data):
            Vx = data[("gas", "velocity_x")]
            Vy = data[("gas", "velocity_y")]
            return np.sqrt(Vx ** 2 + Vy ** 2)

        # Define a derived field: magnitude of magnetic field vector
        def _B_magnitude(field, data):
            Bx = data[("gas", "magnetic_field_x")]
            By = data[("gas", "magnetic_field_y")]
            return np.sqrt(Bx ** 2 + By ** 2)

        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["velocity", "magnetic_field"])

        # Register derived fields in yt dataset
        ds.add_field(("gas", "velocity_field_magnitude"),
                     function=_V_magnitude,
                     units="km/s",
                     sampling_type="cell")

        ds.add_field(("gas", "magnetic_field_magnitude"),
                     function=_B_magnitude,
                     units="gauss",
                     sampling_type="cell")

        # Current simulation time
        time = ds.current_time.to("kyr")
        print(f"Time: {time}")
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting velocity and magnetic field for snapshot {i}...")

        # Set up a figure with two rows for the two panels
        fig = plt.figure(figsize=(10, 14))
        grid = ImageGrid(fig, (0.075, 0.075, 0.85, 0.85),
                         nrows_ncols=(2, 1),
                         axes_pad=0.0,
                         label_mode="all",
                         share_all=False,
                         cbar_mode="each",
                         cbar_size="5%",
                         cbar_pad="1%")

        # ================= Upper Panel: Velocity =================
        slc1 = yt.SlicePlot(ds, "theta", ("gas", "velocity_field_magnitude"))
        slc1.swap_axes()
        slc1.set_cmap(("gas", "velocity_field_magnitude"), "magma")
        slc1.set_zlim(("gas", "velocity_field_magnitude"), Vmin, Vmax)
        slc1.set_figure_size(7)

        # Add velocity vectors (quiver plot)
        slc1.annotate_quiver(("gas", "velocity_y"), ("gas", "velocity_x"), color="white",
                             factor=25, scale=45, normalize=True)

        # Annotate time and label
        st = r"$t=$ " + f"{time:.3f}"
        slc1.annotate_text([0.85, 0.325], st, coord_system="figure", text_args={"color": "black", "fontsize": 2},
                           inset_box_args={'boxstyle': 'square', 'pad': 0.3, 'facecolor': 'white',
                                           'linewidth': 1, 'edgecolor': 'black', 'alpha': 0.5})
        slc1.annotate_text([0.85, 0.55], r"$|\vec{v}|$ km/s", coord_system="figure",
                           text_args={"color": "black", "fontsize": 2},
                           inset_box_args={'boxstyle': 'square', 'pad': 0.3, 'facecolor': 'white',
                                           'linewidth': 1, 'edgecolor': 'black', 'alpha': 0.5})

        # Set plot to corresponding grid location
        plot1 = slc1.plots[("gas", "velocity_field_magnitude")]
        plot1.figure = fig
        plot1.axes = grid[0].axes
        plot1.cax = grid.cbar_axes[0]
        slc1._setup_plots()
        grid.cbar_axes[0].set_ylabel("")

        # ================= Lower Panel: Magnetic Field =================
        slc2 = yt.SlicePlot(ds, "theta", ("gas", "magnetic_field_magnitude"))
        slc2.swap_axes()
        slc2.flip_vertical()  # Flip to match expected orientation
        slc2.set_cmap(("gas", "magnetic_field_magnitude"), "magma")
        slc2.set_zlim(("gas", "magnetic_field_magnitude"), Bmin, Bmax)
        slc2.set_log(("gas", "magnetic_field_magnitude"), True)
        slc2.set_figure_size(7)

        # Add magnetic field vectors
        slc2.annotate_quiver(("gas", "magnetic_field_y"), ("gas", "magnetic_field_x"), color="white",
                             factor=25, scale=45, normalize=True)

        # Annotate label
        slc2.annotate_text([0.15, 0.55], r"$|\vec{B}|$ G", coord_system="figure",
                           text_args={"color": "black", "fontsize": 2},
                           inset_box_args={'boxstyle': 'square', 'pad': 0.3, 'facecolor': 'black',
                                           'linewidth': 1, 'edgecolor': 'black', 'alpha': 0.1})

        # Set plot to corresponding grid location
        plot2 = slc2.plots[("gas", "magnetic_field_magnitude")]
        plot2.figure = fig
        plot2.axes = grid[1].axes
        plot2.cax = grid.cbar_axes[1]
        slc2._setup_plots()
        grid.cbar_axes[1].set_ylabel("")

        # Save the resulting figure to file
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_VBfield2d_{num}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close(fig)

    ###########################################################################
    def get_1ddata_from_2ddata(self, i, length, theta, resolution=64):
        import numpy as np
        from scipy.interpolate import griddata

        # Convert angle to radians
        theta_rad = np.radians(theta)

        # Define the number of points along the line and generate indices
        # For simplicity, we will assume you want to map the length to the grid resolution
        indices = np.linspace(0, resolution - 1, resolution)
        r_indices = indices * np.cos(theta_rad)
        z_indices = indices * np.sin(theta_rad)

        print(r_indices, z_indices)

        line_points = np.vstack((r_indices, z_indices)).T  # shape: (resolution, 2)

        # Load dataset
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["velocity", "magnetic_field", "density"])
        ad = ds.all_data()

        # Get grid cell positions in r and z (indices, not physical coordinates)
        r_index = ad[('index', 'r')].v  # Get r indices
        z_index = ad[('index', 'z')].v  # Get z indices

        # Get the density field directly from the grid data
        density = ad['density'].v

        # Flatten the grid data for easier processing
        pos_indices = np.vstack([r_index.ravel(), z_index.ravel()]).T  # shape: (n_cells, 2)

        print(pos_indices)
        exit(1)

        # Interpolate density along the line using the grid indices
        interpolated_density = griddata(pos_indices, density.ravel(), line_points, method='linear')

        # Calculate distance along the line (Euclidean distance) based on grid indices
        line_distances = np.sqrt(np.diff(r_indices) ** 2 + np.diff(z_indices) ** 2)

        # Compute cumulative distance along the line from origin
        distance = np.cumsum(line_distances)

        # Remove NaNs if interpolation falls outside domain
        mask = ~np.isnan(interpolated_density)
        distance = distance[mask]
        interpolated_density = interpolated_density[mask]

        # Output the distance and interpolated density
        print("Distance along the line (in grid units):", distance)
        print("Interpolated density values:", interpolated_density)

    ###########################################################################
    def plot_Xray_intensity(self, i, jmin, jmax):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["xray_emission"])
        time = (ds.current_time.to("Myr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting x-ray intensity for snapshot {i}...")
        
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
        #prj1 = yt.ProjectionPlot( ds, "y", ("gas", "xray_0.3"),  width=(4.5, "pc"),origin="native", method="integrate", buff_size=(1024, 1024))
        prj1 = yt.ProjectionPlot( ds, "y", ("gas", "xray_0.3"), origin="native", method="integrate", buff_size=(1024, 1024))
        prj1.swap_axes()
        prj1.set_cmap("xray_0.3", "magma")
        prj1.set_zlim(("gas", "xray_0.3"), zmin=(jmin, "erg/cm**2/s/arcmin**2"), zmax=(jmax, "erg/cm**2/s/arcmin**2"))
        prj1.set_log("xray_0.3", True)
        prj1.set_figure_size(7)
        st = r"$t=$ " + f"{time:.3f}" 
        prj1.annotate_text([0.675, 0.265], st, coord_system="figure", text_args={"color":"black", "fontsize":12},
                           inset_box_args={'boxstyle': 'square', 'pad': 0.3, 'facecolor':'white', 'linewidth':1,
                                           'edgecolor':'black', 'alpha': 0.5})
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
        #prj2 = yt.ProjectionPlot( ds, "z", ("gas", "xray_0.3"), width=(4.5, "pc"),origin="native", method="integrate", buff_size=(1024, 1024))
        prj2 = yt.ProjectionPlot( ds, "z", ("gas", "xray_0.3"), origin="native", method="integrate", buff_size=(1024, 1024))
        prj2.set_cmap("xray_0.3", "magma")
        prj2.set_zlim(("gas", "xray_0.3"), zmin=(jmin, "erg/cm**2/s/arcmin**2"), zmax=(jmax, "erg/cm**2/s/arcmin**2"))
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
# i    - index in dataset
# jmin - min of intensity colour scale
# jmax - max of intensity colour scale
# z    - zoom level (integer 2^n) >=1
    def plot_Xray_intensity_XY(self, i, jmin, jmax, z):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = self.get_ds(self.evolution[i], quantities=["xray_emission"])
        time = (ds.current_time.to("yr"))
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting x-ray intensity for snapshot {i}... for $\hat{z}$ projection")
        # XY plane
        #prj2 = yt.ProjectionPlot( ds, "z", ("gas", "xray_0.3"), width=(4.5, "pc"),origin="native", method="integrate", buff_size=(1024, 1024))
        prj2 = yt.ProjectionPlot( ds, "z", ("gas", "xray_0.3"), origin="native", method="integrate", buff_size=(1024, 1024))
        prj2.set_cmap("xray_0.3", "magma")
        prj2.set_zlim(("gas", "xray_0.3"), zmin=(jmin, "erg/cm**2/s/arcmin**2"), zmax=(jmax, "erg/cm**2/s/arcmin**2"))
        prj2.set_log("xray_0.3", True)
        prj2.zoom(z)
        prj2.set_figure_size(7)

        fig = prj2.export_to_mpl_figure((1,1), cbar_mode="single")
        ax = fig.axes[0]
        num = "{:05.4f}".format(time)
        ax.text(0.725, 0.95, num, transform=ax.transAxes, fontsize=16)
        ax.tick_params(which="major", direction="inout", width=1, length=8)
        ax.tick_params(which="minor", direction="in", width=1, length=2)
        
        num = str(i).zfill(5)
        fig.savefig(os.path.join(self.img_dir, f"{self.sim_name}_Proj_Xray003_XY_{num}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

###########################################################################
###########################################################################

