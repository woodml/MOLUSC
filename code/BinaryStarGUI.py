# MOLUSC v.20220321
# Mackenna Wood, UNC Chapel Hill
import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from astroquery.gaia import Gaia
import tkinter as tk
import tkinter.messagebox
from textwrap import wrap
from time import time
import multiprocessing as mp
import sys
import argparse
import gc
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

# Parallel Functions
def calculate_RV_parallel(period, mass_ratio, a, e, cos_i, arg_peri, phase, MJD, calc):
	# Exactly the same as calculate_RV, but with an extra parameter stating whether you need to calculate RV
	# Calculates the RVs for each item when passed arrays of orbital parameters
	# Inputs: Arrays of Period, Mass Ratio, Semi-Major Axis, eccentricity, inclination, arg peri, phase, calculation times
	# Outputs: Velocity Semi-Amplitude (km/s), RVs at each time in MJD

	sin_i = np.sin(np.arccos(cos_i))

	n = len(period)
	RV = [[0.0 for i in range(len(MJD))] for j in range(n)]
	a_star = np.multiply(a, np.divide(mass_ratio, np.add(mass_ratio, 1)))
	K = np.multiply(np.divide((2 * np.pi), period),np.divide(np.multiply(a_star, sin_i), np.sqrt((1 - np.square(e)))))  # AU/days
	K = np.multiply(K, 1731.48)  # km/s

	for i in range(n):  # Iterate over companions
		if calc[i]:
			for j in range(0, len(MJD)):  # Iterate over times
				# Find E
				M = 2 * np.pi * MJD[j] / period[i] - phase[i]
				prev_E = 0.0
				current_E = M

				while abs(current_E - prev_E) > 0.00001:
					prev_E = current_E
					current_E = M + e[i] * np.sin(prev_E)

				# Find true anomaly, f
				f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - e[i]) / (1 + e[i])))
				# Find predicted RV
				RV[i][j] = K[i] * (np.sin(arg_peri[i] + f) + e[i] * np.sin(arg_peri[i]))  # km/s

	return K, RV

def zero_point_model(prediction, a):
	# Alters the prediction by some zero point shift, a
	y = prediction + a
	return y

def zero_point_fit_parallel(experimental, errors, predicted):

	A = [0.] * len(predicted)
	for i in range(len(predicted)):
		[a], _ = scipy.optimize.curve_fit(zero_point_model, predicted[i], experimental, sigma=errors)
		A[i] = a
	return A

# Main Classes
class GUI(tk.Frame):
	# variables
	start_run = 0
	quit_run = 0

	# functions
	def __init__(self, parent):
		self.root = tk.Tk()
		self.parent = parent
		self.root.title('MOLUSC')
		self.create_widgets()

	def start(self):
		self.root.mainloop()

	def create_widgets(self):

		self.root.configure(bg='#ECECEC')

		# Create main panels
		Analysis_Options = tk.Frame(self.root, relief=tk.SUNKEN, borderwidth=4, bg='#ECECEC')
		Output = tk.Frame(self.root, width=200, height=200, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
		Stellar_Info = tk.Frame(self.root, width=200, height=200, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
		Advanced_Options = tk.Frame(self.root, width=200, height=200, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
		Buttons = tk.Frame(self.root, bg='#ECECEC', height=50)
		Display = tk.Frame(self.root, borderwidth=2, relief=tk.SUNKEN, bg='#ECECEC')

		# Place main panels
		Analysis_Options.grid(column=0, row=0, sticky='nsew')
		Stellar_Info.grid(column=0, row=1, sticky='ew')
		Output.grid(column=0, row=2,  sticky='ew')
		Advanced_Options.grid(column=1, row=0, columnspan=2, rowspan=2, sticky='nsew')
		Buttons.grid(column=2, row=2, sticky='se')
		Display.grid(column=3, row=0, rowspan=4, sticky='nsew')

		# Analysis Options
		analysis_label = tk.Label(Analysis_Options, text="Analysis Options", bg='#ECECEC')
		analysis_label.grid(column=0, row=0, columnspan=3)

		vcmd1 = (self.root.register(self.validate_resolution), '%P')

		#  toggle row
		checkbox_frame = tk.Frame(Analysis_Options, borderwidth=0, bg='#ECECEC')
		self.__ao_check = tk.BooleanVar()
		self.__rv_check = tk.BooleanVar()
		self.__ruwe_check = tk.BooleanVar()
		self.__gaia_check = tk.BooleanVar()

		ao_box = tk.Checkbutton(checkbox_frame, text="HRI", variable=self.__ao_check, command=self.activate_ao, bg='#ECECEC')
		rv_box = tk.Checkbutton(checkbox_frame, text="RV", variable=self.__rv_check, command=self.activate_rv, bg='#ECECEC')
		ruwe_box = tk.Checkbutton(checkbox_frame, text='RUWE', variable=self.__ruwe_check, bg='#ECECEC')
		gaia_box = tk.Checkbutton(checkbox_frame, text='Gaia', variable=self.__gaia_check, bg='#ECECEC')

		# sub frame for AO, necessary from multiple AO files
		self.ao_frame = tk.Frame(Analysis_Options, borderwidth=0)

		# AO variables and widgets go into the separate AO frame so that I can add multiple rows of them if needed
		self.ao_rows = 1
		self.ao_file_boxes = []
		self.ao_filter_menus = []
		# label column
		rv_label = tk.Label(Analysis_Options, text='RV File:', bg='#ECECEC')
		ao_label = tk.Label(self.ao_frame, text='Contrast File:', bg='#ECECEC')
		# entry column
		self.__ao_file = tk.StringVar()
		self.__rv_file = tk.StringVar()
		self.__ao_file_box = tk.Entry(self.ao_frame, textvariable=self.__ao_file, state=tk.DISABLED, width=25, disabledbackground='#ECECEC', highlightthickness=0, bd=1)
		self.ao_file_boxes.append(self.__ao_file_box) # adds the first file box to the list
		self.__rv_file_box = tk.Entry(Analysis_Options, textvariable=self.__rv_file, state=tk.DISABLED, width=25, disabledbackground='#ECECEC', highlightthickness=0, bd=1)
		# AO filter menu
		self.__filter_str = tk.StringVar()
		self.__filter_menu = tk.OptionMenu(self.ao_frame, self.__filter_str, 'Filter', 'J', 'H', 'K', 'G','Bp','Rp', 'R', 'I', 'L', 'LL', 'M')
		self.__filter_str.set('Filter')
		self.__filter_menu.config(state=tk.DISABLED, highlightthickness=0)
		self.ao_filter_menus.append(self.__filter_str)
		# AO add button
		self.__add_button = tk.Button(self.ao_frame, text='+', command=self.add_ao, state=tk.DISABLED, width=1)
		# RV resolution Box
		resolution_label = tk.Label(Analysis_Options, text=' Resolution:', bg='#ECECEC')
		self.__resolution = 0.
		self.__resolution_box = tk.Entry(Analysis_Options, validate='focusout', vcmd=vcmd1, state=tk.DISABLED, width=6, disabledbackground='#ECECEC', highlightthickness=0, bd=1)

		# Grid
		# checkbox row
		checkbox_frame.grid(column=0, row=1, columnspan=4, pady=10)
		ao_box.grid(column=0, row=0)
		rv_box.grid(column=1, row=0)
		ruwe_box.grid(column=3, row=0)
		gaia_box.grid(column=4, row=0)
		#  ao
		self.ao_frame.grid(column=0, columnspan=4, row=2)
		ao_label.grid(column=0, row=1, sticky='e')
		self.__ao_file_box.grid(column=1, row=1)
		self.__filter_menu.grid(column=2, row=1, padx=2)
		self.__add_button.grid(column=3, row=1)
		#  rv row
		rv_label.grid(column=0, row=3, sticky='w')
		self.__rv_file_box.grid(column=1, row=3, pady=15)
		resolution_label.grid(column=2, row=3, pady=15)
		self.__resolution_box.grid(column=3, row=3, pady=15)

		# Stellar Info/Basic Options
		stellar_label = tk.Label(Stellar_Info, text='Star Info', bg='#ECECEC')
		stellar_label.grid(column=0, row=0, columnspan=4)

		#  label columns
		number_label = tk.Label(Stellar_Info, text='Generated Companions*:', bg='#ECECEC')
		ra_label = tk.Label(Stellar_Info, text='RA (hms)*:', bg='#ECECEC')
		dec_label = tk.Label(Stellar_Info, text='DEC (dms)*:', bg='#ECECEC')
		mass_label = tk.Label(Stellar_Info, text='Star Mass (M_sun)*:', bg='#ECECEC')
		age_label = tk.Label(Stellar_Info, text='Star Age (Gyr):', bg='#ECECEC')
		added_jitter_label = tk.Label(Stellar_Info, text='Added Jitter (m/s):', bg='#ECECEC')
		rv_floor_label = tk.Label(Stellar_Info, text='RV Floor (m/s):', bg='#ECECEC')

		#  entry columns
		self.__ra_str = ''
		self.__dec_str = ''
		self.__mass = 0
		self.__age = 5
		self.__num_generated = 0
		self.__added_jitter = 20
		self.__rv_floor = 20
		vcmd1 = (self.root.register(self.validate_ra), '%P')
		vcmd2 = (self.root.register(self.validate_dec), '%P')
		vcmd3 = (self.root.register(self.validate_num), '%P')
		vcmd4 = (self.root.register(self.validate_mass), '%P')
		vcmd5 = (self.root.register(self.validate_age), '%P')
		vcmd6 = (self.root.register(self.validate_jitter), '%P')
		vcmd7 = (self.root.register(self.validate_floor), '%P')

		self.number_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd3, highlightthickness=0, bd=1)
		self.ra_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd1, highlightthickness=0, bd=1)
		self.dec_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd2, highlightthickness=0, bd=1)
		self.mass_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd4, highlightthickness=0, bd=1)
		self.age_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd5, highlightthickness=0, bd=1)
		self.age_box.insert(-1, '5')
		self.age_box.config(fg='gray')
		self.added_jitter_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd6, highlightthickness=0, bd=1)
		self.added_jitter_box.insert(-1, '20')
		self.added_jitter_box.config(fg='gray')
		self.rv_floor_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd7, highlightthickness=0, bd=1)
		self.rv_floor_box.insert(-1, '20')
		self.rv_floor_box.config(fg='gray')


		#  help button column
		age_help = tk.Button(Stellar_Info, text='?', command=self.help_age, height=1, width=2, highlightthickness=0, bd=3)
		ra_help = tk.Button(Stellar_Info, text='?', command=self.help_coord, height=1, width=2, highlightthickness=0, bd=3)
		dec_help = tk.Button(Stellar_Info, text='?', command=self.help_coord, height=1, width=2, highlightthickness=0, bd=3)
		added_jitter_help = tk.Button(Stellar_Info, text='?', command=self.help_added_jitter, height=1, width=2, highlightthickness=0, bd=3)
		rv_floor_help = tk.Button(Stellar_Info, text='?', command=self.help_floor, height=1, width=2, highlightthickness=0, bd=3)

		# Grid
		Stellar_Info.grid_columnconfigure(0, weight=1)
		Stellar_Info.grid_columnconfigure(2, weight=1)
		Stellar_Info.grid_rowconfigure(0, weight=1)
		Stellar_Info.grid_rowconfigure(7, weight=1)
		# labels
		number_label.grid(column=0, row=1, sticky='e')
		ra_label.grid(column=0, row=2, sticky='e')
		dec_label.grid(column=0, row=3, sticky='e')
		mass_label.grid(column=0, row=4, sticky='e')
		age_label.grid(column=0, row=5, sticky='e')
		added_jitter_label.grid(column=0, row=6, sticky='e')
		rv_floor_label.grid(column=0, row=7, sticky='e')
		# boxes
		self.number_box.grid(column=1, row=1)
		self.ra_box.grid(column=1, row=2)
		self.dec_box.grid(column=1, row=3)
		self.mass_box.grid(column=1, row=4)
		self.age_box.grid(column=1, row=5)
		self.added_jitter_box.grid(column=1, row=6)
		self.rv_floor_box.grid(column=1, row=7)
		# buttons
		ra_help.grid(column=2, row=2, sticky='w')
		dec_help.grid(column=2, row=3, sticky='w')
		age_help.grid(column=2, row=5, sticky='w')
		added_jitter_help.grid(column=2, row=6, sticky='w')
		rv_floor_help.grid(column=2, row=7, sticky='w')

		# Output
		self.__prefix = tk.StringVar()
		self.__extra = tk.IntVar()
		self.__all_out = tk.IntVar()
		output_label = tk.Label(Output, text='Output', bg='#ECECEC')
		prefix_label = tk.Label(Output, text='File Prefix*:', bg='#ECECEC')
		self.prefix_box = tk.Entry(Output, textvariable=self.__prefix, width=20, highlightthickness=0, bd=1)
		extra_box = tk.Checkbutton(Output, text='Extra Output', variable=self.__extra, pady=9, bg='#ECECEC')
		all_out_box = tk.Checkbutton(Output, text='Write Out All', variable=self.__all_out, bg='#ECECEC')

		Output.grid_columnconfigure(0, weight=1)
		Output.grid_columnconfigure(2, weight=1)
		Output.grid_rowconfigure(0, weight=1)
		output_label.grid(column=0, row=0, columnspan=3)
		prefix_label.grid(column=0, row=1, sticky='e')
		self.prefix_box.grid(column=1, row=1, columnspan=2, sticky='w')
		extra_box.grid(column=0, row=3, sticky='e')
		all_out_box.grid(column=2, row=3)

		# Advanced Options
		self.__P_fixed = None
		self.__P_min = None
		self.__P_max = None
		self.__i_fixed = None
		self.__i_min = None
		self.__i_max = None
		self.__e_fixed = None
		self.__e_min = None
		self.__e_max = None
		self.__arg_peri_fixed = None
		self.__arg_peri_min = None
		self.__arg_peri_max = None
		self.__m_fixed = None
		self.__m_min = None
		self.__m_max = None
		self.__a_fixed = None
		self.__a_min = None
		self.__a_max = None
		self.__phi_fixed = None
		self.__phi_min = None
		self.__phi_max = None
		self.__pd_mu = 5.03
		self.__pd_sig = 2.28
		self.__q_exp = 0.0
		self.__gaia_limit = tk.IntVar()
		self.__gaia_limit.set(18)

		advanced_label = tk.Label(Advanced_Options, text='Advanced Options', bg='#ECECEC')

		# Frames
		grid_frame = tk.Frame(Advanced_Options, borderwidth=0, bg='#ECECEC')
		dist_frame = tk.Frame(Advanced_Options, borderwidth=0, bg='#ECECEC')

		#  Grid Labels
		fixed_label = tk.Label(grid_frame, text='fixed', bg='#ECECEC', width=7)
		min_label = tk.Label(grid_frame, text='min', bg='#ECECEC', width=7)
		max_label = tk.Label(grid_frame, text='max', bg='#ECECEC', width=7)
		period_label = tk.Label(grid_frame, text='Period (days)', bg='#ECECEC')
		inc_label = tk.Label(grid_frame, text='cos(i)', bg='#ECECEC')
		ecc_label = tk.Label(grid_frame, text='Eccentricity', bg='#ECECEC')
		arg_peri_label = tk.Label(grid_frame, text='Arg Periapsis (radians)', bg='#ECECEC')
		mass_ratio_label = tk.Label(grid_frame, text='Mass Ratio (M2/M1)', bg='#ECECEC')
		a_label = tk.Label(grid_frame, text='Semi-Major Axis (AU)', bg='#ECECEC')
		phase_label = tk.Label(grid_frame, text='Phase (radians)', bg='#ECECEC')
		# Split Label
		split_label = tk.Label(Advanced_Options, text='- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -', justify='center', bg='#ECECEC')
		# Distribution Labels
		period_dist_label = tk.Label(dist_frame, text='Period Distribution', justify='center', bg='#ECECEC')
		mu_label = tk.Label(dist_frame, text=u'\u03BC:', bg='#ECECEC')
		sig_label = tk.Label(dist_frame, text=u'\u03C3:', bg='#ECECEC')
		mass_exp_label = tk.Label(dist_frame, text='Mass Ratio Distribution', justify='center', bg='#ECECEC')
		exp_label = tk.Label(dist_frame, text=u'\u03B3:', bg='#ECECEC')
		#  Validation Command registration
		vcmd_p = (self.root.register(self.validate_P), '%P', '%W')
		vcmd_i = (self.root.register(self.validate_i), '%P', '%W')
		vcmd_e = (self.root.register(self.validate_e), '%P', '%W')
		vcmd_w = (self.root.register(self.validate_arg_peri), '%P', '%W')
		vcmd_m = (self.root.register(self.validate_m), '%P', '%W')
		vcmd_a = (self.root.register(self.validate_a), '%P', '%W')
		vcmd_phi = (self.root.register(self.validate_phi), '%P', '%W')
		vcmd_pd_mu = (self.root.register(self.validate_pd_mu), '%P')
		vcmd_pd_sig = (self.root.register(self.validate_pd_sig), '%P')
		vcmd_mass_exp = (self.root.register(self.validate_mass_exp), '%P')
		#  Grid Entries
		entry_width = 7
		self.P1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_p, highlightthickness=0, bd=1)
		self.P2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_p, highlightthickness=0, bd=1)
		self.P3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_p, highlightthickness=0, bd=1)
		self.i1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_i, highlightthickness=0, bd=1)
		self.i2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_i, highlightthickness=0, bd=1)
		self.i3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_i, highlightthickness=0, bd=1)
		self.e1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_e, highlightthickness=0, bd=1)
		self.e2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_e, highlightthickness=0, bd=1)
		self.e3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_e, highlightthickness=0, bd=1)
		self.w1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_w, highlightthickness=0, bd=1)
		self.w2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_w, highlightthickness=0, bd=1)
		self.w3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_w, highlightthickness=0, bd=1)
		self.m1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_m, highlightthickness=0, bd=1)
		self.m2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_m, highlightthickness=0, bd=1)
		self.m3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_m, highlightthickness=0, bd=1)
		self.a1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_a, highlightthickness=0, bd=1)
		self.a2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_a, highlightthickness=0, bd=1)
		self.a3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_a, highlightthickness=0, bd=1)
		self.phi1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_phi, highlightthickness=0, bd=1)
		self.phi2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_phi, highlightthickness=0, bd=1)
		self.phi3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_phi, highlightthickness=0, bd=1)
		# Distribution Entries
		self.PD_mu_box = tk.Entry(dist_frame, width=entry_width, validate='focusout', vcmd=vcmd_pd_mu, highlightthickness=0, bd=1)
		self.PD_mu_box.insert(-1, '5.03')
		self.PD_mu_box.config(fg='gray')
		self.PD_sig_box = tk.Entry(dist_frame, width=entry_width, validate='focusout', vcmd=vcmd_pd_sig, highlightthickness=0, bd=1)
		self.PD_sig_box.insert(-1, '2.28')
		self.PD_sig_box.config(fg='gray')
		self.mass_exp_box = tk.Entry(dist_frame, width=entry_width, validate='focusout', vcmd=vcmd_mass_exp, highlightthickness=0, bd=1)
		self.mass_exp_box.insert(-1, '0.0')
		self.mass_exp_box.config(fg='gray')
		# Gaia Limit
		gaia_label = tk.Label(dist_frame, text='Gaia Completeness Limit:', bg='#ECECEC')
		gaia_18 = tk.Radiobutton(dist_frame, text="18th", variable=self.__gaia_limit, value=18)
		gaia_20 = tk.Radiobutton(dist_frame, text="20th", variable=self.__gaia_limit, value=20)

		#  Help buttons
		help_button = tk.Button(Advanced_Options, text='?', command=self.help_advanced, bg='#ECECEC', width=3, highlightthickness=0, relief=tk.RAISED, bd=3)
		period_help = tk.Button(dist_frame, text='?', command=self.help_period, height=1, width=2, highlightthickness=0, bd=3)
		mass_help = tk.Button(dist_frame, text='?', command=self.help_mass, height=1, width=2, highlightthickness=0, bd=3)

		#  Full Grid
		advanced_label.grid(column=0, columnspan=3, row=0)
		help_button.grid(column=2, row=0)
		grid_frame.grid(column=0, row=1, columnspan=3)
		split_label.grid(column=0, row=2, columnspan=3)
		dist_frame.grid(column=0, row=3, columnspan=3)
		#  limit grid
		fixed_label.grid(column=2, row=1)
		min_label.grid(column=3, row=1)
		max_label.grid(column=4, row=1)
		period_label.grid(column=0, row=2)
		inc_label.grid(column=0, row=3)
		ecc_label.grid(column=0, row=4)
		arg_peri_label.grid(column=0, row=5)
		mass_ratio_label.grid(column=0, row=6)
		a_label.grid(column=0, row=7)
		phase_label.grid(column=0, row=8)
		self.P1_box.grid(column=2, row=2)
		self.P2_box.grid(column=3, row=2)
		self.P3_box.grid(column=4, row=2)
		self.i1_box.grid(column=2, row=3)
		self.i2_box.grid(column=3, row=3)
		self.i3_box.grid(column=4, row=3)
		self.e1_box.grid(column=2, row=4)
		self.e2_box.grid(column=3, row=4)
		self.e3_box.grid(column=4, row=4)
		self.w1_box.grid(column=2, row=5)
		self.w2_box.grid(column=3, row=5)
		self.w3_box.grid(column=4, row=5)
		self.m1_box.grid(column=2, row=6)
		self.m2_box.grid(column=3, row=6)
		self.m3_box.grid(column=4, row=6)
		self.a1_box.grid(column=2, row=7)
		self.a2_box.grid(column=3, row=7)
		self.a3_box.grid(column=4, row=7)
		self.phi1_box.grid(column=2, row=8)
		self.phi2_box.grid(column=3, row=8)
		self.phi3_box.grid(column=4, row=8)
		# distribution
		period_dist_label.grid(column=0, row=0)
		mass_exp_label.grid(column=0, row=1)
		mu_label.grid(column=1, row=0)
		self.PD_mu_box.grid(column=2, row=0)
		sig_label.grid(column=3, row=0)
		self.PD_sig_box.grid(column=4, row=0)
		period_help.grid(column=5, row=0)
		exp_label.grid(column=1,row=1)
		self.mass_exp_box.grid(column=2, row=1)
		mass_help.grid(column=5, row=1)
		# gaia
		gaia_label.grid(column=0, row=2, columnspan=2)
		gaia_18.grid(column=2, row=2)
		gaia_20.grid(column=4, row=2)

		# Buttons
		self.run_button = tk.Button(Buttons, text='Run', command=self.run_code, width=6, height=2, highlightthickness=0, bd=3, bg='#ECECEC')
		close_button = tk.Button(Buttons, text='Close', command=self.close, width=8, height=2, highlightthickness=0, bd=3, bg='#ECECEC')

		self.run_button.grid(column=0, row=0)
		close_button.grid(column=1, row=0, padx=4)

		# Display
		self.display_box_text = tk.StringVar()
		self.status_box_text = tk.StringVar()
		self.status_box_text.set('Status: Waiting for Input')

		self.display_box = tk.Text(Display, height=29, width=40, background='white')
		self.status_box = tk.Label(Display, textvariable=self.status_box_text, width=33, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
		scrollbar = tk.Scrollbar(Display, orient=tk.VERTICAL, command=self.display_box.yview, width=15)

		self. display_box.configure(yscrollcommand=scrollbar.set)

		self.status_box.grid(column=0, row=0, columnspan=2)
		self.display_box.grid(column=0, row=1)
		scrollbar.grid(column=1, row=1)

		# end create_widgets

	def activate_ao(self):
		if self.__ao_check.get() == 1:  # whenever checked
			self.__ao_file_box.config(state=tk.NORMAL)
			self.__filter_menu.config(state=tk.NORMAL)
			self.__add_button.config(state=tk.NORMAL)
		elif self.__ao_check.get() == 0:  # whenever unchecked
			self.__ao_file_box.config(state=tk.DISABLED)
			self.__filter_menu.config(state=tk.DISABLED)
			self.__add_button.config(state=tk.DISABLED)

	def activate_rv(self):
		if self.__rv_check.get() == 1:  # whenever checked
			self.__rv_file_box.config(state=tk.NORMAL)
			self.__resolution_box.config(state=tk.NORMAL)
		elif self.__rv_check.get() == 0:  # whenever unchecked
			self.__rv_file_box.config(state=tk.DISABLED)
			self.__resolution_box.config(state=tk.DISABLED)

	def add_ao(self):
		# Adds another row with entry & filter for mutliple AO files
		self.ao_rows += 1
		ao_file = tk.StringVar()
		filter_str = tk.StringVar()
		# Create Widgits
		label = tk.Label(self.ao_frame, text=('Contrast File ' + str(self.ao_rows) +  ':'), bg='#ECECEC')
		file_box = tk.Entry(self.ao_frame, textvariable=ao_file, width=25, highlightthickness=0, bd=1)
		filter_menu = tk.OptionMenu(self.ao_frame, self.__filter_str, 'Filter', 'J', 'H', 'K', 'G','Bp','Rp', 'R', 'I', 'L', 'LL', 'M')
		filter_str.set('Filter')
		# Place Widgets
		label.grid(column=0, row=self.ao_rows)
		file_box.grid(column=1, row=self.ao_rows)
		filter_menu.grid(column=2, row=self.ao_rows)
		# Add Widgets to list
		self.ao_file_boxes.append(file_box)
		self.ao_filter_menus.append(filter_str)

	def check_completeness(self):
		# When the user presses run need to check if all required fields are filled.
		# file names, coordinates, mass, prefix, num generated all need to be filled

		allow_run = True
		message = 'The following problems were detected:\n'

		if not self.__ao_check.get() and not self.__rv_check.get() and not self.__ruwe_check.get() and not self.__gaia_check.get():
			message = message + '- No analysis is selected. Please select at least one type of analysis.\n'
			allow_run = False
		# Analysis Options
		if self.__ao_check.get():
			# The AO test requires an input file and a filter selection
			if not  self.__ao_file.get():
				self.__ao_file_box.config(bg='lightcoral')
				allow_run = False
				message = message +'- No file containing AO data has been given.\n'
			if self.__filter_str.get() == 'Filter':
				allow_run = False
				message = message + '- No filter is selected. Please choose the filter that the AO data is in.\n'
		if self.__rv_check.get():
			# The RV test requires an input file and a resolution
			if not self.__rv_file.get():
				self.__rv_file_box.config(bg='lightcoral')
				allow_run = False
				message = message +'- No file containing RV data has been given.\n'
			if self.__resolution == 0.:
				self.__resolution_box.config(bg='lightcoral')
				allow_run = False
				message = message +'- No RV Resolution provided.\n'
		# Coordinates
		if not self.__ra_str or self.__ra_str == '00h00m00.00s' or not self.validate_ra(self.__ra_str):
			self.ra_box.config(bg='lightcoral')
			allow_run = False
		if not self.__dec_str or self.__dec_str == '00d00m00.0s' or not self.validate_dec(self.__dec_str):
			self.dec_box.config(bg='lightcoral')
			allow_run=False
		# Mass, Age, Number Generated
		if self.__mass == 0 or not self.validate_mass(self.__mass):
			self.mass_box.config(bg='lightcoral')
			allow_run=False
		if self.__num_generated == 0 or not self.validate_num(self.__num_generated):
			self.number_box.config(bg='lightcoral')
			allow_run = False
		# Prefix
		if self.__prefix.get() == '':
			self.prefix_box.config(bg='lightcoral')
			allow_run = False

		if not allow_run:
			message = message + '\nA required field is missing or invalid. Please check any fields marked in red and try again.\n'
			tk.messagebox.showinfo('Run Message', message)
			return -1
		else:
			return 0

	def gui_print(self, new_text):
		# need to change the text in label to include this message on the end
		text_width = self.display_box.cget('width')
		if new_text == 'clc':
			self.display_box.delete('1.0', 'end')
		else:
			wrap_new_text = wrap(new_text, text_width)
			new_message = ''
			for i in range(0, len(wrap_new_text)):
				new_message += wrap_new_text[i] + '\n'
			self.display_box.insert('end', new_message)
		self.root.update()
		return

	def update_status(self, new_status):
		self.status_box_text.set(('Status: ' + new_status))
		self.gui_print(('\nStatus: '+ new_status))
		print(('\nStatus: '+ new_status))
		self.root.update()

	# Analysis Options validation functions

	def validate_resolution(self, new_text):
		# Check that number generated is an integer greater than zero
		if not new_text:
			# box cleared
			self.__resolution = 0.
			self.__resolution_box.config(bg='white')
			return True
		try:
			n = float(new_text)
			if n > 0.:
				self.__resolution = n
				self.__resolution_box.config(bg='white')
				return True
			else:
				self.__resolution_box.config(bg='lightcoral')
				return False
		except:
			self.__resolution_box.config(bg='lightcoral')
			return False

	# Stellar Info validation functions

	def validate_ra(self, new_text):
		# Check that RA makes sense
		if not new_text:
			# box cleared
			self.__ra_str= '00h00m00.00s'
			self.ra_box.config(bg='white')
			return True
		try:
			hours = new_text[0:2]
			minutes = new_text[3:5]
			seconds = new_text[6:-2]
			if new_text[2] == 'h' and new_text[5] == 'm' and new_text[-1] == 's':
				hours = int(hours)
				minutes = int(minutes)
				seconds = float(seconds)
				if hours < 24 and minutes < 60 and seconds < 60:
					self.__ra_str = new_text
					self.ra_box.config(bg='white')
					return True
				else:
					self.__ra_str = new_text
					self.ra_box.config(bg='lightcoral')
					return False
			else:
				self.__ra_str = new_text
				self.ra_box.config(bg='lightcoral')
				return False
		except:
			self.__ra_str = new_text
			self.ra_box.config(bg='lightcoral')
			return False

	def validate_dec(self, new_text):
		# Check that DEC makes sense
		if not new_text:
			# box cleared
			self.__dec_str = '00d00m00.0s'
			self.dec_box.config(bg='white')
			return True
		try:
			sign = new_text[0]
			degrees = new_text[1:3]
			minutes = new_text[4:6]
			seconds = new_text[7:-2]
			if new_text[3] == 'd' and new_text[6] == 'm' and new_text[-1] == 's':
				degrees = int(degrees)
				minutes = int(minutes)
				seconds = float(seconds)
				if degrees <= 90 and minutes < 60 and seconds < 60:
					if sign == '+' or sign == '-':
						self.__dec_str = new_text
						self.dec_box.config(bg='white')
						return True
					else:
						self.dec_box.config(bg='lightcoral')
						return False
				else:
					self.__dec_str = new_text
					self.dec_box.config(bg='lightcoral')
					return False
			else:
				self.__dec_str = new_text
				self.dec_box.config(bg='lightcoral')
				return False
		except:
			self.__dec_str = new_text
			self.dec_box.config(bg='lightcoral')
			return False

	def validate_num(self, new_text):
		# Check that number generated is an integer greater than zero
		if not new_text:
			# box cleared
			self.__num_generated = 0
			self.number_box.config(bg='white')
			return True
		try:
			n = int(new_text)
			if n > 0:
				self.__num_generated = n
				self.number_box.config(bg='white')
				return True
			else:
				self.number_box.config(bg='lightcoral')
				return False
		except:
			self.number_box.config(bg='lightcoral')
			return False

	def validate_mass(self, new_text):
		# Check that mass is a number greater than zero
		if not new_text:
			# box cleared
			self.__mass = 0
			self.mass_box.config(bg='white')
			return True
		try:
			m = float(new_text)
			if m > 0:
				self.__mass = m
				self.mass_box.config(bg='white')
				return True
			else:
				self.mass_box.config(bg='lightcoral')
				return False
		except:
			self.mass_box.config(bg='lightcoral')
			return False

	def validate_age(self, new_text):
		# Check that age is a number greater than zero
		if not new_text:
			# box cleared
			self.__age = 5
			self.age_box.insert(-1, '5')
			self.age_box.config(bg='white', fg='gray')
			return True
		try:
			a = float(new_text)
			if a > 0:
				self.__age = a
				self.age_box.config(bg='white', fg='black')
				return True
			else:
				self.age_box.config(bg='lightcoral', fg='black')
				return False
		except:
			self.age_box.config(bg='lightcoral', fg='black')
			return False

	def validate_jitter(self, new_text):
		# Check that jitter is a number greater than zero
		if not new_text:
			# box cleared
			self.__added_jitter = 20
			self.added_jitter_box.config(bg='white', fg='gray')
			self.added_jitter_box.insert(-1, '20')
			return True
		try:
			j = float(new_text)
			if j >= 0:
				self.__added_jitter = j
				self.added_jitter_box.config(bg='white', fg='black')
				return True
			else:
				self.added_jitter_box.config(bg='lightcoral', fg='black')
				return False
		except:
			self.added_jitter_box.config(bg='lightcoral', fg='black')
			return False

	def validate_floor(self, new_text):
		# Check that jitter is a number greater than zero
		if not new_text:
			# box cleared
			self.__rv_floor = 20
			self.rv_floor_box.config(bg='white', fg='gray')
			self.rv_floor_box.insert(-1, '20')
			return True
		try:
			j = float(new_text)
			if j >= 0:
				self.__rv_floor = j
				self.rv_floor_box.config(bg='white', fg='black')
				return True
			else:
				self.rv_floor_box.config(bg='lightcoral', fg='black')
				return False
		except:
			self.rv_floor_box.config(bg='lightcoral', fg='black')
			return False

	# Advanced Options validation functions
	def validate_P(self, new_text, box_name):
		# Period must be greater than, or equal to 0.1, if accepted, the other period boxes should be grayed out
		if not new_text:
			# Text box cleared
			if box_name.endswith('entry'):
				self.__P_fixed = None
				self.P1_box.config(bg='white')
				self.P2_box.config(state='normal')
				self.P3_box.config(state='normal')
			elif box_name.endswith('entry2'):
				self.__P_min = None
				self.P2_box.config(bg='white')
				if self.__P_max is None:
					self.P1_box.config(state='normal')
			elif box_name.endswith('entry3'):
				self.__P_max = None
				self.P3_box.config(bg='white')
				if self.__P_min is None:
					self.P1_box.config(state='normal')
			return True
		try:
			p = float(new_text)
			if p >= 0.1:  # If valid set the desired limit to the given value, and disable other boxes as needed
				if box_name.endswith('entry'):
					self.__P_fixed = p
					self.P1_box.config(bg='white')
					self.P2_box.config(state='disabled')
					self.P3_box.config(state='disabled')
				elif box_name.endswith('entry2'):
					self.__P_min = p
					self.P2_box.config(bg='white')
					self.P1_box.config(state='disabled')
				elif box_name.endswith('entry3'):
					self.__P_max = p
					self.P3_box.config(bg='white')
					self.P1_box.config(state='disabled')
				return True
			else:  # If invalid set the desired limit to its default value, color and disable other boxes as needed
				if box_name.endswith('entry'):
					self.P1_box.config(bg='lightcoral')
					self.P2_box.config(state='disabled')
					self.P3_box.config(state='disabled')
				elif box_name.endswith('entry2'):
					self.__P_min = 0.1
					self.P2_box.config(bg='lightcoral')
					self.P1_box.config(state='disabled')
				elif box_name.endswith('entry3'):
					self.__P_max = float('inf')
					self.P3_box.config(bg='lightcoral')
					self.P1_box.config(state='disabled')
				return False
		except:  # If very invalid set the desired limit to its default value, color and disable other boxes as needed
			if box_name.endswith('entry'):
				self.P1_box.config(bg='lightcoral')
				self.P2_box.config(state='disabled')
				self.P3_box.config(state='disabled')
			elif box_name.endswith('entry2'):
				self.__P_min = 0.1
				self.P2_box.config(bg='lightcoral')
				self.P1_box.config(state='disabled')
			elif box_name.endswith('entry3'):
				self.__P_max = float('inf')
				self.P3_box.config(bg='lightcoral')
				self.P1_box.config(state='disabled')
			return False

	def validate_i(self, new_text, box_name):
		# sin(i) must be greater than -1 and less than 1, if accepted, the other boxes should be grayed out
		# a limit value of 'transit' can also be accepted
		if not new_text:
			# box cleared
			if box_name.endswith('entry4'):
				self.__i_fixed = None
				self.i1_box.config(bg='white')
				self.i2_box.config(state='normal')
				self.i3_box.config(state='normal')
			elif box_name.endswith('entry5'):
				self.__i_min = None
				self.i2_box.config(bg='white')
				if self.__i_max is None:
					self.i1_box.config(state='normal')
			elif box_name.endswith('entry6'):
				self.__i_max = None
				self.i3_box.config(bg='white')
				if self.__i_min is None:
					self.i1_box.config(state='normal')
			return True
		try:
			i = float(new_text)
			if -1 <= i <= 1:
				if box_name.endswith('entry4'):
					self.__i_fixed = i
					self.i1_box.config(bg='white')
					self.i2_box.config(state='disabled')
					self.i3_box.config(state='disabled')
				elif box_name.endswith('entry5'):
					self.__i_min = i
					self.i2_box.config(bg='white')
					self.i1_box.config(state='disabled')
				elif box_name.endswith('entry6'):
					self.__i_max = i
					self.i3_box.config(bg='white')
					self.i1_box.config(state='disabled')
				return True
			else:
				if box_name.endswith('entry4'):
					self.i1_box.config(bg='lightcoral')
					self.i2_box.config(state='disabled')
					self.i3_box.config(state='disabled')
				elif box_name.endswith('entry5'):
					self.__i_min = -1.
					self.i2_box.config(bg='lightcoral')
					self.i1_box.config(state='disabled')
				elif box_name.endswith('entry6'):
					self.__i_max = 1.
					self.i3_box.config(bg='lightcoral')
					self.i1_box.config(state='disabled')
				return False
		except:
			if new_text == 'transit':
				if box_name.endswith('entry4'):
					self.__i_fixed = new_text
					self.i1_box.config(bg='white')
					self.i2_box.config(state='disabled')
					self.i3_box.config(state='disabled')
				elif box_name.endswith('entry5'):
					self.__i_min = new_text
					self.i2_box.config(bg='white')
					self.i1_box.config(state='disabled')
				elif box_name.endswith('entry6'):
					self.__i_max = new_text
					self.i3_box.config(bg='white')
					self.i1_box.config(state='disabled')
				return True
			else:
				if box_name.endswith('entry4'):
					self.i1_box.config(bg='lightcoral')
					self.i2_box.config(state='disabled')
					self.i3_box.config(state='disabled')
				elif box_name.endswith('entry5'):
					self.__i_min = -1
					self.i2_box.config(bg='lightcoral')
					self.i1_box.config(state='disabled')
				elif box_name.endswith('entry6'):
					self.__i_max = 1
					self.i3_box.config(bg='lightcoral')
					self.i1_box.config(state='disabled')
				return False

	def validate_e(self, new_text, box_name):
		# e must be greater than 0 and less than 1, if accepted, the other boxes should be grayed out
		if not new_text:
			# Text box cleared
			if box_name.endswith('entry7'):
				self.__e_fixed = None
				self.e1_box.config(bg='white')
				self.e2_box.config(state='normal')
				self.e3_box.config(state='normal')
			elif box_name.endswith('entry8'):
				self.__e_min = None
				self.e2_box.config(bg='white')
				if self.__e_max is None:
					self.e1_box.config(state='normal')
			elif box_name.endswith('entry9'):
				self.__e_max = None
				self.e3_box.config(bg='white')
				if self.__e_min is None:
					self.e1_box.config(state='normal')
			return True
		try:
			e = float(new_text)
			if 0 <= e <= 1:  # If valid set the desired limit to the given value, and disable other boxes as needed
				if box_name.endswith('entry7'):
					self.__e_fixed = e
					self.e1_box.config(bg='white')
					self.e2_box.config(state='disabled')
					self.e3_box.config(state='disabled')
				elif box_name.endswith('entry8'):
					self.__e_min = e
					self.e2_box.config(bg='white')
					self.e1_box.config(state='disabled')
				elif box_name.endswith('entry9'):
					self.__e_max = e
					self.e3_box.config(bg='white')
					self.e1_box.config(state='disabled')
				return True
			else:  # If invalid set the desired limit to its default value, color and disable other boxes as needed
				if box_name.endswith('entry7'):
					self.e1_box.config(bg='lightcoral')
					self.e2_box.config(state='disabled')
					self.e3_box.config(state='disabled')
				elif box_name.endswith('entry8'):
					self.__e_min = 0.
					self.e2_box.config(bg='lightcoral')
					self.e1_box.config(state='disabled')
				elif box_name.endswith('entry9'):
					self.__e_max = 1.
					self.e3_box.config(bg='lightcoral')
					self.e1_box.config(state='disabled')
				return False
		except: # If invalid set the desired limit to its default value, color and disable other boxes as needed
			if box_name.endswith('entry7'):
				self.e1_box.config(bg='lightcoral')
				self.e2_box.config(state='disabled')
				self.e3_box.config(state='disabled')
			elif box_name.endswith('entry8'):
				self.__e_min = 0.
				self.e2_box.config(bg='lightcoral')
				self.e1_box.config(state='disabled')
			elif box_name.endswith('entry9'):
				self.__e_max = 1.
				self.e3_box.config(bg='lightcoral')
				self.e1_box.config(state='disabled')
			return False

	def validate_arg_peri(self, new_text, box_name):
		# arg peri must be greater than 0 and less than pi, if accepted, the other boxes should be grayed out
		if not new_text:
			# box cleared
			if box_name.endswith('entry10'):
				self.__arg_peri_fixed = None
				self.w1_box.config(bg='white')
				self.w2_box.config(state='normal')
				self.w3_box.config(state='normal')
			elif box_name.endswith('entry11'):
				self.__arg_peri_min = None
				self.w2_box.config(bg='white')
				if self.__arg_peri_max is None:
					self.w1_box.config(state='normal')
			elif box_name.endswith('entry12'):
				self.__arg_peri_max = None
				self.w3_box.config(bg='white')
				if self.__arg_peri_min is None:
					self.w1_box.config(state='normal')
			return True
		try:
			w = float(new_text)
			if 0 <= w <= np.pi:
				if box_name.endswith('entry10'):
					self.__arg_peri_fixed = w
					self.w1_box.config(bg='white')
					self.w2_box.config(state='disabled')
					self.w3_box.config(state='disabled')
				elif box_name.endswith('entry11'):
					self.__arg_peri_min = w
					self.w2_box.config(bg='white')
					self.w1_box.config(state='disabled')
				elif box_name.endswith('entry12'):
					self.__arg_peri_max = w
					self.w3_box.config(bg='white')
					self.w1_box.config(state='disabled')
				return True
			else:
				if box_name.endswith('entry10'):
					self.w1_box.config(bg='lightcoral')
					self.w2_box.config(state='disabled')
					self.w3_box.config(state='disabled')
				elif box_name.endswith('entry11'):
					self.__arg_peri_min = 0.
					self.w2_box.config(bg='lightcoral')
					self.w1_box.config(state='disabled')
				elif box_name.endswith('entry12'):
					self.__arg_peri_max = np.pi
					self.w3_box.config(bg='lightcoral')
					self.w1_box.config(state='disabled')
				return False
		except:
			if box_name.endswith('entry10'):
				self.w1_box.config(bg='lightcoral')
				self.w2_box.config(state='disabled')
				self.w3_box.config(state='disabled')
			elif box_name.endswith('entry11'):
				self.__arg_peri_min = 0.
				self.w2_box.config(bg='lightcoral')
				self.w1_box.config(state='disabled')
			elif box_name.endswith('entry12'):
				self.__arg_peri_max = np.pi
				self.w3_box.config(bg='lightcoral')
				self.w1_box.config(state='disabled')
			return False

	def validate_m(self, new_text, box_name):
		# m must be greater than 0 and less than 1, if accepted, the other boxes should be grayed out
		if not new_text:
			# box cleared
			if box_name.endswith('entry13'):
				self.__m_fixed = None
				self.m1_box.config(bg='white')
				self.m2_box.config(state='normal')
				self.m3_box.config(state='normal')
			elif box_name.endswith('entry14'):
				self.__m_min = None
				self.m2_box.config(bg='white')
				if self.__m_max is None:
					self.m1_box.config(state='normal')
			elif box_name.endswith('entry15'):
				self.__m_max = None
				self.m3_box.config(bg='white')
				if self.__m_min is None:
					self.m1_box.config(state='normal')
			return True
		try:
			m = float(new_text)
			if 0 <= m <= 1:
				if box_name.endswith('entry13'):
					self.__m_fixed = m
					self.m1_box.config(bg='white')
					self.m2_box.config(state='disabled')
					self.m3_box.config(state='disabled')
				elif box_name.endswith('entry14'):
					self.__m_min = m
					self.m2_box.config(bg='white')
					self.m1_box.config(state='disabled')
				elif box_name.endswith('entry15'):
					self.__m_max = m
					self.m3_box.config(bg='white')
					self.m1_box.config(state='disabled')
				return True
			else:
				if box_name.endswith('entry13'):
					self.m1_box.config(bg='lightcoral')
					self.m2_box.config(state='disabled')
					self.m3_box.config(state='disabled')
				elif box_name.endswith('entry14'):
					self.__m_min = 0.
					self.m2_box.config(bg='lightcoral')
					self.m1_box.config(state='disabled')
				elif box_name.endswith('entry15'):
					self.__m_max = 1.
					self.m3_box.config(bg='lightcoral')
					self.m1_box.config(state='disabled')
				return False
		except:
			if box_name.endswith('entry13'):
				self.m1_box.config(bg='lightcoral')
				self.m2_box.config(state='disabled')
				self.m3_box.config(state='disabled')
			elif box_name.endswith('entry14'):
				self.__m_min = 0.
				self.m2_box.config(bg='lightcoral')
				self.m1_box.config(state='disabled')
			elif box_name.endswith('entry15'):
				self.__m_max = 1.
				self.m3_box.config(bg='lightcoral')
				self.m1_box.config(state='disabled')
			return False

	def validate_a(self, new_text, box_name):
		# a must be greater than 0, if accepted, the other boxes should be grayed out
		if not new_text:
			# box cleared
			if box_name.endswith('entry16'):
				self.__a_fixed = None
				self.a1_box.config(bg='white')
				self.a2_box.config(state='normal')
				self.a3_box.config(state='normal')
			elif box_name.endswith('entry17'):
				self.__a_min = None
				self.a2_box.config(bg='white')
				if self.__a_max is None:
					self.a1_box.config(state='normal')
			elif box_name.endswith('entry18'):
				self.__a_max = None
				self.a3_box.config(bg='white')
				if self.__a_min is None:
					self.a1_box.config(state='normal')
			return True
		try:
			a = float(new_text)
			if a >= 0:
				if box_name.endswith('entry16'):
					self.__a_fixed = a
					self.a1_box.config(bg='white')
					self.a2_box.config(state='disabled')
					self.a3_box.config(state='disabled')
				elif box_name.endswith('entry17'):
					self.__a_min = a
					self.a2_box.config(bg='white')
					self.a1_box.config(state='disabled')
				elif box_name.endswith('entry18'):
					self.__a_max = a
					self.a3_box.config(bg='white')
					self.a1_box.config(state='disabled')
				return True
			else:
				if box_name.endswith('entry16'):
					self.a1_box.config(bg='lightcoral')
					self.a2_box.config(state='disabled')
					self.a3_box.config(state='disabled')
				elif box_name.endswith('entry17'):
					self.__a_min = 0.
					self.a2_box.config(bg='lightcoral')
					self.a1_box.config(state='disabled')
				elif box_name.endswith('entry18'):
					self.__a_max = float('inf')
					self.a3_box.config(bg='lightcoral')
					self.a1_box.config(state='disabled')
				return False
		except:
			if box_name.endswith('entry16'):
				self.a1_box.config(bg='lightcoral')
				self.a2_box.config(state='disabled')
				self.a3_box.config(state='disabled')
			elif box_name.endswith('entry17'):
				self.__a_min = 0.
				self.a2_box.config(bg='lightcoral')
				self.a1_box.config(state='disabled')
			elif box_name.endswith('entry18'):
				self.__a_max = float('inf')
				self.a3_box.config(bg='lightcoral')
				self.a1_box.config(state='disabled')
			return False

	def validate_phi(self, new_text, box_name):
		# phi must be greater than 0 and less than 2pi, if accepted, the other boxes should be grayed out
		if not new_text:
			# box cleared
			if box_name.endswith('entry19'):
				self.__phi_fixed = None
				self.phi1_box.config(bg='white')
				self.phi2_box.config(state='normal')
				self.phi3_box.config(state='normal')
			elif box_name.endswith('entry20'):
				self.__phi_min = None
				self.phi2_box.config(bg='white')
				if self.__phi_max is None:
					self.phi1_box.config(state='normal')
			elif box_name.endswith('entry21'):
				self.__phi_max = None
				self.phi3_box.config(bg='white')
				if self.__phi_min is None:
					self.phi1_box.config(state='normal')
			return True
		try:
			phi = float(new_text)
			if 0 <= phi <= 2*np.pi:
				if box_name.endswith('entry19'):
					self.__phi_fixed = phi
					self.phi1_box.config(bg='white')
					self.phi2_box.config(state='disabled')
					self.phi3_box.config(state='disabled')
				elif box_name.endswith('entry20'):
					self.__phi_min = phi
					self.phi2_box.config(bg='white')
					self.phi1_box.config(state='disabled')
				elif box_name.endswith('entry21'):
					self.__phi_max = phi
					self.phi3_box.config(bg='white')
					self.phi1_box.config(state='disabled')
				return True
			else:
				if box_name.endswith('entry19'):
					self.__phi_fixed = phi
					self.phi1_box.config(bg='lightcoral')
					self.phi2_box.config(state='disabled')
					self.phi3_box.config(state='disabled')
				elif box_name.endswith('entry20'):
					self.__phi_min = phi
					self.phi2_box.config(bg='lightcoral')
					self.phi1_box.config(state='disabled')
				elif box_name.endswith('entry21'):
					self.__phi_max = phi
					self.phi3_box.config(bg='lightcoral')
					self.phi1_box.config(state='disabled')
				return False
		except:
			if box_name.endswith('entry19'):
				self.phi1_box.config(bg='lightcoral')
				self.phi2_box.config(state='disabled')
				self.phi3_box.config(state='disabled')
			elif box_name.endswith('entry20'):
				self.__phi_min = 0.
				self.phi2_box.config(bg='lightcoral')
				self.phi1_box.config(state='disabled')
			elif box_name.endswith('entry21'):
				self.__phi_max = 2*np.pi
				self.phi3_box.config(bg='lightcoral')
				self.phi1_box.config(state='disabled')
			return False

	def validate_pd_mu(self, new_text):
		# Check that mu is a number greater than zero
		if not new_text:
			# box cleared
			self.__pd_mu = 5.03
			self.PD_mu_box.insert(-1, '5.03')
			self.PD_mu_box.config(bg='white', fg='gray')
			return True
		try:
			j = float(new_text)
			if j >= 0:
				self.__pd_mu = j
				self.PD_mu_box.config(bg='white', fg='black')
				return True
			else:
				self.PD_mu_box.config(bg='lightcoral', fg='black')
				return False
		except:
			self.PD_mu_box.config(bg='lightcoral', fg='black')
			return False

	def validate_pd_sig(self, new_text):
		# Check that sigma is a number greater than zero
		if not new_text:
			# box cleared
			self.__pd_sig = 2.28
			self.PD_sig_box.config(bg='white', fg='gray')
			self.PD_sig_box.insert(-1, '2.28')
			return True
		try:
			j = float(new_text)
			if j >= 0:
				self.__pd_sig = j
				self.PD_sig_box.config(bg='white', fg='black')
				return True
			else:
				self.PD_sig_box.config(bg='lightcoral', fg='black')
				return False
		except:
			self.PD_sig_box.config(bg='lightcoral', fg='black')
			return False

	def validate_mass_exp(self, new_text):
		# Check that gamma is a number
		if not new_text:
			# box cleared
			self.__q_exp = 0.
			self.mass_exp_box.config(bg='white', fg='gray')
			self.mass_exp_box.insert(-1, '0.0')
			return True
		try:
			j = float(new_text)
			self.__q_exp = j
			self.mass_exp_box.config(bg='white', fg='black')
			return True
		except:
			self.mass_exp_box.config(bg='lightcoral', fg='black')
			return False

	# Help button functions
	@ staticmethod
	def help_age():
		message = """Optional. Estimate of the age of the star. This will effect the modeled magnitude of the stars. 
		If none is entered an age of 5 Gyr will be assumed."""
		tk.messagebox.showinfo('Star Age Help', message)

	@staticmethod
	def help_coord():
		message = 'The right ascension (RA) or declination (Dec) of the star. Use the format 00h00m00.00s for RA and +00d00m00.0s for Dec'
		tk.messagebox.showinfo('Coordinates Help', message)

	@staticmethod
	def help_added_jitter():
		message = 'Optional. A jitter term to be added in quadrature to the measurement error.'
		tk.messagebox.showinfo('Added Jitter Help', message)

	@staticmethod
	def help_floor():
		message = """Optional. A lower limit on the radial velocity semi-amplitude. All generated companions with semi-amplitudes lower than this floor will not be rejected by the RV test. See documentation for more detail"""
		tk.messagebox.showinfo('RV Floor Help', message)

	@staticmethod
	def help_advanced():
		message = 'Alter possible limits of orbit parameters for generated companions. All values must be numbers, ' \
				  'except for limits on cos(i), any one of which may be "transit", to limit to only transiting companions'
		tk.messagebox.showinfo('Advanced Help', message)

	@staticmethod
	def help_period():
		message = 'The default period distribution is a Log-Normal distribution with \u03BC=5.03 and \u03C3=2.28. '\
				  'To use a different log normal distribution enter the mean and standard deviation here.'
		tk.messagebox.showinfo('Period Distribution Help', message)

	@staticmethod
	def help_mass():
		message = 'The default mass ratio distribution is a uniform distribution, a.k.a., a power-law distribution'\
				'with \u03B3=0.0. To use a different power law distribution enter the desired exponent here.'
		tk.messagebox.showinfo('Mass Ratio Distribution Help', message)

	# Run and Close Functions
	def run_code(self):
		self.run_button.focus()
		go = self.check_completeness()
		if go == 0:
			self.parent.run()

	def close(self):
		self.root.quit()

	# Access Functions
	def get_ao_filename(self):
		if self.__ao_check.get():
			file_list = [x.get() for x in self.ao_file_boxes]
			return file_list
		else:
			return ['']

	def get_rv_filename(self):
		if self.__rv_check.get()==1:
			return self.__rv_file.get()
		else:
			return ''

	def get_ruwe(self):
		return self.__ruwe_check.get()

	def get_gaia(self):
		return self.__gaia_check.get()

	def get_filter(self):
		filter_list = [x.get() for x in self.ao_filter_menus]
		return filter_list

	def get_resolution(self):
		return self.__resolution

	def get_prefix(self):
		return self.__prefix.get()

	def get_extra(self):
		if self.__extra.get() == 1:
			return True
		else:
			return False

	def get_all_out_bool(self):
		if self.__all_out.get() == 1:
			return True
		else:
			return False

	def get_ra(self):
		return self.__ra_str

	def get_dec(self):
		return self.__dec_str

	def get_age(self):
		return self.__age

	def get_mass(self):
		return self.__mass

	def get_num_generated(self):
		return self.__num_generated

	def get_added_jitter(self):
		return self.__added_jitter

	def get_rv_floor(self):
		return self.__rv_floor

	def get_limits(self):
		# puts all the limit information together in a list
		limits = [self.__P_fixed, self.__P_min, self.__P_max, self.__i_fixed, self.__i_min, self.__i_max, self.__e_fixed,
				  self.__e_min, self.__e_max, self.__arg_peri_fixed, self.__arg_peri_min, self.__arg_peri_max,
				  self.__m_fixed, self.__m_min, self.__m_max, self.__a_fixed, self.__a_min, self.__a_max, self.__phi_fixed,
				  self.__phi_min, self.__phi_max]
		return limits

	def get_p_dist(self):
		return self.__pd_mu, self.__pd_sig

	def get_q_dist(self):
		return self.__q_exp

	def get_gaia_limit(self):
		return self.__gaia_limit.get()


class Application:
	# Input
	input_args = []
	using_gui = True
	# Class variables
	ao_filename = ['']
	filter = ''
	rv_filename = ''
	resolution = 50000
	prefix = ''
	ruwe_check = False
	gaia_check = False
	# Star parameters
	star_ra = '00h00m00.00s'
	star_dec = '00d00m00.0s'
	star_mass = 0
	star_age = 5
	added_jitter = 20
	rv_floor = 20
	# Code Parameters
	num_generated = 0
	limits = []
	pd_mu = 5.03
	pd_sig = 2.28
	q_exp = 0.0
	gaia_limit = 18.
	extra_output = False
	all_output = False
	# Reject lists
	ao_reject_list = []
	rv_reject_list = []
	jitter_reject_list = []
	ruwe_reject_list = []
	gaia_reject_list = []

	# Class functions
	def __init__(self, input_args):
		self.input_args = input_args

	def start(self):
		# Start has to handle the input arguments
		if len(self.input_args) == 1:
			# run with the gui
			self.using_gui = True
			self.gui = GUI(self)
			self.gui.start()
		else:
			# run from the arguments only
			self.using_gui = False
			if self.parse_input():
				self.run()

	def run(self):
		gc.collect()
		if self.using_gui:
			# Clear display
			self.gui.gui_print('clc')
			# Update status
			self.gui.update_status('Running')
			# Get user inputs from GUI
			self.get_inputs()

		# Print out star info  # TESTING
		self.print_out(('Star RA: ' + self.star_ra))
		self.print_out(('Star DEC: ' + self.star_dec))
		self.print_out(('Star Mass: ' + str(self.star_mass)))

		# Generate Companions
		self.print_out('Generating Companions..')
		comps = Companions(self.num_generated, self.limits, self.star_mass, self.pd_mu, self.pd_sig, self.q_exp)
		failure = self.error_check(comps.generate())
		if failure: return
		self.print_out('Companions Generated')
		self.print_out(('Number of companions: ' + str(self.num_generated)))

		# Decide which parts to run, and run them
		#    AO
		if self.ao_filename[0]:
			ao_reject_lists = []
			for i in range(len(self.ao_filename)):
				if self.ao_filename[i]:
					self.print_out(('Analyzing contrast curve in ' + self.ao_filename[i]))
					ao = AO(self.ao_filename[i], comps, self.star_mass, self.star_age, self.star_ra, self.star_dec, self.filter[i])
					# Determine distance
					failure = self.error_check(ao.get_distance(self.star_ra, self.star_dec))
					if failure: return
					if self.extra_output: self.print_out(('Calculated distance to star: %.0f pc' % (ao.star_distance*4.84e-6)))
					# Read contrast file
					failure = self.error_check(ao.read_contrast())
					if failure: return
					if self.extra_output:
						self.print_out('AO Contrast Loaded')
					# Perform test
					result = ao.analyze()
					failure = self.error_check(result)
					if failure: return
					ao_reject_lists.append(result)
					if self.extra_output:
						self.print_out('Star Model Mag: %.2f' %(ao.star_model_mag))
						self.print_out(('AO Low Mass Limit: %.3f' %(ao.low_mass_limit)))
			# Combine reject lists from all AO tests into one
			self.ao_reject_list = np.logical_or.reduce(ao_reject_lists)
		else:
			self.ao_reject_list = np.array([False]*self.num_generated)

		#   RV and Jitter
		if self.rv_filename:
			# Run RV (without Jitter)
			self.print_out('Analyzing RV...')
			rv = RV(self.rv_filename, self.resolution, comps, self.star_mass, self.star_age, added_jitter=self.added_jitter, rv_floor=self.rv_floor, extra_output=self.extra_output)
			# Read in the RV file
			failure = self.error_check(rv.read_in_rv())
			if failure: return
			if self.extra_output: self.print_out('RV Measurements Loaded')
			# Run analysis
			self.rv_reject_list = rv.analyze_rv()
			self.jitter_reject_list = np.array([False]*self.num_generated)
		else:
			self.rv_reject_list = np.array([False]*self.num_generated)
			self.jitter_reject_list = np.array([False] * self.num_generated)

		#   RUWE
		if self.ruwe_check:
			self.print_out('Analyzing RUWE...')
			ruwe = RUWE(self.star_ra, self.star_dec, self.star_age, self.star_mass, comps)
			# Read in RUWE distribution and Normalization tables
			failure = self.error_check(ruwe.read_dist())
			if failure: return
			# Get Gaia information
			failure = self.error_check(ruwe.get_gaia_info())
			if failure: return
			# Perform Test
			self.ruwe_reject_list = ruwe.analyze()
			if self.extra_output:
				self.print_out(('The star has ln(ruwe) of %f.' % (ruwe.ln_ruwe)))
		else:
			self.ruwe_reject_list = np.array([False]*self.num_generated)

		#   Gaia Contrast
		if self.gaia_check:
			self.print_out('Analyzing Gaia Contrast...')
			#todo improve gaia contrast
			gaia = AO('gaia_contrast.txt', comps, self.star_mass, self.star_age, self.star_ra, self.star_dec, 'G', gaia=True)
			# Determine distance
			failure = self.error_check(gaia.get_distance(self.star_ra, self.star_dec))
			if failure: return
			if self.extra_output: self.print_out(('Calculated distance to star: %.0f pc' % (gaia.star_distance*4.84e-6)))
			# Read contrast file
			failure = self.error_check(gaia.read_contrast())
			if failure: return
			if self.extra_output: self.print_out('Gaia Contrast Loaded')
			# Analyze
			self.gaia_reject_list = gaia.analyze_gaia(self.gaia_limit)
		else:
			self.gaia_reject_list = np.array([False]*self.num_generated)

		# Check successes
		w = [True if x == -1 else False for x in self.ao_reject_list]
		x = [True if x == -1 else False for x in self.rv_reject_list]
		y = [True if x == -1 else False for x in self.jitter_reject_list]
		z = [True if x == -1 else False for x in self.ruwe_reject_list]
		if all(w) or all(x) or all(y) or all(z):
			self.gui.update_status('Finished - Unsuccessful')
			if all(w):
				self.gui.gui_print('AO Problem')
			if all(x):
				self.gui.gui_print('RV Problem')
			if all(y):
				self.gui.gui_print('Jitter problem')
			if all(z):
				self.gui.gui_print('RUWE problem')
			return

		# Put together reject lists and output
		keep = np.invert(self.ao_reject_list)*np.invert(self.rv_reject_list)*np.invert(self.jitter_reject_list)*np.invert(self.ruwe_reject_list)*np.invert(self.gaia_reject_list)
		num_kept = sum([1 for x in keep if x])
		self.print_out(('Total number of surviving binaries: ' + str(num_kept)))

		# Write out files
		#   Write out the survivors file
		cols = ['mass ratio', 'period(days)', 'semi-major axis(AU)', 'cos_i', 'eccentricity', 'arg periastron', 'phase']
		keep_table = np.vstack((comps.get_mass_ratio()[keep], comps.get_P()[keep], comps.get_a()[keep],
								comps.get_cos_i()[keep], comps.get_ecc()[keep], comps.get_arg_peri()[keep],
								comps.get_phase()[keep]))
		if self.ao_filename[0]:
			cols = cols + ['Projected Separation(AU)','Model Contrast']
			keep_table = np.vstack((keep_table, np.array(ao.pro_sep)[keep], np.array(ao.model_contrast)[keep]))
		elif self.gaia_check:
			cols = cols + ['Projected Separation(AU)', 'DeltaG']
			keep_table = np.vstack((keep_table, np.array(gaia.pro_sep)[keep], np.array(gaia.model_contrast)[keep]))
		if self.ruwe_check:
			if 'Projected Separation(AU)' in cols and 'DeltaG' not in cols:
				cols = cols + ['DeltaG','Predicted RUWE']
				keep_table = np.vstack((keep_table, np.array(ruwe.delta_g)[keep], np.array(ruwe.predicted_ruwe)[keep]))
			elif 'Projected Separation(AU)' in cols and 'DeltaG' in cols:
				cols = cols + ['Predicted RUWE']
				keep_table = np.vstack((keep_table, np.array(ruwe.predicted_ruwe)[keep]))
			else:
				cols = cols + ['Projected Separation(AU)','DeltaG','Predicted RUWE']
				keep_table = np.vstack((keep_table, np.array(ruwe.projected_sep)[keep], np.array(ruwe.delta_g)[keep], np.array(ruwe.predicted_ruwe)[keep]))
		if self.rv_filename:
			# # Write out RV calculations, makes v. large files
			# cols = cols + ['rv'+str(i) for i in range(0, len(rv.MJD))]
			# keep_table = np.vstack((keep_table, np.transpose(np.array(rv.predicted_RV)[keep])))
			cols = cols + ['RV Amplitude','Binary Type']
			keep_table = np.vstack((keep_table, np.array(rv.amp)[keep], np.array(rv.b_type)[keep]))
		keep_table = np.transpose(keep_table)
		ascii.write(keep_table, (self.prefix + "_kept.csv"), format='csv', names=cols, overwrite=True)
		self.print_out(('Surviving binary parameters saved to: ' + self.prefix + '_kept.csv'))
		#  Write out the input file
		if self.all_output:
			cols = ['mass ratio', 'period(days)', 'semi-major axis(AU)', 'cos_i', 'eccentricity', 'arg periastron', 'phase']
			all_table = np.vstack((comps.get_mass_ratio(), comps.get_P(), comps.get_a(), comps.get_cos_i(),
								   comps.get_ecc(), comps.get_arg_peri(), comps.get_phase()))
			if self.ao_filename[0]:
				cols = cols + ['Projected Semparation(AU)', 'Model Contrast'] + [('AO Rejected ' + str(i+1)) for i in range(len(ao_reject_lists))]
				all_table = np.vstack((all_table, ao.pro_sep, ao.model_contrast, ao_reject_lists))
			if self.rv_filename:
				# Write out RV calculations
				cols = cols + ['RV Amplitude','Binary Type','RV Rejected']
				all_table = np.vstack((all_table, rv.amp, rv.b_type, self.rv_reject_list))
				# Uncomment to  include RV Calcualations in output file (makes it obnoxiously large)
				# cols = cols + ['rv' + str(i) for i in range(0, len(rv.MJD))]
				# all_table = np.vstack((all_table, np.transpose(np.array(rv.predicted_RV))))
			if self.ruwe_check:
				if 'Projected Separation(AU)' not in cols:
					cols = cols + ['Projected Separation(AU)','DeltaG','Predicted RUWE','RUWE Rejected','RUWE Rejection Prob']
					all_table = np.vstack((all_table, ruwe.projected_sep, ruwe.delta_g, ruwe.predicted_ruwe, self.ruwe_reject_list, ruwe.rejection_prob))
				else:
					cols = cols + ['DeltaG','Predicted RUWE','RUWE Rejected','RUWE Rejection Prob']
					all_table = np.vstack((all_table, ruwe.delta_g, ruwe.predicted_ruwe, self.ruwe_reject_list, ruwe.rejection_prob))
			if self.gaia_check:
				if 'Model Contrast' in cols:
					# The AO test has been run, projected sep and contrast already included
					cols = cols + ['Gaia Rejected']
					all_table = np.vstack((all_table, self.gaia_reject_list))
				elif 'DeltaG' in cols or 'Projected Separation(AU)'in cols:
					# RUWE test has been run, projected sep and contrast already included
					cols = cols + ['Gaia Rejected']
					all_table = np.vstack((all_table, self.gaia_reject_list))
				else:
					cols = cols + ['Projected Separation(AU)', 'DeltaG', 'Gaia Rejected']
					all_table = np.vstack((all_table, gaia.pro_sep, gaia.model_contrast, self.gaia_reject_list))

			# Add column containing full rejection info
			cols = cols + ['Full Rejected']
			all_table = np.vstack((all_table, np.invert(keep)))

			all_table = np.transpose(all_table)
			ascii.write(all_table, (self.prefix + "_all.csv"), format='csv', names=cols, overwrite=True)
			self.print_out(('Generated binary parameters saved to: ' + self.prefix + '_all.csv'))

		if self.using_gui: self.gui.update_status('Finished - Successful')
		self.restore_defaults()
		return

	def get_inputs(self):
		# Get variables from the GUI
		#   Analysis Options
		self.ao_filename = self.gui.get_ao_filename()
		self.rv_filename = self.gui.get_rv_filename()
		self.ruwe_check = self.gui.get_ruwe()
		self.gaia_check = self.gui.get_gaia()
		self.filter = self.gui.get_filter()
		self.resolution = self.gui.get_resolution()
		#   Output Options
		self.prefix = self.gui.get_prefix()
		self.extra_output = self.gui.get_extra()
		self.all_output = self.gui.get_all_out_bool()
		#   Star Information
		self.star_ra = self.gui.get_ra()
		self.star_dec = self.gui.get_dec()
		self.star_age = self.gui.get_age()
		self.star_mass = self.gui.get_mass()
		self.num_generated = self.gui.get_num_generated()
		self.added_jitter = self.gui.get_added_jitter()
		self.rv_floor = self.gui.get_rv_floor()
		#   Limits
		self.limits = self.gui.get_limits()
		# Period distribution
		self.pd_mu, self.pd_sig = self.gui.get_p_dist()
		# Mass Ratio distribution
		self.q_exp = self.gui.get_q_dist()
		# Gaia Completeness
		self.gaia_limit = float(self.gui.get_gaia_limit())
		return

	def error_check(self, error_code):
		# If the return value is a list, array, zero or none, no error was found and the process completed successfully
		if isinstance(error_code, list) or isinstance(error_code, np.ndarray):
			return False
		if error_code == 0 or error_code is None:
			return False
		# todo: Add more clear commenting
		# todo: Add error message for cases of too young or too old in AO
		elif error_code == -1:
			self.print_out('ERROR: File Not Found\nCheck filename and try again')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsucessful')
			return True
		elif error_code == -2:
			self.print_out('ERROR: RV file in incorrect format')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -11:
			self.print_out('Period, mass ratio and semi-major axis cannot all be fixed at the same time.')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -12:
			self.print_out('Incompatible limits on period and semi-major axis.')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -21:
			# Contrast file not found
			self.print_out('ERROR: Contrast file Not Found\nCheck filename and try again')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsucessful')
			return True
		elif error_code == -22:
			self.print_out('ERROR: AO file in incorrect format')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -23:
			# AO, mass higher than model grid
			self.print_out('ERROR: The mass of the primary star is larger than the model grid for that age can handle')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -24:
			# AO, mass lower than model grid
			self.print_out('ERROR: The mass of the primary star is smaller than the model grid for that age can handle')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -31:
			# RV file not found
			self.print_out('ERROR: RV file Not Found\nCheck filename and try again')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsucessful')
			return True
		elif error_code == -32:
			self.print_out('ERROR: RV file in incorrect format')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -41:
			# RUWE Distribution file not found
			self.print_out("""ERROR: RUWE Normalization file not found\nThe file should be named table_u0_g_col.txt 
			and located in the folder the code is being run from.""")
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsucessful')
			return True
		elif error_code == -42:
			# RUWE Norm file wrong format
			self.print_out("""ERROR: RUWE Normalization file is wrong format.""")
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsucessful')
			return True
		elif error_code == -43:
			# RUWE Distribution file not found
			self.print_out("""ERROR: RUWE distribution file not found\nThe file should be named RuweTableGP.txt 
			and located in the folder the code is being run from.""")
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsucessful')
			return True
		elif error_code == -44:
			# RUWE Dist file wrong format
			self.print_out("""ERROR: RUWE distribution file is wrong format.""")
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsucessful')
			return True
		elif error_code == -51:
			self.print_out('Unable to find source matching coordinates in Gaia')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		elif error_code == -52:
			# This error code indicates a warning should be printed, but does not necessitate failure
			self.print_out('WARNING: More than one source within 10 arcseconds of given coordinates. Closest one used.')
			return False
		elif error_code == -53:
			self.print_out('Unable to calculate RUWE for given star. The magnitude or color is outside of bounds')
			if self.using_gui:
				self.gui.update_status('Finished - Unsuccessful')
			else:
				self.print_out('Finished - Unsuccessful')
			return True
		return

	def parse_input(self):
		# Read in command line inputs
		# Create Parser
		my_parser = argparse.ArgumentParser(description='Find limits on the orbital parameters of unseen companions')
		# Add arguments
		#  Required arguments (Positional)
		my_parser.add_argument('prefix', help='The prefix to be used on the output files')
		my_parser.add_argument('ra', help='The RA of the star in hms format', metavar='ra[hms]')
		my_parser.add_argument('dec', help='The Dec of the star in dms format', metavar='dec[dms]')
		my_parser.add_argument('n', help='The number of companions to generate', type=int)
		my_parser.add_argument('mass', help='The mass of the primary in solar masses', type=float, metavar='mass[M_sun]')
		# Optional arguments
		my_parser.add_argument('--age', help='The age of the star in Gyr', type=float, required=False, default=5,
							   metavar='AGE[Gyr]')
		my_parser.add_argument('--jitter', help='The RV jitter to be added in quadrature to the error in m/s',
							   required=False, default=20, metavar='JITTER[m/s]')
		my_parser.add_argument('--rv_floor', help='The lowest RV semi-amplitude which can be rejected in m/s',
							   required=False, default=20, metavar='RV FLOOR[m/s]')
		# Analysis Options
		my_parser.add_argument('--rv', help='The path to the file containing the RV data', required=False,
							   metavar='RV_PATH', default='')
		my_parser.add_argument('--resolution', help='The spectral resolution of the RV data', required=False,
							   metavar='RV_PATH', default=50000)
		my_parser.add_argument('--ao', help='The path to the file containing the AO data', required=False,
							   metavar='AO_PATH', default='')
		my_parser.add_argument('--filter', help='The filter in which the AO data was taken', required=False,
							   choices=['J','K','H','G', 'Bp','Rp','R','I','L','LL','M'])
		my_parser.add_argument('--ruwe', action='store_true', help='Apply the RUWE test')
		my_parser.add_argument('--gaia', action='store_true', help='Apply the GAIA contrast test')
		# Other options
		my_parser.add_argument('-v','--verbose', action='store_true', help='Turn on extra output')
		my_parser.add_argument('-a', '--all', action='store_true', help='Write out all generated companions')
		my_parser.add_argument('--transit', action='store_true', help='Turn on transit limits')
		# Run parser
		args = my_parser.parse_args()
		# Input the inputs
		#  Analysis Options
		self.rv_filename = args.rv
		self.resolution = float(args.resolution)
		self.ao_filename = [args.ao]
		self.filter = args.filter
		self.ruwe_check = args.ruwe
		self.gaia_check = args.gaia
		# Output Options
		self.prefix = args.prefix
		self.extra_output = args.verbose
		self.all_output = args.all
		# Stellar Info
		self.num_generated = args.n
		self.star_ra = args.ra
		self.star_dec = args.dec
		self.star_mass = args.mass
		self.star_age = args.age
		self.added_jitter = float(args.jitter)
		self.rv_floor = args.rv_floor
		# Limits
		self.limits = [None]*21
		if args.transit:
			self.limits[3] = 'transit'

		if not self.ao_filename[0] and not self.rv_filename and not self.ruwe_check and not self.gaia_check and not self.star_jitter:
			print('At least one analysis type needs to be chosen')
			return False
		if self.ao_filename[0] and not self.filter:
			print('AO given without filter')
			return False
		if self.star_jitter != -1 and not self.rv_filename:
			print('RV File needs to be given for stellar jitter test')
			return False

		return True

	def print_out(self, message):
		if self.using_gui:
			self.gui.gui_print(message)
		else:
			print(message)
		return

	def restore_defaults(self):
		# Input
		self.input_args = []
		self.using_gui = True
		# Class variables
		self.ao_filename = ''
		self.filter = ''
		self.rv_filename = ''
		self.prefix = ''
		self.ruwe_check = False
		self.gaia_check = False
		# Star parameters
		self.star_ra = '00h00m00.00s'
		self.star_dec = '00d00m00.0s'
		self.star_mass = 0
		self.star_age = 5
		self.added_jitter = 20
		self.rv_floor = 20
		# Code Parameters
		self.num_generated = 0
		self.limits = []
		self.extra_output = False
		self.all_output = False
		# Reject lists
		self.ao_reject_list = []
		self.rv_reject_list = []
		self.jitter_reject_list = []
		self.ruwe_reject_list = []
		self.gaia_reject_list = []
		return


class Companions:
	# class variables
	__P = []  # days
	__mass_ratio = []  # fraction
	__cos_i = []  # unitless
	__a = []  # AU
	__ecc = []  # unitless
	__arg_peri = []  # radians
	__phase = []  # radians
	num_generated = 0

	# class functions
	def __init__(self, num_generated, limits, star_mass, pd_mu, pd_sig, q_exp):
		self.num_generated = num_generated
		self.limits = limits
		self.star_mass = star_mass
		self.mu_log_P = pd_mu
		self.sig_log_P = pd_sig
		self.q_exp = q_exp

	def generate(self):

		np.random.seed()
		# Period (days), Mass Ratio and Semi-Major Axis (AU)
		G = 39.478  # Gravitational constant in AU^3/years^2*M_solar
		P_fixed = self.limits[0]
		P_lower = self.limits[1]
		P_upper = self.limits[2]
		mass_fixed = self.limits[12]
		mass_lower = self.limits[13]
		mass_upper = self.limits[14]
		a_fixed = self.limits[15]
		a_lower = self.limits[16]
		a_upper = self.limits[17]

		print('Generation Mass Ratio Exponent: ', self.q_exp)
		dq = 1e-4  # coarseness

		# Default lower limit on P
		if P_lower is None or P_lower < 0.1:
			P_lower = 0.1

		if a_fixed is None and a_lower is None and a_upper is None:
			# Generate Period independently of mass ratio and separation
			if P_fixed is not None:
				self.__P = np.array([P_fixed] * self.num_generated)
			else:
				if P_upper is None:
					P_upper = float('inf')
				self.__P = np.array([-1.0] * self.num_generated)  # Initializing
				for i in range(0, self.num_generated):
					while self.__P[i] < P_lower or self.__P[i] > P_upper:
						log_P = np.random.normal(self.mu_log_P, self.sig_log_P)
						self.__P[i] = 10 ** log_P
			# Now Generate mass ratio
			if mass_fixed is not None:
				self.__mass_ratio = np.array([mass_fixed] * self.num_generated)
			elif self.q_exp == 0.0:
				# This generates all the  way to q=0 for a uniform distribution
				if mass_lower is None:
					mass_lower = 0.
				if mass_upper is None:
					mass_upper = 1.
				self.__mass_ratio = np.random.uniform(mass_lower, mass_upper, self.num_generated)
			else:
				# This applies a lower mass limit for the power law distribution
				if mass_lower is None:
					mass_lower = 0.01
				if mass_upper is None:
					mass_upper = 1.
				q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
				p = np.power(q_range, self.q_exp)    # probabilities
				p = [x / np.sum(p) for x in p]  # normalized probabilities
				self.__mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
			# Now calculate Semi-Major Axis using P and m
			self.__a = ((self.__P / 365)**2 * G * self.star_mass*(1 + self.__mass_ratio)/(4 * np.pi ** 2))**(1/3)  # AU
			# Simple Cases Done
		elif a_fixed is not None and mass_fixed is not None and P_fixed is not None:
			# Ideally the GUI will not allow you to try to fix all three of them, but for now I will put it here
			return -11
		elif a_fixed is not None:
			if mass_fixed is not None:
				# Set values of a and mass to the given values, calculate P
				self.__a = np.array([a_fixed] * self.num_generated)
				self.__mass_ratio = np.array([mass_fixed] * self.num_generated)
				self.__P = np.sqrt( (4 * np.pi ** 2 * self.__a ** 3) / (G * self.star_mass * (1 + self.__mass_ratio))) * 365
			elif P_fixed is not None:
				# Set values of a and P to the given values, calculate mass
				self.__a = np.array([a_fixed] * self.num_generated)
				self.__P = np.array([P_fixed] * self.num_generated)
				self.__mass_ratio = np.array(4*np.pi**2 * self.__a**3 / (G * self.star_mass * (self.__P / 365)** 2) - 1)
				if self.__mass_ratio[0] > 1:
					# The fixed values of period and semi-major axis require a mass ratio greater than one
					return -12
			else:
				# Only a is fixed.
				self.__a = np.array([a_fixed] * self.num_generated)
				if P_lower is None and P_upper is None:
					# Generate mass ratio as normal and then calculate P
					# Choose to generate mass ratio instead of P, to ensure that it stays within limits
					if mass_lower is None:
						mass_lower = 0.01
					if mass_upper is None:
						mass_upper = 1.
					q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
					p = np.power(q_range, self.q_exp)  # probabilities
					p = [x / np.sum(p) for x in p]  # normalized probabilities
					self.__mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
					self.__P = np.sqrt(
						(4 * np.pi ** 2 * self.__a ** 3) / (G * self.star_mass * (1 + self.__mass_ratio))) * 365
				else:
					# Generate P within limits, and then calculate mass ratio
					if P_lower is None:
						P_lower = 0
					if P_upper is None:
						P_upper = float('inf')
					self.__P = np.array([-1.0] * self.num_generated)  # Initializing
					for i in range(0, self.num_generated):
						while self.__P[i] < P_lower or self.__P[i] > P_upper:
							log_P = np.random.normal(self.mu_log_P, self.sig_log_P)
							self.__P[i] = 10 ** log_P
					self.__mass_ratio = np.array(4 * np.pi ** 2 * self.__a ** 3 / (G * self.star_mass * (self.__P / 365) ** 2) - 1)
		elif a_lower is not None or a_upper is not None:
			if P_fixed is not None:
				self.__P = np.array([P_fixed] * self.num_generated)
				# Generate limits on mass based on limits on a
				if mass_lower is None:
					mass_lower = 0.
				if mass_upper is None:
					mass_upper = 1.
				# The limits will now need to be a list because the mass ratio for a given separation depends on the P
				mass_lower_list = np.array([0.] * self.num_generated)  # Initialization
				mass_upper_list = np.array([float('inf')] * self.num_generated)  # Initialization

				if a_lower is not None:
					mass_lower_list = np.array( 4 * np.pi ** 2 * a_lower ** 3 / (G * self.star_mass * (self.__P / 365) ** 2) - 1)
				if a_upper is not None:
					mass_upper_list = np.array(
						4 * np.pi ** 2 * a_upper ** 3 / (G * self.star_mass * (self.__P / 365) ** 2) - 1)

				mass_lower_list = [max(x, mass_lower) for x in mass_lower_list]
				mass_upper_list = [min(x, mass_upper) for x in mass_upper_list]

				self.__mass_ratio = np.array([0.] * self.num_generated)
				for i in range(0, self.num_generated):
					q_range = np.arange(mass_lower, 1 + dq, dq)  # range of possible q
					p = np.power(q_range, self.q_exp)  # probabilities
					p = [x / np.sum(p) for x in p]  # normalized probabilities
					self.__mass_ratio[i] = np.random.choice(q_range, p=p, size=1)
				# Calculate a
				self.__a = ((self.__P / 365) ** 2 * G * self.star_mass * (1 + self.__mass_ratio) / (
							4 * np.pi ** 2)) ** (1 / 3)  # AU
			else:
				# Generate Mass
				if mass_fixed is not None:
					self.__mass_ratio = np.array([mass_fixed] * self.num_generated)
				elif self.q_exp == 0.0:
					# This generates all the  way to q=0 for a uniform distribution
					if mass_lower is None:
						mass_lower = 0.
					if mass_upper is None:
						mass_upper = 1.
					self.__mass_ratio = np.random.uniform(mass_lower, mass_upper, self.num_generated)
				else:
					# This applies a lower mass limit for the power law distribution
					if mass_lower is None:
						mass_lower = 0.01
					if mass_upper is None:
						mass_upper = 1.
					q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
					p = np.power(q_range, self.q_exp)  # probabilities
					p = [x / np.sum(p) for x in p]  # normalized probabilities
					self.__mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
				# Generate limits on period, based on limits on a
				if P_lower is None:
					P_lower = 0.
				if P_upper is None:
					P_upper = float('inf')

				# The limits will now need to be a list because the period for a given separation depends on the mass
				P_lower_list = np.array([0.] * self.num_generated)  # Initialization
				P_upper_list = np.array([float('inf')] * self.num_generated)  # Initialization

				if a_lower is not None:
					P_lower_list = np.sqrt(
						(4 * np.pi ** 2 * a_lower ** 3) / (G * self.star_mass * (1 + self.__mass_ratio))) * 365
				if a_upper is not None:
					P_upper_list = np.sqrt(
						(4 * np.pi ** 2 * a_upper ** 3) / (G * self.star_mass * (1 + self.__mass_ratio))) * 365

				P_lower_list = [max(x, P_lower) for x in P_lower_list]
				P_upper_list = [min(x, P_upper) for x in P_upper_list]

				# Generate period between limits
				self.__P = np.array([-1.0] * self.num_generated)  # Initializing
				for i in range(0, self.num_generated):
					if P_lower_list[i] > P_upper_list[i]:
						return -11
					while self.__P[i] < P_lower_list[i] or self.__P[i] > P_upper_list[i]:
						log_P = np.random.normal(self.mu_log_P, self.sig_log_P)
						self.__P[i] = 10 ** log_P
				# Calculate a
				self.__a = ((self.__P / 365)**2 * G * self.star_mass * (1 + self.__mass_ratio) / (
							4 * np.pi ** 2))**(1/3)  # AU
		else:
			print("I really don't think you should be able to get here ever...")

		# Inclination
		cos_i_fixed = self.limits[3]
		cos_i_lower = self.limits[4]
		cos_i_upper = self.limits[5]

		# todo figure out what to do about stellar radius in this calculation
		if cos_i_lower == 'transit' or cos_i_upper == 'transit' or cos_i_fixed == 'transit':
			self.__cos_i = np.zeros(self.num_generated)
			b_limit = 1.01
			R_star = 1.  # R_sun
			r_sun_to_au = 0.00465
			for i in range(0, self.num_generated):
				# Determine the upper limit on cos(i), based on semi-major axis, stellar radius, and impact parameter
				cos_i_upper = np.min([(b_limit * r_sun_to_au * R_star)/(self.__a[i]), 1.0])
				#  Generate
				self.__cos_i[i] = np.random.uniform(0., cos_i_upper)
		elif cos_i_fixed is not None:
			self.__cos_i = np.array([cos_i_fixed] * self.num_generated)
		else:
			if cos_i_lower is None:
				cos_i_lower = 0.
			if cos_i_upper is None:
				cos_i_upper = 1.
			self.__cos_i = np.random.uniform(cos_i_lower, cos_i_upper, self.num_generated)  # uniform distribution between 0 and 1

		# Pericenter Phase (radians)
		phase_fixed = self.limits[18]
		phase_lower = self.limits[19]
		phase_upper = self.limits[20]

		# if phase_lower is None and phase_upper is None and phase_fixed is None:
		#     self.__phase = np.random.uniform(0, 2*np.pi, self.num_generated)
		if phase_fixed is not None:
			self.__phase = np.array([phase_fixed] * self.num_generated)
		else:
			if phase_lower is None:
				phase_lower = 0.
			if phase_upper is None:
				phase_upper = 2. * np.pi
			self.__phase = np.random.uniform(phase_lower, phase_upper, self.num_generated)

		# Eccentricity
		e_fixed = self.limits[6]
		e_lower = self.limits[7]
		e_upper = self.limits[8]

		log_P = np.log10(self.__P)

		if e_fixed is not None:
			self.__ecc = np.array([e_fixed] * self.num_generated)
		else:
			# Choose e from within given limits
			a, b, c, d = [0.148, 0.001, 0.042, 0.128]  # parameters from fitting
			dx = 1e-4
			if e_lower is None or e_lower < dx:
				e_lower = dx
			if e_upper is None or e_upper > 0.9999:
				e_upper = 0.9999

			ecc = np.zeros(self.num_generated)
			for i in range(self.num_generated):
				x = np.arange(e_lower, e_upper + dx, dx)  # range of eccentricities
				if log_P[i] < 3:
					# Calculate the mean and std deviation at this period
					e_mu = a * log_P[i] + b
					e_sig = c * log_P[i] + d
					# Construct the probability distribution
					#    probability distribution is Gaussian between 0 and 1, and zero outside those bounds
					y = self.gauss(x, e_mu, e_sig)
					y = np.divide(y, np.sum(y))
					# Generate eccentricity
					ecc[i] = np.random.choice(x, p=y)
				else:
					ecc[i] = np.random.uniform(e_lower, e_upper)
			self.__ecc = ecc

		# Argument of Periapsis (radians)
		arg_fixed = self.limits[9]
		arg_lower = self.limits[10]
		arg_upper = self.limits[11]

		if arg_fixed is not None:
			self.__arg_peri = np.array([arg_fixed] * self.num_generated)
		else:
			if arg_lower is None:
				arg_lower = 0.
			if arg_upper is None:
				arg_upper = np.pi
			self.__arg_peri = np.random.uniform(arg_lower, arg_upper, self.num_generated)

		return

	# Distribution functions
	def gauss(self, x, *p):
		mu, sigma = p
		return np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

	# Accessor functions
	def get_num(self):
		return self.num_generated

	def get_P(self):
		return self.__P

	def get_mass_ratio(self):
		return self.__mass_ratio

	def get_cos_i(self):
		return self.__cos_i

	def get_a(self):
		return self.__a

	def get_ecc(self):
		return self.__ecc

	def get_arg_peri(self):
		return self.__arg_peri

	def get_phase(self):
		return self.__phase

	def get_all(self):
		return np.vstack([self.__P,self.__mass_ratio,self.__cos_i,self.__a,self.__ecc,self.__arg_peri,self.__phase])


class AO:
	# class variables
	ao_filename = ''
	reject_list = []
	age_model = {}  # Table containing stellar model for stars of this age
	star_distance = 0
	star_mass = 0
	mass_ratio = []
	a = []
	a_type = ''
	obs_date = None

	# Input file must be in mas and mags

	def __init__(self, filename, companions, star_mass, star_age, star_ra, star_dec, filter, gaia=False):
		self.ao_filename = filename
		self.mass_ratio = companions.get_mass_ratio()
		self.a = companions.get_a()
		self.companions = companions
		self.star_mass = star_mass
		# Get the appropriate model and year
		self.age_model = self.load_stellar_model(filter, star_age)
		if gaia:
			self.a_type = 'gaia'

	def analyze(self):
		t = time()
		num_generated = len(self.mass_ratio)

		# Unpack companions' orbital parameters
		period = self.companions.get_P()
		e = self.companions.get_ecc()
		arg_peri = self.companions.get_arg_peri()
		phase = self.companions.get_phase()
		cos_i = self.companions.get_cos_i()
		if self.obs_date:
			T_0 = self.obs_date
		else:
			T_0 = 2457388.5  # epoch 2016.0 in JD

		# Read in the contrast
		contrast = self.contrast
		a_type = self.a_type

		#todo add error handling for cases of too young and too old

		# Determine low and high mass limits
		low_mass_limit = self.age_model['M/Ms'][0]
		high_mass_limit = self.age_model['M/Ms'][-1]

		if self.star_mass > high_mass_limit:
			return -23
		elif self.star_mass < low_mass_limit:
			return -24

		t1 = time()

		# Find model mag of primary star
		star_model_mag = self.find_mag(self.star_mass, self.age_model)
		self.star_model_mag = star_model_mag
		print('Star Model Mag', self.star_model_mag)

		# Get masses of companion stars
		cmp_mass = self.star_mass * self.mass_ratio  # companion mass in solar masses
		cmp_mass = [round(i, 3) for i in cmp_mass]

		# Get companion star magnitudes, assign infinite magnitude if below lowest modeled mass
		f_mag = scipy.interpolate.interp1d(self.age_model['M/Ms'], self.age_model['Mag'], kind='cubic', fill_value='extrapolate')
		cmp_model_mag = [f_mag(x) if x >= low_mass_limit else float('Inf') for x in cmp_mass]

		t2 = time()

		# Find Delta Mag
		model_contrast = [x - star_model_mag for x in cmp_model_mag]

		# Calculate projected separation for each generated companion
		pro_sep = [0.0 for i in range(num_generated)]
		for i in range(num_generated):
			# 1. Calculate mean anomaly
			M = 2 * np.pi * T_0 / period[i] - phase[i]
			# 2. Calculate eccentric anomaly iteratively
			prev_E = 0.0
			current_E = M
			while abs(current_E - prev_E) > 0.00001:
				prev_E = current_E
				current_E = M + e[i] * np.sin(prev_E)
			# 3. Calculate true anomaly
			f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - e[i]) / (1 + e[i])))
			# 4. Calculate projected separation in AU
			alpha = f + arg_peri[i]
			sqt = np.sqrt(np.sin(alpha)**2+np.cos(alpha)**2 * cos_i[i]**2)
			pro_sep[i] = self.a[i] * (1-e[i]**2)/(1+e[i]*np.cos(f))*sqt

		four_arc = round(self.star_distance * 0.0000193906, 1)  # 4" in AU at distance of primary

		if a_type == 'hard limit':
			# Find Delta Mag limits, and/or recovery fraction
			# Interpolate linearly between given points to get the estimated contrast limit
			# Reject or accept the hypothetical binary based on the "hard limit" of the experimental contrast
			# 100% of binaries with a contrast less than the experimental contrast are rejected, 100% of those with
			# a greater contrast cannot be rejected
			contrast.sort('Sep (AU)')
			max_sep = contrast['Sep (AU)'][-1]
			min_sep = contrast['Sep (AU)'][0]
			f_con = scipy.interpolate.interp1d(contrast['Sep (AU)'], contrast['Contrast'], kind='linear', fill_value='extrapolate')
			contrast_limit = [f_con(pro_sep[i]) if min_sep < pro_sep[i] < max_sep else 0. for i in range(num_generated)]

			# Compare the model_contrast to the experimental_delta_K
			# If the model contrast is less than the experimental contrast it would have been seen and can be rejected
			self.reject_list = [False] * num_generated
			for i in range(0, num_generated):
				if model_contrast[i] < contrast_limit[i]:
					self.reject_list[i] = True
				else:
					self.reject_list[i] = False

		elif a_type == 'gradient':
			# Find Delta Mag limits, and/or recovery fraction
			# Interpolate linearly between given points to get the estimated contrast limit
			# Reject or accept the hypothetical binary based on a gradient of recovery rates for separation and contrast
			recovery_rate = [0.]*num_generated

			# Get column names and recovery rates
			column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
			column_names = contrast.colnames[1:]

			f = scipy.interpolate.interp2d(contrast['Sep (AU)'], column_rates, [contrast[x] for x in column_names])

			for i in range(num_generated):
				column_names = contrast.colnames  # reset the list of column names
				# Interpolate
				if pro_sep[i] < contrast['Sep (AU)'][0]:  # closer than lowest limit, recovery rate = 0
					recovery_rate[i] = 0.
					continue
				elif pro_sep[i] > contrast['Sep (AU)'][-1]:  # further than farthest limit, recovery rate = 1
					recovery_rate[i] = 1.
					continue
				else:
					new_row = Table(rows=[[float(f(pro_sep[i], x)) for x in column_rates]], names=column_names[1:])
					new_row['Sep (AU)'] = pro_sep[i]
					# Determine which recovery rates the magnitude falls between, and assign it the lower one
					if model_contrast[i] < new_row[column_names[-1]]:  # Greater than largest recovery rate
						recovery_rate[i] = column_rates[-1]
						continue
					elif model_contrast[i] > new_row[column_names[1]]:  # Less than smallest recovery rate
						recovery_rate[i] = 0.
						continue
					else:
						for j in range(1, len(column_names)):
							if new_row[column_names[j]][0] < model_contrast[i]:
								recovery_rate[i] = column_rates[j-2]
								break

			# Make Reject list
			random = np.random.uniform(0, 1, num_generated)
			self.reject_list = [True if random[i] < recovery_rate[i] else False for i in range(0, num_generated)]

		# Write out information for display or output files
		self.model_contrast = model_contrast
		self.pro_sep = pro_sep
		self.low_mass_limit = low_mass_limit

		return np.array(self.reject_list)

	def analyze_gaia(self, gaia_limit):
		# Unpack companions' orbital parameters
		period = self.companions.get_P()
		e = self.companions.get_ecc()
		arg_peri = self.companions.get_arg_peri()
		phase = self.companions.get_phase()
		cos_i = self.companions.get_cos_i()
		num_generated = len(self.mass_ratio)

		# Set date
		T_0 = 2457388.5  # epoch 2016.0 in JD

		# Read in the contrast
		contrast = self.contrast
		a_type = self.a_type

		#todo add error handling for cases of too young and too old

		# Determine low and high mass limits
		low_mass_limit = self.age_model['M/Ms'][0]
		high_mass_limit = self.age_model['M/Ms'][-1]

		if self.star_mass > high_mass_limit:
			return -23
		elif self.star_mass < low_mass_limit:
			return -24

		# Find model mag of primary star
		star_model_mag = self.find_mag(self.star_mass, self.age_model)
		self.star_model_mag = star_model_mag
		print('Star Model Mag', self.star_model_mag)  # TESTING

		# Get masses of companion stars
		cmp_mass = self.star_mass * self.mass_ratio  # companion mass in solar masses
		cmp_mass = [round(i, 3) for i in cmp_mass]

		# Get companion star magnitudes, assign infinite magnitude if below lowest modeled mass
		f_mag = scipy.interpolate.interp1d(self.age_model['M/Ms'], self.age_model['Mag'], kind='cubic', fill_value='extrapolate')
		cmp_model_mag = [f_mag(x) if x >= low_mass_limit else float('Inf') for x in cmp_mass]

		# Find Delta Mag
		model_contrast = [x - star_model_mag for x in cmp_model_mag]

		# Calculate projected separation for each generated companion
		pro_sep = [0.0 for i in range(num_generated)]
		for i in range(num_generated):
			# 1. Calculate mean anomaly
			M = 2 * np.pi * T_0 / period[i] - phase[i]
			# 2. Calculate eccentric anomaly iteratively
			prev_E = 0.0
			current_E = M
			while abs(current_E - prev_E) > 0.00001:
				prev_E = current_E
				current_E = M + e[i] * np.sin(prev_E)
			# 3. Calculate true anomaly
			f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - e[i]) / (1 + e[i])))
			# 4. Calculate projected separation in AU
			alpha = f + arg_peri[i]
			sqt = np.sqrt(np.sin(alpha)**2+np.cos(alpha)**2 * cos_i[i]**2)
			pro_sep[i] = self.a[i] * (1-e[i]**2)/(1+e[i]*np.cos(f))*sqt

		#  Determine Gaia completeness detection limits
		four_arc = round(self.star_distance * 0.0000193906, 1)  # 4" in AU at distance of primary
		completness_absolute = gaia_limit - 5 * np.log10(self.star_distance / 2062650)  # apparent converted to absolute
		completeness_mag = round(completness_absolute - self.star_model_mag, 2)  # delta mag between the primary and gaia's detection limit

		for i in range(len(contrast)):
			for j in range(1, len(contrast.colnames)):
				# Adjusting to the gaia completeness mag
				if contrast[i][j] > completeness_mag:
					contrast[i][j] = completeness_mag

		# Comment this section if you want nearest neighbor limits
		contrast['100%'] = [0. if x < four_arc else completeness_mag for x in contrast['Sep (AU)']]
		contrast.add_row(([four_arc] + [completeness_mag] * (len(contrast.colnames) - 1)))
		contrast.add_row(([self.star_distance] + [completeness_mag] * (len(contrast.colnames) - 1)))
		contrast.sort('Sep (AU)')

		#  end neighbor-less segment

		# Uncomment this section if you want nearest neighbor limits
		# if self.nearest_neighbor_dist < four_arc:
		#     # If the nearest neighbor is less than 4 arcseconds away, I need to not add a 4" row, and truncate
		#     # the existing rows to a maximum of the nearest neighbor distance
		#     # first, find out what the interpolated limits are at the distance of the nearest neighbor
		#     column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
		#     column_names = contrast.colnames[1:]
		#
		#     f_neighbor = scipy.interpolate.interp2d(contrast['Sep (AU)'], column_rates, [contrast[x] for x in column_names])
		#
		#     l = [round(float(f_neighbor(self.nearest_neighbor_dist, x)),2) for x in column_rates]
		#     contrast.add_row(([self.nearest_neighbor_dist]+l))
		#
		#     # Sort by separation
		#     contrast.sort('Sep (AU)')
		#     # Remove all rows after the nearest neighbor. There probably won't be any, since I'm assuming no
		#     # one is going to give me a contrast curve that has another star in it, but the code is here anyways
		#     ind = list(contrast['Sep (AU)']).index(self.nearest_neighbor_dist)
		#     contrast.remove_rows(slice(ind+1, len(contrast)))
		#
		# else:
		#     # Add two rows at the bottom, reaching from 4" to the nearest neighbor
		#     contrast['100%'] = [0.0]*len(contrast)
		#     contrast.add_row(([four_arc]+[completeness_mag]*(len(contrast.colnames)-1)))
		#     contrast.add_row(([self.nearest_neighbor_dist]+[completeness_mag]*(len(contrast.colnames)-1)))
		# # end nearest neighbor segment

		# Get column names and recovery rates
		column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
		column_names = contrast.colnames[1:]
		recovery_rate = [0.] * num_generated

		f = scipy.interpolate.interp2d(contrast['Sep (AU)'], column_rates, [contrast[x] for x in column_names])

		for i in range(num_generated):
			column_names = contrast.colnames  # reset the list of column names
			# Interpolate
			if pro_sep[i] < contrast['Sep (AU)'][0]:  # closer than lowest limit, recovery rate = 0
				recovery_rate[i] = 0.
				continue
			elif pro_sep[i] > contrast['Sep (AU)'][-1]:  # further than farthest limit, recovery rate = 1
				recovery_rate[i] = 1.
				continue
			else:
				new_row = Table(rows=[[float(f(pro_sep[i], x)) for x in column_rates]], names=column_names[1:])
				new_row['Sep (AU)'] = pro_sep[i]
				# Determine which recovery rates the magnitude falls between, and assign it the lower one
				if model_contrast[i] < new_row[column_names[-1]]:  # Greater than largest recovery rate
					recovery_rate[i] = column_rates[-1]
					continue
				elif model_contrast[i] > new_row[column_names[1]]:  # Less than smallest recovery rate
					recovery_rate[i] = 0.
					continue
				else:
					for j in range(1, len(column_names)):
						if new_row[column_names[j]][0] < model_contrast[i]:
							recovery_rate[i] = column_rates[j - 2]
							break

		# Make Reject list
		random = np.random.uniform(0, 1, num_generated)
		self.reject_list = [True if random[i] < recovery_rate[i] else False for i in range(0, num_generated)]

		# Write out information for display or output files
		self.model_contrast = model_contrast
		self.pro_sep = pro_sep
		self.low_mass_limit = low_mass_limit

		return np.array(self.reject_list)

	def find_mag(self, star_mass, t):
		f_mag = scipy.interpolate.interp1d(t['M/Ms'], t['Mag'], kind='cubic', fill_value='extrapolate')
		mag = f_mag(star_mass)
		return mag

	def get_distance(self, star_RA, star_DEC):
		coordinate = SkyCoord(star_RA, star_DEC, frame='icrs')
		width = u.Quantity(10, u.arcsecond)
		height = u.Quantity(10, u.arcsecond)
		job_str = ("SELECT TOP 10 DISTANCE(POINT('ICRS', %f, %f), POINT('ICRS', ra, dec)) AS dist, * FROM gaiaedr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', %f, %f),CIRCLE('ICRS', ra, dec, 0.08333333)) ORDER BY dist ASC)" % (coordinate.ra.degree, coordinate.dec.degree, coordinate.ra.degree, coordinate.dec.degree))
		job = Gaia.launch_job(job_str)
		gaia_info = job.get_results()
		if gaia_info:
			if len(gaia_info) > 1:
				# Multiple possible sources. Sort by search distance and take the closest one. Print warning.
				# Nearest neighbor distance is set to distance between chosen source and it's closest neighbor in the search region
				gaia_info.sort('dist')
				parallax = gaia_info['parallax'][0]
				picked_coord = SkyCoord(gaia_info['ra'][0], gaia_info['dec'][0], unit='degree')
				# Convert star parallax to distance in AU: d[AU] = 1/p["] * 206265 AU/parsec
				star_distance = 1 / (parallax / 1000) * 206265
				self.star_distance = star_distance

				# Find nearest neighbor
				#    calculate distance between the picked star and all other stars in search radius
				gaia_info['separation'] = [picked_coord.separation(SkyCoord(x['ra'], x['dec'], unit='degree')).mas for x in gaia_info]
				#    sort the search results by separation from picked star (the picked star should have a separation of zero)
				gaia_info.sort('separation')
				#    choose the closest one to the picked star
				nearest_neighbor_dist = gaia_info['separation'][1]  # distance to n.n. in mas
				#   convert from mas to AU
				self.nearest_neighbor_dist = round(self.star_distance * np.tan(np.radians(nearest_neighbor_dist/(3.6e6))), 1)  # distance to n.n. in AU
				return -52
			else:
				# Only once source in search area. Nearest neighbor distance is set to search width
				parallax = gaia_info['parallax'][0]
				# Convert star parallax to distance in AU: d[AU] = 1/p[as] * 206265AU/parsec
				star_distance = 1 / (parallax / 1000) * 206265
				self.star_distance = star_distance

				# Set nearest neighbor distance to search distance
				nearest_neighbor_dist = width.to('mas').value
				#   convert from mas to AU
				self.nearest_neighbor_dist = round(self.star_distance * np.tan(np.radians(nearest_neighbor_dist/(3.6e6))), 1)  # distance to n.n. in AU
				return 0
		else:
			return -51

	def load_stellar_model(self, filter, star_age):
		print('LOADING STELLAR MODEL...')
		# Read in file containing stellar model with the filter needed
		if filter == 'J' or filter == 'H' or filter == 'K':  # 2MASS filters
			model_chart = {}
			BHAC_file = 'BHAC15_2MASS.txt'
			with open(BHAC_file, 'r') as content_file:
				content = content_file.read()
			tables = content.split(
				sep='\n-----------------------------------------------------------------------------------------------\n')
			tables = [x for x in tables if len(x) > 1]

			for table in tables:
				n = table.find('M')
				time_segment = table[0:n]
				table_segment = table[n:]
				age = float(time_segment[time_segment.find('=') + 1:])
				year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')
				if filter == 'J':
					year_chart = year_chart['M/Ms', 'Mj']
					year_chart.rename_column('Mj', 'Mag')
				elif filter == 'H':
					year_chart = year_chart['M/Ms', 'Mh']
					year_chart.rename_column('Mh', 'Mag')
				elif filter == 'K':
					year_chart = year_chart['M/Ms', 'Mk']
					year_chart.rename_column('Mk', 'Mag')
				model_chart[age] = year_chart
		elif filter == 'R' or filter == 'I':
			model_chart = {}
			BHAC_file = 'BHAC15_CFHT.txt'
			with open(BHAC_file, 'r') as content_file:
				content = content_file.read()
			tables = content.split(
				sep='\n----------------------------------------------------------------------------------------------------------------------------------------\n')
			tables = [x for x in tables if len(x) > 1]
			for table in tables:
				n = table.find('M')
				time_segment = table[0:n]
				table_segment = table[n:]
				age = float(time_segment[time_segment.find('=') + 1:])
				year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')
				if filter == 'R':
					year_chart = year_chart['M/Ms', 'R']
					year_chart.rename_column('R', 'Mag')
				elif filter == 'I':
					year_chart = year_chart['M/Ms', 'I']
					year_chart.rename_column('I', 'Mag')
				model_chart[age] = year_chart
		elif filter == 'G' or filter == 'Rp' or filter ==  'Bp':
			model_chart = {}
			BHAC_file = 'BHAC15_GAIA.txt'
			with open(BHAC_file, 'r') as content_file:
				content = content_file.read()
			tables = content.split(
				sep='\n----------------------------------------------------------------------------------------------------------------------------------------\n')
			tables = [x for x in tables if len(x) > 1]
			for table in tables:
				n = table.find('M')
				time_segment = table[0:n]
				table_segment = table[n:]
				age = float(time_segment[time_segment.find('=') + 1:])
				year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')
				if filter == 'G':
					year_chart = year_chart['M/Ms', 'G']
					year_chart.rename_column('G', 'Mag')
				elif filter == 'Rp':
					year_chart = year_chart['M/Ms', 'G_RP']
					year_chart.rename_column('G_RP', 'Mag')
				elif filter == 'Bp':
					year_chart = year_chart['M/Ms', 'G_BP']
					year_chart.rename_column('G_BP', 'Mag')
				model_chart[age] = year_chart
		elif filter == 'L' or filter == 'LL' or filter == 'M':
			print('MODEL:  CIT2')
			model_chart = {}
			BHAC_file = 'BHAC15_CIT2.txt'
			with open(BHAC_file, 'r') as content_file:
				content = content_file.read()
			content = content[content.find('\n\n\n'):]  # Cut off the intro material
			tables = content.split(
				sep='\n!-----------------------------------------------------------------------------------------------\n\n')
			tables = [x for x in tables if len(x) > 1]
			for table in tables:
				n1 = table.find('t (Gyr)')
				n2 = table.find('M')
				time_segment = table[n1:table[n1:].find('!') + n1]
				table_segment = table[n2:]

				age = float(time_segment[time_segment.find('=') + 1:])
				year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')

				if filter == 'L':
					year_chart = year_chart['M/Ms', 'Ml']
					year_chart.rename_column('Ml', 'Mag')
				elif filter == 'LL':
					year_chart = year_chart['M/Ms', 'Mll']
					year_chart.rename_column('Mll', 'Mag')
				elif filter == 'M':
					year_chart = year_chart['M/Ms', 'Mm']
					year_chart.rename_column('Mm', 'Mag')
				model_chart[age] = year_chart
		# Set up interpolation
		ages = np.array(list(model_chart.keys()))
		#  check if age is modeled, if it is simply take that chart, if not interpolate within the model to get it
		if star_age in ages:
			return model_chart[star_age]
		#  find ages above and below the desired age
		diff = [star_age - x for x in ages]
		low_age = ages[np.argmin([x if x > 0 else float('inf') for x in diff])]
		high_age = ages[np.argmin([abs(x) if x < 0 else float('inf') for x in diff])]
		#  get masses
		new_model = Table()
		new_model['M/Ms'] = model_chart[high_age]['M/Ms']
		for y in model_chart[low_age].colnames[1:]:
			y_list = []
			for m in list(new_model['M/Ms']):
				# interpolate
				low_i = list(model_chart[low_age]['M/Ms']).index(m)
				high_i = list(model_chart[high_age]['M/Ms']).index(m)
				xs = [low_age, high_age]
				ys = [model_chart[low_age][y][low_i], model_chart[high_age][y][high_i]]
				f = scipy.interpolate.interp1d(xs, ys, kind='linear')
				y_list.append(f(star_age))
			#  add to table
			new_model[y] = y_list
		return new_model

	def read_contrast(self):
		contrast = {}
		# Each line in the file should have the separation in mas and the contrast for different recovery rates,
		# separated by whitespace
		# If a date is included it should be the first line of the file, starting with a # and followed by the JD date
		try:
			contrast_in = open(self.ao_filename, 'r')
		except FileNotFoundError:
			return -21
		try:
			f = open(self.ao_filename, 'r')
			read_data = f.read()
			if read_data.startswith('#'):
				# Remove the date from the first part
				split = read_data.find('\n')
				# Get the date
				self.obs_date = float(read_data[1:split].lstrip())
				# Get the contrast
				contrast = read_data[split:]
				contrast_table = Table.read(contrast, format='ascii')
			else:
				# No date given, just get the contrast
				contrast_table = Table.read(read_data, format='ascii', delimiter=' ', fast_reader=False)
			# Rename columns, assuming they are in the correct format of separation, magnitude
			contrast_table.rename_column(list(contrast_table.columns)[0], 'Sep')
			# Convert separation from mas to AU, order columns correctly
			contrast_table['Sep (AU)'] = [round(self.star_distance * np.tan(np.radians(x/(3.6e6))), 1) for x in contrast_table['Sep']]
			order = ['Sep (AU)'] + list(contrast_table.columns)[1:-1]
			contrast_table = contrast_table[order]
		except TypeError:
			return -22

		# If there is only one contrast given (hard limit)
		if len(contrast_table.columns) == 2:
			self.a_type = 'hard limit'
			contrast_table.rename_column(list(contrast_table.columns)[1], 'Contrast')
		# Gradient contrast given
		# The first line should be a header, and should list the recovery rates used
		else:
			if self.a_type != 'gaia':
				self.a_type = 'gradient'

		self.contrast = contrast_table

		return 0


class RV:
	# Needs from application: filename, star mass, all generated values (P, e, arg_peri, cos_i, mass_ratio, a)
	# Needs to give back reject list
	# class variables
	rv_filename = ''
	added_jitter = 0.02 # km/s
	rv_floor = 0.02  # km/s
	rv_reject_list = []
	jitter_reject_list = []
	MJD = []
	experimental_RV = []
	measurement_error = []
	predicted_RV = []
	jitter = []

	# class functions
	def __init__(self, filename, resolution, companions, mass, age, added_jitter=0, rv_floor=20, extra_output=True):
		self.restore_defaults()
		self.rv_filename = filename
		self.resolution = resolution
		self.companions = companions
		self.star_mass = mass
		self.added_jitter = added_jitter/1000  # convert to km/s
		self.rv_floor = rv_floor/1000  # km/s
		self.extra_output = extra_output
		self.model = self.load_stellar_model('G', age)

	def analyze_rv(self):
		# Calculate predicted RV
		#     Using orbital mechanics equations from Perryman 2011, and solving for E(anomaly) numerically
		#     For time comparison I want to generate times with the same range as the times in MJD, and then calculate
		#     the predicted RV at that time, which I can then compare to the experimental values using least squares
		num_generated = self.companions.get_num()
		# Unpack parameters
		period = self.companions.get_P()
		mass_ratio = self.companions.get_mass_ratio()
		a = self.companions.get_a()
		e = self.companions.get_ecc()
		cos_i = self.companions.get_cos_i()
		arg_peri = self.companions.get_arg_peri()
		phase = self.companions.get_phase()

		# Determine velocity limit
		delta_v = 2.998e5 / self.resolution  # m/s

		t = time()

		# Determine contrast
		cmp_mass = np.multiply(self.star_mass, mass_ratio)  # companion mass in solar masses

		f_mag = scipy.interpolate.interp1d(self.model['M/Ms'], self.model['Mag'], kind='cubic', fill_value='extrapolate')

		prim_model_mag = f_mag(self.star_mass)
		cmp_model_mag = [f_mag(x) if x >= self.model['M/Ms'][0] else float('inf') for x in cmp_mass]  # companion mags, assign infinite magnitude if below lowest modeled mass
		contrast = np.subtract(cmp_model_mag, prim_model_mag)

		# Determine luminosity
		f_lum = scipy.interpolate.interp1d(self.model['M/Ms'], self.model['L/Ls'], kind='cubic', fill_value='extrapolate')
		prim_lum = np.power(10, f_lum(self.star_mass))  # primary luminosity
		cmp_lum = [np.power(10, f_lum(x)) if x >= self.model['M/Ms'][0] else 0. for x in cmp_mass]  # companion luminosity, 0 if below limit

		# Choose to run either parallelized or non-parallelized based on the number of generated companions
		if num_generated < 50000:  # Serial
			self.predicted_RV = [np.zeros(len(self.MJD)) for x in range(num_generated)] # pre-allocate

			# calculate RV curves
			prim_K, prim_rv = self.calculate_RV(period, mass_ratio, a, e, cos_i, arg_peri, phase, self.MJD)
			cmp_K, cmp_rv = self.calculate_RV(period, np.divide(1, mass_ratio), a, e, cos_i, arg_peri, phase, self.MJD)
			cmp_rv = np.multiply(-1, cmp_rv)

			max_delta_rv = np.max(np.absolute(np.subtract(prim_rv, cmp_rv)), axis=1)

			# Determine the overall predicted RV
			for i in range(num_generated):
				if contrast[i] > 5:
					# SB1
					self.predicted_RV[i] = prim_rv[i]
				elif contrast[i] <= 5:
					# SB2, looks like SB1. RV is weighted average
					rv = np.average([prim_rv[i], cmp_rv[i]], axis=0, weights=[prim_lum, cmp_lum[i]])
					self.predicted_RV[i] = rv

			for i in range(0, num_generated):
				# Fit the zero point
				[zero_point], pcov = scipy.optimize.curve_fit(self.zero_point_model, self.predicted_RV[i], self.experimental_RV, sigma=self.measurement_error)
				# Shift all predicted values by the zero_point
				self.predicted_RV[i] += zero_point

			# Compare experimental and predicted RVs
			amp = [np.ptp(self.predicted_RV[i]) for i in range(num_generated)]
			chi_squared = [sum(np.divide(np.square(np.subtract(self.experimental_RV, self.predicted_RV[i])),
						np.add(np.square(self.measurement_error), self.added_jitter ** 2))) for i in range(num_generated)]
			# The degrees of freedom is equal to (N-1)+1, for the number of data points and the applied velocity shift
			prob = [stats.chi2.cdf(chi_squared[i], len(self.MJD)-1) for i in range(0, num_generated)]

		else:  # Parallel
			# Determine cpu count
			cpu_count = mp.cpu_count()-1

			contrast_check = [True if x <= 5 else False for x in contrast]

			if __name__ == "__main__":
				# Split the parameters into chunks to pass to workers
				divisor = int(np.ceil(min(num_generated / cpu_count, 200000)))
				n_divisor = int(np.ceil(num_generated / divisor))

				period = [period[i:i + divisor] for i in range(0, num_generated, divisor)]
				mass_ratio = [mass_ratio[i:i + divisor] for i in range(0, num_generated, divisor)]
				a = [a[i:i + divisor] for i in range(0, num_generated, divisor)]
				e = [e[i:i + divisor] for i in range(0, num_generated, divisor)]
				cos_i = [cos_i[i:i + divisor] for i in range(0, num_generated, divisor)]
				arg_peri = [arg_peri[i:i + divisor] for i in range(0, num_generated, divisor)]
				phase = [phase[i:i + divisor] for i in range(0, num_generated, divisor)]
				contrast_check = [contrast_check[i:i+divisor] for i  in range(0,  num_generated, divisor)]

				# Move into parallel processing
				#  Create Processes
				pool = mp.Pool(cpu_count)

				# Use Pool to calculate RVs
				prim_results = pool.starmap(self.calculate_RV, [(period[j], mass_ratio[j], a[j], e[j], cos_i[j], arg_peri[j], phase[j], self.MJD) for j in range(n_divisor)])
				cmp_results = pool.starmap(calculate_RV_parallel, [(period[j], np.divide(1, mass_ratio[j]), a[j],
						e[j], cos_i[j], arg_peri[j], phase[j], self.MJD, contrast_check[j]) for j in range(n_divisor)])

				# Concatenate Results
				prim_K = np.hstack([prim_results[i][0] for i in range(int(np.ceil(num_generated / divisor)))])
				prim_rv = np.vstack([prim_results[i][1] for i in range(int(np.ceil(num_generated / divisor)))])
				cmp_K = np.hstack([cmp_results[i][0] for i in range(int(np.ceil(num_generated / divisor)))])
				cmp_rv = np.vstack([cmp_results[i][1] for i in range(int(np.ceil(num_generated / divisor)))])
				cmp_rv = np.multiply(-1, cmp_rv)

				max_delta_rv = np.max(np.absolute(np.subtract(prim_rv, cmp_rv)), axis=1)

				# Determine the overall predicted RV
				self.predicted_RV = [np.zeros(len(self.MJD)) for x in range(num_generated)]  # pre-allocate

				for i in range(num_generated):
					if contrast[i] > 5:
						# SB1
						self.predicted_RV[i] = prim_rv[i]
					elif contrast[i] <= 5:
						# SB2, looks like SB1. RV is weighted average
						rv = np.average([prim_rv[i], cmp_rv[i]], axis=0, weights=[prim_lum, cmp_lum[i]])
						self.predicted_RV[i] = rv

				# Use Pool to calculate zero point
				split_RV  = [self.predicted_RV[i:i+divisor] for i  in range(0,  num_generated, divisor)]
				zero_points = pool.starmap(zero_point_fit_parallel, [(self.experimental_RV, self.measurement_error, split_RV[j]) for j in range(n_divisor)])
				zero_points = np.concatenate(zero_points, axis=0)
				pool.close()

				# Shift all by zero point
				self.predicted_RV = [np.add(self.predicted_RV[i], zero_points[i]) for i in range(num_generated)]

				# Compare experimental and predicted RVs
				amp = [np.ptp(self.predicted_RV[i]) for i in range(num_generated)]
				chi_squared = [sum(np.divide(np.square(np.subtract(self.experimental_RV, self.predicted_RV[i])),
							   np.add(np.square(self.measurement_error), self.added_jitter**2))) for i in range(num_generated)]
				prob = [stats.chi2.cdf(chi_squared[i], len(self.MJD)-1) for i in range(0, num_generated)]
			# End Parallelized

		# Reject things with a rejection probability greater than 0.997, corresponding to 3 sigma
		rv_fit_reject = np.array([True if np.random.rand() < x else False for x in prob])
		# Check amplitude and resolution
		above_amplitude = np.array([True if abs(x) > self.rv_floor else False for x in amp])
		self.amp = amp
		visible_sb2 = np.array([True if contrast[i] < 5 and max_delta_rv[i] > delta_v else False for i in range(num_generated)])
		self.b_type = ['              ']*num_generated
		for i in range(num_generated):
			if contrast[i] > 5:
				self.b_type[i] = 'SB1'
			elif contrast[i] <=5 and max_delta_rv[i] < delta_v:
				self.b_type[i] = 'Unresolved SB2'
			elif contrast[i] <=5 and max_delta_rv[i] > delta_v:
				self.b_type[i] = 'Resolved SB2'

		self.rv_reject_list = np.array([True if (rv_fit_reject[i] and above_amplitude[i]) or visible_sb2[i] else False for i in range(num_generated)])
		return self.rv_reject_list

	def calculate_jitter(self, predicted, experimental, error, threshold=0.00001):
		# Calculates the stellar jitter term necessary to produce a chi_squared value of 1
		# Inputs: Measured RV, Predicted RV, Measurement Error
		# Outputs: Jitter
		a = [experimental[i]-predicted[i] for i in range(0, len(experimental))]
		b = np.std(a)**2
		c = np.median(error)**2
		initial_guess = np.sqrt(b+c)
		prev_jitter = 0
		new_jitter = initial_guess
		itr = 0
		# Using Newton's method to iterate
		while abs(new_jitter - prev_jitter) > threshold and itr < 10000:
			# Change prev_jitter
			prev_jitter = new_jitter
			# Calculate f(j) and f'(j)
			chi_square = sum([(experimental[i] - predicted[i]) ** 2 / (error[i]**2 + prev_jitter**2) for i in range(0, len(error))])
			chi_square_prime = sum([-(experimental[i] - predicted[i])**2 / (error[i]**2 + prev_jitter**2)**2 * 2 * prev_jitter for i in range(0, len(error))])
			# Calculate new jitter
			new_jitter = prev_jitter - (chi_square - 1) / chi_square_prime
			itr += 1
		if abs(new_jitter - prev_jitter) < threshold:
			return new_jitter
		else:
			return -1

	@ staticmethod
	def calculate_RV(period, mass_ratio, a, e, cos_i, arg_peri, phase, MJD):
		# Calculates the RVs for each item when passed arrays of orbital parameters
		# Inputs: Arrays of Period, Mass Ratio, Semi-Major Axis, eccentricity, inclination, arg peri, phase, calculation times
		# Outputs: Velocity Semi-Amplitude (km/s), RVs at each time in MJD

		sin_i = np.sin(np.arccos(cos_i))

		n = len(period)
		RV = [[0.0 for i in range(len(MJD))] for j in range(n)]
		a_star = np.multiply(a, np.divide(mass_ratio, np.add(mass_ratio, 1)))
		K = np.multiply(np.divide((2 * np.pi), period),np.divide(np.multiply(a_star, sin_i), np.sqrt((1 - np.square(e)))))  # AU/days
		K = np.multiply(K, 1731.48)  # km/s

		for i in range(n):  # Iterate over companions
			for j in range(0, len(MJD)):  # Iterate over times
				# Find E
				M = 2 * np.pi * MJD[j] / period[i] - phase[i]
				prev_E = 0.0
				current_E = M

				while abs(current_E - prev_E) > 0.00001:
					prev_E = current_E
					current_E = M + e[i] * np.sin(prev_E)

				# Find true anomaly, f
				f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - e[i]) / (1 + e[i])))
				# Find predicted RV
				RV[i][j] = K[i] * (np.sin(arg_peri[i] + f) + e[i] * np.sin(arg_peri[i]))  # km/s

		return K, RV

	def read_in_rv(self):
		try:
			rv_in = open(self.rv_filename, 'r')
		except FileNotFoundError:
			return -31

		try:
			rv = Table.read(self.rv_filename, format='ascii', delimiter=' ', fast_reader=False)
			col_names = list(rv.columns)
			assert len(col_names) == 3
		except:
			try:
				rv = Table.read(self.rv_filename, format='ascii', delimiter='\t', fast_reader=False)
			except Exception as e:
				return -32
		# rename the columns so they match what i want
		col_names = list(rv.columns)
		if len(col_names) != 3:
			return -32
		rv.rename_column(col_names[0],'JD')
		rv.rename_column(col_names[1],'RV')
		rv.rename_column(col_names[2],'RVerr')

		rv.sort('JD')

		t_0 = 2400000.5
		self.MJD = rv['JD']-t_0
		self.experimental_RV = rv['RV']
		self.measurement_error = rv['RVerr']
		return 0

	def restore_defaults(self):
		self.rv_filename = ''
		self.added_jitter = 0.02  # km/s
		self.rv_floor = 0.02  # km/s
		self.rv_reject_list = []
		self.jitter_reject_list = []
		self.MJD = []
		self.experimental_RV = []
		self.measurement_error = []
		self.predicted_RV = []
		self.jitter = []
		gc.collect()

	def zero_point_model(self, prediction, a):
		# Alters the prediction by some zero point shift, a
		y = prediction + a
		return y

	def load_stellar_model(self, filter, star_age):
		# Read in file containing stellar model with the filter needed
		# RV always loads G filter
		if filter == 'J' or filter == 'H' or filter == 'K':  # 2MASS filters
			model_chart = {}
			BHAC_file = 'BHAC15_2MASS.txt'
			with open(BHAC_file, 'r') as content_file:
				content = content_file.read()
			tables = content.split(
				sep='\n-----------------------------------------------------------------------------------------------\n')
			tables = [x for x in tables if len(x) > 1]

			for table in tables:
				n = table.find('M')
				time_segment = table[0:n]
				table_segment = table[n:]
				age = float(time_segment[time_segment.find('=') + 1:])
				year_chart = Table.read(table_segment, format='ascii', fast_reader=False)
				if filter == 'J':
					year_chart = year_chart['M/Ms', 'Mj']
					year_chart.rename_column('Mj', 'Mag')
				elif filter == 'H':
					year_chart = year_chart['M/Ms', 'Mh']
					year_chart.rename_column('Mh', 'Mag')
				elif filter == 'K':
					year_chart = year_chart['M/Ms', 'Mk']
					year_chart.rename_column('Mk', 'Mag')
				model_chart[age] = year_chart

		elif filter == 'G' or filter == 'R' or filter == 'I':
			model_chart = {}
			BHAC_file = 'BHAC15_CFHT.txt'
			with open(BHAC_file, 'r') as content_file:
				content = content_file.read()
			tables = content.split(
				sep='\n----------------------------------------------------------------------------------------------------------------------------------------\n')
			tables = [x for x in tables if len(x) > 1]

			for table in tables:
				n = table.find('M')
				time_segment = table[0:n]
				table_segment = table[n:]
				age = float(time_segment[time_segment.find('=') + 1:])
				year_chart = Table.read(table_segment, format='ascii', fast_reader=False)
				if filter == 'G':
					year_chart = year_chart['M/Ms', 'G', 'L/Ls']
					year_chart.rename_column('G', 'Mag')
				elif filter == 'R':
					year_chart = year_chart['M/Ms', 'R', 'L/Ls']
					year_chart.rename_column('R', 'Mag')
				elif filter == 'I':
					year_chart = year_chart['M/Ms', 'I', 'L/Ls']
					year_chart.rename_column('I', 'Mag')
				model_chart[age] = year_chart

		ages = np.array(list(model_chart.keys()))
		#  Check if age is modeled, and if it is simply return that table
		if star_age in ages:
			return model_chart[star_age]
		# If the age is not included in the models, linearly interpolate the parameters from included ages
		#  Find ages above and below the desired age
		diff = [star_age - x for x in ages]
		low_age = ages[np.argmin([x if x > 0 else float('inf') for x in diff])]
		high_age = ages[np.argmin([abs(x) if x < 0 else float('inf') for x in diff])]
		young_chart = model_chart[low_age]
		old_chart = model_chart[high_age]

		#  Get masses
		common_mass = np.intersect1d(young_chart['M/Ms'], old_chart['M/Ms'])

		young_chart = young_chart[[x in common_mass for x in  young_chart['M/Ms']]]
		old_chart = old_chart[[x in common_mass for x in old_chart['M/Ms']]]

		new_model = Table()
		new_model['M/Ms'] = common_mass

		#  Interpolate
		for col in model_chart[low_age].colnames[1:]:
			col_list = []
			for i in range(len(common_mass)):
				# interpolate the new value of the parameter for each mass
				xs = [low_age, high_age]
				ys = [young_chart[col][i], old_chart[col][i]]
				f = scipy.interpolate.interp1d(xs, ys, kind='linear')
				col_list.append(f(star_age))
			# add to table
			new_model[col] = col_list

		return new_model


class RUWE:
	# class variables
	reject_list = []
	star_ra = ''
	star_dec = ''
	star_age = 5
	star_mass = 1
	num_generated = 0
	u0_dictionary = {}
	ln_ruwe = 0.
	binary_prob = 0.
	projected_sep = []

	# class functions
	def __init__(self, star_ra, star_dec, age, mass, companions):
		self.star_ra = star_ra
		self.star_dec = star_dec
		self.star_age = age
		self.star_mass = mass
		self.num_generated = companions.get_num()
		self.a = companions.get_a()
		self.period = companions.get_P()
		self.e = companions.get_ecc()
		self.phase = companions.get_phase()
		self.arg_peri = companions.get_arg_peri()
		self.cos_i = companions.get_cos_i()
		self.mass_ratio = companions.get_mass_ratio()

	def analyze(self):
		# a. Calculate projected separation
		pro_sep = [0.0 for i in range(self.num_generated)]
		T_0 = 2457388.5  # epoch 2016.0 in JD, corresponding to Gaia eDR3

		for i in range(self.num_generated):
			# 1. Calculate mean anomaly
			M = 2 * np.pi * T_0 / self.period[i] - self.phase[i]
			# 2. Calculate eccentric anomaly iteratively
			prev_E = 0.0
			current_E = M
			while abs(current_E - prev_E) > 0.00001:
				prev_E = current_E
				current_E = M + self.e[i] * np.sin(prev_E)
			# 3. Calculate true anomaly
			f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - self.e[i]) / (1 + self.e[i])))
			# 4. Calculate projected separation in AU
			# Added an absolute value around the cos(f) in the separation calculation
			alpha = f + self.arg_peri[i]
			sqtmp = np.sqrt(np.sin(alpha)**2+np.cos(alpha)**2 * self.cos_i[i]**2)
			pro_sep[i] = self.a[i] * (1-self.e[i]**2)/(1+self.e[i]*np.cos(f))*sqtmp  # AU
		self.projected_sep = pro_sep

		# b. Calculate distance
		star_distance = 1 / (self.parallax)  # kpc

		# c. Get Barraffe models, and set up the  interpolation to get magnitude
		model = self.load_stellar_model(self.star_age)
		f_mag =  scipy.interpolate.interp1d(model['M/Ms'], model['G'], fill_value='extrapolate')

		# d. Calculate contrast
		primary_mag = f_mag(self.star_mass)
		companion_mag = np.array([f_mag(self.star_mass*x) for x in self.mass_ratio])
		delta_g = np.subtract(companion_mag, primary_mag)
		self.delta_g = delta_g

		# e. Get ruwe
		ruwe = np.exp(self.ln_ruwe)
		log_ruwe = np.log10(ruwe)

		# f. Get predicted RUWE
		#    convert from mas to AU
		star_distance = star_distance * 2.063e+8 # AU
		self.ruwe_dist['Sep(AU)'] = [star_distance * np.tan(np.radians(x/(3.6e6))) for x in np.power(10, self.ruwe_dist['log(sep)']) ]

		#  2D interpolation functions for ruwe and sigma_ruwe
		x_edges = np.unique(np.array(self.ruwe_dist['Sep(AU)']))
		y_edges = np.unique(np.array(self.ruwe_dist['DeltaG']))
		z = np.reshape(np.array(self.ruwe_dist['log(RUWE)']), [len(y_edges), len(x_edges)])
		z_sigma = np.reshape(np.array(self.ruwe_dist['sigma_log(RUWE)']), [len(y_edges), len(x_edges)])

		f_ruwe = scipy.interpolate.interp2d(x_edges, y_edges, z)
		f_sigma = scipy.interpolate.interp2d(x_edges, y_edges, z_sigma)

		pred_log_ruwe = np.concatenate([f_ruwe(self.projected_sep[i], delta_g[i]) for i in range(self.num_generated)])
		self.predicted_ruwe = 10**pred_log_ruwe
		pred_sigma = np.concatenate([f_sigma(self.projected_sep[i], delta_g[i]) for i in range(self.num_generated)])

		# g. Determine rejection probabilities
		rejection_prob = [stats.halfnorm.cdf(10**pred_log_ruwe[i], loc=10**log_ruwe, scale=10**pred_sigma[i]) for i in range(self.num_generated)]

		# never reject something where the observed ruwe is higher than the predicted (halfnorm)
		# never reject something that is outside the RUWE Distribution grid
		rejection_prob = [0.0 if (not -.1 < delta_g[i] < 7.1) or
						 (not np.min(self.ruwe_dist['Sep(AU)']) < self.projected_sep[i] < np.max(self.ruwe_dist['Sep(AU)']))
						 else rejection_prob[i] for i in range(self.num_generated)]
		self.rejection_prob = rejection_prob

		rejection = np.random.rand(self.num_generated)

		reject_list = [True if rejection[i] < rejection_prob[i] else False for i in range(self.num_generated)]

		return reject_list

	def load_stellar_model(self, star_age):
		# Read in file containing stellar model
		# todo Interpolate to get the chart for the exact age or binned age or something that doesnt rely on the age being in the chart
		model_chart = {}
		BHAC_file = 'BHAC15_CFHT.txt'
		with open(BHAC_file, 'r') as content_file:
			content = content_file.read()
		tables = content.split(
			sep='\n----------------------------------------------------------------------------------------------------------------------------------------\n')
		tables = [x for x in tables if len(x) > 1]

		for table in tables:
			n = table.find('M')
			time_segment = table[0:n]
			table_segment = table[n:]
			age = float(time_segment[time_segment.find('=') + 1:])
			year_chart = Table.read(table_segment, format='ascii', fast_reader=False)
			year_chart['Age(Gyr)'] = age
			model_chart[age] = year_chart

		# Set up interpolation
		ages = np.array(list(model_chart.keys()))
		#  check if age is modeled
		if star_age in ages:
			return model_chart[star_age]
		#  find ages above and below the desired age
		diff = [star_age-x for x in ages]
		low_age = ages[np.argmin([x if x > 0 else float('inf') for x in diff])]
		high_age = ages[np.argmin([abs(x) if x < 0 else float('inf') for x in diff])]
		#  get masses
		new_model = Table()
		new_model['M/Ms'] = model_chart[high_age]['M/Ms']
		for y in model_chart[low_age].colnames[1:]:
			y_list = []
			for m in list(new_model['M/Ms']):
				# interpolate
				low_i = list(model_chart[low_age]['M/Ms']).index(m)
				high_i = list(model_chart[high_age]['M/Ms']).index(m)
				xs = [low_age, high_age]
				ys = [model_chart[low_age][y][low_i], model_chart[high_age][y][high_i]]
				f = scipy.interpolate.interp1d(xs, ys, kind='linear')
				y_list.append(f(star_age))
			#  add to table
			new_model[y] = y_list

		return new_model

	def get_gaia_info(self):
		#   Query Gaia for parallax, g_mag, color, astrometric_chi2 and n_good_obs, calculate RUWE
		coordinate = SkyCoord(self.star_ra, self.star_dec, frame='icrs')
		width = u.Quantity(10, u.arcsecond)
		height = u.Quantity(10, u.arcsecond)
		job_str = ("SELECT TOP 10 DISTANCE(POINT('ICRS', %f, %f), POINT('ICRS', ra, dec)) AS dist, * FROM gaiaedr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', %f, %f),CIRCLE('ICRS', ra, dec, 0.08333333)) ORDER BY dist ASC)" % (
			coordinate.ra.degree, coordinate.dec.degree, coordinate.ra.degree, coordinate.dec.degree))
		job = Gaia.launch_job(job_str)
		gaia_info = job.get_results()
		if gaia_info and len(gaia_info) >= 1:
			# Only one star fits the coordinates, all is well
			self.gmag = gaia_info['phot_g_mean_mag'][0]
			self.color = gaia_info['bp_rp'][0]
			self.n_good_obs = gaia_info['astrometric_n_good_obs_al'][0]
			self.astrometric_chi2 = gaia_info['astrometric_chi2_al'][0]
			self.parallax = gaia_info['parallax'][0]
			self.parallax_error = gaia_info['parallax_error'][0]
			self.ln_ruwe = np.log(gaia_info['ruwe'][0])

			#    Check if magnitude, color and Gaia solution are valid for calculating RUWE
			if 3.6 <= self.gmag <= 21. and -1 <= self.color <= 10 and gaia_info['astrometric_params_solved'][0] == 31:
				# All is well
				return
			else:
				# The magnitude or color is outside of bounds
				return -53
		else:
			# No Gaia results
			return -51

	def read_dist(self):
		file_name = 'RuweTableGP.txt'
		t = Table.read(file_name, format='ascii', delimiter=' ')

		self.ruwe_dist = t
		return


if __name__ == '__main__':
	app = Application(sys.argv)
	app.start()
