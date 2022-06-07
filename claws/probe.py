#!/usr/bin/env python
"""Abstract class defining a probe.
"""

# Import modules
import re
import numpy as np
import matplotlib.pyplot as plt

from claws.helpers import indent
from claws.custom_exceptions import InputError
from claws.claws import output_options, \
    get_unit_factors
import claws.opendrift_wrapper as wrapper
import claws.postpro as postpro

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "hystrath@gmail.com"
__status__ = "Production"


# ---------------------------------------------------------------------------- #
# Classes 
# ---------------------------------------------------------------------------- #

class Probe(object):
    def __init__(self, Station):
        # __station_obj: Station object; Monitoring station location
        self.__station_obj = Station
        # __is_inside_domain: bool; Whether or not a probe is indide the
        # 'concentration domain'. default is False/outside
        self.__is_inside_domain = False
        # __lonbin, __latbin: integer; bin indices longitudinally and
        # latitudinally. default is -1 (not set)
        self.__lonbin = -1
        self.__latbin = -1
        
        self._sanitize()
        
    def __str__(self, indent_lvl=1):
        return """{1}\n{0}Location = {2}""".format(indent(indent_lvl),
            self.name(), self.__station_obj.__str__(2))
    
    def Station(self):
        return self.__station_obj
        
    def _sanitize(self):
        pass
                
    def name(self):
        return ' '.join(re.findall('([A-Z][a-z]+)', type(self).__name__))
        
    def check_validity(self, corners, concentration):
        """Determine whether or not a probe is inside the 'concentration
        domain'
        
        Arguments:
            corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
                extent of the domain/histogram
                
            concentration: list of Xarray histograms; as generated by Opendrift
                get_histogram function from the data stored in the .nc outfile
                and the get_quadtree_histograms function   
        """
        nlonbins, nlatbins = np.shape(concentration[0])[1:3]
        lonext = abs(corners[1] - corners[0])
        latext = abs(corners[3] - corners[2])
        
        # Shorthands
        lon = self.__station_obj.lon()
        lat = self.__station_obj.lat()
        
        # Extract probe longitude and latitude bin indices
        lonbin = int((lon - corners[0])/lonext*nlonbins)
        latbin = int((lat - corners[2])/latext*nlatbins)
        
        # If any of the following assertions fail, the probe is outside the
        # 'concentration' domain
        try:
            assert(lonbin >= 0)
            assert(latbin >= 0)
            assert(lonbin < nlonbins)
            assert(latbin < nlatbins)
            
            self.__lonbin = lonbin
            self.__latbin = latbin
            self.__is_inside_domain = True
        except AssertionError:
            print("{}\n\tis outside the concentration domain and will " \
                      "not be considered".format(self.__str__()))
            
    def record_history(self, corners, concentration, quadtree_obj,
                       quadtree_conc_lvl):
        """Store concentration history for all root and leaf bins, for all
        output times
        
        Arguments:
            corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
                extent of the domain/histogram
            
            concentration: list of Xarray histograms; as generated by Opendrift
                get_histogram function from the data stored in the .nc outfile
                and the get_quadtree_histograms function
                
            quadtree_obj: Quadtree object as implemented in
                claws.quadtree
                
            quadtree_conc_lvl: nested python list; output of the
                get_quadtree_histograms function. Indicates the number of
                sub-histograms per output time bins and their depth level.
        """
        ndt = np.shape(concentration[0])[0]
        ntb = len(output_options["time_bins"])
        
        self.__history = np.zeros(shape=(ndt))
        self.__quadtree_history = np.zeros(shape=(ntb))
        
        # Skip probe if outside the 'concentration domain'
        if not self.__is_inside_domain: return
        
        # Store concentration time series
        history = np.array(wrapper.get_series_probe(concentration[0],
                              self.__lonbin, self.__latbin))
        try:
            if quadtree_obj.is_active():
                # Initialise quadtree history array to history array taken at
                # output times
                self.__quadtree_history = history[output_options["time_bins"]]
                # Search if there are any negative concentration values denoting
                # that the probe bin presents a sub-histogram
                if np.shape(np.where(history < 0.0))[1] > 0:
                    # Loop over all time bins where a negative concentration was
                    # found for the probe root bin
                    for time_pos in np.where(history < 0.0)[0]:
                        tb_pos = np.where(time_pos ==
                                          output_options["time_bins"])[0][0]
                        offset = np.sum([len(quadtree_conc_lvl[t]) for t in
                                            range(tb_pos+1)])
                        
                        # Loop over all subhistograms created for this time bin
                        for c in np.arange(offset, offset + 
                                len(quadtree_conc_lvl[tb_pos+1])):
                            # Determine is the probe is located in this
                            # sub-histogram. If so, return a finite value
                            subconc = wrapper.is_point_in_subhisto_bbox((
                                self.__lonbin, self.__latbin), concentration[c])
                            if np.isfinite(subconc):
                                # If the finite value returned is positive, the
                                # 'leaf' quadtree level has been reached for 
                                # this probe root bin. Otherwise, keep searching
                                if subconc >= 0.0:
                                    self.__quadtree_history[tb_pos] = subconc
                                    history[time_pos] = np.abs(history[time_pos])
                                    break
                        
                        # If the concentration remains negative at the end of
                        # the search, a problem occurred            
                        if history[time_pos] < 0:
                            raise ValueError
        except ValueError:
            print("Sub-histogram for probe {} not found".format(self.__str__))
            exit()
        finally:
            self.__history = history
            
    def plot_concentration(self, wf, time, time_last_treatment, chemical_obj,
                           pixelsize_m, dz, quadtree_obj, pidx,
                           last_treatment_label='Last treatment',
                           ylabel='Concentration', filename='probe'):
        """Plot concentration time series at this station for a given chemical
        
        Arguments:
            wf: string; working folder
            
            time: numpy array; time series
            
            time_last_treatment: float; time of latest treatment over all
                seeding locations (units: output_options["time_units"])
                
            chemical_obj: ChemicalSubstance object as implemented in
                claws.ChemicalSubstance
            
            pixelsize_m: float; size of the bins along the longitudinal and
                latitudinal directions, in meters 
                
            dz: float; depth bin size in meters.       
            
            quadtree_obj: Quadtree object as implemented in
                claws.quadtree 
            
            pidx: integer; probe index in the listf of probe objects
            
            last_treatment_label: string; message to print on the graph at time
                of last treatment. Default is 'Last treatment'
            
            ylabel: string; label of the ordinate axis without units. default is
                'Concentration'
                
            filename: string; plot filename without the path nor the file
                extension. default is 'probes'
        """
        # Skip probe if outside the 'concentration domain'
        if not self.__is_inside_domain: return
        
        # Conversion factors from program units to user-defined output units
        cf, lf, tf = get_unit_factors()
        
        # Plot figure
        fig, ax = plt.subplots(1)
        psize_str = ('%f'%pixelsize_m).rstrip('0').rstrip('.')
        dz_str = ('%f'%dz).rstrip('0').rstrip('.')
        plt.plot(time, self.__history*cf, color='k', linestyle='-',
                 label='{0} m x {0} m x {1} m'.format(psize_str, dz_str), lw=1)
        if quadtree_obj.is_active():         
            plt.scatter(output_options["time_bins"]*tf,
                        self.__quadtree_history*cf, color='blue', marker='+',
                        label='quadtree (max depth: {})'.format(
                            quadtree_obj.max_quadtree_depth())) 
                    
        if np.isfinite(chemical_obj.eqs_3hour()):
            plt.axhline(y=chemical_obj.eqs_3hour()*cf, color='green',
                        linestyle='dotted', label='3-hour EQS')
        if np.isfinite(chemical_obj.eqs_72hour()):
            plt.axhline(y=chemical_obj.eqs_72hour()*cf, color='green',
                        linestyle='--', label='72-hour EQS')
        if np.isfinite(chemical_obj.mac_72hour()):
            plt.axhline(y=chemical_obj.mac_72hour()*cf, color='red',
                        linestyle='--', label='72-hour MAC')
        
        ax.set_yscale('log')
        
        # Adjust y-axis bounds
        min_val = -np.inf
        if np.shape(np.nonzero(self.__history))[1] > 0:
            min_val = np.min(self.__history[np.nonzero(self.__history)]*cf)
        postpro.round_logscale_yaxis(ax, set_min=min_val)
        
        postpro.print_last_treatment(plt, ax, time_last_treatment,
                                     last_treatment_label)
        plt.tick_params(axis='y', which='minor')
        ax.yaxis.set_minor_formatter('')
        plt.xlabel('Time {}'.format(output_options["time_units"]))
        plt.ylabel('{} {}'.format(ylabel,
                                  output_options["concentration_units"]))
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                   mode="expand", borderaxespad=0, ncol=2)
        plt.savefig(wf + filename + '_' + str(pidx+1) + '.png')
        plt.close(fig)
