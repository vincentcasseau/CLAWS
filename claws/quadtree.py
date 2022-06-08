#!/usr/bin/env python
"""Class to define a quadtree histogram structure.

Use the derived class NoQuadtree
"""

# Import modules
import re
import numpy as np

from claws.custom_exceptions import InputError
from claws.helpers import *
from claws.claws import output_options
import claws.opendrift_wrapper as wrapper

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"


# ---------------------------------------------------------------------------- #
# Variables 
# ---------------------------------------------------------------------------- #

# _quadtree_max_tolerated_depth: integer; maximum quadtree depth permitted.
# Quadtree.__max_depth cannot exceed this value. default is 4
_quadtree_max_tolerated_depth = 4


# ---------------------------------------------------------------------------- #
# Classes 
# ---------------------------------------------------------------------------- #

class Quadtree(object):
    def __init__(self, is_active=None, max_quadtree_depth=0,
                 min_particles_per_bin=1000, variable_root_bin_width=False,
                 concentration_target=1.0, leaf_bin_width=np.nan,
                 variable_leaf_bin_width=False, leaf_to_seeding_area_ratio=5.0,
                 input_conc_units='ng/L'):
        # is_active: bool; Whether or not to create a quadtree from the root
        # histogram. default is None (ie, not set)
        self.__is_active = is_active
        # max_depth: integer; Maximum quadtree depth. default is 0
        # (ie, no quadtree)
        self.__max_depth = max_quadtree_depth
        # min_particles_per_bin: integer; Number of particles per bin that
        # triggers the creation of a child histogram within the bounds of a
        # parent histogram. default is 1000
        self.__min_particles_per_bin = min_particles_per_bin
        # variable_root_bin_width: bool; Whether or not the root bin width can
        # be edited automatically using a target concentration. default is False
        self.__variable_root_bin_width = variable_root_bin_width
        # concentration_target: float; Concentration target in ng/L to recompute
        # the root bin width. default is 1 ng/L
        self.__concentration_target = concentration_target
        # leaf_bin_width: float; Set target value for the leaf bin width
        # (meters). default is NaN
        self.__leaf_bin_width = leaf_bin_width,
        # variable_leaf_bin_width: bool; Whether or not the leaf bin width can
        # be edited automatically using the ratio of the leaf bin area to the
        # seeding area. default is False
        self.__variable_leaf_bin_width = variable_leaf_bin_width
        # leaf_to_seeding_area_ratio: float; ratio of the leaf bin area to the
        # seeding area of a reference seeding location. default is 5.0
        self.__leaf_to_seeding_area_ratio = leaf_to_seeding_area_ratio
        # method: integer; method used to create the quadtree structure
        self.__method = None
        
        self._sanitize(input_conc_units)
        
    def __str__(self):
        if not __is_active:
            print("Quadtree structure\n\tActive = False")
        else:    
            print("""Quadtree structure\n\tActive = True\n\Max depth = {}"\
                "\n\tMin particles per bin = {}""".format(self.__max_depth,
                self.__min_particles_per_bin))
                
    def is_active(self):
        return self.__is_active
      
    def max_quadtree_depth(self):
        return self.__max_depth
      
    def min_particles_per_bin(self):
        return self.__min_particles_per_bin
      
    def variable_root_bin_width(self):
        return self.__variable_root_bin_width
      
    def concentration_target(self):
        return self.__concentration_target
      
    def leaf_bin_width(self):
        return self.__leaf_bin_width
      
    def variable_leaf_bin_width(self):
        return self.__variable_leaf_bin_width
      
    def leaf_to_seeding_area_ratio(self):
        return self.__leaf_to_seeding_area_ratio
        
    def _sanitize(self, input_conc_units):
        "Sanitise quadtree inputs"
        
        # Determine if a quadtree structure is used
        if self.__is_active is not None:
            # "is_active" used to decide whether a quadtree structure exists
            assert(type(self.__is_active) is bool)
            if not self.__is_active:
                self.__max_depth = 0
        else:
            # "max_depth" used to decide whether a quadtree structure exists
            assert(type(self.__max_depth) is int)
            if self.__max_depth < 0:
                raise InputError(self.__max_depth,
                    'Quadtree max depth should be a positive integer')
            elif self.__max_depth == 0:
                self.__is_active = False
            else:
                self.__is_active = True
                
        if not self.__is_active: return
        
        # method 0: fix the root bin width defined by the user, and use the 
        #           minimum number of particles per bin and maximum quadtree
        #           depth for bin subdivisions
        # method 1: fix the root bin width defined by the user, derive the 
        #           maximum quadtree depth from the leaf bin width, and use the
        #           minimum number of particles per bin and maximum quadtree
        #           depth for bin subdivisions (as in method 0)
        # method 2: recompute the root bin width using the concentration target,
        #           and use the minimum number of particles per bin and maximum
        #           quadtree depth for bin subdivisions (as in method 0)
        # method 3: recompute the root bin width using the concentration target
        #           (as in method 2), derive the maximum quadtree depth from the
        #           leaf bin width, and use the minimum number of particles per
        #           bin and maximum quadtree depth for bin subdivisions (as in
        #           method 1)
        # method 4: recompute the root bin width using the concentration target
        #           (as in method 2), derive the leaf bin width from the 
        #           leaf-to-seeding area ratio, the maximum quadtree depth from
        #           the leaf bin width. Fianlly, use the minimum number of
        #           particles per bin and maximum quadtree depth for bin
        #           subdivisions (as in method 1)
        
        self.__method = 0
        
        assert(type(self.__variable_root_bin_width) is bool)
        if self.__variable_root_bin_width:
            assert(type(self.__variable_leaf_bin_width) is bool)
            if self.__variable_leaf_bin_width:
                self.__method = 4
            elif np.isfinite(self.__leaf_bin_width):
                assert(type(self.__leaf_bin_width) is float)
                if self.__leaf_bin_width <= 0.0:
                    raise InputError(self.__leaf_bin_width,
                        'Leaf bin width should be a positive number')
                self.__method = 3
            else:
                self.__method = 2
        elif np.isfinite(self.__leaf_bin_width):
            assert(type(self.__leaf_bin_width) is float)
            if self.__leaf_bin_width <= 0.0:
                raise InputError(self.__leaf_bin_width,
                    'Leaf bin width should be a positive number')
            self.__method = 1
        
        # Concentration target stored in ng/L
        self.__concentration_target = convert_conc_to_prog_units(
            self.__concentration_target, input_conc_units, self.name())
        
    def name(self):
        return ' '.join(re.findall('([A-Z][a-z]+)', type(self).__name__))
        
    def update(self, pixelsize_m, bin_depth, mass_particle, seeding_radius):
        """Apply quadtree method defined in the sanitise function and return the
        new pixelsize_m if necessary

        Arguments:
            pixelsize_m: float; size of the bins along the longitudinal and
                latitudinal directions, in meters
                
            bin_depth: float; depth bin size in meters
            
            mass_particle: float; mass of a single particle, in nanograms
            
            seeding_radius: float; representative seeding radius in meters
        """
        pixelsize_m = float(pixelsize_m)
        
        # Don't do anything if no quadtree structure
        if not self.__is_active:
            return pixelsize_m
            
        if self.__method in [2, 3, 4]:
            # Recompute root bin width, pixelsize_m, using the concentration
            # target
            # voxel volume = mass particle / target concentration
            # root bin width = sqrt(mass particle / 
            #                                (target concentration * bin depth))
            pixelsize_m = np.sqrt(mass_particle/(self.__concentration_target
                                      *bin_depth))
            pixelsize_m = round(pixelsize_m, 3)
            
        # Create an internal copy of pixelsize_m
        self.__root_bin_width = pixelsize_m

        if self.__method in [1, 3, 4]:
            if self.__method == 4:
                # Compute the leaf bin width from the leaf to seeding area ratio
                if type(seeding_radius) in [list, np.ndarray]:
                    seeding_radius = seeding_radius[0]
                assert(type(seeding_radius) in [int, float])    
                seeding_area = np.pi*seeding_radius**2
                self.__leaf_bin_width = np.sqrt(
                    self.__leaf_to_seeding_area_ratio*seeding_area)
            leaf_pixelsize_m = self.__leaf_bin_width + 1.e-6
            
            # Compute the maximum depth from the root and leaf bin widths
            try:
                qd_depth_float = np.log2(pixelsize_m/leaf_pixelsize_m)
                self.__max_depth = np.ceil(qd_depth_float).astype('int')
                self.__max_depth = min(self.__max_depth, 
                                                _quadtree_max_tolerated_depth)
                                                
                pixelsize_m = leaf_pixelsize_m * 2**self.__max_depth
                pixelsize_m = round(pixelsize_m, 3)
                self.__root_bin_width = pixelsize_m                            
            except ZeroDivisionError:
                self.__max_depth = 0
            finally:
                if self.__max_depth == 0:
                    self.__is_active = False
                    
        if self.__method > 1:
            print(self.__str__)
            print("Please review the quadtree structure info ...")
        return pixelsize_m
        
    def get_histograms(self, nsimulations, outfile, root_histo,
                       root_pixelsize_m, **kwargs):
        """Return multidimensional histograms binning particles into voxels,
        irrespective of their seeding source
        Each bin of an input root histogram is scanned and if the local number
        of particles exceed the defined minimum number of particles per bin (*),
        'min_particles_per_bin', another histogram is created with the
        following properties:
            - its extent is that of the parent histogram bin
            - the size of its bins is half of the parent histogram
        The search is stopped when no more bins satifify the condition (*) or
        if the search depth is larger than the maximum search depth allowed,
        'max_quadtree_depth'    
        
        Arguments:
            nsimulations: integer; number of ensembles
            
            outfile: string; full path to the Opendrift .nc outfile(s), minus
                its/their .nc extension
            
            root_histo: Xarray; background histogram as generated by Opendrift
                get_histogram function from the data stored in the .nc outfile
                
            root_pixelsize_m: float; size of the root histogram bins along the
                longitudinal and latitudinal directions, in meters
                
            **kwargs:    
                status: integer; status of the particles to consider. default
                     is 0, that is active particles only. If status is None, all
                     particles are considered irrespective of their status
                    
                z0: float; elevation in meters at the top of the depth bin.
                    default is 0
                
                dz: float; depth bin size in meters. default is 5
                
                marker_index: integer; marker index to identify different
                    chemical substances. default: None (summing over all seeding
                    sources)
        """
        time_bins = output_options["time_bins"]
        
        # Put the root histogram (level 0) in the list of histograms
        histos_lvl = [[0]]
        histos = [root_histo]
        
        if self.__max_depth == 0:
            return [histos_lvl, histos] 
        
        for t in range(len(time_bins)):
            print("Creating sub-histograms for time bin", t+1, "/", len(time_bins))
            # Initialise loop variables
            if nsimulations > 1:
                of = outfile + '_1.nc'
            else:
                of = outfile + '.nc'
            
            current_level = 0
            histos_lvl_ptr = [0]
            histos_lvl.append([])
            cmpt = 0
            
            while current_level < self.__max_depth:
            
                if current_level == 0:
                    h_idx_beg = 0
                else:
                    h_idx_beg = histos_lvl_ptr[current_level-1] + 1
                
                # Loop through all histograms in current 'quadtree level'
                for hi in range(h_idx_beg, histos_lvl_ptr[current_level] + 1):
                    h = histos[hi]
                    lenh = len(np.shape(h))
                    if lenh == 3:
                        # dims are (time / lon_bin / lat_bin)
                        nlonbins, nlatbins = np.shape(h)[1:3]
                    else:
                        # dims are (lon_bin / lat_bin)
                        nlonbins, nlatbins = np.shape(h)[:2]
                        
                    h_dlon = h.lon_bin[1].item() - h.lon_bin[0].item()
                    h_dlat = h.lat_bin[1].item() - h.lat_bin[0].item()
                    
                    # Scan through all bins
                    for lati in range(nlatbins):
                        for lonj in range(nlonbins):
                            # Check if there are enough particles in bin (lonj,lati)
                            if lenh == 3:
                                if h[time_bins[t],lonj,lati] < self.__min_particles_per_bin:
                                    # not enough particles in root histogram bin,
                                    # cycle
                                    continue
                            else:
                                if h[lonj,lati] < self.__min_particles_per_bin:
                                    # not enough particles in sub-histogram bin,
                                    # cycle
                                    continue
                            
                            # Create child histogram
                            cmpt += 1
                            # Set pixelsize_m to None to split a pixel into
                            # four even pixels using linspace
                            pixelsize_m = None
                            corners = [h.lon_bin[lonj].item(), h_dlon/2.,
                                       h.lat_bin[lati].item(), h_dlat/2.]
                            new_histo = wrapper.get_histogram(
                                outfile=of, pixelsize_m=pixelsize_m,
                                corners=corners, time_bin=time_bins[t],
                                **kwargs)
                             
                            # Verify that the sub-histogram creation was successful
                            assert((h.lon_bin[lonj].item() - h_dlon/4.) - new_histo.lon_bin[0].item() < 1e-10)
                            assert((h.lat_bin[lati].item() - h_dlat/4.) - new_histo.lat_bin[0].item() < 1e-10)
                            assert((h.lon_bin[lonj].item() + h_dlon/4.) - new_histo.lon_bin[1].item() < 1e-10)
                            assert((h.lat_bin[lati].item() + h_dlat/4.) - new_histo.lat_bin[1].item() < 1e-10)

                            # Average histogram values over all ensembles and
                            # seeding locations
                            if nsimulations > 1:
                                new_histo = average_ensembles(
                                    h=new_histo,
                                    nsimulations=nsimulations,
                                    outfile=outfile,
                                    pixelsize_m=pixelsize_m,
                                    corners=corners,
                                    time_bin=time_bins[t], **kwargs)                          
                                                      
                            histos.append(new_histo)
                            histos_lvl[-1].append(current_level+1)
                            # Flip the sign of the histogram bins that present
                            # sub-histograms to keep track of them without
                            # storing additional data
                            # For example, this is useful for filtering out
                            # certain bins that present sub-histograms when
                            # calculating the area > EQS
                            if lenh == 3:
                                histos[hi][time_bins[t],lonj,lati] *= -1.
                            else:
                                histos[hi][lonj,lati] *= -1.
                            
                # If the histogram counter is equal to histos_lvl_ptr[-1], no new
                # histograms have been created, stop the search for this time bin
                if histos_lvl_ptr[-1] == cmpt:
                    break
                else:
                    histos_lvl_ptr.append(cmpt)
                    # Increment level counter
                    current_level += 1
            
        # Return the list storing the 'quadtree level' of each histogram and the
        # list of histograms
        return [histos_lvl, histos]
