#!/usr/bin/env python
"""Defines functions looping over all probe objects.
"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "hystrath@gmail.com"
__status__ = "Production"


# ---------------------------------------------------------------------------- #
# Functions 
# ---------------------------------------------------------------------------- #

def check_probes_validity(probes, corners, concentration):
    """Determine whether or not probes are inside the 'concentration
    domain'
    
    Arguments:
        probes: list of Probe objects as implemented in
            claws.Probe
        
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            extent of the domain/histogram
            
        concentration: list of Xarray histograms; as generated by Opendrift
            get_histogram function from the data stored in the .nc outfile
            and the get_quadtree_histograms function   
    """
    for P in probes:
        P.check_validity(corners, concentration)
        
def record_probes_history(probes, corners, concentration, quadtree_obj,
                          quadtree_conc_lvl):
    """Store concentration history of all probes for all root and leaf bins,
    and for all output times
    
    Arguments:
        probes: list of Probe objects as implemented in
            claws.Probe
            
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
    for P in probes:
        P.record_history(corners, concentration, quadtree_obj,
                         quadtree_conc_lvl)
        
def plot_probes_concentration(probes, wf, time, time_last_treatment,
                              chemical_obj, pixelsize_m, dz, quadtree_obj,
                              last_treatment_label='Last treatment',
                              ylabel='Concentration', filename='probe'):
    """Plot concentration time series at this station for a given chemical
        
    Arguments:
        probes: list of Probe objects as implemented in
            claws.Probe
            
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
        
        last_treatment_label: string; message to print on the graph at time
            of last treatment. Default is 'Last treatment'
        
        ylabel: string; label of the ordinate axis without units. default is
            'Concentration'
            
        filename: string; plot filename without the path nor the file
            extension. default is 'probes'
        """
    for pidx, P in enumerate(probes):
        P.plot_concentration(wf, time, time_last_treatment, chemical_obj,
                             pixelsize_m, dz, quadtree_obj, pidx,
                             last_treatment_label, ylabel, filename)
