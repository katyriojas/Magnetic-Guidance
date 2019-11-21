MATLAB folder for RA-L Experiments 
Last Updated: 11/21/19
    0) Preop plan is generated with MagneticGuidance.m

Scripts for generating post-experiment plots (these are called in plot scripts - don't need to run independently)
    1) LoadRALData_Manual.m
    2) LoadRALData_Phantom_Robotic.m
    3) LoadRALData_Cadaver_Robotic.m
Next we plot and compare forces for different methods/trials
*all of these scripts have default to call scripts 1-3, so should not need to run those
    4) MagneticGuidanceCompareForces_UG_G_phantom.m (will take a while to run)
        - Plots Fmag vs. LID all trials (UG-EA,UG-MEA,G)
        - Plots Fmag vs. AID all trials (UG-EA,UG-MEA,G)
    5) MagneticGuidanceCompareForces_UG_G_cadaver.m
        - Plots Fmag vs. LID all 6 trials (3 UG and 3 G)
            - also plots the trim points generated in LoadRALData_Cadaver_Robotic)
    6) MagneticGuidanceCompareForces_manual.m
        - Plots Fmag vs. time all trials as subplots
More plotting scripts/data manipulation (these default to use the saved matrices)
    7) MagneticGuidance_Cadaver_Avgs.m
        - Plots average cadaver data
    8) MagneticGuidanceCompareForces_BinnedResults.m
        - TODO: Trim robotic data?
    9) MagneticGuidanceColorbars.m
        - Generates Individual trial Fmag vs. AID for phantom  
    10) CDF plot
        - TODO: change phantom data to all stop within 10 degrees of end (need to trim robotic data probably?)
    
