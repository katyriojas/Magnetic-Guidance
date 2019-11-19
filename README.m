Preop plan is generated with MagneticGuidance.m

Scripts for generating post-experiment plots:
First we load in all of the data:
    1) LoadRALData_Manual.m
    2) LoadRALData_Phantom_Robotic.m
    3) LoadRALData_Cadaver_Robotic.m
Next we plot and compare forces for different methods/trials:
    4) MagneticGuidanceCompareForces_UG_G_phantom.m
        - Plots Fmag vs. LID all trials (UG-EA,UG-MEA,G)
        - Plots Fmag vs. AID all trials (UG-EA,UG-MEA,G)
    5) MagneticGuidanceCompareForces_UG_G_cadaver.m
        - Plots Fmag vs. LID all 6 trials (3 UG and 3 G)
            - also plots the trim points generated in LoadRALData_Cadaver_Robotic)
    6) MagneticGuidance_Cadaver_Avgs.m
        - Plots average cadaver data
    7) MagneticGuidanceCompareForces_manual.m
        - Plots Fmag vs. time all trials as subplots
    8) MagneticGuidanceCompareForces_BinnedResults.m (have to run this before running colorbars)
    9) MagneticGuidanceColorbars.m
        - Generates Individual trial Fmag vs. AID for phantom    
