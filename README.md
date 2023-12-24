[![Citation Badge](https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.3758%2Fs13428-016-0822-1&color=blue)](https://scholar.google.com/citations?view_op=view_citation&citation_for_view=uRUYoVgAAAAJ:KlAtU1dfN6UC)
[![DOI](https://zenodo.org/badge/DOI/10.3758/s13428-016-0822-1.svg)](https://doi.org/10.3758/s13428-016-0822-1)
# I2MC - MATLAB implementation

## About I2MC
The I2MC algorithm was designed to accomplish fixation classification in data across a wide range of noise levels and when periods of data loss may occur.

Cite as:
[Hessels, R.S., Niehorster, D.C., Kemner, C., & Hooge, I.T.C. (2017). Noise-robust fixation detection in eye-movement data - Identification by 2-means clustering (I2MC). Behavior Research Methods, 49(5): 1802--1823. doi: 10.3758/s13428-016-0822-1](https://link.springer.com/article/10.3758/s13428-016-0822-1)

For more information, questions, or to check whether we have updated to a better version, e-mail: royhessels@gmail.com / dcnieho@gmail.com. I2MC is available from www.github.com/royhessels/I2MC

A Python implementation of the I2MC algorithm is available at https://github.com/dcnieho/I2MC_Python, please ensure to read the readme before using it.

Most parts of the I2MC algorithm are licensed under the Creative Commons Attribution 4.0 (CC BY 4.0) license. Some functions are under MIT license, and some may be under other licenses.

## How to use
Quick start guide for adopting this script for your own data:
1. Build an import function specific for your data (see [importTobiiTX300.m](/functions/import/importTobiiTX300.m) for an example). 

2. Change line 141 of [I2MC.m](/I2MC.m#L141) to use your new import function. The format should be:
    `data.time` for the timestamp
    
    `data.left.X` & `data.left.Y` for left gaze coordinates
    
    `data.right.X` & `data.right.Y` for right gaze coordinates
    
    `data.average.X` & `data.average.Y` for average of right and left gaze coordinates
    
    You may provide coordinates from both eyes, only the left, only the right, or only the average. 
    Gaze coordinates should be in pixels, timestamps should be in milliseconds

3. Adjust the variables in the "necessary variables" section to match your data

4. Adjust the variables in the "optional variables" section. For data with sampling frequencies > 120Hz, the defaults should be a good starting point. For 120Hz data, we suggest to start with opt.downsamples set to [2 3 5], and opt.chebyOrder to 7. For data with lower frequencies, we suggest to start with  opt.downsampFilter set to 0, opt.downsamples to [2 3], and opt.steptime to 0.

5. Run the algorithm

**Note**: _Signal Processing Toolbox_ is required for the default downsampling procedure. If not available, set `opt.downsampFilter` to 0. This will use a different downsampling procedure.

Tested on MATLAB R2012a, R2014b, R2016a, R2017a, & R2019b
