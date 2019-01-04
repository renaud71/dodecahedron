# dodecahedron
central projection of a spherical map on a dodecahedron, matlab script
author: R Toussaint, Institut de Physique du Globe de Strasbourg, CNRS / University of Strasbourg, Dec 2018
[![DOI](https://zenodo.org/badge/163993570.svg)](https://zenodo.org/badge/latestdoi/163993570)

Matlab script
the script is called dodemap.m
it reads a .png file containing the image to map: color (RGB) as function of longitue and latitude, mapped in a 2x1 standard image representation
(each column is a meridian, each line is a parallel, latitudes vary linearly as lines from top to bottom, longitudes linearly as column numbers from left to rightexport_fig_out_l5.png

%colorcarte=imread('imagesourcered2.jpg'); % opening the image where color is stored as function ot longitude, latitude (image with horizontal pixels twice larger than vertical pixel size, 360 degrees versus 180)
in this case, 'export_fig_out_l5.png'
a tomography map of the Earth mantle, at 530 km under the earth surface - the colors correspond to the percentage of variation of the S waves.
Structures corresponding to subducting plate slabs, hot spots, accretion of oceanic cruts, disappeared oceans under the crust, as Farallon and Thetys, 
are visible. data: C. Zaroli, Geoph. Journal International, 2016, DOI: 10.1093/gji/ggw315
Any spherical map can be used as input, only the image dimension have to be adjusted in the script.

The output dodecahedron is exported as export_fig_out_result.png

It can be built by printing it on paper and using scisors and tape.

The output part, for saving purposes, uses the library  https://fr.mathworks.com/matlabcentral/fileexchange/23629-export_fig
(or https://github.com/altmany/export_fig)
that should be present in the directory of the script, copied in a subdirectory altmany-export_fig-5b3965b/

This part at the end of the script can also be commented, and the figures saved directly, commenting the last 7 lines of the script.

