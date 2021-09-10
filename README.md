# G305 Paper 1
All the analysis done on G305 - Paper 1

Link to the paper : https://www.aanda.org/articles/aa/full_html/2021/06/aa40205-20/aa40205-20.html

The ".class" files are written in Gildas-CLASS package which is a state-of-the-art software package used
for analyzing (sub-)millimeter radioastronomical data. It is mainly written in FORTRAN-90/95 with a few
parts in C/C++ (https://www.iram.fr/IRAMFR/GILDAS/).

A short explanation of each file:
- redall_lasmagal.class : Reduction pipeline for the raw data
- noise-map.class : calculate channel-wise noise map from a spectral cube
- moment-0-map.class : calculate integrated intensity map from a spectral cube
- moment-1-2-map.class : calculate 1st and 2nd moment maps from a spectral cube
- spectrum_sum.class : calculate the average spectrum from a given spectral cube
- new_channel_maps.py : extract and plot channel-wise maps from a spectral cube
- do_lte.class : calculate the channel-wise excitation temperature and optical depth of a spectral cube
- column_density_map.py : calculate a column density data cube given the excitation temp. and optical depth cubes.
- *plot*.py : plotting with python. Names are mostly self-explanatory
- Ionization_Length_vs_Time.py : modelling the radius of ionizaiton for a spherically expanding shell from an OB star.
- excitation_vs_8um.py :  All analyses pertaining to section 6 of the paper.
- velocity_centroid_pdf_stat.py : calculate centroid velocity pdf and its moments (sec-7.1 of the paper)
- feedback_with_random_sqrt_pixels.py : analysis for section 7.2
- feedback_intervals_with_random_sqrt_pixels.py : analysis for section 7.3
