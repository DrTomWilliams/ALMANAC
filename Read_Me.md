Author Dr Thomas Williams, Aberystwyth University

This repository contains the data associated with the method paper publication submitted to Space Weather (DOI to be added upon publication).

In the code tree are the ALMANAC IDL source files. To run ALMANAC please see cme.pro. The custom IDL files developed as part of ALMANAC are also included here in collaboration with Prof Huw Morgan, Aberystwyth University.

In addition to these files one will require the IDL coyote library and a SolarSoft installation on IDL. The repositories I have installed are listed below but note that not all of these may be required:

aia hmi eve stereo chianti eis sot xrt hessi hic iris swap lyra azam cactus cmes corimp mjastereo nlfff nrl s3drs sbrowser spex sunspice xray wispr smei spice stix cds eit lasco mdi sumer uvcs impact plastic secchi ssc trace helio ontology vso bcs hxt sxt ucon wbs

The list of complete compiled functions/procedures to get ALMANAC running is listed below:

SIZEARR
INT2STR_HUW
INT2STR
GET_SHIFT
ONE2N
ARRAY_INDICES
REGION_SIZE_FILTER
PAD_2D
TOTALSTRING
GET_REV_IND
UNPAD_2D
PAD_3D
UNPAD_3D
FITSHEAD2WCS
NTRIM
BYTE2INT
GET_FITS_PAR
GET_FITS_TIME
GET_FITS_ROLL
GET_FITS_CDELT
GET_FITS_CEN
COMP_FITS_CRVAL
WCS_FIND_KEYWORD
WCS_SIMPLE
VALID_WCS
WCS_FIND_TIME
WCS_FIND_POSITION
PB0R
ANYTIM2JD
JULDAY
SUN_POS
WCS_RSUN
WCS_AU
TIM2CARR
GET_SUN
TIM2JD
RECPOL
WCS_FIND_SPECTRUM
WCS_FIND_DISTORTION
WCS_GET_COORD
WCS_PROJ_TAB
WCS_PROJ_TAN
WCS_CONVERT_FROM_COORD
WCS_CONV_FIND_ANG_UNITS
WCS_CONV_HPC_HCC
WCS_CONV_FIND_DSUN
WCS_PARSE_UNITS
WCS_PARSE_UNITS_BASE
WCS_CONV_HCC_HG
WCS_CONV_FIND_HG_ANGLES
MINMAX
MOVIEMAKER
MGN
PLOT_IMAGE
SETSCALE
GET_VIEWPORT
EXPAND_TV
GET_IM_KEYWORD
TVSELECT
IM_KEYWORD_SET
BSCALE
WHERE_MISSING
BYTSCLI
STORE_TV_SCALE
TVUNSELECT
CGOPLOT
CGPLOT
CGSETCOLORSTATE
CGGETCOLORSTATE
SETDEFAULTVALUE
CGCHECKFORSYMBOLS
CGDEFAULTCOLOR
COLORSAREIDENTICAL
CGDEFCHARSIZE
CGCOLOR
CGCOLOR24
JPEG
TOM_FFMPEG
