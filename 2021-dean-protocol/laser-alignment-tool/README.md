# LaserAlignmentTool
## Overview
Originally developed by the Gustafsson Lab, this  alignment tool is designed for the coaxial alignment of optical elements with a laser diode.  

Additional information is available in Sarah Abrahamsson's manuscript - MultiFocus Polarization Microscope (MF- PolScope) for 3D polarization imaging of up to 25 focal planes simultaneously.  https://www.osapublishing.org/DirectPDFAccess/80AF3C88-EC2F-A4F9-3AD496CAE0FFFA2A_313825/oe-23-6-7734.pdf - Appendix E nicely describes how to alignn the laser alignment tool, which involves retroreflecting the laser off of a series of co-planar optical surfaces.

This alignment device consists of a hollow brass cylinder (we have also manufactured one using aluminum without any problems) with multiple set screws for holding and aligning a cylindrical laser diode, which fits within the hollow interior.  

## Potential Useful Modifications
The original design called for 0-80 set screws, and an RMS external thread.  However, We find that the 0-80 set screws are bit small, and larger set screws are convenient.  Furthermore, we recommend that you change the external threads to accomodate your optomechanics.  Student and academic versions of Autodesk Inventor are available for free (https://www.autodesk.com/education/free-software/inventor-professional), and there are numerous tutorials online on how to use the software (e.g., https://www.youtube.com/watch?v=SsLkAokkeR8).

### Common thread types:
* 1.035"-40 - SM1 thread, ThorLabs
* M25 or M32 Nikon Instruments
* RMS - Olympus
* M26-32 Mitotoyo
* M27 - Zeiss

## Laser Diode Modules
We have used laser diodes from Newport (https://www.newport.com/f/laser-diode-modules-cw), and we recommend that you choose a wavelength that is appropriate for your application.  For example, if you are planning on mainly imaging GFP, then select a wavelength of ~500 nm, which is intermediate to both the excitation (488 nm) and emission (509 nm) maxima.  We also caution you to use a reasonable laser power.  You will also need a power supply (LPMS-8-110, LPMS-5 220, or LPMS-5-110).  

## Machine Shops
For individuals without a machine shop, we recommend that you contact a third-party rapid prototyping CNC company.  Potential examples include:
* Protolabs - https://www.protolabs.com
* Prismier - https://prismier.com/service/prototype-metal-cnc-machining-turning
* Xometry - https://www.xometry.com/rapid-prototyping-service
