# MSc-Medical-Physics-Project-Imperial
Super-resolution reconstruction algorithm for MRI based on COMBINE with an example 2D dataset from a spherical phantom on a 3T Siemens Magnetom Verio (MATLAB .mat file)

## Structure
* 'SR_COMBINE': source code
* 'external': external source code, see Acknowledgements
* 'image_282_data': example 2D dataset from a spherical phantom on a 3T Siemens Magnetom Verio (MATLAB .mat file)

## Example data
The sequence parameters were as follows: TR/TE=10/5ms; N<sub>PE</sub> x N<sub>FE</sub> = 50x248; FOV = 250x250mm<sup>2</sup>; slice thickness = 10mm; α=0.4°; 1000 dummy TRs (10s of steady-state preparation); 5 separate images at equidistant phase increments (i.e. 0°,72°,144°,...), with the unbalanced gradient along the PE direction.

## Acknowledgements
The ifft2c and fft2c functions were written by Michael Lustig for ESPIRiT. All use and distribution rights are as described in the original code.
This code is based on Lally PJ, Matthews PM, Bangerter NK. Unbalanced SSFP for super-resolution in MRI. Magn Reson Med. 2020;00:1–13. https://doi.org/10.1002/mrm.28593
