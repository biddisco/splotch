/*
 * Copyright (c) 2002-2005 Max-Planck-Society
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef PLANCK_CONSTANTS_H
#define PLANCK_CONSTANTS_H

#include <cmath>

// mathematical constants

const double pi=3.141592653589793238462643383279502884197;
const double twopi=2.0*pi;
const double inv_twopi=1.0/twopi;
const double fourpi=4.0*pi;
const double halfpi=pi/2.0;
const double inv_halfpi=1.0/halfpi;
const double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;

const double ln2 = std::log(2.);
const double ln10 = std::log(10.);

const double onethird=1.0/3.0;
const double twothird=2.0/3.0;
const double fourthird=4.0/3.0;

const double degr2rad=pi/180.0;
const double rad2degr=180.0/pi;

/// calculate the FWHM of a Gauss curve from its sigma.
const double sigma2fwhm=2.3548200450309493; // sqrt(8*log(2.))
const double fwhm2sigma=1/sigma2fwhm;

// physical constants

const double Jansky2SI=1.0e-26;
const double SI2Jansky=1.0e+26;

//! light speed in m/s
const double speedOfLight=2.99792458e8;

//! Boltzmann's constant in J/K
const double kBoltzmann=1.380658e-23;

//! Stefan-Boltzmann constant in W/m^2/K^4
const double sigmaStefanBoltzmann=5.67051e-8;

//! Planck's constant in J s
const double hPlanck=6.6260755e-34;

//! Astronomical unit in m
const double astronomicalUnit=1.49597893e11;

//! Solar constant in W/m^2
const double solarConstant=1368.0;

//! Tropical year in s
const double tropicalYear=3.15569259747e7;

//! Average CMB temperature in K
const double tcmb = 2.726;

//! Colatitude of the solar system motion relative to CMB
//! (ecliptical coordinates)
const double solsysdir_ecl_theta = 1.7678013480275747;

//! Longitude of the solar system motion relative to CMB
//! (ecliptical coordinates)
const double solsysdir_ecl_phi = 3.0039153062803194;

//! Speed of the solar system motion relative to CMB in m/s
const double solsysspeed = 371000.0;


// technical constants

//! Healpix value representing "undefined"
const double Healpix_undef=-1.6375e30;

#endif
