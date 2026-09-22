function [varargout] = atomic_constants(name)
% ATOMIC_CONSTANTS Returns the values of various physical constants in Hartree units
%
% Supported constants include:
%    "c" - speed of light
%    "ε0" - vacuum permittivity
%    "μ0" - vacuum permeability
%    "kB" - Boltzman's constant
%    "G" - gravitational constant
%    "ħ" - Reduced Plank constant
%    "e" - electron charge
%    "μB" - Bohr magneton
%
% SEE ALSO: ATOMIC_UNITS
arguments (Repeating)
	name {mustBeMember(name,["c","ε0","μ0","kB","G","ħ","e","μB"])}
end

% base constants for Hartree units
ee = 1.602176634e-19; % electron charge in coulombs
me = 9.109383702e-31; % electron mass in kg
a0 = 5.291772109e-11; % bohr radius in meters
hb = 1.054571818e-34; % reduced plank constant in joule-seconds
Eh = 4.359744722e-18; % hartrees(energy) in joules

% physical constants in SI units
G  = 6.6743015000e-11; % gravitational constant: (m^3 / kg * s^2)
c  = 2.9979245800e+08; % speed of light: (m/s)
e0 = 8.8541878188e-12; % vacuum permittivity: (C^2 * s^2 / kg * m^3)
u0 = 1.2566370613e-06; % vacuum permeability: (kg * m / C^2)
kB = 1.3806490000e-23; % boltzman constant: (J/K)

% base SI units given in hartrees:
% meters, seconds, kilograms, joules, coulombs, and kelvin
[m,s,kg,J,C,K] = deal(1/a0,Eh/hb,1/me,1/Eh,1/ee,1);

for ii = 1:nargin
	switch name{ii}
		case "c"
			val = c  * (m/s);
		case "ε0"
			val = e0 * (C^2 * s^2)/(kg * m^3); %#ok
			val = 1/(4*pi); % use the more precise equivalent
		case "μ0"
			val = u0 * (kg * m / C^2);
		case "kB"
			val = kB * (J/K);
		case "G"
			val = G * (m^3 / (kg * s^2));
		case "ħ"
			val = hb * (J*s);
		case "e"
			val = ee * (C);
		case "μB"
			val = hb*ee/(2*me);
	end
	varargout{ii} = val; %#ok
end
% DONE
end