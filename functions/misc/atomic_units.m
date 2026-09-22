function [varargout] = atomic_units(name)
% ATOMIC_UNITS Returns the SI-unit equivalent in Hartree atomic units.
%
% [VAL1,VAL2,...] = ATOMIC_UNITS(NAME1,NAME2,...) returns conversion
% factors for multiple units. Each input NAME must be a string specifying
% a supported SI unit, optionally preceded by an SI prefix.
%
% Supported SI units:
%    m - meter
%    s - second
%    J - joule
%    eV - electron volt
%    g - gram
%    N - newton
%    V - volt
%    T - tesla
%
% Supported SI prefixes:
%    P - peta (10^15)
%    T - tera (10^12)
%    G - giga (10^9)
%    M - mega (10^6)
%    k - kilo (10^3)
%    c - centi (10^-2)
%    m - milli (10^-3)
%    u - micro (10^-6)
%    n - nano (10^-9)
%    p - pico (10^-12)
%    f - femto (10^-15)
%
% Examples
% =========
% Write the Bohr radius in angstroms:
%    disp("1 Bohr = "+1/(1e-10*atomic_units("m"))+" A") 
%
% Write the Hartree unit of energy in eV:
%    disp("1 Hartree = "+1/atomic_units("eV")+" eV") 
%
% Write a femtosecond in terms of Hartree-seconds:
%    disp("1 femtosecond = "+atomic_units("fs")+" Hartree-seconds")
%
% Get conversion factors for several units at once:
%    [nm,fs,meV] = atomic_units("nm","fs","meV");
%
%
% SEE ALSO: ATOMIC_CONSTANTS
arguments (Repeating)
	name string
end
% base constants for Hartree units
ee = 1.602176634e-19; % electron charge in coulombs
me = 9.109383702e-31; % electron mass in kg
a0 = 5.291772109e-11; % bohr radius in meters
hb = 1.054571818e-34; % reduced plank constant in joule-seconds
Eh = 4.359744722e-18; % hartrees(energy) in joules

% conversion factors into hartrees
meter  = 1/a0;         % one meter
second = Eh/hb;        % one second
Joule  = 1/Eh;         % one Joule
gram   = (0.001)/me;   % one gram
Newton = a0/Eh;        % one Newton
Volt   = ee/Eh;        % one Volt
Tesla  = a0^2 * ee/hb; % one Tesla

% unit labels and values
unitName = ["m", "s", "J", "eV", "g", "N", "V", "T"];
unitVals = [meter, second, Joule, Volt, gram, Newton, Volt, Tesla];

% prefix labels and values
prefixName = ["P", "T", "G", "M", "k", "", "c", "m", "u", "n", "p", "f"];
prefixVals = [1e15, 1e12, 1e9, 1e6, 1e3, 1, 1e-2, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15];

% array containing all possible combinations of prefixes and units
validUnits = (prefixName(:) + unitName);

% return unit conversion factors
for nn = 1:nargin
	[row,col] = find(name{nn} == validUnits);
	
	if isequal(numel(row),numel(col),1)
		varargout{nn} = prefixVals(row) * unitVals(col); %#ok
	else
		% input is not a valid unit
		error("ERROR: unit requested at position "+nn+" is invalid")	
	end
end
% DONE
end