clear; close all; clc

altitude = (500:1:3000)*1e3;
latitude = 0;
longitude = 0;
year = 2000;
dayOfYear = 4;
UTseconds = 1;


[T rho] = atmosnrlmsise00(altitude, latitude, longitude, year, dayOfYear, UTseconds);
semilogy(altitude,rho(:,6))