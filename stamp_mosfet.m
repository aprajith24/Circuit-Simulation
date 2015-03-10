function [new_M] = stamp_mosfet(old_M,D);
%STAMP_MOSFET : Stamps mosfets into the MNA matrix
% syntax : [new_M]=stamp_mosfet(old_M,D)
%
% new_M,old_M are self-explanatory
% D is the data vector corresponding to the conductance or resistance.
