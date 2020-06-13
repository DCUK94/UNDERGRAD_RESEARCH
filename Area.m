  function [ delta ] = Area(cells, c, nb_c )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global LAMDA_AREA  ;

e_area=0;
if c>0
    e_area = e_area + LAMDA_AREA*(1 - 2*cells.area(c) + 2*cells.target_area(c));
    %-
end
if nb_c>0
    e_area = e_area + LAMDA_AREA*(1 + 2*cells.area(nb_c) - 2*cells.target_area(nb_c));
    %+
end
delta = e_area;

end


