function [ delta, delta_c_perimeter, delta_nb_c_perimeter  ] = AdhesionAndPerimeter( cells, mat , c, nb_c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
e_adhesion = 0; 
new_e_adhesion = 0;
e_perimeter = 0;


before_latt_c = zeros(3,3);
before_latt_nb_c = zeros(3,3);
before_mat = zeros(5,5);

after_latt_nb_c = zeros(3,3);
after_latt_c = zeros(3,3);
after_mat = zeros(5,5);


global J_CC J_CM LAMDA_PERIMETER ;

nhood = mat(2:4,2:4);
if c>0 && nb_c>0
   for i = 1:3
        for j = 1:3
    if nhood(i,j) == c
        before_latt_c(i,j) = 1;
    else
        before_latt_c(i,j) = 0;
    end
    if nhood(i,j) == nb_c
        before_latt_nb_c(i,j) = 1;
    else
        before_latt_nb_c(i,j) = 0;
    end
    
        end
   end
   
   
   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(before_latt_c, lut);
   before_c_perimeter = sum(out(:));
   
   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(before_latt_nb_c, lut);
   before_nb_c_perimeter = sum(out(:));
   
   
end  
 if c>0 && nb_c == 0
   for i = 1:3
        for j = 1:3
    if nhood(i,j) == c
        before_latt_c(i,j) = 1;
    else
        before_latt_c(i,j) = 0;
    end
    
        end
   end 
   

   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(before_latt_c, lut);
   before_c_perimeter = sum(out(:));
    %-
end
if nb_c>0 && c==0
   for i = 1:3
        for j = 1:3
    if nhood(i,j) == nb_c
        before_latt_nb_c(i,j) = 1;
    else
        before_latt_nb_c(i,j) = 0;
    end
        end
   end
   
   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(before_latt_nb_c, lut);
   before_nb_c_perimeter = sum(out(:));
   
    %+
end

dist = [     4, sqrt(5), 2, sqrt(5),       4;
       sqrt(5), sqrt(2), 1, sqrt(2), sqrt(5);
             2,       1, 1,       1,       2;
       sqrt(5), sqrt(2), 1, sqrt(2), sqrt(5);
             4, sqrt(5), 2, sqrt(5),       4;] ;  
      
if c>0     
for i =1:5
    for j = 1:5
if(mat(i,j)==c)      
before_mat(i,j) = 0;  
elseif(mat(i,j)==0) 
before_mat(i,j) = J_CM/dist(i,j) ; 
else
    before_mat(i,j) = J_CC/dist(i,j);
end
    end
end
end

if c==0
   for i = 1:5
       for j = 1:5       
if(mat(i,j)==0) 
before_mat(i,j) = 0 ; 
else
    before_mat(i,j) = J_CM/dist(i,j) ; 
       end
   end
   end
end

 e_adhesion = sum(sum(before_mat));

mat(3,3) = nb_c;
nhood = mat(2:4,2:4);

if nb_c > 0    
for i =1:5
    for j = 1:5
if(mat(i,j)==nb_c)      
after_mat(i,j) = 0;  
elseif(mat(i,j)==0) 
after_mat(i,j) = J_CM/dist(i,j) ; 
else
    after_mat(i,j) = J_CC/dist(i,j);
end
    end
end
end

if nb_c==0
   for i = 1:5
       for j = 1:5       
if(mat(i,j)==0) 
after_mat(i,j) = 0 ; 
else
    after_mat(i,j) = J_CM/dist(i,j) ; 
       end
   end
   end
end

    new_e_adhesion = sum(sum(before_mat));


delta1  = new_e_adhesion - e_adhesion;




if c>0 && nb_c>0
   
   for i = 1:3
        for j = 1:3
    if nhood(i,j) == c
        after_latt_c(i,j) = 1;
    else
        after_latt_c(i,j) = 0;
    end
    if nhood(i,j) == nb_c
        after_latt_nb_c(i,j) = 1;
    else
        after_latt_nb_c(i,j) = 0;
    end
    
        end
   end
%    after_latt_c
   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(after_latt_c, lut);
   after_c_perimeter = sum(out(:));
   
   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(after_latt_nb_c, lut);
   after_nb_c_perimeter = sum(out(:));
    %-
end

if c>0 && nb_c == 0
   
   for i = 1:3
        for j = 1:3
    if nhood(i,j) == c
        after_latt_c(i,j) = 1;
    else
        after_latt_c(i,j) = 0;
    end
        end
   end
%    after_latt_c
   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(after_latt_c, lut);
   after_c_perimeter = sum(out(:));
end



if nb_c>0 && c == 0
 
   for i = 1:3
        for j = 1:3
    if nhood(i,j) == nb_c
        after_latt_nb_c(i,j) = 1;
    else
        after_latt_nb_c(i,j) = 0;
    end
        end
   end
   
   lut = makelut(@checkForPixelEdges, 3);
   out = bwlookup(after_latt_nb_c, lut);
   after_nb_c_perimeter = sum(out(:));
   
    %+
end


if c>0
 delta_c_perimeter = after_c_perimeter -  before_c_perimeter;
 e_perimeter = e_perimeter + LAMDA_PERIMETER*( (cells.perimeter(c) + delta_c_perimeter - cells.target_perimeter(c))^2 - (cells.perimeter(c) - cells.target_perimeter(c))^2 );
else
   delta_c_perimeter = 0; 
end
if nb_c>0
 delta_nb_c_perimeter = after_nb_c_perimeter -  before_nb_c_perimeter;  
 e_perimeter = e_perimeter + LAMDA_PERIMETER*( (cells.perimeter(nb_c) + delta_nb_c_perimeter - cells.target_perimeter(nb_c))^2 - (cells.perimeter(nb_c) - cells.target_perimeter(nb_c))^2) ;
else
   delta_nb_c_perimeter = 0;
end 
 delta2 = e_perimeter;



delta = delta1 + delta2;
end