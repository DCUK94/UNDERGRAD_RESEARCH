clear all
clc;

global J_CC J_CM  LAMDA_AREA LAMDA_PERIMETER XMAX YMAX

XMAX = 256;
YMAX = 256;

MCS = 500;

STEP_MAX=XMAX*YMAX*16*MCS;

init_diameter_y = 8;
init_diameter_x = 8;

density = 0.9;

LAMDA_AREA=1; 
LAMDA_PERIMETER=0.2;

TARGET_CELL_SIZE = pi*64;
TARGET_CELL_PERIMETER = pi*16;

J_CC = 0.1; 
J_CM = 0.01;

x = floor(sqrt((XMAX*YMAX*density)/(init_diameter_y*init_diameter_y)));

TEMPERATURE = 2.0;

X_CELL_NUMB=x; 
Y_CELL_NUMB=x; 

mat  =  1:X_CELL_NUMB*Y_CELL_NUMB;
mat  =  reshape(mat,Y_CELL_NUMB,X_CELL_NUMB);
mat  =  repelem(mat,init_diameter_y,init_diameter_x);
[len_row, len_col] = size(mat);

latt = zeros(YMAX,XMAX);
latt(1:len_row ,1:len_col) = mat;
latt = circshift(latt,[round((YMAX-len_row)/2) round((XMAX-len_col)/2)]);


cell_numb=X_CELL_NUMB*Y_CELL_NUMB;



initial_area = sum(sum(latt==1));

initial_perimeter = 2*init_diameter_y + 2*init_diameter_x;

nhd = zeros(5,5);

for i=1:cell_numb
    
    cells.area(i,1)= initial_area;
    cells.perimeter(i,1) = initial_perimeter;
end

cells.target_area(1:cell_numb,1)=TARGET_CELL_SIZE;
cells.target_perimeter(1:cell_numb,1)=TARGET_CELL_PERIMETER;


for step = 1:STEP_MAX
    
    
    rnd_x = randi(XMAX);
    rnd_y = randi(YMAX);
    c = latt(rnd_y,rnd_x);
    
    nb_indx=randi(8);
    switch nb_indx
        case 1
            nb_x=-1; nb_y=-1;
        case 2
            nb_x=-1; nb_y=0;
        case 3
            nb_x=-1; nb_y=1;
        case 4
            nb_x=0; nb_y=-1;
        case 5
            nb_x=0; nb_y=1;
        case 6
            nb_x=1; nb_y=-1;
        case 7
            nb_x=1; nb_y=0;
        case 8
            nb_x=1; nb_y=1;
    end
    nb_c = [rnd_y+nb_y,rnd_x+nb_x];
    
     %Periodic Boundary Conditions
   if nb_c(2) < 1 %Left
       nb_c(2) = nb_c(2) + XMAX;
   end
   if nb_c(2) > XMAX %Right
       nb_c(2) = nb_c(2) - XMAX;
   end
   if nb_c(1) < 1 %Top
       nb_c(1) = nb_c(1) + YMAX;
   end    
   if nb_c(1) > YMAX %Bottom
       nb_c(1) = nb_c(1) - YMAX;
   end
      
      nb_c = latt(nb_c(1),nb_c(2));
     
      
      if c ~= nb_c
        % Calculation of energy difference and transition probability
        
        k =1;
        l = 1;
  for i = rnd_y - 2:rnd_y + 2
  for j =rnd_x - 2:rnd_x + 2
   if j < 1 %Left
       j = j + XMAX;
   end
   if j > XMAX %Right
       j = j - XMAX;
   end
   if i < 1 %Top
       i = i + YMAX;
   end    
   if i > YMAX %Bottom
       i = i - YMAX;
   end
   nhd(k,l) = latt(i,j);
   l = l + 1;
    end
    l = 1;
    k = k + 1;
  end
  
  
  
        [AdeAndPeri_energy , delta_c, delta_nb_c]   = AdhesionAndPerimeter(cells, nhd , c, nb_c);
        e_all =  Area(cells, c, nb_c ) + AdeAndPeri_energy ;
         
        if e_all >= 0
            prob = exp(-e_all/TEMPERATURE);
        elseif e_all < 0
            prob = 1.0;
        end
        
        % Update of the state
        if prob >= rand
            latt(rnd_y,rnd_x) = nb_c;
            if c > 0
                cells.area(c) = cells.area(c)-1;
                cells.perimeter(c) = cells.perimeter(c) + delta_c;
            end
            if nb_c > 0
                cells.area(nb_c) = cells.area(nb_c)+1;
                cells.perimeter(nb_c) = cells.perimeter(nb_c) + delta_nb_c;
            end
        end   
    end
    
  
    remain_steps = STEP_MAX - step
   
end

save('density_0.9.mat');

