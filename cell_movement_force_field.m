% % % % % % % % % % % % % % % % % CELL MOTILITY PROCESS
clear all
clc;
load('density_06_64.mat');
global J_CC J_CM  LAMDA_AREA LAMDA_PERIMETER XMAX YMAX F
COM = zeros(XMAX,YMAX);
tau = 10;
MCS = 2000;
del_psi = pi/3;
F = 10;
shi = 0.95;
STEP_MAX=XMAX*YMAX*16*MCS;

BaseName='density06_';
ii = 1;

for i = 1:cell_numb
   for j = 1:XMAX
       for k = 1:YMAX
           if latt(j,k) == i
               COM(j,k) = 1;
           else
               COM(j,k) = 0;
           end   
       end
   end
    [x_com,y_com] = CenterOfMass(COM);
    
    cells.comx(i,1) = x_com;
    cells.comy(i,1) = y_com; 
    cells.x_com(i,1) =  cells.comx(i,1);
    cells.y_com(i,1) =  cells.comy(i,1);
end

for i = 1:cell_numb
       cells.psi(i,1) = 2*pi*rand;
end


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
        [ delta_force ,c_x_com, c_y_com, nb_c_x_com, nb_c_y_com] = force( cells, latt, c, nb_c, rnd_y,rnd_x);
        e_all =  Area(cells, c, nb_c ) + AdeAndPeri_energy + delta_force ;
         
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
                cells.comx(c,1) = c_x_com;
                cells.comy(c,1) = c_y_com;
            end
            if nb_c > 0
                cells.area(nb_c) = cells.area(nb_c)+1;
                cells.perimeter(nb_c) = cells.perimeter(nb_c) + delta_nb_c;
                cells.comx(nb_c,1) = nb_c_x_com;
                cells.comy(nb_c,1) = nb_c_y_com;
            end
        end   
    end
    
  
%     remain_steps = STEP_MAX - step
    
    if mod(step,XMAX*YMAX*16)==0
        
        idx = step/(XMAX*YMAX*16);
        disp(idx);

for i = 1:cell_numb
cells.x_com(i,idx + 1) =  cells.comx(i,1);
cells.y_com(i,idx + 1) =  cells.comy(i,1);
end
    end
    
    
    if mod(step,tau*XMAX*YMAX*16)==0
        epsilon = normrnd(0,1);
        for i = 1:cell_numb
        cells.psi(i,1) =  cells.psi(i,1)*shi + del_psi*epsilon;
        end
    end
    
     if mod(step,400*XMAX*YMAX*16)==0
        FileName=[BaseName,num2str(ii)];
        save(FileName);
        ii = ii + 1;
     end
    
end

save('final_presentT_06_64_anim.mat');

