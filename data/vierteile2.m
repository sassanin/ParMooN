clear
% Anzahl der Punkte pro Kante
N0 = 2;
N1 = 2;
N2 = 3;
N3 = 5;
N4 = 4;
rotations = 16;

% Kanta 0, z = 140, -50 <=x <= 50
x0 = 0;
y0 = 140;
h = 50/(N0-1);
for i=1:N0
    x(i) = -(i-1)*h;
%    if i== N0/2
%        x(i) = -24;
%    end
%    if i== N0/2 +1 
%        x(i) = -29;
%    end
    y(i) = 140;
end
points = N0;


% Kante 1, x=-50, 19.39 <= y <= 140
x0 = -50;
y0 = 140;
h = (140-19.377)/(N1-1);
for i=2:N1
    x(points+(i-1)) = -50;
    y(points+(i-1)) = 140 - (i-1)*h;
end
points = points + N1 -1;

% berechne Endpunkt der Kante 2
a = 40;
b = -(100-sqrt(90*90-40*40));
c = -100;
r1 = 10;
r2 = 100;
t = r1*r1-r2*r2+c*c-a*a-b*b

p(1) = 1 + (c-b)^2/(a*a)
p(2) = 2*c+t*(c-b)/(a*a)
p(3) = c*c-r2*r2+t*t/(4*a*a)

q = roots(p);

N2y = (q(1)+q(2))/2;
N2x = r2*r2-(N2y+c)^2;
N2x = -sqrt(abs(N2x));

% Winkel f"ur den Kreisbogen
beta = atan((100+b)/40);
% Mittelpunkt des Kreisbogens
N2x = -40;
N2y = 100-sqrt(90*90-40*40)

% Kante 2
h = beta/(N2-1);
for i=2:N2
    phi = pi + (i-1)*h;
    x(points+(i-1)) = r1 * cos(phi) + N2x;
    y(points+(i-1)) = r1 * sin(phi) + N2y;
end

points = points + N2 - 1;

% Kante 3
alpha1=asin(2.5*sqrt(2)/100.);
alpha = pi/2-beta-alpha1;
N3x = 0;
N3y = 100;


h = alpha/(N3-2);

  % x(points+(N3-4)) =-36.664;
   % y(points+(N3-4)) = 6.9636;
   % x(points+(N3-3)) =-26.257;
   % y(points+(N3-3)) = 3.5086;
   %  x(points+(N3-2)) =-3.005;
   % y(points+(N3-2)) = 0.051;
   % x(points+(N3-1)) = 0.;
   % y(points+(N3-1)) = 0.;
 
for i=2:N3
    phi = pi+beta + (i-1)*h;
    x(points+(i-1)) = r2 * cos(phi) + N3x;
    y(points+(i-1)) = r2 * sin(phi) + N3y;

   %if i == N3-1
    %  x(points+(i-1)) =-3.005 ;
   % y(points+(i-1)) = 0.051;
  % else
     if i == N3
    x(points+(i-1)) =0.0 ;
    y(points+(i-1)) = 0.0;
   
    end
end
points = points + N3 - 1;

% Ebenendarstellung f"ur den Deckel, ben"otigt in PRM-Datei
% Punkt, Richtungsvektor, Normalenvektor (aussen)
% Ebenengleichung z=140
% pcl_NO => plane counter
pcl_NO = 1;
point(pcl_NO,1:3) = [0 0 140];
richtung(pcl_NO,1:3) = [1 0 0];
normale(pcl_NO,1:4) = [0 0 1 pcl_NO];

polygon(pcl_NO,1:16) =  [1 9 16 23 30 37 44 51 58 65 72 79 86 93 100 107];  
pcl_NO = pcl_NO + 1;
% compute the points in space

gamma = 2*pi/rotations;

% z"ahler
index = 1;
offset = 1;
offset1 = 1;
offset2 = 1;
offset3 = 1;
offset4 = 1;
offset5 = 1;
offset6 = 1;
% rotiere die Kurve
for i=1:rotations
    % Winkel    
    phi = (i-1) * gamma;
    % Schleife "uber die Punkte
    for j=1:points
        % der Punkt (0,0,140) nur einmal
        if i>1 && j==1
            continue;
        end
        % der Punkt (0,0,0) nur einmal        
        if i>1 && j==points
            continue;
        end
        % berechne die Punkte im Raum
        x3D(index,1) = index-1;
        x3D(index,2) = -x(j) * cos(phi);
        x3D(index,3) = -x(j) * sin(phi);
        x3D(index,4) = y(j);
        index = index + 1;    
        % berechne die Vektoren und den Index f"ur die PRM-Datei
        % obere Ebene
        if j<= N0
            % Fl"achenindex f"ur diesen Punkt
            x3D(index-1,5) = 1;
            continue;
        end
	%%%%
	%%%% CYLINDER
	%%%%
        % Seitenfl"achen, die Ebene ist die gleiche f"ur alle
        % Fl"achen dieser Rotation
        if j< N0+N1 && i>1
            % im ersten Punkt wird die Ebenengleichung berechnet
            % nehme daf"ur den berechneten Punkt, den der auch an der oberen
            % Ebene liegt und einen Punkt der
            % vorherigen Rotation, das geht erst f"ur i>1
            if j==N0+1
                punkt(1,1:3) = x3D(index-2,2:4);
                punkt(2,1:3) = x3D(index-1,2:4);
                punkt(3,1:3) = x3D(index-2-points+offset,2:4);
                punkt(4,1:3) = x3D(index-1-points+offset,2:4);
                poly(1,1) = x3D(index-2,1);
                poly(1,2) = x3D(index-1,1);
                poly(1,3) = x3D(index-2-points+offset,1);
                poly(1,4) = x3D(index-1-points+offset,1);
                % von nun an werden immer zwei Punkte nicht ber"ucksichtg
                % in der Liste (0,0,0), (0,0,140)
                offset = 2;
                point(pcl_NO,1:3) = punkt(1,1:3);
                polygon(pcl_NO,1:4) = poly(1,1:4);
                % nehme Einheitsvektor
                richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
                richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
                % berechne die Normale
                % zun"achst zweiten Richtungsvektor
                dir2 = punkt(3,1:3) - punkt(1,1:3);
                dir2 = dir2/norm(dir2);
                % nun Formel f"ur Vektorproduct
		c = cross(dir2, richtung(pcl_NO,1:3));
		normale(pcl_NO,1:3) = c(1:3);
		normale(pcl_NO, 4) = pcl_NO;
                pcl_NO = pcl_NO + 1;
            end
            x3D(index-1,5) = pcl_NO-1;
        end
%%%%
	%%%% TORUS 1
	%%%%
        % Seitenfl"achen, die Ebene ist die gleiche f"ur alle
        % Fl"achen dieser Rotation
        if j< N0+N1+N2+N3 && i>1
            % im ersten Punkt wird die Ebenengleichung berechnet
            % nehme daf"ur den berechneten Punkt, den der auch an der oberen
            % Ebene liegt und einen Punkt der
            % vorherigen Rotation, das geht erst f"ur i>1
            if j==N0+N1
              offset = 1;
                punkt(1,1:3) = x3D(index-2,2:4);
                punkt(2,1:3) = x3D(index-1,2:4);
                punkt(3,1:3) = x3D(index-2-points+offset1,2:4);
                punkt(4,1:3) = x3D(index-1-points+offset1,2:4);
                poly(2,1) = x3D(index-2,1);
                poly(2,2) = x3D(index-1,1);
                poly(2,3) = x3D(index-2-points+offset1,1);
                poly(2,4) = x3D(index-1-points+offset1,1);
                % von nun an werden immer zwei Punkte nicht ber"ucksichtg
                % in der Liste (0,0,0), (0,0,140)
               offset1 = 2;
                point(pcl_NO,1:3) = punkt(1,1:3);
                polygon(pcl_NO,1:4) = poly(2,1:4);
                % nehme Einheitsvektor
                richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
                richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
                % berechne die Normale
                % zun"achst zweiten Richtungsvektor
                dir2 = punkt(3,1:3) - punkt(1,1:3);
                dir2 = dir2/norm(dir2);
                % nun Formel f"ur Vektorproduct
		c = cross(dir2, richtung(pcl_NO,1:3));
		normale(pcl_NO,1:3) = c(1:3);
		normale(pcl_NO, 4) = pcl_NO;
                pcl_NO = pcl_NO + 1;
            end
            x3D(index-1,5) = pcl_NO-1;
	end	
	%%%%
	%%%% TORUS 2
	%%%%
        % Seitenfl"achen, die Ebene ist die gleiche f"ur alle
        % Fl"achen dieser Rotation
        if j< N0+N1+N2+N3 && i>1
            % im ersten Punkt wird die Ebenengleichung berechnet
            % nehme daf"ur den berechneten Punkt, den der auch an der oberen
            % Ebene liegt und einen Punkt der
            % vorherigen Rotation, das geht erst f"ur i>1
            if j==N0+N1+1
                punkt(1,1:3) = x3D(index-2,2:4);
                punkt(2,1:3) = x3D(index-1,2:4);
                punkt(3,1:3) = x3D(index-2-points+offset2,2:4);
                 punkt(4,1:3) = x3D(index-1-points+offset2,2:4);
                poly(3,1) = x3D(index-2,1);
                poly(3,2) = x3D(index-1,1);
                poly(3,3) = x3D(index-2-points+offset2,1);
                poly(3,4) = x3D(index-1-points+offset2,1);
                % von nun an werden immer zwei Punkte nicht ber"ucksichtg
                % in der Liste (0,0,0), (0,0,140)
                offset2 = 2;
                point(pcl_NO,1:3) = punkt(1,1:3);
                polygon(pcl_NO,1:4) = poly(3,1:4);
                % nehme Einheitsvektor
                richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
                richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
                % berechne die Normale
                % zun"achst zweiten Richtungsvektor
                dir2 = punkt(3,1:3) - punkt(1,1:3);
                dir2 = dir2/norm(dir2);
                % nun Formel f"ur Vektorproduct
		c = cross(dir2, richtung(pcl_NO,1:3));
		normale(pcl_NO,1:3) = c(1:3);
		normale(pcl_NO, 4) = pcl_NO;
                pcl_NO = pcl_NO + 1;
            end
            x3D(index-1,5) = pcl_NO-1;
	end

		%%%%
	%%%% CALOTE 1
	%%%%
        % Seitenfl"achen, die Ebene ist die gleiche f"ur alle
        % Fl"achen dieser Rotation
        if j< N0+N1+N2+N3 && i>1
            % im ersten Punkt wird die Ebenengleichung berechnet
            % nehme daf"ur den berechneten Punkt, den der auch an der oberen
            % Ebene liegt und einen Punkt der
            % vorherigen Rotation, das geht erst f"ur i>1
            if j==N0+N1+2
                punkt(1,1:3) = x3D(index-2,2:4);
                punkt(2,1:3) = x3D(index-1,2:4);
                punkt(3,1:3) = x3D(index-2-points+offset3,2:4);
                punkt(4,1:3) = x3D(index-1-points+offset3,2:4);
                poly(4,1) = x3D(index-2,1);
                poly(4,2) = x3D(index-1,1);
                poly(4,3) = x3D(index-2-points+offset3,1);
                poly(4,4) = x3D(index-1-points+offset3,1);
               offset3 = 2;
                %n nun an werden immer zwei Punkte nicht ber"ucksichtg
                % in der Liste (0,0,0), (0,0,140)
                offset = 2;
                point(pcl_NO,1:3) = punkt(1,1:3);
                 polygon(pcl_NO,1:4) = poly(4,1:4);
                % nehme Einheitsvektor
                richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
                richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
                % berechne die Normale
                % zun"achst zweiten Richtungsvektor
                dir2 = punkt(3,1:3) - punkt(1,1:3);
                dir2 = dir2/norm(dir2);
                % nun Formel f"ur Vektorproduct
		c = cross(dir2, richtung(pcl_NO,1:3));
		normale(pcl_NO,1:3) = c(1:3);
		normale(pcl_NO, 4) = pcl_NO;
                pcl_NO = pcl_NO + 1;
            end
            x3D(index-1,5) = pcl_NO-1;
	end

	%%%%
	%%%% CALOTE 2
	%%%%
        % Seitenfl"achen, die Ebene ist die gleiche f"ur alle
        % Fl"achen dieser Rotation
        if j< N0+N1+N2+N3 && i>1
            % im ersten Punkt wird die Ebenengleichung berechnet
            % nehme daf"ur den berechneten Punkt, den der auch an
            % der oberen
            % Ebene liegt und einen Punkt der
            % vorherigen Rotation, das geht erst f"ur i>1
            if j==N0+N1+3
                punkt(1,1:3) = x3D(index-2,2:4);
                punkt(2,1:3) = x3D(index-1,2:4);
                punkt(3,1:3) = x3D(index-2-points+offset4,2:4);
                punkt(4,1:3) = x3D(index-1-points+offset4,2:4);
                poly(5,1) = x3D(index-2,1);
                poly(5,2) = x3D(index-1,1);
                poly(5,3) = x3D(index-2-points+offset4,1);
                poly(5,4) = x3D(index-1-points+offset4,1);
                offset4 = 2;
                % von nun an werden immer zwei Punkte nicht ber"ucksichtg
                % in der Liste (0,0,0), (0,0,140)
               
                point(pcl_NO,1:3) = punkt(1,1:3);
                polygon(pcl_NO,1:4) = poly(5,1:4);
                % nehme Einheitsvektor
                richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
                richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
                % berechne die Normale
                % zun"achst zweiten Richtungsvektor
                dir2 = punkt(3,1:3) - punkt(1,1:3);
                dir2 = dir2/norm(dir2);
                % nun Formel f"ur Vektorproduct
		c = cross(dir2, richtung(pcl_NO,1:3));
		normale(pcl_NO,1:3) = c(1:3);
		normale(pcl_NO, 4) = pcl_NO;
                pcl_NO = pcl_NO + 1;
            end
            x3D(index-1,5) = pcl_NO-1;
	end
        
  	%%%%
	%%%% CALOTE 3
	%%%%
        % Seitenfl"achen, die Ebene ist die gleiche f"ur alle
        % Fl"achen dieser Rotation
        if j< N0+N1+N2+N3 && i>1
            % im ersten Punkt wird die Ebenengleichung berechnet
            % nehme daf"ur den berechneten Punkt, den der auch an
            % der oberen
            % Ebene liegt und einen Punkt der
            % vorherigen Rotation, das geht erst f"ur i>1
            if j==N0+N1+4
                punkt(1,1:3) = x3D(index-2,2:4);
                punkt(2,1:3) = x3D(index-1,2:4);
                punkt(3,1:3) = x3D(index-2-points+offset6,2:4);
                punkt(4,1:3) = x3D(index-1-points+offset6,2:4);
                poly(13,1) = x3D(index-2,1);
                poly(13,2) = x3D(index-1,1);
                poly(13,3) = x3D(index-2-points+offset6,1);
                poly(13,4) = x3D(index-1-points+offset6,1);
                offset6 = 2;
                % von nun an werden immer zwei Punkte nicht ber"ucksichtg
                % in der Liste (0,0,0), (0,0,140)
               
                point(pcl_NO,1:3) = punkt(1,1:3);
                polygon(pcl_NO,1:4) = poly(13,1:4);
                % nehme Einheitsvektor
                richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
                richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
                % berechne die Normale
                % zun"achst zweiten Richtungsvektor
                dir2 = punkt(3,1:3) - punkt(1,1:3);
                dir2 = dir2/norm(dir2);
                % nun Formel f"ur Vektorproduct
		c = cross(dir2, richtung(pcl_NO,1:3));
		normale(pcl_NO,1:3) = c(1:3);
		normale(pcl_NO, 4) = pcl_NO;
                pcl_NO = pcl_NO + 1;
            end
            x3D(index-1,5) = pcl_NO-1;
	end      
%%%%
	%%%% CALOTE 4
	%%%%
        % Seitenfl"achen, die Ebene ist die gleiche f"ur alle
        % Fl"achen dieser Rotation
        if j< N0+N1+N2+N3 && i>1 
            % im ersten Punkt wird die Ebenengleichung berechnet
            % nehme daf"ur den berechneten Punkt, den der auch an
            % der oberen
            % Ebene liegt und einen Punkt der
            % vorherigen Rotation, das geht erst f"ur i>1
            if j==N0+N1+4
                punkt(1,1:3) = x3D(index-1,2:4);
                punkt(2,1:3) = x3D(9,2:4);
                punkt(3,1:3) = x3D(index-1-points+offset5,2:4);
                %punkt(4,1:3) = x3D(index-1-points+offset5,2:4);
                poly(6,1) = x3D(index-1,1);
                poly(6,2) = x3D(9,1);
                poly(6,3) = x3D(index-1-points+offset5,1);
                %poly(4,4) = x3D(index-1-points+offset,1);
                
                
                % von nun an werden immer zwei Punkte nicht ber"ucksichtg
                % in der Liste (0,0,0), (0,0,140)
                offset5 = 2;
                point(pcl_NO,1:3) = punkt(1,1:3);
                polygon(pcl_NO,1:3) = poly(6,1:3);
                % nehme Einheitsvektor
                richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
                richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
                % berechne die Normale
                % zun"achst zweiten Richtungsvektor
                dir2 = punkt(3,1:3) - punkt(1,1:3);
                dir2 = dir2/norm(dir2);
                % nun Formel f"ur Vektorproduct
		c = cross(dir2, richtung(pcl_NO,1:3));
		normale(pcl_NO,1:3) = c(1:3);
		normale(pcl_NO, 4) = pcl_NO;
                pcl_NO = pcl_NO + 1;
            end
            x3D(index-1,5) = pcl_NO-1;
	end

    end
end
pcl=pcl_NO;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Cylinder
%%%
% Die Ebenen zwischen der letzten und ersten Reihe m"ussen noch beschrieben
% werden
for j=N0+1:points
% Seitenfl"achen, die Ebene ist die gleiche f"ur alle
% Fl"achen dieser Rotation
   if j< N0+N1
   % im ersten Punkt wird die Ebenengleichung berechnet
   % nehme daf"ur den berechneten Punkt, den der auch an der oberen
   % Ebene liegt und einen Punkt der
   % vorherigen Rotation, das geht erst f"ur i>1
   index0 = j;
   if j==N0+1
       index1 = N0 + (rotations-1)*(points-2)+1;
       punkt(1,1:3) = x3D(index0-1,2:4);
       punkt(2,1:3) = x3D(index0,2:4);
       punkt(3,1:3) = x3D(index1,2:4);
       punkt(4,1:3) = x3D(index1,2:4);
       poly(7,1) = x3D(index0-1,1);
       poly(7,2) = x3D(index0,1);
       poly(7,3) = x3D(index1,1);
       poly(7,4) = x3D(index1+1,1);
       
       point(pcl_NO,1:3) = punkt(1,1:3);
       polygon(pcl_NO,1:4) = poly(7,1:4)
       % nehme Einheitsvektor
       richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
       richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
       % berechne die Normale
       % zun"achst zweiten Richtungsvektor           
       dir2 = punkt(3,1:3) - punkt(1,1:3);
       dir2 = dir2/norm(dir2);
       % nun Formel f"ur Vektorproduct
       c = cross(dir2, richtung(pcl_NO, 1:3));
       normale(pcl_NO,1:3) = c(1:3);
       normale(pcl_NO,4) = pcl_NO;
       pcl_NO = pcl_NO + 1;
   end
   x3D(index0,5) = pcl_NO-1; 
   end
   %%%
   %%% Torus 1
   %%%
% Die Ebenen zwischen der letzten und ersten Reihe m"ussen noch beschrieben
% werden

% Seitenfl"achen, die Ebene ist die gleiche f"ur alle
% Fl"achen dieser Rotation
   if j< N0+N1+N2
   % im ersten Punkt wird die Ebenengleichung berechnet
   % nehme daf"ur den berechneten Punkt, den der auch an der oberen
   % Ebene liegt und einen Punkt der
   % vorherigen Rotation, das geht erst f"ur i>1
   index0 = j;
   if j==N0+N1
       index1 = N0 + N1 - 1+ (rotations-1)*(points-2)+1;
       punkt(1,1:3) = x3D(index0-1,2:4);
       punkt(2,1:3) = x3D(index0,2:4);
       punkt(3,1:3) = x3D(index1,2:4);
       punkt(4,1:3) = x3D(index1+1,2:4);
       poly(8,1) = x3D(index0-1,1);
       poly(8,2) = x3D(index0,1);
       poly(8,3) = x3D(index1,1);
       poly(8,4) = x3D(index1+1,1);
       
       point(pcl_NO,1:3) = punkt(1,1:3);
       polygon(pcl_NO,1:4) = poly(8,1:4)
       % nehme Einheitsvektor
       richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
       richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
       % berechne die Normale
       % zun"achst zweiten Richtungsvektor           
       dir2 = punkt(3,1:3) - punkt(1,1:3);
       dir2 = dir2/norm(dir2);
       % nun Formel f"ur Vektorproduct
       c = cross(dir2, richtung(pcl_NO, 1:3));
       normale(pcl_NO,1:3) = c(1:3);
       normale(pcl_NO,4) = pcl_NO;
       pcl_NO = pcl_NO + 1;
   end
   x3D(index0,5) = pcl_NO-1; 
   end
   %%%
   %%% Torus 2
   %%%
% Die Ebenen zwischen der letzten und ersten Reihe m"ussen noch beschrieben
% werden

% Seitenfl"achen, die Ebene ist die gleiche f"ur alle
% Fl"achen dieser Rotation
   if j< N0+N1+N2
   % im ersten Punkt wird die Ebenengleichung berechnet
   % nehme daf"ur den berechneten Punkt, den der auch an der oberen
   % Ebene liegt und einen Punkt der
   % vorherigen Rotation, das geht erst f"ur i>1
   index0 = j;
   if j==N0+N1+1
       index1 = N0 + N1 + (rotations-1)*(points-2)+1;
       punkt(1,1:3) = x3D(index0-1,2:4);
       punkt(2,1:3) = x3D(index0,2:4);
       punkt(3,1:3) = x3D(index1,2:4);
       poly(9,1) = x3D(index0-1,1);
       poly(9,2) = x3D(index0,1);
       poly(9,3) = x3D(index1,1);
       poly(9,4) = x3D(index1+1,1);
       
       point(pcl_NO,1:3) = punkt(1,1:3);
       polygon(pcl_NO,1:4) = poly(9,1:4)
       % nehme Einheitsvektor
       richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
       richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
       % berechne die Normale
       % zun"achst zweiten Richtungsvektor           
       dir2 = punkt(3,1:3) - punkt(1,1:3);
       dir2 = dir2/norm(dir2);
       % nun Formel f"ur Vektorproduct
       c = cross(dir2, richtung(pcl_NO, 1:3));
       normale(pcl_NO,1:3) = c(1:3);
       normale(pcl_NO,4) = pcl_NO;
       pcl_NO = pcl_NO + 1;
   end
   x3D(index0,5) = pcl_NO-1; 
   end
   %%%
   %%% Calote 1
   %%%
% Die Ebenen zwischen der letzten und ersten Reihe m"ussen noch beschrieben
% werden

% Seitenfl"achen, die Ebene ist die gleiche f"ur alle
% Fl"achen dieser Rotation
   if j< N0+N1+N2+N3
   % im ersten Punkt wird die Ebenengleichung berechnet
   % nehme daf"ur den berechneten Punkt, den der auch an der oberen
   % Ebene liegt und einen Punkt der
   % vorherigen Rotation, das geht erst f"ur i>1
   index0 = j;
   if j==N0+N1+2
       index1 = N0 + N1 + 1 + (rotations-1)*(points-2)+1;
       punkt(1,1:3) = x3D(index0-1,2:4);
       punkt(2,1:3) = x3D(index0,2:4);
       punkt(3,1:3) = x3D(index1,2:4);
       poly(10,1) = x3D(index0-1,1);
       poly(10,2) = x3D(index0,1);
       poly(10,3) = x3D(index1,1);
       poly(10,4) = x3D(index1+1,1);
       
       
       point(pcl_NO,1:3) = punkt(1,1:3);
       % nehme Einheitsvektor
       polygon(pcl_NO,1:4) = poly(10,1:4)
       richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
       richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
       % berechne die Normale
       % zun"achst zweiten Richtungsvektor           
       dir2 = punkt(3,1:3) - punkt(1,1:3);
       dir2 = dir2/norm(dir2);
       % nun Formel f"ur Vektorproduct
       c = cross(dir2, richtung(pcl_NO, 1:3));
       normale(pcl_NO,1:3) = c(1:3);
       normale(pcl_NO,4) = pcl_NO;
       pcl_NO = pcl_NO + 1;
   end
   x3D(index0,5) = pcl_NO-1; 
   end
   %%%
   %%% Calote 2
   %%%
% Die Ebenen zwischen der letzten und ersten Reihe m"ussen noch beschrieben
% werden

% Seitenfl"achen, die Ebene ist die gleiche f"ur alle
% Fl"achen dieser Rotation
   if j< N0+N1+N2+N3
   % im ersten Punkt wird die Ebenengleichung berechnet
   % nehme daf"ur den berechneten Punkt, den der auch an der oberen
   % Ebene liegt und einen Punkt der
   % vorherigen Rotation, das geht erst f"ur i>1
   index0 = j;
   if j==N0 + N1 + 3
       index1 = N0 + N1 + 2 + (rotations-1)*(points-2)+1;
       punkt(1,1:3) = x3D(index0-1,2:4);
       punkt(2,1:3) = x3D(index0,2:4);
       punkt(3,1:3) = x3D(index1,2:4);
       poly(11,1) = x3D(index0-1,1);
       poly(11,2) = x3D(index0,1);
       poly(11,3) = x3D(index1,1);
       poly(11,4) = x3D(index1+1,1);
       
       point(pcl_NO,1:3) = punkt(1,1:3);
       polygon(pcl_NO,1:4) = poly(11,1:4)
       % nehme Einheitsvektor
       richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
       richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
       % berechne die Normale
       % zun"achst zweiten Richtungsvektor           
       dir2 = punkt(3,1:3) - punkt(1,1:3);
       dir2 = dir2/norm(dir2);
       % nun Formel f"ur Vektorproduct
       c = cross(dir2, richtung(pcl_NO, 1:3));
       normale(pcl_NO,1:3) = c(1:3);
       normale(pcl_NO,4) = pcl_NO;
       pcl_NO = pcl_NO + 1;
   end
   x3D(index0,5) = pcl_NO-1; 
   end
   %%%
   %%% Calote 3
   %%%
% Die Ebenen zwischen der letzten und ersten Reihe m"ussen noch beschrieben
% werden

% Seitenfl"achen, die Ebene ist die gleiche f"ur alle
% Fl"achen dieser Rotation
   if j< N0+N1+N2+N3
   % im ersten Punkt wird die Ebenengleichung berechnet
   % nehme daf"ur den berechneten Punkt, den der auch an der oberen
   % Ebene liegt und einen Punkt der
   % vorherigen Rotation, das geht erst f"ur i>1
   index0 = j;
   if j==N0 + N1 + 4
       index1 = N0 + N1 + 3 + (rotations-1)*(points-2)+1;
       punkt(1,1:3) = x3D(index0-1,2:4);
       punkt(2,1:3) = x3D(index0,2:4);
       punkt(3,1:3) = x3D(index1,2:4);
       poly(14,1) = x3D(index0-1,1);
       poly(14,2) = x3D(index0,1);
       poly(14,3) = x3D(index1,1);
       poly(14,4) = x3D(index1+1,1);
       
       point(pcl_NO,1:3) = punkt(1,1:3);
       polygon(pcl_NO,1:4) = poly(14,1:4)
       % nehme Einheitsvektor
       richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
       richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
       % berechne die Normale
       % zun"achst zweiten Richtungsvektor           
       dir2 = punkt(3,1:3) - punkt(1,1:3);
       dir2 = dir2/norm(dir2);
       % nun Formel f"ur Vektorproduct
       c = cross(dir2, richtung(pcl_NO, 1:3));
       normale(pcl_NO,1:3) = c(1:3);
       normale(pcl_NO,4) = pcl_NO;
       pcl_NO = pcl_NO + 1;
   end
   x3D(index0,5) = pcl_NO-1; 
   end   
   
   %%%
   %%% Calote 4
   %%%
% Die Ebenen zwischen der letzten und ersten Reihe m"ussen noch beschrieben
% werden

% Seitenfl"achen, die Ebene ist die gleiche f"ur alle
% Fl"achen dieser Rotation
   if j< N0+N1+N2+N3
   % im ersten Punkt wird die Ebenengleichung berechnet
   % nehme daf"ur den berechneten Punkt, den der auch an der oberen
   % Ebene liegt und einen Punkt der
   % vorherigen Rotation, das geht erst f"ur i>1
   index0 = j;
   if j==N0+N1+4
       index1 = N0 + N1 + 3 + (rotations-1)*(points-2)+1;
       punkt(1,1:3) = x3D(index0,2:4);
       punkt(2,1:3) = x3D(index0+1,2:4);
       punkt(3,1:3) = x3D(index1+1,2:4);
       poly(12,1) = x3D(index0,1);
       poly(12,2) = x3D(index0+1,1);
       poly(12,3) = x3D(index1+1,1);
       
       point(pcl_NO,1:3) = punkt(1,1:3);
       polygon(pcl_NO,1:3) = poly(12,1:3)
       % nehme Einheitsvektor
       richtung(pcl_NO,1:3) = punkt(2,1:3) - punkt(1,1:3);
       richtung(pcl_NO,1:3) = richtung(pcl_NO,1:3)/norm(richtung(pcl_NO,1:3));
       % berechne die Normale
       % zun"achst zweiten Richtungsvektor           
       dir2 = punkt(3,1:3) - punkt(1,1:3);
       dir2 = dir2/norm(dir2);
       % nun Formel f"ur Vektorproduct
       c = cross(dir2, richtung(pcl_NO, 1:3));
       normale(pcl_NO,1:3) = c(1:3);
       normale(pcl_NO,4) = pcl_NO;
       pcl_NO = pcl_NO + 1;
   end
   x3D(index0,5) = pcl_NO-1; 
   end



   
end











%x3D
%plot3(x3D(:,2),x3D(:,3),x3D(:,4),'-+');

% schreibe die Daten weg
fid=fopen('vierteile1.poly','w');

s = sprintf('%d %d %d %d \n\n',index-2,3,0,1);
  fprintf(fid,'%s',s);

for i=2:index-1
  s = sprintf('%d %g %g %g %d\n',x3D(i,1),x3D(i,2),x3D(i,3),x3D(i,4),x3D(i,5));
  fprintf(fid,'%s',s);
end
%%%%%%%%% %%


s = sprintf('\n %d %d \n',pcl_NO-1,1);
  fprintf(fid,'%s',s);
for i=1:pcl_NO-1
  s = sprintf('%d %d %d \n',1,0,i);
  fprintf(fid,'%s',s);
  if polygon(i,4)==0
        s = sprintf('%d %d %d %d \n', 3, polygon(i,1),polygon(i,2),polygon(i,3));
        fprintf(fid,'%s',s);
    elseif polygon(i,5)==0
      
        s = sprintf('%d %d %d %d %d\n' ,4, polygon(i,3),polygon(i,1),polygon(i,2),polygon(i,4));
        fprintf(fid,'%s',s);
  else
    
     s = sprintf('%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n' ,16, polygon(i,1),polygon(i,2),polygon(i,3),polygon(i,4), polygon(i,5),polygon(i,6),polygon(i,7),polygon(i,8), polygon(i,9),polygon(i,10),polygon(i,11),polygon(i,12), polygon(i,13),polygon(i,14),polygon(i,15),polygon(i,16));
        fprintf(fid,'%s',s);
end
end
s = sprintf('\n %d \n',0);
  fprintf(fid,'%s',s);
s = sprintf('\n %d \n',0);
  fprintf(fid,'%s',s);

fclose(fid);

%%%%%%%%%%
fid=fopen('vierteile1.prm','w');
s = sprintf('%s\n','NBCT');
fprintf(fid,'%s',s);
s = sprintf('%d\n',1);
fprintf(fid,'%s',s);
s = sprintf('%s\n','IBCT');
fprintf(fid,'%s',s);
s = sprintf('%d\n',1);
fprintf(fid,'%s',s);
s = sprintf('%s\n','NCOMP');
fprintf(fid,'%s',s);
s = sprintf('%d\n',pcl_NO-1);
fprintf(fid,'%s',s);
s = sprintf('%s %s %s\n','ITYP','NSPLINE','NPAR');
fprintf(fid,'%s',s);
for i=1:pcl_NO-1
    s = sprintf(' %d\t\t%d\t\t%d\n',10,0,3);
 fprintf(fid,'%s',s);
end   
s = sprintf('%s\n','PARAMETERS');
fprintf(fid,'%s',s);

for i=1:pcl_NO-1
  s = sprintf('%g %g %g \n',point(i,1),point(i,2),point(i,3));
  fprintf(fid,'%s',s);
  s = sprintf('%g %g %g\n',richtung(i,1),richtung(i,2),richtung(i,3));
  fprintf(fid,'%s',s);
  s = sprintf('%g %g %g\n',normale(i,1),normale(i,2),normale(i,3));
  fprintf(fid,'%s',s);
end


fclose(fid);
