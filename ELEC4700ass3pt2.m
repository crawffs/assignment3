%Assignment 3 part 2
%Mike Crawford 100952432
x = 30;
y = 20;

conductOut = 1;
conductBox = 01e-2;
conduct = zeros(x,y);
G = sparse (x*y, x*y);
V = zeros(1, x*y);
Voltagex = 0.1;

for i = 1:x
    for j = 1:y
        %sets boundaries of boxes
        if (i > (0.3*x) || i < (0.6*x)) && (j > (0.6*y) || j < (0.3*y))
            conduct(i,j) = conductBox;
        else
            conduct(i,j) = conductOut;
        end
        
    end
end

for i = 1:x
    for j = 1:y
        n = j + (i - 1)*y;
        nx1 = j + ((i-1) - 1)*y;
        nx2 = j + ((i+1) - 1)*y;
        ny1 = (j-1) + (i - 1)*y;
        ny2 = (j+1) + (i - 1)*y;
       
        if i == 1
            V(n) = Voltagex;
            G(n,n) = 1;
        elseif i == x
            V(n) = 0;
            G(n,n) =1;
        elseif j == 1 
            Cx1 = (conduct(i,j) + conduct(i-1,j))/2;
            Cx2 = (conduct(i,j) + conduct(i+1,j))/2;
        
            Cy2 = (conduct(i,j) + conduct(i,j+1))/2;
            G(n,n) = -Cx1 - Cx2 - Cy2;
            G(n, nx1) = Cx1;
            G(n,nx2) = Cx2;
            G(n, ny2) = Cy2;
        elseif j == y
             Cx1 = (conduct(i,j) + conduct(i-1,j))/2;
             Cx2 = (conduct(i,j) + conduct(i+1,j))/2;
             Cy1 = (conduct(i,j) + conduct(i,j-1))/2;
            
            G(n,n) = -Cx1 - Cx2 - Cy1;
            G(n,nx1) = Cx1;
            G(n,nx2) = Cx2;
            G(n,ny1) = Cy1;
        else
             Cx1 = (conduct(i,j) + conduct(i-1,j))/2;
            Cx2 = (conduct(i,j) + conduct(i+1,j))/2;
            Cy1 = (conduct(i,j) + conduct(i,j-1))/2;
            Cy2 = (conduct(i,j) + conduct(i,j+1))/2;
            G(n,n) = -Cx1 - Cx2 - Cy1 - Cy2;
            G(n,nx1) = Cx1;
            G(n,nx2) = Cx2;
            G(n,ny1) = Cy1;
            G(n,ny2) = Cy2;
        end
    end
end
S = G\V';
surface = zeros(x,y);
for i = 1:x
    for j = 1:y
        n = j + (i - 1)*y;
        surface(i,j) = S(n);
    end
end
figure(1)
surf(surface)
title('Surface plot of V(x,y)')

[Ex, Ey] = gradient(-surface);
figure(2)
quiver(Ex, Ey)
title('Electric field quiver plot')

