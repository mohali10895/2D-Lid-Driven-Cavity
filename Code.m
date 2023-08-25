clc
clear  
close all

%% Inputs
U = 5;
L = 0.2;
no_points = 81;
nu = 0.001;                     
Re=U*L/nu;

%% SETUP 2D GRID
delta_x_y=L/(no_points-1);  
X=zeros(no_points);
Y=zeros(no_points);
for j=1:no_points
    if j<no_points
        X(:,j+1)=j*delta_x_y;
    end
end
for i=1:no_points
    if i<no_points
        Y(i+1,:)=i*delta_x_y;
    end
end
Y=flip(Y);

%% Variables and Counters
omega_new = zeros(no_points,no_points);  
ebsia = omega_new;  
omega_old = omega_new;  
u = omega_new;  
v = omega_new;
im = 1:no_points-2;  
ip = 3:no_points;  
jm = 1:no_points-2;  
jp = 3:no_points;
i = 2:no_points-1;
j = 2:no_points-1;

%% Solver
dt = 0.0001;  
maxIt = 1000000;  
maxerr = 1e-7; % Time Step; Max iter; Max error
for iter = 1:maxIt 
    omega_old = omega_new;
    
    %%% BOUNDARY CONDITIONS
    omega_new(1:no_points,no_points) = -2*ebsia(1:no_points,no_points-1)/(delta_x_y^2) - U*2/delta_x_y; % Top
    omega_new(1:no_points,1)  = -2*ebsia(1:no_points,2)/(delta_x_y^2);                                  % Bottom
    omega_new(1,1:no_points)  = -2*ebsia(2,1:no_points)/(delta_x_y^2);                                  % Left
    omega_new(no_points,1:no_points) = -2*ebsia(no_points-1,1:no_points)/(delta_x_y^2);                 % Right
    
    %%%VORTICITY TRANSPORT EQUATION
    omega_new(i,j) = omega_old(i,j) + ...
            ((-1*(ebsia(i,jp)-ebsia(i,jm))/(2*delta_x_y)).* ((omega_old(ip,i)-omega_old(im,j))/(2*delta_x_y))+...
            ((ebsia(ip,j)-ebsia(im,j))/(2*delta_x_y)).*((omega_old(i,jp)-omega_old(i,jm))/(2*delta_x_y))+...
             nu*(omega_old(ip,j)+omega_old(im,j)-4*omega_old(i,j)+omega_old(i,jp)+omega_old(i,jm))/(delta_x_y^2))*dt;
         
    %%% EQUATION FOR STREAM FUNCTION     
    ebsia(i,j) = (omega_new(i,j)*delta_x_y^2 + ebsia(ip,j) + ebsia(i,jp) + ebsia(i,jm) + ebsia(im,j))/4;
    
    %%% CHECK FOR CONVERGENCE
    if iter > 10
        error = max(max(omega_new - omega_old));
        if error < maxerr
            break;
        end
    end
    if(rem(iter, 1000)) == 0
       figure(1);
       semilogy(iter, error, '-ko')
       hold on
       xlabel('Iterations')
       ylabel('Residual Error')
       grid on
    end
end
 
%% VELOCITY FROM STREAM
u(2:no_points-1,no_points) = U;
for n = 2:no_points-1
    for m = 2:no_points-1
        u(n,m)=(ebsia(n,m+1)-ebsia(n,m-1))/(2*delta_x_y);  
        v(n,m)=(-ebsia(n+1,m)+ebsia(n-1,m))/(2*delta_x_y);
    end
end


%% Figures
% Grid Generation
figure
plot(X,Y,"k")
hold on
plot(Y,X,"k")
xlabel('X (m)')
ylabel('Y (m)')
title('Grid Generation')

% Streamline pattern for the flow inside the cavity.
figure 
x = X(1,:);  
y = flip(Y(:,1));
N = 1000;  
xstart = max(x)*rand(N,1);  
ystart = max(y)*rand(N,1);
[x,y] = meshgrid(x,y);
delta_x_y=streamline(x,y,u',v',xstart,ystart,[0.1, 200]);
title('Streamline pattern for the flow inside the cavity')
xlabel('X (m)')
ylabel('Y (m)')
grid on
set(delta_x_y,'color','k')


u=flip(u');
v=flip(v');
x = X(1,:);  
y = flip(Y(:,1));
omega=flip(omega_new');
mid=(no_points+1)/2;
u_mid_ver=u(:,mid);
u_mid_hor=u(mid,:);
v_mid_ver=v(:,mid);
v_mid_hor=v(mid,:);
omega_mid_ver=omega(:,mid);
omega_mid_hor=omega(mid,:);


%u, v, and ? along the vertical centerline of the cavity.
figure
plot(u_mid_ver,y)
xlabel('U (m/s)')
ylabel('Y (m)')
title('U along the vertical centerline of the cavity')
grid on

figure
plot(v_mid_ver,y)
xlabel('V (m/s)')
ylabel('Y (m)')
title('V along the vertical centerline of the cavity')
grid on

figure
plot(omega_mid_ver,y)
xlabel('{\omega} (s^{-1})')
ylabel('Y (m)')
title('{\omega} along the vertical centerline of the cavity')
grid on

%u, v, and ? along the horizontal centerline of the cavity.
figure
plot(x,u_mid_hor)
xlabel('X (m)')
ylabel('U (m/s)')
title('U along the horizontal centerline of the cavity')
grid on

figure
plot(x,v_mid_hor)
xlabel('X (m)')
ylabel('V (m/s)')
title('V along the horizontal centerline of the cavity')
grid on

figure
plot(x,omega_mid_hor)
xlabel('X (m)')
ylabel('{\omega} (s^{-1})')
title('{\omega} along the horizontal centerline of the cavity')
grid on

%Vorticity contours in the cavity.
figure
contourf(X,Y,omega,'LineStyle', 'none')
xlabel('X (m)')
ylabel('Y (m)')
title('Vorticity contours in the cavity')
colorbar

%Vorticity along the walls
omega_bottom=omega(no_points,:);
omega_left=omega(:,1);
omega_right=omega(:,no_points);
omega_bottom_mat=zeros(no_points);
omega_left_mat=zeros(no_points);
omega_right_mat=zeros(no_points);
for k=1:no_points
    omega_bottom_mat(k,:)=omega_bottom;
    omega_left_mat(:,k)=omega_left;
    omega_right_mat(:,k)=omega_right;
end

figure
surf(X,Y,omega_bottom_mat)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('{\omega} (s^{-1})')
title('Vorticity along the bottom wall')

figure
surf(X,Y,omega_left_mat)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('{\omega} (s^{-1})')
title('Vorticity along the left wall')

figure
surf(X,Y,omega_right_mat)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('{\omega} (s^{-1})')
title('Vorticity along the right wall')


for o=1:13
    fig=figure(o);
    filename = sprintf('fig%1d.png', o);
    saveas(gcf,filename);
end
