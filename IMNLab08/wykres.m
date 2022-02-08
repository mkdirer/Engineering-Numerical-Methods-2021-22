clc;clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%vx%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
data = readmatrix('vx.dat');
dataM = data.';
h = pcolor(dataM);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
caxis([-5 45]);
title('$v_{x}$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%vy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
data = readmatrix('vy.dat');
dataM = data.';
h = pcolor(dataM);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
caxis([-20 20]);
title('$v_{y}$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%c, xsr%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = readtable('c_xsr.dat');

figure;

plot(data.Var1(1:10500), data.Var2(1:10500), 'r',  'LineWidth', 4);
hold on;
plot(data.Var1(1:10500), data.Var3(1:10500), 'm',  'LineWidth', 4);
hold on;
plot(data.Var1(10501:21000), data.Var2(10501:21000), 'b',  'LineWidth', 4);
hold on;
plot(data.Var1(10501:21000), data.Var3(10501:21000), 'g',  'LineWidth', 4);
legend({'$c(D=0)$', '$x_{sr}(D=0)$', '$c(D=0.1)$', '$x_{sr}(D=0.1)$'},...
    'Location','best','Orientation','vertical','FontSize', 20, 'FontWeight', 'bold','Interpreter','latex');
title('$c, x_{sr}$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');
% %%%%%%%%%%%%%%%%%%%%%%%%%%mapa_0.0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = readmatrix('mapa_0.dat');
data2000 = zeros(401, 91);
data4000 = zeros(401, 91);
data6000 = zeros(401, 91);
data8000 = zeros(401, 91);
data10000 = zeros(401, 91);
for i = 1:401
    for j = 1:91
        data2000(i,j) = data(i,j);
        data4000(i,j) = data(i+401,j);
        data6000(i,j) = data(i+802,j);
        data8000(i,j) = data(i+1203,j);
        data10000(i,j) = data(i + 1604,j);        
    end
end
dataM1 = data2000.';
dataM2 = data4000.';
dataM3 = data6000.';
dataM4 = data8000.';
dataM5 = data10000.';

figure;
subplot(3,2,1);
h = pcolor(dataM1);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
caxis([0 20]);
title('$D = 0.0$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,2);
h = pcolor(dataM2);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
caxis([0 20]);
title('$D = 0.0$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,3);
h = pcolor(dataM3);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
caxis([0 20]);
title('$D = 0.0$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,4);
h = pcolor(dataM4);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
caxis([0 20]);
title('$D = 0.0$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,5);
h = pcolor(dataM5);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
caxis([0 20]);
title('$D = 0.0$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');


% %%%%%%%%%%%%%%%%%%%%%%%%%%mapa_0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = readmatrix('mapa_0.1.dat');
data2000 = zeros(401, 91);
data4000 = zeros(401, 91);
data6000 = zeros(401, 91);
data8000 = zeros(401, 91);
data10000 = zeros(401, 91);
for i = 1:401
    for j = 1:91
        data2000(i,j) = data(i,j);
        data4000(i,j) = data(i+401,j);
        data6000(i,j) = data(i+802,j);
        data8000(i,j) = data(i+1203,j);
        data10000(i,j) = data(i + 1604,j);        
    end
end
dataM1 = data2000.';
dataM2 = data4000.';
dataM3 = data6000.';
dataM4 = data8000.';
dataM5 = data10000.';
figure;
subplot(3,2,1);
h = pcolor(dataM1);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
title('$D = 0.1$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,2);
h = pcolor(dataM2);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
title('$D = 0.1$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,3);
h = pcolor(dataM3);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
title('$D = 0.1$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,4);
h = pcolor(dataM4);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
title('$D = 0.1$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');

subplot(3,2,5);
h = pcolor(dataM5);
set(h, 'EdgeColor', 'none');
colormap(gca,'jet');
colorbar
title('$D = 0.1$','Interpreter','latex');
xlabel('x','Interpreter','latex'); 
ylabel('y','Interpreter','latex');




