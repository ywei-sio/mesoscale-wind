% An example of solving wind stress fields from divergence and curl
% Referring to:
% Wei et al (2019), Mesoscale wind stress-SST coupled perturbations in the Kuroshio Extension, Progress in Oceanography, 172, 108¨C123.
clear;clc;close all
nx=562;ny=242; nt=90;
indexKE(nx,ny)=0;
indexNum(nx,ny)=0;

lon=99.875+(0:561)*0.25;
lat=-0.125+(0:241)*0.25;

r=6371;% radius of earth [km]
pi=3.1415926;
delta_y(1:nx,1:ny)=r*0.25*pi/180; % [km]

delta_x(1:nx,1:ny)=0;
for i=1:nx
    delta_x(i,1:length(lat))=r*cos(lat/180*pi)*0.25*pi/180; % [km]
end

load indexKE
N=sum(sum(indexKE));%

num=0;
M=0;%
for i=1:nx
    for j=1:ny
        if indexKE(i,j)
            num=num+1;
            indexNum(i,j)=num;
            if indexKE(i-1,j)&&indexKE(i+1,j)&&indexKE(i,j-1)&&indexKE(i,j+1)
                M=M+1;
            end
        end
    end
end

%% ---------grid points---------------------------
h1=figure('Position',[10,50,1000 400]);
set(gca, 'color', 'w','XAxisLocation','bottom');
m_proj('miller','longitudes',[140 190], 'latitudes', [32 45]);hold on;
m_grid('tickdir', 'out','xtick',140:10:190,'fontsize',12,'ytick', 30:5:45,'fontsize',12,'box','on');
m_coast('patch','k','edgecolor',[1 1 1]*0);
hold on;
for i=1:nx
    for j=1:ny
        if indexKE(i,j)
            if indexKE(i-1,j)&&indexKE(i+1,j)&&indexKE(i,j-1)&&indexKE(i,j+1)
                m_plot(lon(i),lat(j),'o','MarkerEdgeColor','w','MarkerFaceColor','b','MarkerSize',3);
            else
                m_plot(lon(i),lat(j),'o','MarkerEdgeColor','w','MarkerFaceColor','r','MarkerSize',3);
            end
        end
    end
end

%% weight function
% can be used to ensure the wind stress variation in the boundaries smooth
w(1:nx,1:ny)=0;
lon1=190;lat1=45;
lon2=190;lat2=32;
lon3=140;lat3=32;
lon4=140;lat4=45;

%---line equation: (x-x1)/(x2-x1)=(y-y1)/(y2-y1)-------
%---line equation: (x-x1)*(y2-y1)-(y-y1)*(x2-x1)=0-------
s=0.5;
for i=1:nx
    for j=1:ny
        x=lon(i);y=lat(j);
        %         line1
        x1=lon1;x2=lon2;y1=lat1;y2=lat2;
        d1=( (x-x1)*(y2-y1)-(y-y1)*(x2-x1) )/sqrt( (y2-y1)^2+(x2-x1)^2);
        %         line2
        x1=lon2;x2=lon3;y1=lat2;y2=lat3;
        d2=( (x-x1)*(y2-y1)-(y-y1)*(x2-x1) )/sqrt( (y2-y1)^2+(x2-x1)^2);
        %         line3
        x1=lon3;x2=lon4;y1=lat3;y2=lat4;
        d3=( (x-x1)*(y2-y1)-(y-y1)*(x2-x1) )/sqrt( (y2-y1)^2+(x2-x1)^2);
        %         line4
        x1=lon4;x2=lon1;y1=lat4;y2=lat1;
        d4=( (x-x1)*(y2-y1)-(y-y1)*(x2-x1) )/sqrt( (y2-y1)^2+(x2-x1)^2);
        
        if (d1>0)&&(d2>0)&&(d3>0)&&(d4>0)
            if indexKE(i,j)
                w(i,j)=(1-exp(-d1^2/s))*(1-exp(-d2^2/s))*(1-exp(-d3^2/s/20))*(1-exp(-d4^2/s));
            end
        end
    end
end

h2=figure('Position',[10,50,1000 400]);
m_proj('miller','longitudes',[140 190], 'latitudes', [32 45]);hold on;
m_grid('tickdir', 'out','xtick',140:10:190,'fontsize',12,'ytick', 30:5:45,'fontsize',12,'box','on');
h=m_pcolor(lon,lat,w');set(h,'LineStyle','none');
m_coast('patch','w','edgecolor',[1 1 1]*0);
colorbar;

%%
load KE_Wind_Jan2003
%--------A*y=c--------------
%------- y=(taux, tauy)-------------
%------- c=(div, curl)-------------
%--------y=A_inv*c---------------
% unit of div and curl: N/m^2 per 10 000km
% unit of taux and tauy: N/m^2
% unit of matrix A: 10 1000km

A(1:2*M,1:2*N)=0;
%----------construct matrix A-----------------
num=0;
for i=1:nx
    for j=1:ny
        if indexKE(i,j)
            
            if indexKE(i-1,j)&&indexKE(i+1,j)&&indexKE(i,j-1)&&indexKE(i,j+1)
                num=num+1;
                %------------------------div---------
                % [taux(i,j)-taux(i-1,j)]/delta_x(i,j)+[tauy(i,j)-tauy(i,j-1)]/delta_y(i,j)
                if indexNum(i,j)>0
                    A(num,indexNum(i,j))=1/delta_x(i,j)*10000;
                end
                
                if indexNum(i-1,j)>0
                    A(num,indexNum(i-1,j))=-1/delta_x(i-1,j)*10000;
                end
                
                if indexNum(i,j)>0
                    A(num,N+indexNum(i,j))=1/delta_y(i,j)*10000;
                end
                
                if indexNum(i,j-1)>0
                    A(num,N+indexNum(i,j-1))=-1/delta_y(i,j-1)*10000;
                end
                
                %------------------------curl---------
                % dv/dx-du/dy
                % [tauy(i,j)-tauy(i-1,j)]/delta_x(i,j)-[taux(i,j)-taux(i,j-1)]/delta_y(i,j)
                if indexNum(i,j)>0
                    A(num+M,N+indexNum(i,j))=1/delta_x(i,j)*10000;
                end
                
                if indexNum(i-1,j)>0
                    A(M+num,N+indexNum(i-1,j))=-1/delta_x(i-1,j)*10000;
                end
                
                if indexNum(i,j)>0
                    A(num+M,indexNum(i,j))=-1/delta_y(i,j)*10000;
                end
                
                if indexNum(i,j-1)>0
                    A(num+M,indexNum(i,j-1))=1/delta_y(i,j-1)*10000;
                end
                
            end
        end
    end
end

clear delta_x delta_y
alpha=10;
A_inv=A'*inv(A*A'+alpha*(eye(2*M)));

%%
%-------------- c--------------------
c_div(1:2*M,1:1)=0;% divergence only
c_curl(1:2*M,1:1)=0;% curl only
num=0;
for i=1:nx
    for j=1:ny
        if indexKE(i,j)
            if indexKE(i-1,j)&&indexKE(i+1,j)&&indexKE(i,j-1)&&indexKE(i,j+1)
                num=num+1;
                c_div(num)=wind_div(i,j);
                c_div(num+M)=wind_curl(i,j)*0;
                
                c_curl(num)=wind_div(i,j)*0;
                c_curl(num+M)=wind_curl(i,j);
            end
        end
    end
end

%-----solve the wind stress field-------------------
y_div=A_inv*c_div;
y_curl=A_inv*c_curl;

taux_div(1:nx,1:ny)=0;
tauy_div(1:nx,1:ny)=0;

taux_curl(1:nx,1:ny)=0;
tauy_curl(1:nx,1:ny)=0;

num=0;
for i=1:nx
    for j=1:ny
        if indexKE(i,j)
            num=num+1;
            taux_div(i,j)=y_div(num);
            tauy_div(i,j)=y_div(num+N);
            
            taux_curl(i,j)=y_curl(num);
            tauy_curl(i,j)=y_curl(num+N);
        end
    end
end

%% --------plot----------------------
map_type=1;
LonTick=120:10:200;
LatTick=30:5:50;
delta=2;
scale=0.0015;

h3=figure('position',[100,50,1000,900]);
%
subplot('position',[.1 .7 .8 .25]);box on;hold on;
set(gca, 'color', 'w','XAxisLocation','bottom');
m_proj('miller','longitudes',[140 190], 'latitudes', [32 45]);
m_grid('tickdir', 'out','xtick',LonTick,'ytick', LatTick,'box','on','linestyle','none','fontsize',12);
m_coast('patch',[.7 .7 .7],'edgecolor',[0 0 0],'linewidth',.5);
grid off
for i=1:nx/delta
    for j=1:ny/delta
        if indexKE(i*delta,j*delta)
            fun_vector(lon(i*delta),lat(j*delta),taux_obs(i*delta,j*delta)/scale,tauy_obs(delta*i,delta*j)/scale,[1 1 1]*0,0.5,map_type);hold on
        end
    end
end
fun_vector(142.0,43.5,0.1/scale,0/scale,[1 1 1]*0,1,map_type);
m_text(141.5,42.9,'0.1 Nm^{-2}','fontweight','bold','fontsize',12);
m_text(141.5,33,'(a)','fontweight','bold','fontsize',12);
%
subplot('position',[.1 .4 .8 .25]);box on; hold on
set(gca, 'color', 'w','XAxisLocation','bottom');
m_proj('miller','longitudes',[140 190], 'latitudes', [32 45]);
m_grid('tickdir', 'out','xtick',LonTick,'ytick',LatTick,'box','on','linestyle','none','fontsize',12);
m_coast('patch',[.7 .7 .7],'edgecolor',[0 0 0],'linewidth',.5);
grid off
for i=1:nx/delta
    for j=1:ny/delta
        if indexKE(i*delta,j*delta)
            fun_vector(lon(i*delta),lat(j*delta),taux_div(i*delta,j*delta)/scale,tauy_div(i*delta,j*delta)/scale,[1 1 1]*0,0.5,map_type);
        end
    end
end
fun_vector(142.0,43.5,0.1/scale,0/scale,[1 1 1]*0,1,map_type);
m_text(141.5,42.9,'0.1 Nm^{-2}','fontweight','bold','fontsize',12);
m_text(141.5,33,'(b)','fontweight','bold','fontsize',12);
%
subplot('position',[.1 .1 .8 .25]);box on;hold on
set(gca, 'color', 'w','XAxisLocation','bottom');
m_proj('miller','longitudes',[140 190], 'latitudes', [32 45]);hold on;
m_grid('tickdir', 'out','xtick',LonTick,'ytick', LatTick,'box','on','linestyle','none','fontsize',12);
m_coast('patch',[.7 .7 .7],'edgecolor',[0 0 0],'linewidth',.5);
grid off
for i=1:nx/delta
    for j=1:ny/delta
        if indexKE(i*delta,j*delta)
            fun_vector(lon(i*delta),lat(j*delta),(taux_curl(i*delta,j*delta))/scale,(tauy_curl(delta*i,delta*j))/scale,[1 1 1]*0,0.5,map_type);
        end
    end
end
fun_vector(142.0,43.5,0.1/scale,0/scale,[1 1 1]*0,1,map_type);
m_text(141.5,42.9,'0.1 Nm^{-2}','fontweight','bold','fontsize',12);
m_text(141.5,33,'(c)','fontweight','bold','fontsize',12);
