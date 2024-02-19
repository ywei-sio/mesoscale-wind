function fun_vector(lon,lat,U,V,s,lw,map_type)
% for vector plot
% scalfac=0.012;
scalfac=0.036;
HEADA=5*pi/180;
HEADL=.75;
vel=U+sqrt(-1)*V;
ig=isfinite(vel);
x=lon(ig);
y=lat(ig);
vel=vel(ig);
igg=find(vel~=0);
x=x(igg);
y=y(igg);
vel=vel(igg);
z=x(:)+sqrt(-1)*y(:);
vel=vel(:)*scalfac;
r=vel*HEADL;
wr1=r*exp(+sqrt(-1)*HEADA);
wr2=r*exp(-sqrt(-1)*HEADA);

w2plot=ones(length(z),2);
w2plot(:,1)=z;w2plot(:,2)=z+vel;

w1plot=ones(length(z),3);
w1plot(:,1)=z+wr1;w1plot(:,2)=z+vel;w1plot(:,3)=z+wr2;

if map_type==1
    for i=1:length(z)
        m_plot(real(w2plot(i,:)),imag(w2plot(i,:)),'color',s,'linewidth',lw);
        m_plot(real(w1plot(i,:)),imag(w1plot(i,:)),'color',s,'linewidth',lw);
    end
else
    for i=1:length(z)
        plot(real(w2plot(i,:)),imag(w2plot(i,:)),'color',s,'linewidth',lw);
        plot(real(w1plot(i,:)),imag(w1plot(i,:)),'color',s,'linewidth',lw);
    end
end

