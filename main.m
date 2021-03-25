clear;
clc;
global gamma;
global p_L;
global p_R;
global u_L;
global u_R;
global rho_L;
global rho_R;
gamma=1.4;

% Sod 激波管初值条件
rho_L=1.0;
  u_L=0.0;
  p_L=1.0;
rho_R=0.125;
  u_R=0.0;
  p_R=0.1;

% 时间、空间范围离散
T=0.5;
dT=1.e-2;
t=dT:dT:T;

L=1;
dL=1.e-3;
x=-L:dL:L;

p  =zeros(length(t),length(x));
u  =zeros(length(t),length(x));
rho=zeros(length(t),length(x));

% 赋初始值
p(1,x<=0)=p_L;
p(1, x>0)=p_R;
u(1,x<=0)=u_L;
u(1, x>0)=u_R;
rho(1,x<=0)=rho_L;
rho(1, x>0)=rho_R;

% 求解接触间断的压强、速度
[p_,~]=solve_F(u_L-u_R,p_L,p_R,rho_L,rho_R);
u_=u_L-f(p_,p_L,rho_L);

% 求解接触间断的密度,行波速度
xs_R=0;
xr_R1=0;
xr_R2=0;
xs_L=0;
xr_L1=0;
xr_L2=0;
x_gap=0;
for i=1:length(t)
    
    if p_>p_R
        % 右行波为激波
        rhoR=(p_R-p_)*rho_R/((p_R-p_)+(u_R-u_)^2*rho_R);
        vs_R=(rhoR*u_-rho_R*u_R)/(rhoR-rho_R);
        xs_R=xs_R+vs_R*dT;
        p(i,x>=xs_R)=p_R;
        u(i,x>=xs_R)=u_R;
        rho(i,x>=xs_R)=rho_R;
        x_inner_R=xs_R;
    elseif p_<p_R
        % 右行波为膨胀波
        rhoR=rho_R*((p_/p_R)^(1/gamma));
        vr_R1=u_R+sqrt(gamma*p_R/rho_R);
        xr_R1=xr_R1+vr_R1*dT;
        vr_R2=u_+sqrt(gamma*p_/rhoR);
        xr_R2=xr_R2+vr_R2*dT;
        xr_R_head=max([xr_R1,xr_R2]);
        xr_R_tail=min([xr_R1,xr_R2]);
        p(i,x>=xr_R_head)=p_R;
        u(i,x>=xr_R_head)=u_R;
        rho(i,x>=xr_R_head)=rho_R;
        x_inner_R=min([xr_R1,xr_R2]);
        
        if (xr_R_head-xr_R_tail)~=0
            rare=and(x<xr_R_head,x>xr_R_tail);
            c=(gamma-1)/(gamma+1)*(-u_R+x(rare)/t(i))+2/(gamma+1)*sqrt(gamma*p_R/rho_R);
            u(i,rare)=-c+x(rare)/t(i);
            p(i,rare)=p_R*(c/sqrt(gamma*p_R/rho_R)).^(2*gamma/(gamma-1));
            rho(i,rare)=gamma*p(i,rare)./(c.^2);
        end
    end

    if p_>p_L
        % 左行波为激波
        rhoL=(p_L-p_)*rho_L/((p_L-p_)+(u_L-u_)^2*rho_L);
        vs_L=(rhoL*u_-rho_L*u_L)/(rhoL-rho_L);
        xs_L=xs_L+vs_L*dT;
        p(i,x<=xs_L)=p_L;
        u(i,x<=xs_L)=u_L;
        rho(i,x<=xs_L)=rho_L;
        x_inner_L=xs_L;
    elseif p_<p_L
        % 左行波为膨胀波
        rhoL=rho_L*((p_/p_L)^(1/gamma));
        vr_L1=u_L-sqrt(gamma*p_L/rho_L);
        vr_L2=u_-sqrt(gamma*p_/rhoL);
        xr_L2=xr_L2+vr_L2*dT;
        xr_L1=xr_L1+vr_L1*dT;
        xr_L_head=min([xr_L1,xr_L2]);
        xr_L_tail=max([xr_L1,xr_L2]);
        p(i,x<=xr_L_head)=p_L;
        u(i,x<=xr_L_head)=u_L;
        rho(i,x<=xr_L_head)=rho_L;
        x_inner_L=max([xr_L1,xr_L2]);
        
        if (xr_L_head-xr_L_tail)~=0
            rare=and(x>xr_L_head,x<xr_L_tail);
            c=(gamma-1)/(gamma+1)*(u_L-x(rare)/t(i))+2/(gamma+1)*sqrt(gamma*p_L/rho_L);
            u(i,rare)=c+x(rare)/t(i);
            p(i,rare)=p_L*(c/sqrt(gamma*p_L/rho_L)).^(2*gamma/(gamma-1));
            rho(i,rare)=gamma*p(i,rare)./(c.^2);
        end
    end
    
    x_gap=x_gap+u_*dT;
    
    p(i,and(x>x_inner_L,x<x_inner_R))=p_;
    u(i,and(x>x_inner_L,x<x_inner_R))=u_;
    rho(i,and(x>x_gap,x<x_inner_R))=rhoR;
    rho(i,and(x<x_gap,x>x_inner_L))=rhoL;
    
end

for i=1:length(t)
    figure(1);
    color=[54,100,139]/255;
    x1=min(x);
    x2=max(x);
    y1=min(rho(1,:))-0.1*(max(rho(1,:))-min(rho(1,:)));
    y2=max(rho(1,:))+0.1*(max(rho(1,:))-min(rho(1,:)));
    plot(x,rho(i,:),'.-','color',color);
    axis([x1,x2,y1,y2]);
    legend('$\rho$','interpreter','latex','fontsize',13);
    xlabel('$x$','interpreter','latex','fontsize',13);
    ylabel('$\rho$','interpreter','latex','fontsize',13);
    title(strcat('$t=',num2str(t(i)),'$'),'interpreter','latex','fontsize',13);
    M_rho(i)=getframe;
end

for i=1:length(t)
    figure(2);
    color=[238,201,0]/255;
    x1=min(x);
    x2=max(x);
    y1=min(u(1,:))-0.1*(max(u(1,:))-min(u(1,:)));
    y2=max(u(1,:))+0.1*(max(u(1,:))-min(u(1,:)));
    plot(x,u(i,:),'.-','color',color);
    axis([x1,x2,y1,y2]);
    legend('$u$','interpreter','latex','fontsize',13);
    xlabel('$x$','interpreter','latex','fontsize',13);
    ylabel('$u$','interpreter','latex','fontsize',13);
    title(strcat('$t=',num2str(t(i)),'$'),'interpreter','latex','fontsize',13);
    M_u(i)=getframe;
end
    
for i=1:length(t)
    figure(3);
    color=[255,48,48]/255;
    x1=min(x);
    x2=max(x);
    y1=min(p(1,:))-0.1*(max(p(1,:))-min(p(1,:)));
    y2=max(p(1,:))+0.1*(max(p(1,:))-min(p(1,:)));
    plot(x,p(i,:),'.-','color',color);
    axis([x1,x2,y1,y2]);
    legend('$p$','interpreter','latex','fontsize',13);
    xlabel('$x$','interpreter','latex','fontsize',13);
    ylabel('$p$','interpreter','latex','fontsize',13);
    title(strcat('$t=',num2str(t(i)),'$'),'interpreter','latex','fontsize',13);
    M_p(i)=getframe;
end

%movie(M);


