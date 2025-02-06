function myelinated_FHN_HFS()
clc
clear all
close all
% modelling pulse propagation in myealinated FHN neuron. As a model we use
% FHN equations with high frequency stimulation.


    P.N = 500; % number of nodes
    P.A = 0.5; % stimulation amplitude (normalized to frequency, i.e. A=a/omega)
    P.D = 0.015;% Diffussion coefficient. When A=0 pulse decay and does not travel. 
                % When A is increased by value larger than some threshold
                % pulse can travel
    P.omega = 10;
             
    P.eps   = 0.0008; % FHN parameters
    P.gamma = 0.8;
    P.beta  = 0.7;

    
    A = P.A;
    N  = P.N;
    gamma = P.gamma;
    beta = P.beta;


    % Finding stationary solution and setting initial conditions
    v00=roots([1/3 0 -(1-A^2/2)+1/gamma beta/gamma]);
    [idx,~,~]=find(imag(v00)==0);
    v0=v00(idx)
    w0=(v0+beta)/gamma
    
    xinit=zeros(2*N,1);
    xinit(1:N)=v0;
    xinit(N+1:end)=w0;
    
    % Exciting elements between 25 and 35 for pulse initiation
    xinit(25:35)=v0+2;
    xinit(N+22:N+32)=w0+1;
    tint=0:2:10000;
    
    %integrating system
    [T, X]=ode45(@(t,x)sys_eqns(t,x,P),tint, xinit);

%     %animation
%     maxY=max(max(X));
%     minY=min(min(X));
%     for ii=1:length(T)
%         plot(X(ii,1:N),'r.')
%         ylabel('v_n')
%         xlabel('n')
%         pav=sprintf('t=%f',T(ii));
%         axis([-0.1 N+10 minY-0.1 maxY+0.1])
%         title(pav);
%         drawnow
%         pause(0.001);
%     end

    figure
    plot(X(1,1:N),'b-'), hold on
    plot(X(1000,1:N),'r-')
    plot(X(end,1:N),'g-')
    ylabel('v_n')
    xlabel('n')
    tm2 = sprintf('T=%.2f',T(1000));
    tm3 = sprintf('T=%.2f',T(end));
    legend({'initial',tm2,tm3});
    
    
    figure
    plot(X(end,1:N),'b+','MarkerSize',2,'MarkerFaceColor','b' ), hold on
    plot(X(end,N+1:1:end),'ro','MarkerSize',1.25,'MarkerFaceColor','r')
  % axis([100 225 -2 2])
    xlabel('n')
    legend({'v_n','w_n'})

   figure
   plot(T,X(:,200))
        

end


function dx= sys_eqns(t,x,P)
    
    D = P.D;
    A = P.A;
    N  = P.N;
    eps = P.eps;
    gamma = P.gamma;
    beta = P.beta;
    omega = P.omega;

    dx=zeros(2*N,1);
    dx(1:N)=D*DDx( x(1:N) )-1/3*x(1:N).^3+x(1:N)-x(N+1:end)+A*omega*cos(omega*t);
    dx(N+1:end)= eps*( x(1:N)+beta-gamma*x(N+1:end));
end


function Dx=DDx(x)
    nn=length(x);
    Dx=zeros(nn,1);
    for ii=2:nn-1
        Dx(ii)=(x(ii+1)+x(ii-1)-2*x(ii));
    end
%     Dx(nn)=x(nn-1)-x(nn);
    Dx(1)=x(2)-2*x(1)+x(nn);
    Dx(nn)=x(nn-1)-2*x(nn)+x(1);
end