function myelinated_HFS_avg_comparrison()
clc
clear all
close all
% modelling pulse propagation in myealinated FHN neuron. Comparing results
% between avg and HFS equations

    P.N = 500; % number of nodes
    P.A = 0.5; % stimulation amplitude (normalized to frequency, i.e. A=a/omega)
    P.D = 0.015;% Diffussion coefficient. When A=0 excitation decay and does not travel. 
                % When A is increased by value larger than some threshold,
                % pulse can travel

    P.omega = 10; % the higer frequency, the better coincidance between avg and HFS solutions
             
    P.eps   = 0.0008; % FHN parameters
    P.gamma = 0.8;
    P.beta  = 0.7;

    
    A = P.A;
    N  = P.N;
    gamma = P.gamma;
    beta = P.beta;

    omega = P.omega;

  % Finding stationary solution of averaged eqs and setting initial conditions
    v00=roots([1/3 0 -(1-A^2/2)+1/gamma beta/gamma]);
    [idx,~,~]=find(imag(v00)==0);
    v0=v00(idx)
    w0=(v0+beta)/gamma
    
    xinit=zeros(2*N,1);
    xinit(1:N)=v0;
    xinit(N+1:end)=w0;
    
    xinit(25:35)=v0+2;
    xinit(N+22:N+32)=w0+1;
    
    %integration time is multiple of 2*pi/omega 
    %necessary for matching initial conditions for averaged and HFS systems
    tint = linspace(0,(2*pi/omega)*floor(3000*omega/(2*pi)),6000);

    options=odeset('RelTol', 1.0e-7, 'AbsTol', 1.0e-9); 
    
    %integrating avg system for transient processes
    [T, X]=ode45(@(t,x)sys_eqns_avg(t,x,P),tint, xinit, options);
    
    xinit=X(end,:);
%     xinit(25:35)=xinit(25:35)+2;
%     xinit(N+22:N+32)=xinit(N+22:N+32)+1;
    tint=0:2:4000;
    
 
    
    %integrating averaged system
    [T, X]=ode45(@(t,x)sys_eqns_avg(t,x,P),tint, xinit, options);
    %integrating HFS system
    [T1, X1]=ode45(@(t,x)sys_eqns_HFS(t,x,P),tint, xinit, options);
  
    
    figure
    subplot(211)
    plot(X(end,1:1:N),'b+','MarkerSize',2,'MarkerFaceColor','b' ), hold on
    plot(X(end,N+1:1:end),'ro','MarkerSize',1.25,'MarkerFaceColor','r')
    title('averaged solutions at final time')
    xlabel('n','FontSize',10)
    ylabel('v_n, w_n','FontSize',10)
    set(gca,'LineWidth',0.4,'FontSize',8)

    
    subplot(212)
    
    % remove HFS component from v variable
    %plot(T1,X1(:,150)-A*sin(omega*T),'b','MarkerSize',2), hold on 
    plot(T1,X1(:,150),'b','MarkerSize',2), hold on
    plot(T,X(:,150),'r','MarkerSize',3)
    title('comparisson between averaged and HFS v variable')
    xlabel('t','Fontsize',10)
    ylabel('v_{150}','FontSize',10)
        
    

end


function dx= sys_eqns_avg(t,x,P)
    
    D = P.D;
    A = P.A;
    N  = P.N;
    eps = P.eps;
    gamma = P.gamma;
    beta = P.beta;

    dx=zeros(2*N,1);
    dx(1:N)=D*DDx( x(1:N) )-1/3*x(1:N).^3+x(1:N)*(1-A^2/2)-x(N+1:end);
    dx(N+1:end)= eps*( x(1:N)+beta-gamma*x(N+1:end));
end


function dx= sys_eqns_HFS(t,x,P)
    
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