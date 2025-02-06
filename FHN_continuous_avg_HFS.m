function FHN_continuous_avg_HFS()  
% coninuous Fitzhugh Nagumo model under affect of high frequency
% stimulation. Comparisson between original and averaged equations
    clear all; close all
    % FN parameters
    P.a   = 0.7;
    P.b   = 0.8;
    P.eps = 0.008; 
    
    P.Amp =  0.; % A=a/w stimulation amplitude
    P.w   =  50; % stimulation frequency

    L=800; % neuron length in x direction (dimensionless units)
    P.N = 1600; % number of spacial discretization steps
    P.dx = L/P.N;
    
    w   = P.w; 
    N   = P.N; 
    

    %preparing initial pulse solution
    v0(1:N) = 0; %pradines salygos V kintamajam
    v0(N+1:2*N) = 0; %pradines salygos W kintamajam
    tint = [0 (2*pi/w)*floor(300*w/(2*pi))];
    [t,vv] = ode45(@(t,x)FHN_avg(t,x,P),tint,v0');
    v0 = vv(end,:);
    v0(1:100) = v0(1)+3;
    [t,vv] = ode45(@(t,x)FHN_avg(t,x,P),tint,v0');  
    v0 = vv(end,:);

    % array of stimulation amplitudes
    AA=[0 0.6  1  1.13 ];

    
    for ii=1:length(AA)
        P.Amp = AA(ii);
        fprintf('A = %.2f\n',P.Amp);
        
        tint = [0 (2*pi/w)*floor(250*w/(2*pi))];
        [t,v] = ode45(@(t,x)FHN_HFS(t,x,P),tint,v0'); %integravimas
        [t,vv] = ode45(@(t,x)FHN_avg(t,x,P),tint,v0');
        
        tint = linspace(0, (2*pi/w)*floor(1200*w/(2*pi)), 1200 );
        [t,v] = ode45(@(t,x)FHN_HFS(t,x,P),tint,v(end,:)); %integravimas
        [t,vv] = ode45(@(t,x)FHN_avg(t,x,P),tint,vv(end,:));
    
%     %   plot pulse propagation
%         figure
%         imagesc(1:N,t,vv(:,1:N))
%         xlabel('n')
%         ylabel('t')
%         pav = sprintf("A = %f",P.Amp);
%         title(pav)
        
        figure (22)
        subplt=410+ii;
        subplot(subplt);
        
        plot(t,v(:,1200)) % under HFS stimulation
        hold on
        plot(t,vv(:,1200),'r-') % averaged solution
        axis([0 1200 -2.5 2.5])
        ylabel('v');
        drawnow;
    end
    xlabel('t')
    



end

% Fitzhugh Nagum under HFS
function vdot = FHN_HFS(t,v,P)
    a   = P.a  ; 
    b   = P.b  ;
    eps = P.eps;
    Amp = P.Amp;
    N   = P.N;
    dx  = P.dx;
    w   = P.w; 
    vdot =zeros(2*N,1);
    vdot(1:N) = v(1:N).*(1-v(1:N).*v(1:N)/3)- v(N+1:end) + DDx(v(1:N),dx,P)'+Amp*w*cos(w*t);
    vdot(N+1:2*N) = eps*(v(1:N)+a-b*v(N+1:end));
end


% Fitzhugh Nagum averaged equations
function vdot = FHN_avg(~,v,P)
    a   = P.a  ; 
    b   = P.b  ;
    eps = P.eps;
    Amp = P.Amp;
    N   = P.N;
    dx  = P.dx;
    vdot =zeros(2*N,1);
    vdot(1:N) = v(1:N).*(1-Amp*Amp/2-v(1:N).*v(1:N)/3)- v(N+1:end) + DDx(v(1:N),dx,P)';
    vdot(N+1:2*N) = eps*(v(1:N)+a-b*v(N+1:end));
end

%Calculation of the spacial derivative
function V=DDx(v,dx,P)
    N = P.N;
    V=zeros(1,N);
    for j=2:N-1
        V(j)=(v(j-1)-2*v(j)+v(j+1))/dx^2;
    end
    %periodic boundary conditions
    %V(1)=(v(N)-2*v(1)+v(2))/dx^2; V(N)=(v(N-1)-2*v(N)+v(1))/dx^2;
    % Free boundaries
    V(1)=-2*v(1)+v(2); V(N)=v(N-1)-2*V(N);
end

