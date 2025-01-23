close all;
clear all;
clc;

L=100; % Length of the domain
N=[70 100 110 120 130 140 150 200 230 250 270 300 330 350 370 400 450 500 700 ]; % Number of points in 1D space
nt=[1 2 3 4 5 6 7 8 10 15 20 250 30 35 40 50 60 70 80 90 100]; % Number of points in time
c=340; %Speed of wave
T=0.001; %Time travelled by the wave

%% Exact solution
dy=L/(N(end));
y=0:dy:L;
uexact=sinc(y-L/2 - c*T)/4;
[a,b]=max(uexact);
uexact(1:round(b/2))=0;
uexact(end-round(b/2):end)=0;
uexact=uexact';

erl=zeros(length(N),length(nt));
alpha=zeros(length(N),length(nt));

%% Numerical simulation using FEM
for(p=1:numel(N))

    dx=L/(N(p)); %Space discretization

    x=0:dx:L;
    u0=zeros(numel(x),1); % Initial condition
    u0 = sinc(x-L/2)/4;
    u0(1:L/(4*dx))=0;
    u0(end-(L/(4*dx)):end)=0;
    [uamp,uind]=max(u0);
    u0l=u0';

    K=zeros(numel(x)); % Stiffness matrix
    K(1,1)=-1/2;
    K(1,2)=1/2;
    M=zeros(numel(x)); %Mass matrix
    M(1,1)=2;
    M(1,2)=1;
    for(i=2:N(p))
        K(i,i-1)=-1/2;
        K(i,i)=0;
        K(i,i+1)=1/2;

        M(i,i-1)=1;
        M(i,i)=4;
        M(i,i+1)=1;
    end

    M(numel(x),numel(x)-1)=1;
    M(numel(x),numel(x))=2;
    K(numel(x),numel(x)-1)=-1/2;
    K(numel(x),numel(x))=1/2;
    M=M/6;

    %identity matrix
    idma=zeros(numel(x),1);
    idma=eye(numel(x));



    %% LAX Wendroff Method
    k=M\K; %Global Matrix

    for(j=1:numel(nt)) %space loop
        u0l=u0';

        dt=T/nt(j); % Time discretization
        t=0:dt:T;

        c1=c*dt/dx;     %CFL aka alpha singular value

        k1=c1*k;
        alpha(p,j)=c*(dt/dx); %CFL plugged into matrix
        for(i=1:nt(j)+1)
            u1l=(idma - k1 + 0.5*(c1)^2 * k * k)*u0l; %Laxwendroff eqn
            u0l=u1l;
            %             if(j==10 && p==16) %For wavefield plots
            %                 u3(i,:)=u1l;
            %                 figure(p+5)
            %                 plot(x,u0,x,u1l,"LineWidth",1);
            %                 legend("U0","Un")
            %                 xlabel("Distance in m")
            %                 ylabel("Amplitude")
            %                 xlim([25 75])
            %             end

            hold on;
        end
        utemp=interp1(x,u1l,y); %Interpolation function to compare exact solution with simulation
        erl(p,j)=sum((utemp-uexact')/uexact')^2; %Relative error squared
        %         figure(p) % Plot for space iteration
        %         plot(x,u1l,"linewidth",1)
        %         fontsize(gca,13,'points')
        %         hold on

    end
    %     plot(y,uexact,"LineWidth",1,"LineStyle","--") %Plotting exact solution with space iteration
    %     xlabel("Distance in m")
    %     ylabel("Amplitude")
    %     hold all
end

%% error maps and alpha

figure(p+1);
imagesc(nt,N,erl)
xlabel("Number of points in time")
ylabel("Number of points in space")
title("Error Map")
colorbar

figure(p+2);
imagesc(nt,N,alpha)
xlabel("Number of points in time")
ylabel("Number of points in space")
title("Value of alpha for various discretizations")
colorbar

%% Forward and backward propagation
figure(p+3)
plot(x,u0,x,u1l,"linewidth",1.5)
legend("u0","Numerical simulation of U_(n+1)")
xlabel("Distance in m")
ylabel("Amplitude")

%% Wavefield plots
% figure(p+4)
% fontsize(gca,13,'points')
%
% surf(x,t,u3);
% ylabel("Time");
% xlabel("Distance")
% zlabel("Amplitude")

