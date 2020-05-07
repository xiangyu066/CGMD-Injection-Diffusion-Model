echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: CGMD_Flagellation.m
% - Author: XYZ
% - Created date: December 28, 2016
% - Modified date: May 4, 2020
% - Notes:
%       1.) Using coarse-grained MD simulate the flagellation following the injection-diffusion model.
%       2.) Two particles(2nm) connect with a spring.
%       3.) Fixed loading rate.
% - Version: 2.0.1
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
clear all, warning('off')
disp('Running...')

%% Define units
global nm um sec msec usec pN KT
nm = 1;
um = 1E+3 *(nm);
sec = 1;
msec = 1E-3 *(sec);
usec = 1E-6 *(sec);
pN = 1;
KT = 4.114 *(pN*nm);

%% Define parameters of imaging system
T_loading = 0.035 *(sec);           % loading time
F_loading = 5 *(pN);                % loading force only work on the first particle
DIFF_Const = 1.25E+4 *(nm^2/sec);
dt = 0.1 *(usec);                   % simulation time step
L0 = 56 *(nm);                      % the length of the unfolded flagellin
stiffness = 20 *(pN/nm);            % spring constant
sigma = 2 *(nm);                    % the size of the particle
epsilon = 1 *(KT);                  
dL = 0.47 *(nm);                    % increment after folding flagellin
stopFL = 1500 *(nm);

monitorTime = 0.5 *(sec);
plotSubFlaLen = 300 *(nm);

%% Preallocating 1vavriables and functions
drag_coef = KT/DIFF_Const;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FL = 0 *(nm);
x = [-L0; 0];
%%%%%%%%%%%%%%%% Check %%%%%%%%%%%%%%%%
% stopFL = 200 *(nm);
% monitorTime = dt;
% FL = 168 *(nm);
% x = [-L0,sigma, L0+2*sigma;...
%     0,L0+sigma,2*L0+2*sigma]+24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = ones(size(x));
nFlagellins = size(x,2);            % count for total number of flagellin
Status_loading = 0;                 % determine whether pump can load new flagellin into the channel
Fx = zeros(size(x));                % total force
F_k = zeros(size(x));               % stifness force
F_LJ = zeros(size(x));              % L-J force

tCounts = 0;

%%
outputname = ['CGMD-Fla_T',num2str(T_loading),'_F',num2str(F_loading),'_D',num2str(DIFF_Const)];

vSession = VideoWriter([outputname,'.avi']);
vSession.FrameRate = 100;
open(vSession)

fid = fopen([outputname,'.txt'],'w');
fprintf(fid,'%s\t\t%s\n','Units','nm/sec/pN');
fprintf(fid,'%s\t\t%.3f\n','KT',KT);
fprintf(fid,'%s\t%.3f\n','T_loading',T_loading);
fprintf(fid,'%s\t%.1f\n','F_loading',F_loading);
fprintf(fid,'%s\t\t%.2e\n','dt',dt);
fprintf(fid,'%s\t\t%.2f\n','L0',L0);
fprintf(fid,'%s\t%.1f\n','stiffness',stiffness);
fprintf(fid,'%s\t\t%.2f\n','sigma',sigma);
fprintf(fid,'%s\t\t%.3f\n','epsilon',epsilon);
fprintf(fid,'%s\t\t%.2f\n\n','dL',dL);
fprintf(fid,'%s\t%s\n','Time','FL');

figure(1), set(gcf,'Windowstate','maximized'), 
if (rem(stopFL,plotSubFlaLen)==0)
    nSubs = stopFL/plotSubFlaLen;
else
    nSubs = 1+floor(stopFL/plotSubFlaLen);
end

nplot = 1;
for nSub = 1:2:2*nSubs
    subplot(nSubs,2,nSub),
    yticks([])
    xlim([-L0+plotSubFlaLen*(nplot-1),plotSubFlaLen*nplot]), ylim([0,sigma]), grid on, daspect([1 0.05 1])
    set(gca,'fontsize',16)
    
    if (nplot==1)
        title(['T = 0 [sec] / ','FL = 0 [nm]'], 'fontsize',24,'fontweight','bold')
    elseif (nplot==nSubs)
        xlabel('FL [nm]','fontsize',16,'fontweight','bold')
    end
    
    nplot = nplot+1;
end
subplot(nSubs,2,(2:2:2*nSubs)), grid on
set(gca,'fontsize',16), xlabel('Time [sec]','fontsize',16,'fontweight','bold'), ylabel('FL [nm]','fontsize',16,'fontweight','bold')
ylim([0,stopFL])

%%  Coarse-grained MD simulation
tic
while (FL<stopFL)
    % loading condition
    if (rem(tCounts,round(T_loading/dt))==0)
        if (isempty(x))
            Status_loading = 1;
        else
            if x(1)>sigma
                Status_loading = 1;
            end
        end
    end
    
    if (Status_loading)
        x = [-L0, x(1,:);...
            0, x(2,:)];
        y = [1, y(1,:);...
            1, y(2,:)];
        F_x = zeros(size(x));
        F_k = zeros(size(x));
        F_LJ = zeros(size(x));
        Status_loading = 0;
        nFlagellins = nFlagellins+1;
    end
    
    % the injection-diffusion model
    if (~isempty(x))
        % calculate thermal motion
        dx_thermal = sqrt(2*DIFF_Const*dt)*randn(size(x));
        x = x+dx_thermal;
        
        % calculate stiffness force
        dx = abs(x(2,:)-x(1,:))-L0;
        F_k = [stiffness*dx; -stiffness*dx];
        
        % calculate LJ force
        rij = abs(circshift(x(:),1)-x(:));
        F_LJ_L2R = -4*epsilon*((6*sigma^6)./rij.^7 - (12*sigma^12)./rij.^13);
        rij = circshift(rij,-1);
        F_LJ_R2L = +4*epsilon*((6*sigma^6)./rij.^7 - (12*sigma^12)./rij.^13);
        F_LJ = reshape(F_LJ_L2R+F_LJ_R2L,size(x));

        % calculate total force and determine whether the first particle is suffer from a loading force
        Fx = F_k+F_LJ;
        if (x(1)<sigma/2)
            Fx(1) = Fx(1)+F_loading;
        end
        
        % update position
        Vx = Fx/drag_coef;
        x = x+Vx*dt;
        
        % folded event
        if (x(end-1)>FL)
            x(:,end) = [];
            y(:,end) = [];
            Fx(:,end) = [];
            F_k(:,end) = [];
            F_LJ(:,end) = [];
            FL = FL+dL;
        end
    end
    
    % monitoring
    if (rem(tCounts,round(monitorTime/dt))==0)
%         disp(['T = ',num2str(tCounts*dt), ' [sec] / ','FL = ',num2str(FL),' [nm]'])
        
        % save simulation data
        fprintf(fid,'%.3f\t%.2f\n',tCounts*dt,FL);

        % real-time plot simulation data
        nplot = 1;
        for nSub = 1:2:2*nSubs
            subplot(nSubs,2,nSub); cla(gca), plot(x,y,'-','linewidth',2)
            hold on, line([FL,FL],ylim,'color','r','linewidth',2)
            yticks([])
            if nplot ==1
                xlim([-L0+plotSubFlaLen*(nplot-1),plotSubFlaLen*nplot]), ylim([0,sigma]), grid on, daspect([1 0.05 1])
                title(['T = ',num2str(tCounts*dt), ' [sec] / ','FL = ',num2str(FL),' [nm]'], 'fontweight','bold')
            else
                xlim([plotSubFlaLen*(nplot-1),plotSubFlaLen*nplot]), ylim([0,sigma]), grid on, daspect([1 0.05 1])
            end
            set(gca,'fontsize',16)
            
            nplot=nplot+1;
        end
        subplot(nSubs,2,(2:2:2*nSubs)), hold on, plot(tCounts*dt,FL,'bo'), grid on
        xlim([0,1+floor(tCounts*dt)])
        drawnow
        
        % record animation
        frame = getframe(gcf);
        writeVideo(vSession, frame);
    end
    
    tCounts = tCounts+1;
end
fclose(fid);
close(vSession), toc

%%
disp('Done.')

