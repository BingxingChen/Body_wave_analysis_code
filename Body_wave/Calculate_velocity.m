% reference：John H. Long Jr 《Flapping flexible fishPeriodic and secular body reconfigurations in swimming lamprey, Petromyzon marinus》 
% by Bingxing Chen 20181121
%   x/body_x：the axial direction  y/body_y：the lateral direction (not rotated)
clc
clear
close all

disp('________________________________________________________________')
disp('Calculate the velocity and amplitude of fish body')
% load bd_data.mat
% DF=1.72;%estimated frequency
% timestep = 1/30;%the time of frames
% num_per_cycle=round(1/DF/timestep)+1;
% num_orign=73;%%the initial frame
% num_end=135;%%the last frame
% num_jpgs=num_end-num_orign+1;%%the numbers of frame for the midline_gen
% NumFr=num_jpgs-num_per_cycle;% the number of cycles
% n_points=30;%fish body points

load bd_data_1.2HZ.mat
DF=1.2;timestep = 1/30;num_per_cycle=round(1/DF/timestep)+1;
NumFr=36-num_per_cycle;


% load bd_datadeg632.9.mat
% DF=2.9;timestep = 1/30;num_per_cycle=round(1/DF/timestep)+1;
% num_orign=47;%%要分析图片的开始值
% num_end=117;%%要分析图片的末端值
% num_jpgs=num_end-num_orign+1;%%midline_gen中图片的总数
% NumFr=num_jpgs-num_per_cycle;

% load bd_data_deg452.75.mat
% DF=2.75;timestep = 1/30;num_per_cycle=round(1/DF/timestep)+1;
% num_orign=65;%%要分析图片的开始值
% num_end=140;%%要分析图片的末端值
% num_jpgs=num_end-num_orign+1;%%midline_gen中图片的总数
% NumFr=num_jpgs-num_per_cycle;


% 51 pixel represents 100mm
Fish_length_calibration=51.5*3.9;

U.x=[];
% strt_frm=1
for strt_frm=1:NumFr% different cycles
Data.one.x =body_x(strt_frm:strt_frm+1*num_per_cycle,:)';%% row: body points line:time series
Data.one.y = body_y(strt_frm:strt_frm+1*num_per_cycle,:)';
fishL =[];sz=size(Data.one.x);%%sz 求取两个值 sz(1)行：表示鱼体的位置 体长的数据 有多少个鱼体的点 sz(2)列：时间
for il=1:sz(2)%表示
fishL(il) =sum(sqrt(diff(Data.one.x(:,il)).^2+diff(Data.one.y(:,il)).^2));
end
fishLength=mean(fishL); %the body length  
FittData.one.H=[];FittData.H0=[];


%% obtain the rotation angle through iterations
Data.one.X =Data.one.x;Data.one.Y =Data.one.y;
iteration.flag=0;%
iteration.eps=1e-2;%iteration precision
while(1)
FittData.one.H=[];
for i=1:num_per_cycle+1;%one cycle
% FittData.one.H0=[1 i*timestep (i*timestep)^2/2  sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep) sin(3*DF*2*pi*i*timestep) cos(3*DF*2*pi*i*timestep) sin(4*DF*2*pi*i*timestep) cos(4*DF*2*pi*i*timestep)  ];
FittData.one.H0=[1 i*timestep sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) ];% the preliminary Fourier transforms.. and it just needs the velocity term and fundamental harmonic
FittData.one.H=[FittData.one.H;FittData.one.H0;];% 
end
Data.one.xInverse=Data.one.X';Data.one.yInverse=Data.one.Y';
FittData.one.S_u=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.xInverse);%least square fitting method
FittData.one.S_v=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.yInverse);%least square fitting method
FittData.one.u_model=[];FittData.one.v_model=[];
FittData.one.u_model=FittData.one.H*FittData.one.S_u;FittData.one.v_model=FittData.one.H*FittData.one.S_v;
FittData.one.u_model_Periodic=FittData.one.H(:,3:end)*FittData.one.S_u(3:end,:)+FittData.one.H(:,1)*FittData.one.S_u(1,:);
FittData.one.v_model_Periodic=FittData.one.H(:,3:end)*FittData.one.S_v(3:end,:)+FittData.one.H(:,1)*FittData.one.S_v(1,:);
FittData.P=polyfit(FittData.one.u_model_Periodic,FittData.one.v_model_Periodic,1);
FittData.f =polyval(FittData.P,FittData.one.u_model_Periodic);
% % plot to verify
figure(1)
plot(Data.one.xInverse,Data.one.yInverse,'r',FittData.one.u_model,FittData.one.v_model,'k')%%鱼体波的图 the all fish body;
title('body wave model and actualbody')
axis equal
figure(2)%% the 'eight' of figure
plot(FittData.one.u_model_Periodic(1:end,30)',FittData.one.v_model_Periodic(1:end,30)','r') %转动之前的鱼体波数据
title('tranjectories of body points')
figure(3)
for iiii = 1:size(Data.one.xInverse,1)%%说明这有多少列 也就是  不同时间进行求值
plot(FittData.one.u_model_Periodic(iiii,:)',FittData.one.v_model_Periodic(iiii,:)','r') %转动之前的鱼体波数据
title('periodic components of midlines.... ')
hold on
end
hold on
plot(FittData.one.u_model_Periodic,FittData.f,'k')% the fitting curve
hold off
axis equal
pause(0.5)

if abs(atan(FittData.P(1)))<iteration.eps%%
       iteration.flag=1;
       break;%break the while
end
FittData.theta(strt_frm)=atan(FittData.P(1));
Data.one.X = cos(FittData.theta(strt_frm))*Data.one.x+sin(FittData.theta(strt_frm))*Data.one.y;%%%the axial direction  Rotate raw data
Data.one.Y= -sin(FittData.theta(strt_frm))*Data.one.x+cos(FittData.theta(strt_frm))*Data.one.y;%the axial direction    Rotate raw data进行旋转
end
% after rotating the raw data
%% Fourier fitting in order to calculate the velocity and tail amplitude
FittData.one.H=[];
for i=1:num_per_cycle+1;%one cycle
% %new    
FittData.one.H0=[1 i*timestep sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep) sin(3*DF*2*pi*i*timestep) cos(3*DF*2*pi*i*timestep) sin(4*DF*2*pi*i*timestep) cos(4*DF*2*pi*i*timestep)  ];
% FittData.one.H0=[1 i*timestep sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep)  ];
FittData.one.H=[FittData.one.H;FittData.one.H0;];% 
end
Data.one.xInverse=Data.one.X';Data.one.yInverse=Data.one.Y';
FittData.one.S_u=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.xInverse);%least square fitting method
FittData.one.S_v=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.yInverse);%least square fitting method
FittData.one.u_model=[];FittData.one.v_model=[];
FittData.one.u_model=FittData.one.H*FittData.one.S_u;FittData.one.v_model=FittData.one.H*FittData.one.S_v;
FittData.one.u_model=FittData.one.u_model';
FittData.one.v_model=FittData.one.v_model';

U.calculateFit(:,strt_frm)=FittData.one.S_u(2,:);%%the velocity term of the fitting curve in the axial direction (u).   it is the velocity of fish body
% verify 
U.calculate(:,strt_frm)=30*abs(FittData.one.u_model(:,end)-FittData.one.u_model(:,1))/size(FittData.one.u_model,2);%%鱼头用X数据  end 表示鱼头 体长的数据 记在行  时间记再列
U.head(strt_frm)=30*abs(FittData.one.u_model(end,end)-FittData.one.u_model(end,1))/size(FittData.one.u_model,2);%%the velocity of head
%%but the tail amplitude should be recalculated. because the dundamental amplitude is not the entire amplitude..
A.cycle(strt_frm) =(max(FittData.one.v_model(1,:))-min(FittData.one.v_model(1,:)));%%鱼尾用Y数据  1表示鱼尾

% for i_actual= 1:n_points%%body points
% % U.x(i_actual)=30*(Data.one.x(i_actual,end)-Data.one.x(i_actual,1))/fishLength/ size(Data.one.x,2);%%velocity of x direction
% % U.y(i_actual)=30*(FittData.one.v_model(i_actual,end)-FittData.one.v_model(i_actual,1))/fishLength/ size(FittData.one.v_model,2);%%velocity of y direction
% % the velocity also can use the actual velocity by the calibration for the actual fish length 
% % % % U.one.x(i_actual)=30*(FittData.one.u_model(i_actual,end)-FittData.one.u_model(i_actual,1))/ size(FittData.one.u_model,2);%%x velocities of all body points in one cycle
% % % % U.one.y(i_actual)=30*(FittData.one.v_model(i_actual,end)-FittData.one.v_model(i_actual,1))/ size(FittData.one.v_model,2);%%y velocity of all body points in one cycle
% %record velocities of all body points and all cycles
% U.all.x(i_actual,strt_frm)=U.one.x(i_actual);%
% U.all.y(i_actual,strt_frm)=U.one.y(i_actual);%
% end
%the x/y components of the average velocity for individual body points
% % U.one.mean=[mean(U.one.x(:)) mean(U.one.y(:))];% the value of the composite average velocity in one cycle.
% % U.one.xy=sqrt(U.one.x.^2+U.one.y.^2);
% % RE.one.x=(U.one.x-mean(U.one.x(:)))/norm(U.one.mean);RE.one.y=(U.one.y-mean(U.one.y(:)))/norm(U.one.mean);
% % RE.one.Square=sqrt(RE.one.x.^2+RE.one.y.^2);RE.all.Square(strt_frm,:)=RE.one.Square;
% % % The unsteadliness index in one cycle
% % UI.one=RE.one.Square.^2*norm(U.one.mean)./U.one.xy;UI.all(strt_frm,:)=UI.one;%The unsteadliness index of all cycles
% % %calculate theta by the composite average velocity in one cycle.
% % theta.one.U=atan2(mean(U.one.y(:)),mean(U.one.x(:)));theta.all.U(strt_frm)=theta.one.U;
% % 
end
% UI.mean=mean(UI.all(:));%The means of unsteadliness index
U.headmean=mean(U.head);
U.head_tailmean=mean(mean(U.calculate));
U.mean=mean(mean(U.calculateFit));
Stride_length=U.mean/DF;
% U.meanstandDeviation=sqrt(sum(((U.calculateFit-U.mean)/U.mean).^2)/size(U.calculateFit,2));
U.meanstandDeviation=sqrt(mean(sum(((U.calculateFit-U.mean)/U.mean).^2)/size(U.calculateFit,2)));
A.mean=mean(A.cycle);
A.standDeviation=sqrt(sum(((A.cycle-A.mean)/A.mean).^2)/size(A.cycle,2));
% The unsteadliness index in one cycle
St.mean=DF*A.mean/U.mean;
St.meanhead_tail=DF*A.mean/U.head_tailmean;

disp( '************The calculations*************');
fprintf(1,'U.mean is (the velocity term of the fitting curve in the axial direction (u)): %1.4f\n',U.mean/Fish_length_calibration)%%we use this data
fprintf(1,'The stand Deviation of U.mean of 30 body points in this trial(many cycles)is: %3.4f\n',U.meanstandDeviation)
fprintf(1,'A.mean  is: %3.4f\n',A.mean/Fish_length_calibration)
fprintf(1,'The stand Deviation of A.mean of 30 body points in this trial(many cycles) is: %3.4f\n',A.standDeviation)
fprintf(1,'St.mean is: %3.4f\n',St.mean)
fprintf(1,'Stride length is: %3.4f\n',Stride_length/Fish_length_calibration)
disp( '**********comparision to verify***************');
fprintf(1,'U.headmean (just the velocity of head point) (pixel/s) is: %3.4f\n',U.headmean/Fish_length_calibration)
fprintf(1,'St.meanhead_tail is: %3.4f\n',St.meanhead_tail)
fprintf(1,'U.head_tailmean (pixel/s) is: %3.4f\n',U.head_tailmean/Fish_length_calibration)
data_body=[DF U.mean/Fish_length_calibration A.mean/Fish_length_calibration St.mean  Stride_length/Fish_length_calibration]


