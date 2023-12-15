function [aout1,aout2,aout3,aout4,aout5]=function_sl_R_P_bi_RRMSE(data_in,ind_mod)

%purpose of this function is to compute slope, R. value, p-value (at 95% confidence interval),bias, RRMSE.
%--input: 
%--the only input is the data_in file which is a two column array. The
%first column is used for observations and the second column is used for
%modeled values. The model can take data with NaN values.

% ind_mod: is an index, which if -1 means that the second row of the data_in is residual 
% values. If it is set to 1 then it means the second row of the data_in is
% the modeled values.

%----output:
% aout1: slope value
% 
% aout2: Correlation coefficient R
% 
% aout3: p-value at 95% confidence intetval
% 
% aout4: bias value
% 
% aout 5: RRMSE value---


%%%%%%%%%%%%%%%%%%% load the data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; clc; close all;
% % cd('C:\Users\ravindra\Desktop\MGC_data_analysis\PCA_analysis\function_work')  %--laptop---
% cd('\\HAS-Meade.catnet.arizona.edu\UserFolders\ravindradwivedi\Desktop\PCA_analysis') %---computer lab
% data_in=load('test_data_for_sl_R_p_bi_RRMSE_funct.txt');
% %--format of the above file is:
% %--(1): Obs
% %--(2): Mod
% % example:
% % 0.089423769	0.189423769
% % 0.10440446	0.20440446
% % 0.121738968	0.221738968
% ind_mod=-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sl=0;
R_val=0;
p_val=0;
bi_val=0;
rrmse_val=0;

data_in1=data_in;
[m0,n0]=size(data_in1);

len_ch=0;
for i1=1:m0
   if  (~isnan(data_in1(i1,1))) && (~isnan(data_in1(i1,2)))
     len_ch=len_ch+1;       
   end    
end

data_in2=zeros(len_ch,2);

count_1=0;
for i1=1:m0
   if  (~isnan(data_in1(i1,1))) && (~isnan(data_in1(i1,2)))
       count_1=count_1+1;
       data_in2(count_1,:)=data_in1(i1,:);     
   end    
end

data_in3=zeros(len_ch,3);
%---format of the data_in3 file is:
%---(1)---: Obs
%---(2)---: Mod
%---(3)---: Res (Mod-Obs)

[md,nd]=size(data_in3);

data_in3(:,1)=data_in2(:,1);

if ind_mod==-1  %---means the second column of the input data is residulas
    data_in3(:,2)=data_in2(:,1)+data_in2(:,2);
    data_in3(:,3)=data_in2(:,2);
elseif ind_mod==1  %---means the second column of the input data is modeled values
   data_in3(:,2)= data_in2(:,2);
   data_in3(:,3)=data_in2(:,2)-data_in2(:,1);
end


ab_temp=zeros(size(data_in3,1));
ba_temp=zeros(size(data_in3,1));
      x101=zeros(size(data_in3,1),1);
      y101=zeros(size(data_in3,1),1);
      [ab_temp,ba_temp]=sort(data_in3(:,1)); %----sort in the ascending (default) order by obs values
      
      for k5=1:length(ba_temp)
       x101(k5,1)=data_in3(ba_temp(k5,1),1);%---obs value ---- ascending order-
       y101(k5,1)=data_in3(ba_temp(k5,1),3);%---residuals---- in the order of the ascending order of obs values-
      end
       
      %plot(x101,y101,'o');
      [R_Ci_Q,P_Ci_Q]=corrcoef(x101,y101,'alpha',0.05); %---get the R value (with sign), p value at 95% confidence interval--- 
      [p101,s101] = polyfit(x101,y101,1);
      [yfit101,dy101] = polyconf(p101,x101,s101,'alpha',0.05,'predopt','observation');  %,'observation'
      %slope_C_Q=p101(1,1); %---get the slope of the best fit line---
      R_P_vec=[R_Ci_Q(1,2),P_Ci_Q(1,2)];

      sum_1=0; sum_2=0;
      for i2=1:len_ch
          sum_1=sum_1+(data_in3(i2,2)-data_in3(i2,1));
          sum_2=sum_2+((data_in3(i2,2)-data_in3(i2,1))^2);
      end
      mean_obs=mean(data_in3(:,1)); %---get the mean value of the observed series
      
      
sl_val=p101(1,1); %---get the slope of the best fit line---;
R_val=R_Ci_Q(1,2);
p_val=P_Ci_Q(1,2);
bi_val=sum_1/mean_obs;
rrmse_val=(sqrt(sum_2))/(len_ch*mean_obs);

aout1=sl_val;
aout2=R_val;
aout3=p_val;
aout4=bi_val;
aout5=rrmse_val;

end
