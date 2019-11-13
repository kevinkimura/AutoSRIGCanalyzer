%% This script will filter and average the current
clear;clc;close all

EC = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\Cufoil Data\Constant_PulsedLONG\190818-19_CuFoil165 Const -1650mV\190818-19 CuFoil165 Long test EC data.xlsx']);
len = length(EC.data);
figure
plot(EC.data(:,1),EC.data(:,2))
%%
filtered = zeros(round(len/20),2); %% rounds up, divides length by number of GC samples taken during part, makes a matrix of zeros with two columns
ct = 1; % initializing the counter

for i = 1:len
    % finding faradaic current
    if EC.data(i,2) < -0.5 && EC.data(i,2) > -1.5 %% taking current data that is negative but not huge negative charging spikes. this value will change based on case
        filtered(ct,:) = EC.data(i,:); % adding cathodic Faradaic current data to new matrix
        ct = ct +1;
    else
    end
end
figure
hold on
plot(filtered(:,1),filtered(:,2))
% load('CuFoil168filtered_ECdata.mat')
% 
% [val,idx] = min(abs(filtered(:,1)-ECdata168(ptnum-1).data(end,1))); %finding the end point from previous set
% figure
% hold on
% plot(filtered(idx:end,1),filtered(idx:end,2))
% 
% ECdata168(ptnum).data = filtered(idx:end,:); %changes with the part number
% % ECdata168(ptnum).data = filtered; %changes with the part number
% save('CuFoil168filtered_ECdata.mat','ECdata168')


%% Finding Average Faradaic Current from Previously Filtered Data

% load('CuFoil168filtered_ECdata.mat')
% filteredpt8 = filteredpt8(any(filteredpt8,2),:); %removes zeros
% ECdata = [];
% for a = 1:3
%     ECdata = [ECdata; ECdata168(a).data];
% end
% 
% for a = 4:7
%     ECdata168(a).data(:,1) = ECdata168(a).data(:,1) - ECdata168(a).data(1,1); %start all at zero time
%     ECdata = [ECdata(:,1) ECdata(:,2); ECdata168(a).data(:,1)+ECdata(end,1) ECdata168(a).data(:,2)];
% end
ECdata = EC.data;
% save('CuFoil168TotalECdata.mat','data')
figure
plot(ECdata(:,1), ECdata(:,2))
% AllECdata(:,1) = ECtime;
AllECdata = ECdata;
Area165 = 0.088;
AllECdata(:,1) = AllECdata(:,1) - AllECdata(1,1);  % scaling the time to start at 0
AllECdata(:,1) = AllECdata(:,1)/60; %converting to minutes
AllECdata(:,2) = AllECdata(:,2)/(2*Area165); %half the current (equal pulse), density (mA/cm^2)
GC1strun_starttime = 10; %minutes
GCrunrest_time = 26; %min

GCsections = floor((AllECdata(end,1)-GC1strun_starttime)/GCrunrest_time)+1;
LeftoverTime = round(AllECdata(end,1)) - GC1strun_starttime - (GCsections-1)*GCrunrest_time;

CAavgI = zeros(GCsections,3);   % 1st col time, 2nd col indexes, 3rd col average current for each GC section

%finding the time points for each GC sample
for i = 1:GCsections
    n=GC1strun_starttime + GCrunrest_time*(i-1); %value to find
    [val,idx] = min(abs(AllECdata(:,1)-n)); %finding the closest time point
    CAavgI(i,1) = AllECdata(idx,1); %transfer to the average current matrix
    CAavgI(i,2) = idx;
    if i == 1 % first GC run after ECstart
        CAavgI(i,3) = mean(AllECdata(1:CAavgI(1,2),2)); %calculating average current
    else
        CAavgI(i,3) = mean(AllECdata(CAavgI(i-1,2):CAavgI(i,2),2));
    end
end


figure
hold on
plot(AllECdata(:,1),AllECdata(:,2))
for x = 1:GCsections
    plot(CAavgI(x,1),CAavgI(x,3),'o')
end

CAavgI