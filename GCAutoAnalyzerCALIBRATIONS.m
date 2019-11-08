%% CO2R Autoanalyzer Calibration
% This code will take the GC calibration .ASC files and generate the
% calibration curve values.
%
% Developed by Kevin Kimura (2019) while part of Prof. Tobias Hanrath's Lab
% at Cornell University

%% Setting up the Analysis
% 
% Make sure the file path is updated
clear; clc; close all;

%GC Aug 2019 Calibration 'Hugh'
TCD1_100A = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_100 H2 KK TCD02.ASC']); %
TCD1_100B = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_100 H2 KK TCD03.ASC']); %
TCD1_250A = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_250 H2 KK TCD01.ASC']); %
TCD1_250B = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_250 H2 KK TCD02.ASC']); %
TCD1_250C = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_250 H2 KK TCD03.ASC']); %
TCD1_500A = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_500 H2 KK TCD01.ASC']); %
TCD1_500B = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_500 H2 KK TCD02.ASC']); %
TCD1_500C = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_500 H2 KK TCD04.ASC']); %

FID1_100A = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_100 cal KK FID03.ASC']); %
FID1_100B = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190817 1_100 cal KK FID04.ASC']); %
FID1_250A = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190816 1_250 cal KK FID02.ASC']); %
FID1_500A = importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190816 1_500 cal KK FID01.ASC']); %
FID1_500B= importdata(['C:\Users\Kevin\Documents\Graduate_Work_MAIN\Electrochemistry\Data\GCData_analyzer\Aug2019_GCcalibrations\190816 1_500 cal KK FID04.ASC']); %

time = 0.2/60:0.2/60:16; %last value is the total minutes, 5Hz intervals
H2peak1 = 0.8067; %first H2 peak location
CH4peak1 = 2.513; %first CH4 peak location
COpeak1 = 3.777; %first CO peak location
CH4peak2 = 8.283; %second CH4 peak location
C2H4peak1 = 10.8; %C2H4 peak location

SummaryH2(1,:) = [0 0];
SummaryH2(2,:) = [(H2peak1AreaFinder(H2peak1, time, TCD1_500A.data)+H2peak1AreaFinder(H2peak1, time, TCD1_500B.data)+H2peak1AreaFinder(H2peak1, time, TCD1_500C.data))/3 100];
SummaryH2(3,:) = [(H2peak1AreaFinder(H2peak1, time, TCD1_250A.data)+H2peak1AreaFinder(H2peak1, time, TCD1_250B.data)+H2peak1AreaFinder(H2peak1, time, TCD1_250C.data))/3 200];
SummaryH2(4,:) = [(H2peak1AreaFinder(H2peak1, time, TCD1_100A.data)+H2peak1AreaFinder(H2peak1, time, TCD1_100B.data))/2 500];

    Area = 0:1000;
    PH2 = polyfit(SummaryH2(:,1),SummaryH2(:,2),1);
    H2fit = PH2(1)*Area;
    figure
    hold on
    plot(SummaryH2(:,1),SummaryH2(:,2),'o')
    plot(Area,H2fit)
    title('H_2')
    xlabel('Area')
    ylabel('ppm')

SummaryCH4_1(1,:) = [0 0];
SummaryCH4_1(2,:) = [(CH4peak1AreaFinder(CH4peak1, time, FID1_500A.data)+CH4peak1AreaFinder(CH4peak1, time, FID1_500B.data))/2 20];
SummaryCH4_1(3,:) = [CH4peak1AreaFinder(CH4peak1, time, FID1_250A.data) 40];
SummaryCH4_1(4,:) = [(CH4peak1AreaFinder(CH4peak1, time, FID1_100A.data)+CH4peak1AreaFinder(CH4peak1, time, FID1_100B.data))/2 100];

    PCH4A = polyfit(SummaryCH4_1(:,1),SummaryCH4_1(:,2),1);
    CH4fit = PCH4A(1)*Area;
    figure
    hold on
    plot(SummaryCH4_1(:,1),SummaryCH4_1(:,2),'o')
    plot(Area,CH4fit)
    title('CH_4 1st peak')
    xlabel('Area')
    ylabel('ppm')
    
SummaryCO(1,:) = [0 0];
SummaryCO(2,:) = [(COpeak1AreaFinder(4.407, time, FID1_500A.data)+COpeak1AreaFinder(4.407, time, FID1_500B.data))/2 20];
SummaryCO(3,:) = [COpeak1AreaFinder(4.407, time, FID1_250A.data) 40];
SummaryCO(4,:) = [(COpeak1AreaFinder(4.407, time, FID1_100A.data)+COpeak1AreaFinder(4.407, time, FID1_100B.data))/2 100];

    PCO = polyfit(SummaryCO(:,1),SummaryCO(:,2),1);
    COfit = PCO(1)*Area;
    figure
    hold on
    plot(SummaryCO(:,1),SummaryCO(:,2),'o')
    plot(Area,COfit)
    title('CO')
    xlabel('Area')
    ylabel('ppm')
    
SummaryCH4_2(1,:) = [0 0];
SummaryCH4_2(2,:) = [(CH4peak2AreaFinder(8.293, time, FID1_500A.data, 55, 25, 20, 55)+CH4peak2AreaFinder(8.293, time, FID1_500B.data, 55, 25, 20, 55))/2 20];
SummaryCH4_2(3,:) = [CH4peak2AreaFinder(8.293, time, FID1_250A.data, 55, 25, 20, 55) 40];
SummaryCH4_2(4,:) = [(CH4peak2AreaFinder(8.293, time, FID1_100A.data, 55, 25, 20, 55)+CH4peak2AreaFinder(8.293, time, FID1_100B.data, 55, 25, 20, 55))/2 100];

    PCH4B = polyfit(SummaryCH4_2(:,1),SummaryCH4_2(:,2),1);
    CH4fit = PCH4B(1)*Area;
    figure
    hold on
    plot(SummaryCH4_2(:,1),SummaryCH4_2(:,2),'o')
    plot(Area,CH4fit)
    title('CH_4 2nd peak')
    xlabel('Area')
    ylabel('ppm')
    
SummaryC2H4_1(1,:) = [0 0];
SummaryC2H4_1(2,:) = [(C2H4peak1AreaFinder(C2H4peak1, time, FID1_500A.data, 200, 80, 60, 120)+C2H4peak1AreaFinder(C2H4peak1, time, FID1_500B.data, 200, 80, 60, 120))/2 20];
SummaryC2H4_1(3,:) = [C2H4peak1AreaFinder(C2H4peak1, time, FID1_250A.data, 200, 80, 60, 120) 40];
SummaryC2H4_1(4,:) = [(C2H4peak1AreaFinder(C2H4peak1, time, FID1_100A.data, 200, 80, 60, 120)+C2H4peak1AreaFinder(C2H4peak1, time, FID1_100B.data, 200, 80, 60, 120))/2 100];

    PC2H4 = polyfit(SummaryC2H4_1(:,1),SummaryC2H4_1(:,2),1);
    C2H4fit = PC2H4(1)*Area;
    figure
    hold on
    plot(SummaryC2H4_1(:,1),SummaryC2H4_1(:,2),'o')
    plot(Area,C2H4fit)
    title('C_2H_4')
    xlabel('Area')
    ylabel('ppm')
    
%generating the calibration linear slope values
csvwrite('Autogenerated_Calibration.csv',[PH2(1); PCH4A(1); PCO(1); PCH4B(1); PC2H4(1)])
    
%%
function Area = CH4peak1AreaFinder(CH4peak1location, time, GCdata)
[val,idx] = min(abs(time-CH4peak1location)); %finding the index for Peak position

    %make all data of interest positive
    if min(GCdata((idx-200):(idx+180))) < 0
        shiftup = abs(min(GCdata((idx-200):(idx+180)))); %take absolute value
    else
        shiftup = 0;
    end
    
    if GCdata(idx) < 8000
        Rightin = 80;
        Rightout = 150;
    else
        Rightin = 150;
        Rightout = 200;
    end
    
tpseudo = [time((idx-200):(idx-150)),time((idx+Rightin):(idx+Rightout))];
datapseudo = [GCdata((idx-200):(idx-150)),GCdata((idx+Rightin):(idx+Rightout))];
P = polyfit(tpseudo,datapseudo,1);
yfit = P(1)*tpseudo+P(2);
BaselineArea = trapz(tpseudo, yfit+shiftup);
CH4peak1Area = trapz(time((idx-200):(idx+Rightout)), GCdata((idx-200):(idx+Rightout))+shiftup);
Area = CH4peak1Area - BaselineArea;
    figure %graphing to check fit
    hold all
    area(time((idx-200):(idx+Rightout)), GCdata((idx-200):(idx+Rightout))+shiftup)
    area(tpseudo, yfit+shiftup)
    title(['CH4 peak 1',num2str(Rightin)])
end
%%
function Area = H2peak1AreaFinder(H2peak1location, time, GCdata)
[val,idx]=min(abs(time-H2peak1location)); %finding the index for Peak position
    %make all data of interest positive
    if min(GCdata((idx-30):(idx+50))) < 0
        shiftup = abs(min(GCdata((idx-30):(idx+50))));
    else
        shiftup = 0;
    end
    
    if GCdata(idx) < 8000
        Rightin = 20;
        Rightout = 50;
    else
        Rightin = 30;
        Rightout = 50;
    end
tpseudo = [time((idx-30):(idx-18)),time((idx+Rightin):(idx+Rightout))];
datapseudo = [GCdata((idx-30):(idx-18)),GCdata((idx+Rightin):(idx+Rightout))];
timeH2 = time((idx-30):(idx+Rightout)); %all time
P = polyfit(tpseudo,datapseudo,4); %3rd order polyfit
yfit = P(1)*timeH2.^4+P(2)*timeH2.^3+P(3)*timeH2.^2+P(4)*timeH2+P(5);
BaselineArea = trapz(timeH2, yfit+shiftup);
H2peak1Area = trapz(time((idx-30):(idx+Rightout)), shiftup+GCdata((idx-30):(idx+Rightout)));
Area = H2peak1Area - BaselineArea;
    figure
    hold all
    area(time((idx-30):(idx+Rightout)), shiftup+GCdata((idx-30):(idx+Rightout)))
    area(timeH2, yfit+shiftup)
    title(['H2 peak 1 Rightin',num2str(Rightin)])
end
%%
function Area = COpeak1AreaFinder(COpeak1location, time, GCdata)
[val,idx]=min(abs(time-COpeak1location)); %finding the index for Peak position
    %make all data of interest positive
    if min(GCdata((idx-200):(idx+200))) < 0
        shiftup = abs(min(GCdata((idx-200):(idx+200)))); %take absolute value
    else
        shiftup = 0;
    end
tpseudo = [time((idx-200):(idx-150)),time((idx+150):(idx+200))];
datapseudo = [GCdata((idx-200):(idx-150)),GCdata((idx+150):(idx+200))];
P = polyfit(tpseudo,datapseudo,2);
yfit = P(1)*tpseudo.^2+P(2)*tpseudo+P(3);
BaselineArea = trapz(tpseudo, yfit+shiftup);
COpeak1Area = trapz(time((idx-200):(idx+200)), GCdata((idx-200):(idx+200))+shiftup);
Area = COpeak1Area - BaselineArea;
    figure %graphing to double check fit
    hold all
    area(time((idx-150):(idx+150)), GCdata((idx-150):(idx+150))+shiftup)
    area(tpseudo, yfit+shiftup)
    title('CO peak 1')
end
%%
function Area = CH4peak2AreaFinder(CH4peak2location, time, GCdata, OutLeft, InLeft, InRight, OutRight)
[val,idx]=min(abs(time-CH4peak2location)); %finding the index for Peak position
    %make all data of interest positive
    if min(GCdata((idx-OutLeft):(idx+OutRight))) < 0
        shiftup = abs(min(GCdata((idx-OutLeft):(idx+OutRight)))); %take absolute value
    else
        shiftup = 0;
    end
tpseudo = [time((idx-OutLeft):(idx-InLeft)),time((idx+InRight):(idx+OutRight))]; %time span of interest for peak
datapseudo = [GCdata((idx-OutLeft):(idx-InLeft)),GCdata((idx+InRight):(idx+OutRight))]; %GCdata span of interest for peak
timeCH4 = time((idx-OutLeft):(idx+OutRight)); %all time
P = polyfit(tpseudo,datapseudo,2); % calculating quadratic fit
yfit = P(1)*timeCH4.^2+P(2)*timeCH4+P(3); %quadratic fit for baseline
BaselineArea = trapz(timeCH4, yfit+shiftup); %ares of baseline
CH4peak2Area = trapz(time((idx-OutLeft):(idx+OutRight)), GCdata((idx-OutLeft):(idx+OutRight))+shiftup);
Area = CH4peak2Area - BaselineArea;
    figure
    hold all
    area(time((idx-OutLeft):(idx+OutRight)), GCdata((idx-OutLeft):(idx+OutRight))+shiftup)
    area(timeCH4, yfit+shiftup)
    title('CH4 peak 2')
end
%%
function Area = C2H4peak1AreaFinder(C2H4peak1location, time, GCdata,OutLeft, InLeft, InRight, OutRight)
[val,idx]=min(abs(time-C2H4peak1location)); %finding the index for Peak position
    %make all data of interest positive
    if min(GCdata((idx-OutLeft):(idx+OutRight))) < 0
        shiftup = abs(min(GCdata((idx-OutLeft):(idx+OutRight)))); %take absolute value
    else
        shiftup = 0;
    end
tpseudo = [time((idx-OutLeft):(idx-InLeft)),time((idx+InRight):(idx+OutRight))];
datapseudo = [GCdata((idx-OutLeft):(idx-InLeft)),GCdata((idx+InRight):(idx+OutRight))];
timeC2H4 = time((idx-OutLeft):(idx+OutRight));
% dataC2H4 = GCdata((idx-200):(idx+200));
P = polyfit(tpseudo,datapseudo,2);
yfit = P(1)*timeC2H4.^2+P(2)*timeC2H4+P(3);
BaselineArea = trapz(timeC2H4, yfit+shiftup); %ares of baseline
C2H4peak1Area = trapz(time((idx-OutLeft):(idx+OutRight)), GCdata((idx-OutLeft):(idx+OutRight))+shiftup);
Area = C2H4peak1Area - BaselineArea;
    figure
    hold all
    area(time((idx-OutLeft+20):(idx+OutRight)), GCdata((idx-OutLeft+20):(idx+OutRight))+shiftup)
    area(timeC2H4(20:end), yfit(20:end)+shiftup)
    title('C2H4 peak 1')
end