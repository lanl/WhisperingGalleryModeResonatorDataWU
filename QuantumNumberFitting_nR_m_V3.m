clear all clc

%This is a MATLAB script for assigning quantum numbers to Whispering
%gallery mode (WGM) spectra and determining the radius and index of
%refraction of the environment.
%This code first finds the peaks in the data spectrum and stores their
%locations in a matrix. Next, TE and TM mode locations are predicted using
%analytical equations over a specified range of radii and refractive
%indices of the surrounding media. For each radius-refractive index pair,
%observed peaks are matched to predicted modes based on their proximity. To
%obtain the best fit, a subset of the fits are selected that containn the
%highest number of matched peaks. Next, this subset of fits is searched to
%find the fit with the lowest error (i.e. difference in observed and
%predicted WGM locations). The fit with the lowest error is deemed the best
%fit.It should be noted that it is necessary to first identify fits that
%match the highest number of peaks since each new peak matched inherently
%adds to the amount of error possible. Now that quantum numbers have been
%assignned to peakss, the peak positions are used to fit for the radius and
%refractive inndex of the surrounding ennvironment.
%In version7a of this code, a new feature was added that would itteratively
%do the above proceedure to scan through parameter space and if not enough
%peaks were matched or the error was above a specified amount the search
%area in R-space is shifted by a user specified amount and the process is
%repeated until enough peaks are matched and the error is low enough.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Directions:
%To use this file, one must first make sure several parameters are suitable
%to the data that is to be fit. These include the index of refraction of
%the sphere material (nc),the range of environmental indices of refraction to
%search through (defined by the minimum (ne_min), maximum (ne_max) values
%to search through and the number of points (ne_np)), the range of radii to
%search throughh(defined by the minimum (Rmin)and maximum (Rmax) values to
%search through in micrometers as well as the number of points to search
%through (Rnp)), and finally, the range of angular quantum numbers to
%initially predict (n). Both the radius and refractive index search space
%should be tuned based on the environmental media in question and the
%size/distribution of the spheres employed. The range of angular quantum
%numbers to predict will also need to be adjusted depending on the sphere
%sizes in use with the larger sphere sizes necessitating a larger range. It
%should be noted that for larger numbers of refractive index and radius
%points to scan through the optimization process can take quite a lot of
%time, so minimizing your search space may be necessary to reduce
%computation time.
%Additionally, the amount of peaks that can be found in the data spectrum
%will depend on the peakfinding parameters thresh (basically a threshold
%intensity for finding peaks) and sel (basically defines how sharp and
%intense the feature needs to be relative to its surroundings and the min
%and max of the data set). The peaks that are found or left out can have a
%dramatic influence on the final fit results. These parameters should be
%changed as necessary to find the desired peaks in your data spectrum.
%The parameters mentioned above can all be found in the first few lines of
%code below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%Important Parameters:
nc=1.5597;%Index of refraction of sphere material
eps=nc^2;mu=1;
m_min=1.05;
m_max=1.2;
m_np=60;%Number of Refractive index points to scan through
Deltam=(m_max-m_min)/(m_np-1);%Step size for varying ne
mvect=m_min:Deltam:m_max;%Array containing all values of ne to span
c=299792458;%Speed of light in m/s
Rnmin=6;%Minimum sphere radius in micrometers to span for fittinng.
Rnmax=8;%Maximum sphere radius in micrometers to span for fittinng.
Rnnp=2001;%Number of points to span Radius for fittinng.
DeltaRn=(Rnmax-Rnmin)/(Rnnp-1);%Step size to vary R for fitting.
Rnvect=Rnmin:DeltaRn:Rnmax;%Array containing all values of R to span
RnspaceChange=1;%Specifies the increment that R-space will change by every
%time the later while loop proceeds to run which matches peaks and makes
%sure a certain number of peaks are matched.
Rnfinalshiftwind=0.45;%Specified the range around the first radius that the
%upcoming while loop finds that satisfies termination conditions. The loop
%will commense one more time trying to find the optimal R, ne pair with a
%new Rvector ranging from R_previousbest-Rfinalshiftwind to
%R_previousbest+Rfinalshiftwind.
n=(1:100);%Angular quantum numbers to predict
%First, lets get the data.
[filename, pathname]=uigetfile('*.*','Get Data');
ss=strcat(pathname,filename);
Datai=xlsread(ss); %Contains data to be analyzed

FilePathFileNameFileExt=[pathname filename];
[filepath, name, Ext]=fileparts(FilePathFileNameFileExt);


Datax=Datai(:,1);%Wavelength portion of data
Datay=Datai(:,2);%Intensity portion of data
[Datax,sortIdx]=sort(Datax,'descend');%Rearanges Data x wave in descending
%order.
Datay=Datay(sortIdx);%Rearanges the corresponding y wave in the same order.

figure(1)
plot3(Datax,Datay,1:length(Datax))
view(0,90)
xlabel('Wavelength (nm)')
ylabel('Emission Intensity (a.u.)')
datacursormode on

flag1='n';
while flag1=='n'
    nstart=input('How many starting points do you want to cut?');
    nstart=nstart+1;
    nend=input('How many end points do you want to cut?');
    ndend=length(Datax)-nend;
    Dataxnew=Datax(nstart:ndend);
    Dataynew=Datay(nstart:ndend);
    clf
    plot(Datax,Datay,'r',Dataxnew,Dataynew,'g')
    xlabel('Wavelength (nm)')
    ylabel('Emission Intensity (a.u.)')
    title('Endpoint Removal Area')
    flag1=input('Is it OK? (y/n)','s');
end
Datax=Dataxnew;
Datay=Dataynew;


Datay=Datay./sum(Datay.^2).^0.5;%Normalizes data
%Now lets find the peaks in the data:
selFactor=6;
sel=(max(Datay)-min(Datay))/selFactor;
thresh=0.022;
[peakloc,peakmag]=peakfinder(Datay,sel,thresh);
%Note that peakfinder is an external function

peakloc=Datax(peakloc);
figure(1)
plot(Datax,Datay,'-',peakloc,peakmag,'o')
L=length(peakloc);
flag2=input('Is this ok? (y/n)','s');
while flag2=='n'
    selFactor=input('Specify Selectivity Factor (#)');
    sel=(max(Datay)-min(Datay))/selFactor;
    thresh=input('Specify Threshold Value (#)');
    [peakloc,peakmag]=peakfinder(Datay,sel,thresh);
    peakloc=Datax(peakloc);
    clf
    figure(1)
    plot(Datax,Datay,'-',peakloc,peakmag,'o')
    L=length(peakloc);
    flag2=input('Is this ok? (y/n)','s');
end

%Now we set up some equations to predict WGM positions taking equations
%from [Pang, S.; Beckham, R. E.; Meissner, K. E. Appl. Phys. Lett. 2008,
%92, 2211108].

%Roots of the airy function
airy=[2.338, 4.088, 5.521, 6.787];

wTEqnr=@(q,n,nr,m) (c./(nr*1e-6))*((n+0.5)+2.^(-1/3).*airy(q)*(n+0.5).^...
    (1/3)-(m/sqrt(m^2-1))+(0.3*2^(-2/3))*(airy(q).^2)*(n+0.5).^(-1/3)+...
    ((2^(-1/3))*m^3*(2*1/3-1)/(m^2-1)^(3/2))*airy(q)*(n+0.5).^(-2/3));
wTMqnr=@(q,n,nr,m) (c./(nr*1e-6))*((n+0.5)+2.^(-1/3).*airy(q)*(n+0.5).^...
    (1/3)-(1/m/sqrt(m^2-1))+(0.3*2^(-2/3))*(airy(q).^2)*(n+0.5).^(-1/3)+...
    ((2^(-1/3))*m*(2/3/m^4-1)/(m^2-1)^(3/2))*airy(q)*(n+0.5).^(-2/3));

LamPredTE=@(q,n,nr,m) 2*pi()*c*(1e9)./wTEqnr(q,n,nr,m)';
LamPredTM=@(q,n,nr,m) 2*pi()*c*(1e9)./wTMqnr(q,n,nr,m)';

%Finally, before we evaluate parameter space we need to specify how many
%peaks must match before an acceptable fit is found.
NumberOfPeaks=L
NumbToMatchRequirement=input('How many peaks must be matched to have a suitable fit? (#)');
minAcceptableError=input('What is the minimum acceptable Error? (#)')';

%Now we start our while loop to evaluate chunks of parameter space:
PreStopCondition='n';
StopCondition='n';
while StopCondition=='n'
if PreStopCondition=='y'
    Rnvect=Rnvect(BestMinCols(1))-Rnfinalshiftwind:2*Rnfinalshiftwind/(Rnnp-1):...
        Rnvect(BestMinCols(1))+Rnfinalshiftwind;
end
%Now we are about to evaluate the parameter space. To increase efficiency,
%we will begin by prealocating several cells in advance for speed.

MatchUnique=zeros(length(mvect),length(Rnvect));
MatchUnique=num2cell(MatchUnique);%Prealocating for speed
TEq1_peaks=zeros(length(mvect),length(Rnvect));
TEq1_peaks=num2cell(TEq1_peaks);%Prealocating for speed
TMq1_peaks=zeros(length(mvect),length(Rnvect));
TMq1_peaks=num2cell(TMq1_peaks);%Prealocating for speed
%The below values are being prealocated for speed but will contain the
%Error, the number of TEq1 matches, number of TMq1 matches, and the number
%of total matches of each value of the radius tested.
ErrorFunc=zeros(length(mvect),length(Rnvect));%Prealocating for speed
LTEq1=zeros(length(mvect),length(Rnvect));%Prealocating for speed
LTMq1=zeros(length(mvect),length(Rnvect));%Prealocating for speed
LtotMatch=zeros(length(mvect),length(Rnvect));%Prealocating for speed

%Now we search the parameter space of m annd ns x R
for ii=1:m_np
for i=1:Rnnp
%Set up for loop to run through all desired values of R and ne
%First, predict the TE/TM mode locations
TEq1pred=LamPredTE(1,n,Rnvect(i),mvect(ii));
TMq1pred=LamPredTM(1,n,Rnvect(i),mvect(ii));
%Combine predictions and differentiate them with a 1 (TE modes) or 2 (TM
%modes)
modetype=vertcat(ones(size(TEq1pred)),2*ones(size(TMq1pred)));
modeN=vertcat(n',n');
ModePredFull=vertcat(TEq1pred,TMq1pred);
[ModePredFull,sortIdx]=sort(ModePredFull,'descend');
modetype=modetype(sortIdx);
modeN=modeN(sortIdx);
ModePredFull=horzcat(ModePredFull,modetype,modeN);
%Now we have an array where the first collumn containss each of the
%predicted wavelengths for WGM resonances, the second collumn contains a 1
%or 2 denoting mode type, i.e. 1 for TE and 2 for TM, and the third collum
%contains the n angular quantum number.

%Now we will begin matching observed and predicted peaks.
Tolerance=2;%Tolerance in nm for matching peaks.
Match=[0,0,0,0,0,0];
Matchpart=[0,0,0,0,0,0];
for k=1:L
    dk=abs(peakloc(k)-ModePredFull(:,1));%Column of absolute differences
    %beteen the kth observed peak and each of the predicted peaks.
    if k==1
    dsub=dk;
    else
    dsub=horzcat(dsub,dk); %Building the difference matrix
    end
    dmin=find(dsub(:,k)==min(dsub(:,k)),1,'first');
    if dsub(dmin,k) <= Tolerance
    IndexofBestMatch=dmin;
    elseif dsub(dmin,k) > Tolerance
    IndexofBestMatch=[];
    end
    if IndexofBestMatch >= 1
    Matchpart=[peakloc(k),ModePredFull(IndexofBestMatch,1),...
        ModePredFull(IndexofBestMatch,2),ModePredFull(IndexofBestMatch,3)...
        ,ModePredFull(IndexofBestMatch,1)-peakloc(k),peakmag(k)];
        if Match==[0,0,0,0,0,0]
            Match=Matchpart;
        elseif Match ~=[0,0,0,0,0,0]
            Match=vertcat(Match,Matchpart);
        end
    end
    
end
%Now that we have matched the predicted and observed WGM resonnances, we
%must ensure that the matches are unique, i.e. that there are not multiple
%peaks in the observed spectrum that are matched to a single predicted
%peak.
MatchUnique{ii,i}=Match;
%Step 1: Find repeats in Match(:,2)
X=Match(:,2);

if length(X)>1 %Adresses bug in case of 1 or fewer matches
uniqueX=unique(X);%Finds unique peak matches in the second column of the
%peak matching list.
if length(uniqueX) > 1
countofX=hist(X,uniqueX);%Makes histogram with number of bins equal to the
%number of unique peak matches.
IndexToRepeatedValues=(countofX ~=1);
RepeatedValues=uniqueX(IndexToRepeatedValues);
%Step 2: Determine best matches for repeated values in Match(:,2).
elseif length(uniqueX) == 1
    RepeatedValues=X(find(X == uniqueX));
end
ForFinalRemoval=0;
for k=1:length(RepeatedValues)
    IndexofRepeatsLeft=find(X==RepeatedValues(k));
    IndexofBetterMatch=find(abs(Match(IndexofRepeatsLeft,5))==...
        min(abs(Match(IndexofRepeatsLeft,5))));
    IndexofRepeatsLeft(IndexofBetterMatch)=[];
    %Removes the index of the peak that we will want to keep because it is
    %the best possible match.
    if k==1
    ForFinalRemoval=IndexofRepeatsLeft;
    else
    ForFinalRemoval=vertcat(ForFinalRemoval,IndexofRepeatsLeft);
    end
end
%Step 3: Remove least matching duplicates
if ForFinalRemoval~=0
MatchUnique{ii,i}(ForFinalRemoval,:)=[];
end
end
%Step 4: Sort peaks by peak type
TEq1_MatchIndexes=find(MatchUnique{ii,i}(:,3)==1);
TEq1_peaks{ii,i}=MatchUnique{ii,i}(TEq1_MatchIndexes,:);
TMq1_MatchIndexes=find(MatchUnique{ii,i}(:,3)==2);
TMq1_peaks{ii,i}=MatchUnique{ii,i}(TMq1_MatchIndexes,:);

ErrorFunc(ii,i)=norm(MatchUnique{ii,i}(:,1)-MatchUnique{ii,i}(:,2));


LTEq1(ii,i)=length(TEq1_peaks{ii,i}(:,1));
LTMq1(ii,i)=length(TMq1_peaks{ii,i}(:,1));
LtotMatch(ii,i)=LTEq1(ii,i)+LTMq1(ii,i);
end
end

figure(2)
surf(Rnvect,mvect,ErrorFunc)
grid off
view(0,90)
xlabel('n_s x Radius (RIU x um)')
ylabel('Contrast Ratio')
zlabel('Error')
title('Error')
shading interp
colorbar
figure(3)
surf(Rnvect,mvect,LtotMatch)
shading interp
view(0,90)
xlabel('n_s x Radius (RIU x um)')
ylabel('Contrast Ratio')
zlabel('Number of Matched Peaks')
title('Number of Matched Peaks')
colorbar
%Now that the parameter space has been assessed, we need to find the
%optimal R-ne pair.

%First, we find the coordinates that contain the highest number of matched
%peaks
idx=find(LtotMatch==max(LtotMatch,[],'all'));
%Now we convert our findings to linear indices:

%Now we use these indices of the locations of the highest matched peaks and
%look through the error function at these same locations to figure out
%where the lowest error of this set of positions is.
%[BestMinRows,BestMinCols]=mink(ErrorFunc(idx),3);
[bestmins,bestminidx]=mink(ErrorFunc(idx),3);
BestMinLocs=idx(bestminidx);
[BestMinRows,BestMinCols]=ind2sub(size(ErrorFunc),BestMinLocs);

%Plot the best solution found
figure(4)
plot3(Datax,Datay,1:length(Datax),'b-',TEq1_peaks{BestMinLocs(1)}(:,2),...
    TEq1_peaks{BestMinLocs(1)}(:,6),TEq1_peaks{BestMinLocs(1)}(:,4),'rx',...
    TMq1_peaks{BestMinLocs(1)}(:,2),TMq1_peaks{BestMinLocs(1)}(:,6),...
    TMq1_peaks{BestMinLocs(1)}(:,4),'ko')
view(0,90)
xlabel('Wavelength (nm)')
ylabel('Intensity (a.u.)')
title('1st best Fit')
legend('Data','TE Modes','TM Modes')


Rnvect(BestMinCols)%Displays Radii acquired from top 3 fits
mvect(BestMinRows)%Displays the ne acquired from top 3 fits
idx=BestMinLocs(1);%Gives index of best fit.
if PreStopCondition=='y'
    StopCondition='y';
    break
end
if max(LtotMatch,[],'all')>= NumbToMatchRequirement && ...
        ErrorFunc(idx) < minAcceptableError
    PreStopCondition='y';
else
    Rnvect=Rnvect+RnspaceChange;

end
end
%Now that we have a good approximated fit, we can use the quantum numbers
%acquired above to determine both the radius and ratio of the index of
%refraction of the sphere and the surrounding media.

figure(5)
obj_fun=@(params) norm(LamPredTE(1,TEq1_peaks{idx}(:,4),params(1),...
    params(2))'-TEq1_peaks{idx}(:,1))+norm(LamPredTM(1,TMq1_peaks{idx}...
    (:,4),params(1),params(2))'-TMq1_peaks{idx}(:,1));

options=optimset('PlotFcns',@optimplotfval);
%Note, this plotting step doesn't always work but it does seem to still
%optimize the parameters
sol=fminsearch(obj_fun,[Rnvect(BestMinCols(1)),mvect(BestMinRows(1))],options)
error=obj_fun(sol)%The error of the final fit
ContrastRatio=sol(2)% Index of refraction of surrounding medium

%Now we catalogue the TM and TE peaks. The first column contains the
%observed peak locations, the second contains the best predicted locations,
%the third contains the angular quantum numbers, the fourth contains the
%difference in wavelength, and the fifth contains the magnitude or hight of
%each peak in the normalized spectrum.
TMpeaks=[TMq1_peaks{idx}(:,1), LamPredTM(1,TMq1_peaks{idx}(:,4),...
    sol(1),sol(2))',TMq1_peaks{idx}(:,4),TMq1_peaks{idx}(:,5),...
    TMq1_peaks{idx}(:,6)];
TEpeaks=[TEq1_peaks{idx}(:,1), LamPredTE(1,TEq1_peaks{idx}(:,4),...
    sol(1),sol(2))',TEq1_peaks{idx}(:,4),TEq1_peaks{idx}(:,5),...
    TEq1_peaks{idx}(:,6)];


%Now we plot the predicted peak locations of the final result.
figure(5)
clf
plot3(Datax,Datay,1:length(Datax),'b-',TEpeaks(:,2),TEpeaks(:,5),...
    TEpeaks(:,3),'rx',TMpeaks(:,2),TMpeaks(:,5),TMpeaks(:,3),'ko')
view(0,90)
xlabel('Wavelength (nm)')
ylabel('Intensity (a.u.)')
title('Optimized Fit of WGM Spectrum')
legend('Data','TE Modes','TM Modes')


%Now we use two forloops to add peak labels
for i=1:length(TMpeaks(:,1))
    numbtext=num2str(TMpeaks(i,3));
    PeakQN={strcat('TM','_',numbtext(1),'_',numbtext(2))};
    text(TMpeaks(i,1),TMpeaks(i,5),PeakQN)
end
for i=1:length(TEpeaks(:,1))
    numbtext=num2str(TEpeaks(i,3));
    PeakQN={strcat('TE','_',numbtext(1),'_',numbtext(2))};
    text(TEpeaks(i,1),TEpeaks(i,5),PeakQN)
end

%Finally, we add a textbox with the optimized sphere RI x Radius product
%and the contrast ratio.

str1=strcat( 'R = ',num2str(sol(1)),'RIU x um');
str2=strcat('m = ', num2str(sol(2)));
strfull={str1,str2};
annotation('textbox',[.15 .6 .3 .3], 'String',strfull,'FitBoxToText','on')
f=gcf;

TEforSaving=TEpeaks(:,1:3);
TEdesignator=ones(size(TEpeaks(:,1)));
TEforSaving=horzcat(TEforSaving,TEdesignator);
TMforSaving=TMpeaks(:,1:3);
TMdesignator=2*(ones(size(TMpeaks(:,1))));
TMforSaving=(TMforSaving);
TMforSaving=horzcat(TMforSaving,TMdesignator);
Headers=[string('Peaks_exp') string('Peaks_theo') string('Mode Number') ...
    string('1TE_2TM')];
ForSaving=vertcat(Headers,TEforSaving, TMforSaving);

flag3='n';
flag3=input('Do you want to save fit results? (y/n)','s');
while flag3=='y'
    nametag=[name,'Fit'];
    fullFileName=fullfile(pathname, nametag);
    xlswrite(fullFileName,ForSaving);
    exportgraphics(f,[fullFileName '.png'],'Resolution',600')
    flag3='n';
end
