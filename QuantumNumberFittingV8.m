
%© 2022. Triad National Security, LLC. All rights reserved.
%This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
%National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
%Department of Energy/National Nuclear Security Administration. All rights in the program are
%reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
%Security Administration. The Government is granted for itself and others acting on its behalf a
%nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
%derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
%others to do so.

%This program is open source under the BSD-3 License.
%Redistribution and use in source and binary forms, with or without modification, are permitted
%provided that the following conditions are met:

%1. Redistributions of source code must retain the above copyright notice, this list of conditions and
%the following disclaimer.

%2.Redistributions in binary form must reproduce the above copyright notice, this list of conditions
%and the following disclaimer in the documentation and/or other materials provided with the
%distribution.

%3.Neither the name of the copyright holder nor the names of its contributors may be used to endorse
%or promote products derived from this software without specific prior written permission.

%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
%IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
%CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
%PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
%OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
%WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
%OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.




clear all clc

%This is a MATLAB script for assigning quantum numbers to Whispering
%gallery mode (WGM) spectra and determining the radius and index of
%refraction of the environment.

%This code first uploads the desired data sets to be fit and asks the user
%if they would like to remove some of the endpoint of each dataset. Next,
%the peaks in each dataset are identified using a peakfinding function
%written for MATLAB by Nathanael Yoder ([Nathanael Yoder (2021).
%peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
%https://www.mathworks.com/matlabcentral/fileexhange/25500-peakfinder-x0
%-extrema-includeendpoints-interpolate).MATLAB Central File Exchange]) and
%their locations stored in a matrix. Each spectrum is then fit one at a
%time. The fitting program seeks to identify every peak that was found
%using the above peakfinding function. To do this, TE and TM mode locations
%are predicted using analytical equations over a specified initial range of
%radii and refractive indices of the surrounding media. For each
%radius-refractive index pair, observed peaks are matched to predicted
%modes based on their proximity. To be considered a possible match, the
%predicted peaks must be within a specified tolerance of the observed peak
%(currently set at 2 nm). In the case that multiple predicted peaks are
%matched to a single observed peak, only the predicted peak with the
%closest resonance frequency is kept. The quality of the assignments are
%then judged both by the number of peaks that could be matched and the
%error in the predicted and observed peak locations as defined by the sum
%of the square of the differences in the observed and predicted peak
%locations in nm. The subset of fits that yields the highest  number of
%matched peaks is searched for the lowest error. If this fit could not
%match every observed peak then a new fitting session is initiated in which
%the new search grid occurs over larger sphere radii. This process is
%repeated until a fit can be found that can identify each of the peaks
%found in the experimental spectrum. At this point, a new fit sequence is
%initiated that covers a smaller range of radii around the previous best
%fit and has a higher density of points to search through. The results and
%peak assignments for this fit are then fed into a non-linear solver to
%obtain a refined radius, external refractive index output. The same
%procedure is then repeated for each fo the uploaded spectra and the
%results are saved as they are obtained.

%Additionally, the amount of peaks that can be found in the data spectrum
%will depend on the peakfinding parameters thresh (basically a threshold
%intensity for finding peaks) and sel (basically defines how sharp and
%intense the feature needs to be relative to its surroundings and the min
%and max of the data set). The peaks that are found or left out can have a
%dramatic influence on the final fit results. These parameters should be
%changed as necessary to find the desired peaks in your data spectrum.
%The parameters mentioned above can all be found in the first few lines of
%code below.

%Note: Each dataset must be uploaded individually. The data files should be
%excel files with the first column containing the wavelength information
%and the second column containing the intensity information of the
%spectrum. The top line in each column should contain the name of the
%dataset.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=1;%Indexing for while loop that controls data uploading
NameCell=[];
NameCell=num2cell(NameCell);
PathCell=[];
PathCell=num2cell(PathCell);
DataCell=[];
DataCell=num2cell(DataCell);
PeaklocCell=[];
PeaklocCell=num2cell(PeaklocCell);
PeakmagCell=[];
PeakmagCell=num2cell(PeakmagCell);
Lvect=[];%Vector containing the #of peaks in each dataset to fit for.
flagalpha='y';
while flagalpha == 'y'




%Important Parameters:
nc=1.5497;%1.5597;%1.44091;%Index of refraction of sphere material
ne_min=1.3;%Initial Minimum external RI
ne_max=1.4;%Initial Maximum external RI to search for
ne_np=80;%Number of Refractive index points to scan through
Deltane=(ne_max-ne_min)/(ne_np-1);%Step size for varying ne
ne=ne_min:Deltane:ne_max;%Array containing all values of ne to span
c=299792458;%Speed of light in m/s
Rmin=3.0;%Minimum sphere radius in micrometers to span for fittinng.
Rmax=4.0;%Maximum sphere radius in micrometers to span for fittinng.
Rnp=2001;%Number of points to span Radius for fittinng.
DeltaR=(Rmax-Rmin)/(Rnp-1);%Step size to vary R for fitting.
R=Rmin:DeltaR:Rmax;%Array containing all values of R to span
Rcopy=R;
RspaceChange=1;%Specifies the increment that R-space will change by every
%time the later while loop proceeds to run which matches peaks and makes
%sure a certain number of peaks are matched.
Rfinalshiftwind=0.3;%Specified the range around the first radius that the
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

NameCell{t}=name;
PathCell{t}=pathname;   

Datax=Datai(:,1);%Wavelength portion of data
Datay=Datai(:,2);%Intensity portion of data
[Datax,sortIdx]=sort(Datax,'descend');%Rearanges Data x wave in descending
%order.
Datay=Datay(sortIdx);%Rearanges the corresponding y wave in the same order.

figure(1)
clf
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

DataCell{t}=horzcat(Datax,Datay);%Storing Data

selFactor=4;
sel=(max(Datay)-min(Datay))/selFactor;
thresh=0.022;
[peakloc,peakmag]=peakfinder(Datay,sel,thresh);
%Note that peakfinder is an external function

peakloc=Datax(peakloc);

figure(1)
plot(Datax,Datay,'-',peakloc,peakmag,'o')
L=length(peakloc);
Lvect(t)=L;%Stores the # of peaks used for fitting
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
    Lvect(t)=L;%Stores the # of peaks used for fitting
    flag2=input('Is this ok? (y/n)','s');
end
PeaklocCell{t}=peakloc;
PeakmagCell{t}=peakmag;
t=t+1;
flagalpha=input('Do you want to add another data set? (y/n)','s');
NameOfLast=name
end
t=t-1;

%Now we set up some equations to predict WGM positions taking equations
%from [Pang, S.; Beckham, R. E.; Meissner, K. E. Appl. Phys. Lett. 2008,
%92, 2211108].

%Roots of the airy function
airy=[2.338, 4.088, 5.521, 6.787];

wTEqnr=@(q,n,r,m) (c./(r*1e-6)/nc)*((n+0.5)+2.^(-1/3).*airy(q)*(n+0.5).^...
    (1/3)-(m/sqrt(m^2-1))+(0.3*2^(-2/3))*(airy(q).^2)*(n+0.5).^(-1/3)+...
    ((2^(-1/3))*m^3*(2*1/3-1)/(m^2-1)^(3/2))*airy(q)*(n+0.5).^(-2/3));
wTMqnr=@(q,n,r,m) (c./(r*1e-6)/nc)*((n+0.5)+2.^(-1/3).*airy(q)*(n+0.5).^...
    (1/3)-(1/m/sqrt(m^2-1))+(0.3*2^(-2/3))*(airy(q).^2)*(n+0.5).^(-1/3)+...
    ((2^(-1/3))*m*(2/3/m^4-1)/(m^2-1)^(3/2))*airy(q)*(n+0.5).^(-2/3));

LamPredTE=@(q,n,r,m) 2*pi()*c*(1e9)./wTEqnr(q,n,r,m)';
LamPredTM=@(q,n,r,m) 2*pi()*c*(1e9)./wTMqnr(q,n,r,m)';

%Begin forloop that will analyze each dataset
for Dataidx=1:t
    
%Finally, before we evaluate parameter space we need to specify how many
%peaks must match before an acceptable fit is found.
NumberOfPeaks=Lvect(Dataidx)
NumbToMatchRequirement=NumberOfPeaks;
minAcceptableError=10;

%Now we start our while loop to evaluate chunks of parameter space:
PreStopCondition='n';
StopCondition='n';
while StopCondition=='n'
if PreStopCondition=='y'
    R=R(BestMinCols(1))-Rfinalshiftwind:2*Rfinalshiftwind/(Rnp-1):...
        R(BestMinCols(1))+Rfinalshiftwind;
end
%Now we are about to evaluate the parameter space. To increase efficiency,
%we will begin by prealocating several cells in advance for speed.

MatchUnique=zeros(length(ne),length(R));
MatchUnique=num2cell(MatchUnique);%Prealocating for speed
TEq1_peaks=zeros(length(ne),length(R));
TEq1_peaks=num2cell(TEq1_peaks);%Prealocating for speed
TMq1_peaks=zeros(length(ne),length(R));
TMq1_peaks=num2cell(TMq1_peaks);%Prealocating for speed
%The below values are being prealocated for speed but will contain the
%Error, the number of TEq1 matches, number of TMq1 matches, and the number
%of total matches of each value of the radius tested.
ErrorFunc=zeros(length(ne),length(R));%Prealocating for speed
LTEq1=zeros(length(ne),length(R));%Prealocating for speed
LTMq1=zeros(length(ne),length(R));%Prealocating for speed
LtotMatch=zeros(length(ne),length(R));%Prealocating for speed

%Now we search the parameter space of ne annd R
for ii=1:ne_np
for i=1:Rnp
%Set up for loop to run through all desired values of R and ne
%First, predict the TE/TM mode locations
TEq1pred=LamPredTE(1,n,R(i),nc/ne(ii));
TMq1pred=LamPredTM(1,n,R(i),nc/ne(ii));
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
for k=1:Lvect(Dataidx)
    dk=abs(PeaklocCell{Dataidx}(k)-ModePredFull(:,1));%Column of absolute
    %differences between the kth observed peak and each of the predicted
    %peaks.
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
    Matchpart=[PeaklocCell{Dataidx}(k),ModePredFull(IndexofBestMatch,1),...
        ModePredFull(IndexofBestMatch,2),ModePredFull(IndexofBestMatch,3)...
        ,ModePredFull(IndexofBestMatch,1)-PeaklocCell{Dataidx}(k),...
        PeakmagCell{Dataidx}(k)];
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
surf(R,ne,ErrorFunc)
grid off
view(0,90)
xlabel('Radius (um)')
ylabel('n_e (RIU)')
zlabel('Error')
title('Error')
shading interp
colorbar
figure(3)
surf(R,ne,LtotMatch)
shading interp
view(0,90)
xlabel('Radius (um)')
ylabel('n_e (RIU)')
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
plot3(DataCell{Dataidx}(:,1),DataCell{Dataidx}(:,2),1:length...
    (DataCell{Dataidx}(:,1)),'b-',TEq1_peaks{BestMinLocs(1)}(:,2),...
    TEq1_peaks{BestMinLocs(1)}(:,6),TEq1_peaks{BestMinLocs(1)}(:,4),'rx'...
    ,TMq1_peaks{BestMinLocs(1)}(:,2),TMq1_peaks{BestMinLocs(1)}(:,6),...
    TMq1_peaks{BestMinLocs(1)}(:,4),'ko')
view(0,90)
xlabel('Wavelength (nm)')
ylabel('Intensity (a.u.)')
title('1st best Fit')
legend('Data','TE Modes','TM Modes')


R(BestMinCols)%Displays Radii acquired from top 3 fits
ne(BestMinRows)%Displays the ne acquired from top 3 fits
idx=BestMinLocs(1);%Gives index of best fit.
if PreStopCondition=='y'
    StopCondition='y';
    break
end
if max(LtotMatch,[],'all')>= NumbToMatchRequirement && ...
        ErrorFunc(idx) < minAcceptableError
    PreStopCondition='y';
else
    R=R+RspaceChange;

end
end
%Reinitialize stop conditions for the above while loop
PreStopCondition='n';
StopCondition='n';

%Now that we have a good approximated fit, we can use the quantum numbers
%acquired above to determine both the radius and ratio of the index of
%refraction of the sphere and the surrounding media.

figure(5)
obj_fun=@(params) norm(LamPredTE(1,TEq1_peaks{idx}(:,4),params(1),...
    nc/params(2))'-TEq1_peaks{idx}(:,1))+norm(LamPredTM(1,TMq1_peaks{idx}...
    (:,4),params(1),nc/params(2))'-TMq1_peaks{idx}(:,1));

options=optimset('PlotFcns',@optimplotfval);
%Note, this plotting step doesn't always work but it does seem to still
%optimize the parameters
sol=fminsearch(obj_fun,[R(BestMinCols(1)),ne(BestMinRows(1))],options)
error=obj_fun(sol)%The error of the final fit
n_medium=sol(2)% Index of refraction of surrounding medium

%Now we catalogue the TM and TE peaks. The first column contains the
%observed peak locations, the second contains the best predicted locations,
%the third contains the angular quantum numbers, the fourth contains the
%difference in wavelength, and the fifth contains the magnitude or hight of
%each peak in the normalized spectrum.
TMpeaks=[TMq1_peaks{idx}(:,1), LamPredTM(1,TMq1_peaks{idx}(:,4),...
    sol(1),nc/sol(2))',TMq1_peaks{idx}(:,4),TMq1_peaks{idx}(:,5),...
    TMq1_peaks{idx}(:,6)];
TEpeaks=[TEq1_peaks{idx}(:,1), LamPredTE(1,TEq1_peaks{idx}(:,4),...
    sol(1),nc/sol(2))',TEq1_peaks{idx}(:,4),TEq1_peaks{idx}(:,5),...
    TEq1_peaks{idx}(:,6)];


%Now we plot the predicted peak locations of the final result.
figure(5)
clf
plot3(DataCell{Dataidx}(:,1),DataCell{Dataidx}(:,2),1:length(...
    DataCell{Dataidx}(:,1)),'b-',TEpeaks(:,2),TEpeaks(:,5),...
    TEpeaks(:,3),'rx',TMpeaks(:,2),TMpeaks(:,5),TMpeaks(:,3),'ko')
view(0,90)
xlabel('Wavelength (nm)')
ylabel('Intensity (a.u.)')
title('Optimized Fit of WGM Spectrum')
legend('Data','TE Modes','TM Modes')
ylim([0 0.02+max(DataCell{Dataidx}(:,2))])

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

%Finally, we add a textbox with the optimized radius, contrast ratio, and
%environmental index of refraction.
str1=strcat( 'R = ',num2str(sol(1)),'um');
str2=strcat('m = ', num2str(nc/sol(2)));
str3=strcat('n_e = ', num2str(n_medium));
strfull={str1,str2,str3};
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

%Now we save the results
    %nametag=[NameCell{Dataidx},'Fit'];
    %fullFileName=fullfile(PathCell{Dataidx}, nametag);
    %xlswrite(fullFileName,ForSaving);
    nametag1=[NameCell{Dataidx},'Fit','.xls'];
    fullFileName1=fullfile(PathCell{Dataidx}, nametag1);
   writematrix(ForSaving,fullFileName1,'WriteMode','overwrite') 
    
    nametag2=[NameCell{Dataidx},'Fit'];
    fullFileName2=fullfile(PathCell{Dataidx}, nametag2);
    exportgraphics(f,[fullFileName2 '.png'],'Resolution',600)
    
    R=Rcopy;%Resetting R for the next itteration of the for loop/
end
