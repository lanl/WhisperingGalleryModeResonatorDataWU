clear all clc

%First, we obtain data


[filename_Ext, pathname]=uigetfile('*.txt*','Get Data');
FilePathFileNameFileExt=[pathname filename_Ext];
[filepath, name,ext]=fileparts(FilePathFileNameFileExt);
ss=strcat(pathname,filename_Ext);
data=xlsread(ss);



%Now we need to reorganize the data so that each spectrum is separately
%accessible.

Npix=length(data(:,1));%Number of pixels
Nspectra=length(data(1,:))-1%Gives the number of spectra collected
Specnumb=1:Nspectra;%Vector containing the index number of each spectrum
wavelength=data(:,1);%Contains the wavelength vector
DataCell=zeros(Nspectra,1);
DataCell=num2cell(DataCell);
A_fullDatMat=zeros(Npix,Nspectra); %Prealocating for speed
WavelengthCuttoff=540;%Wavelength in nm to serve as the low end wavelength
%cuttoff
for i=1:Nspectra
   Spec_i=data(:,i+1);
   DataCell{i}=horzcat(wavelength,Spec_i);
   
   Datax=DataCell{i}(:,1);
   Datay=DataCell{i}(:,2);
   [Datax,sortIdx]=sort(Datax,'descend');%Rearanges Data x wave in 
   %descending order.
   Datay=Datay(sortIdx);%Rearanges the corresponding y wave in same order
   DataCell{i}(:,1)=Datax;
   DataCell{i}(:,2)=Datay;
  A_fullDatMat(:,i)=DataCell{i}(:,2);
end
wavelength=sort(wavelength,'descend');%Rearanges wavelengths in descending
%order

%Add a universal artifact remover
weights=ones(size(wavelength));
N=length(wavelength);
D=diff(speye(N),2);
H=D'*D;

flag1='n';
flag1=input('Remove Universal Artifact? (y/n)','s');
while flag1=='y'
    nstart=input('Specify starting point?');
    nend=input('Specify ending point?');
    weights(nstart:nend)=0;
    W=spdiags(weights,0,N,N);
    C=chol(W+H);
    for i=1:Nspectra
        z0=C\(C'\(weights.*DataCell{i}(:,2)));
        DataCell{i}(:,2)=z0;
        A_fullDatMat(:,i)=z0;
    end
    flag1='n';
end


%Since the Witec Instrument tends to yield an intrinsic baseline offset, we
%wil first remove this offset.
for i=1:Nspectra
DataCell{i}(:,2)=DataCell{i}(:,2)-mean(DataCell{i}(Npix-30:Npix,2));
end

%Next, we remove all pixels corresponding to wavelengths below the
%specified wavelength cuttoff.
CutIdx=find(wavelength < WavelengthCuttoff);
wavelength(CutIdx)=[];
A_fullDatMat(CutIdx,:)=[];
for i=1:Nspectra
    DataCell{i}(CutIdx,:)=[];
end


%Now that the data is organized we can start to inspect it. The following
%allows the user to inspect the data at any desired time step from 1 to
%Nspectra
if Nspectra>1
figure(1)
surf(Specnumb,wavelength,A_fullDatMat)
shading interp
xlabel('Spectrum Number')
ylabel('Wavelength (nm)')
zlabel('Emission Intensity')
end

flag2='y';
flag3ref=1;
while flag2=='y'
    Nspectra
    
    flag3=input('Plot Data? 0 for no, 1-Nspectra for desired spectrum     ');
    if flag3==0
        flag2='n';
        flag3=flag3ref;
    end
    figure(2)
    plot(DataCell{flag3}(:,1),DataCell{flag3}(:,2))
    xlabel('Wavelength (nm)')
    ylabel('Intensity')
    
end



%Now we can start to process each collected spectrum. We will start with an
%automatic endpoint truncation and follow with an automatic baseiline
%correction.
EndThreshold=20; %Threshold value for designating an endpoint to a dataset
MinNumbPixels=300; %Specifies the minimum number of pixels required to keep
% a given spectrum. If it is below this value, we get rid of it.

for i=1:Nspectra
nstart=find(DataCell{i}(:,2) >=EndThreshold,1,'first');%Finds first point
%above a threshold value
nend=find(DataCell{i}(:,2) >=EndThreshold,1,'last');%Finds last point above
%a threshold value
DataCell{i}=DataCell{i}(nstart:nend,1:2);%Removes all points outside the
%boundaries set by the above nstart and nend points.
end
%Now we remove the data with insufficiently long spectra
SpecLengths=zeros(Nspectra,1);%Prealocating for speed
for i=1:Nspectra
    SpecLengths(i)=length(DataCell{i}(:,2));
end
SaveIdx=find(SpecLengths>=MinNumbPixels);%Finds all data with at least the
%minimum number of pixels allowed. The next line keeps only these spectra.
DataCell=DataCell(SaveIdx);




%Next, we perform baseline correction on each spectrum that is left using
%the whittaker smoothing algorithm

lambda=10000;
ratio=0.0001;
Threshold=0.0003;
PP=0.001;%Asymmetry Parameter
%Estimate baseline with arPLS

for i=1:length(DataCell)
    y0=DataCell{i}(:,2);
    N=length(y0);
    D=diff(speye(N),2);
    H=lambda*(D)'*D;
    weights=ones(N,1);
    t=1;
   W=spdiags(weights(:,t),0,N,N);
    %Cholesky decomposition
    C=chol(W+H);
    z0=C\(C'\(weights(:,t).*y0));
    zzz=z0;
    t=t+1;
    d=y0-z0;
    %creates new weighting wave based on initial fit.
    for ii=1:N;
        if d(ii)<0
            weights(ii,t)=1-PP;
        else
            weights(ii,t)=PP;
        end
    end
%Iteratively uses Whittaker smoothing and weighting redetermination to
%calculate a baseline. After the while loop is completed the residual
%between the initial data and the new background is calculated then
%normalized and a broad moving average of this wave is calculated. This
%moving average is meant to raise the value of this residual in the
%vicinity of peaks in the data. Finally a new weighting wave is generated
%by using a very low threshold value. Anything below the threshold is
%weighted highly and anything above the threshold is weighted very little.
%Then the data is reinterpolated using the Whittaker method and the new
%weighting wave.
    while true
    W=spdiags(weights(:,t),0,N,N);
    %Cholesky decomposition
    C=chol(W+H);
    z0=C\(C'\(weights(:,t).*y0));
    d=y0-z0;
    t=t+1;
    for iii=1:N
        if d(iii)<0
            weights(iii,t)=1-PP;
        else
            weights(iii,t)=PP;
        end
    end
    QQQ=norm(weights(:,t-1)-weights(:,t))/norm(weights(:,t-1));
        if QQQ< ratio, break; end
        if t>50 break; end
    end
    DataCell{i}(:,2)=(y0-z0)./z0;%Subtracts and scales data by baseline
end

%Next, we attempt to remove spectra that don't have strong WGMs by
%analyzing their derivatives.
lambdaStep=abs(DataCell{1}(1,1)-DataCell{1}(2,1));%x-axis step size in nm
DataDeriv=zeros(size(DataCell));
DataDeriv=num2cell(DataDeriv);%Prealocating for speed
DataDerivMax=zeros(size(DataDeriv));%To store the max amplitude of the 2nd
%derivatives calculated below
zerocount=DataDerivMax;
StoreThresh=1;%derivative magnitude Threshold for keeping data. i.e.
%if the max absolute value of the 2nd derivative is below the storage
%threshold chosen here, then that data set will be discrarded.


%zeroThresh=500;%Number of zeros allowed in the derivative
for i=1:length(DataCell)
   DataDeriv{i}=diff(DataCell{i}(:,2)./max(DataCell{i}(:,2)))./lambdaStep;
   
   DataDeriv{i}=smoothdata(DataDeriv{i},'movmedian',8);
   DataDerivMax(i)=max(abs(DataDeriv{i}));
   zerocount(i)=sum(DataDeriv{i} == 0);
end
SaveIdx=find(DataDerivMax >= StoreThresh); %Finds all the spectra that
%satisfy the storage conditions specified above.
DataCell=DataCell(SaveIdx);%Gets rid of data which don't have large enough
%2nd derivative magnitudes

%zerocount=zerocount(SaveIdx);
%SaveIdx=find(zerocount <= zeroThresh);
%DataCell=DataCell(SaveIdx);%Gets rid of data for which there are too many
%zeros in the 2nd derivative which can happen if the detector is maxed out
%over too much of the spectrum.


%Next, we can decide if we want to visualize any baselie corrected
%spectrum.
flag4='y';
flag2='y';
flag3ref=1;
NumberOfSpectra=length(DataCell);
while flag4=='y'
    flag4=input('Would you like to view spectra? (y/n)','s');
    if flag4=='n'
        flag2='n';
    end
while flag2=='y'
    NumberOfSpectra
    flag3=input...
        ('Plot Data? 0 for no, 1 through #of spectra for desired spectrum     ');
    if flag3==0
        flag2='n';
        flag3=flag3ref;
    end
    figure(2)
    plot(DataCell{flag3}(:,1),DataCell{flag3}(:,2))
    xlabel('Wavelength (nm)')
    ylabel('Intensity')
    
end
flag4='n';
end

%Now that we have baseline corrected data, we need to begin spectral
%comparisons to make sure that we don't have more than one data set on any
%particular sphere. We can start by first making enough duplicate data sets
%such that each one can be compared to every other spectrum. Next, we will
%match the spectral range in each comparisson spectrum.

DataGrid=zeros(length(DataCell));%Prealocating for speed
DataGrid=num2cell(DataGrid);%Converting to cell.

for i =1:length(DataCell)
    DataGrid(:,i)=DataCell;%Prepares a grid of duplicate datasets
end


for i=1:NumberOfSpectra-1
    for ii=i+1:NumberOfSpectra
        if i<ii
        N1=length(DataGrid{ii,i}(:,1));%Length of Data_iixi
        N2=length(DataGrid{i,ii}(:,1));%Length of Data_ixii
        Dat1_min=DataGrid{ii,i}(N1,1);%Min x-value of Data iixi
        Dat1_max=DataGrid{ii,i}(1,1);%Max x-value of data iixi
        Dat2_min=DataGrid{i,ii}(N2,1);%Min x-value of Data ixii
        Dat2_max=DataGrid{i,ii}(1,1);%Max x-value of Data ixii
        
        %Now we set up begining and end truncators
        if Dat1_min < Dat2_min
        Datax=DataGrid{ii,i}(:,1);
        Datay=DataGrid{ii,i}(:,2);
        min_loc=find(Datax==Dat2_min);
        Datax=Datax(1:min_loc);
        Datay=Datay(1:min_loc);
        Dataxy=horzcat(Datax,Datay);
        DataGrid{ii,i}=Dataxy;
        elseif Dat2_min < Dat1_min
        Datax=DataGrid{i,ii}(:,1);
        Datay=DataGrid{i,ii}(:,2);
        min_loc=find(Datax==Dat1_min);
        Datax=Datax(1:min_loc);
        Datay=Datay(1:min_loc);
        Dataxy=horzcat(Datax,Datay);
        DataGrid{i,ii}=Dataxy;
        end
 N1=length(DataGrid{ii,i}(:,1));%Length of Data_iixi
 N2=length(DataGrid{i,ii}(:,1));%Length of Data_ixii
 if Dat1_max > Dat2_max
     Datax=DataGrid{ii,i}(:,1);
     Datay=DataGrid{ii,i}(:,2);
     max_loc=find(Datax==Dat2_max);
     Datax=Datax(max_loc:N1);
     Datay=Datay(max_loc:N1);
     Dataxy=horzcat(Datax,Datay);
     DataGrid{ii,i}=Dataxy;
 elseif Dat2_max > Dat1_max
     Datax=DataGrid{i,ii}(:,1);
     Datay=DataGrid{i,ii}(:,2);
     max_loc=find(Datax==Dat1_max);
     Datax=Datax(max_loc:N2);
     Datay=Datay(max_loc:N2);
     Dataxy=horzcat(Datax,Datay);
     DataGrid{i,ii}=Dataxy;
        end
    end
    end
end

%Now that we have a grid of spectra that have similar wavelength ranges, we
%can begin making some comparisons. We start by calculating the Hit Quality
%Index (HQI).
HQI=zeros(NumberOfSpectra);
for i=1:NumberOfSpectra
    for ii=1:NumberOfSpectra
    HQI(i,ii)=100*(sum(DataGrid{i,ii}(:,2).*DataGrid{ii,i}(:,2)).^2)/...
        (sum(DataGrid{i,ii}(:,2).^2)*sum(DataGrid{ii,i}(:,2).^2));
    end
end
figure(3)
h=heatmap(HQI,'ColorLimits',[0 100],'ColorMap',parula)
h.Title='HQI';
xlabel('SpectrumNumber')
ylabel('SpectrumNumber')

%Now we begin to set up the machinery for peak finding and peak location
%comparisons.

Peaks=zeros(size(DataGrid));%Prealocating for speed.
Peaks=num2cell(Peaks);

for i=1:NumberOfSpectra
    for ii=1:NumberOfSpectra
        sel=(max(DataGrid{i,ii}(:,2))-min(DataGrid{i,ii}(:,2)))/4;
%sel - The amount above surrounding data for a peak to be, identified.
%Larger values mean the algorithm is more selective in finding peaks.
thresh=0.01*norm(DataGrid{i,ii}(:,2));
%thresh - A threshold value which peaks must be larger than to be
%identified.
[peakloc,peakmag]=peakfinder(DataGrid{i,ii}(:,2),sel,thresh,1,false,true);
%Since you're choosing to use quadratic interpolation to find the magnitude
%and location of the peak (Thereby not requireing consistent point spacing
%between each dataset) your output for the peak locations are fractional
%indicies. Thus the standard method of retrieving the desired wavelength by
%using the indicies won't work. We will consequently use the point slope
%form of a line to generate the wavelength positions of each peak from
%their fractional indicies. Note that this method requires a constant
%spacing between points within a given dataset.
slope=DataGrid{i,ii}(2,1)-DataGrid{i,ii}(1,1);
peakloc=DataGrid{i,ii}(1,1)+slope*(peakloc);
Peaks{i,ii}=[peakloc,peakmag]; %Gives Peak locations and magnitudes of each
%data set.
    end
end
%Now that we have peak locations for each spectrum we need to match up the
%corresponding peaks.
D=zeros(size(DataGrid));%Initialize Array
D=num2cell(D);%Convert to cell
Match=zeros(size(D));
Match=num2cell(Match);
Tolerance=0.2;%Tolerancce for assigning matching peaks


for i=1:NumberOfSpectra
for j=1:NumberOfSpectra
if isempty(Peaks{j,i}) == 0 && isempty(Peaks{i,j}) == 0
    %In previous versions of this code there is an error if either
    %Peaks{i,j} or j,i are empty. The above line checks to make sure they
    %are nonempty before proceding.
for k=1:length(Peaks{i,j}(:,1))
 
    dk=abs(Peaks{i,j}(k,1)-Peaks{j,i}(:,1));
    if k == 1
        dsub=dk;
    elseif k~=1
        dsub=horzcat(dsub,dk);
    end
    dmin=find(dsub(:,k)==min(dsub(:,k)),1,'first');
    if dsub(dmin,k) <= Tolerance
        IndexofBestMatch=dmin;
    elseif dsub(dmin,k)> Tolerance
        IndexofBestMatch=[];
    end
    if IndexofBestMatch >=1
        Matchpart=[Peaks{i,j}(k,1),Peaks{i,j}(k,2),...
            Peaks{j,i}(IndexofBestMatch,1),...
            Peaks{j,i}(IndexofBestMatch,2),...
            Peaks{j,i}(IndexofBestMatch,1)-Peaks{i,j}(k,1)];
        if i~=j
            if Match{i,j}==[0]
                Match{i,j}=Matchpart;
            elseif Match{i,j}~=[0]
                Match{i,j}=vertcat(Match{i,j},Matchpart);
            end
        elseif i==j
            Match{i,j}=Peaks{i,j};
        end
    end
 end   
 end
D{i,j}=dsub;
end
end

%Now that we have found matching peaks within a given tolerance, we want to
%make sure that each match is unique. Thus, for the cases where matches are
%not unique, we will only keep the best of the non-unique matches.
MatchUnique=Match;

for i=1:NumberOfSpectra
    for j=1:NumberOfSpectra
        if i~=j
        %Step 1: Find reapeats in Match{i,j}(:,3)
        if Match{i,j}~= 0 %Adresses bug in case of zero matches
            X=Match{i,j}(:,3);
        if length(X)>1 %Addresses bug in case of 1 or fewer matches
        uniqueX=unique(X);%Finds unique peak matches in the third
        %column of the peak matching list.
        countofX=hist(X,length(uniqueX));
        %Makes histogram with number of bins equal to the number of
        %unique peak matches.
        IndexToRepeatedValues=(countofX ~=1);
        RepeatedValues=uniqueX(IndexToRepeatedValues);
        %Step 2: Determine best matches for repeated values in
        %Match{i,j}(:,3)
        if length(RepeatedValues)>=1
        for k=1:length(RepeatedValues)
            IndexofRepeatsLeft=find(X==RepeatedValues(k));
            IndexofBetterMatch=find(abs(Match{i,j}...
                (IndexofRepeatsLeft,5))==min(abs(Match{i,j}...
                (IndexofRepeatsLeft,5))));
            IndexofRepeatsLeft(IndexofBetterMatch)=[];
            %Removes the index of teh peak that we will want to keep
            %because it is the best possible match.
        if k==1
            ForFinalRemoval=IndexofRepeatsLeft;
        else
            ForFinalRemoval=vertcat(ForFinalRemoval,IndexofRepeatsLeft);
        end
        end
        %Step 3: remove Least matching duplicates
        MatchUnique{i,j}(ForFinalRemoval,:)=[];
        end
        end
        end
        end
    end
end

%Now that the peaks are matched, we need to assess the quality of the match
%between each spectrum.
MatchQualityBeta=zeros(size(Match));%Prealocating for speed
for i=1:NumberOfSpectra
    for j=1:NumberOfSpectra
        if i~=j
            if Match{i,j} ~= 0
                MatchQualityBeta(i,j)=100*(((sum(MatchUnique{i,j}(:,1).*...
                    MatchUnique{i,j}(:,3))^2)/(sum(MatchUnique{i,j}...
                    (:,1).^2)*sum(MatchUnique{i,j}(:,3).^2)))-(1-(0.5*...
                    length(MatchUnique{i,j}(:,1))*((1/length(Peaks{i,j}...
                    (:,1)))+(1/length(Peaks{j,i}(:,1)))))));
            end
        elseif i==j
            MatchQualityBeta(i,j)=100;
        end
    end
end


figure(4)
hh=heatmap(MatchQualityBeta,'ColorLimits',[0 100],'ColorMap',parula)
hh.Title='PPI_b_e_t_a';
xlabel('SpectrumNumber')
ylabel('SpectrumNumber')


%Now that we've calculated the HQI and PPI for each spectral pair to assess
%similarities, we need to decide which spectra we want to keep. We can now
%plot several datasets together to compare them to decide which to
%ultimately keep.

flag5='y';
flag6='y';
flag7ref=1;

flag5=input('Would you like to compare spectra (y/n)','s');
while flag5=='y'

    if flag5=='y'
        flag6='y';
    end
while flag6=='y'
    NumberOfSpectra
    flag7=input('Pick data to plot. (0 for none and 1 through # of total spectra for desired spectrum     ');
    if flag7==0
        flag6='n';
        flag7=flag7ref;
        
    end
    figure(5)
    hold on
    plot(DataCell{flag7}(:,1),DataCell{flag7}(:,2))
    xlabel('Wavelength (nm)')
    ylabel('Intensity')
end
flag5=input('Would you like to view another set? (y/n)','s');
clf
end



%Now that we've made comparisons between datasets, we need to select the
%data that we want to save. 


flag8=input('Do you want to save some data? (y/n)','s');
%flag10

while flag8=='y'
    flag9=input('Which spectrum do you whish to save? (#)');
    flag9string=num2str(flag9);
    savename=[name,'_',flag9string,'BC'];
    NameAndData=vertcat([string(savename) string(savename)],...
        DataCell{flag9});
    fullFileName=fullfile(pathname, savename);
    xlswrite(fullFileName,NameAndData);
    flag8=input('Do you want to save another file? (y/n)','s');
end