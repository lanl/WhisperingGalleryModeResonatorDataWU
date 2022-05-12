clear all clc
%This script is designed for spectral clustering. The instructions for how
%to accomplish this were provided by
%[https://changyaochen.github.io/spectral-clustering/]
%First lets uploadand plot the data.
[filename, pathname]=uigetfile('*.*','Get Data');
ss=strcat(pathname,filename);
Data=xlsread(ss); %Contains data to be analyzed

FilePathFileNameFileExt=[pathname filename];
[filepath, name, Ext]=fileparts(FilePathFileNameFileExt);

Rvect=Data(:,1);%Vector containing all Radius values
nevect=Data(:,2);%Vector containing all external R.I. values

%now we Z-scale the data to account for different variance in both
%dimensions

%RvectZScaled=(Rvect-mean(Rvect))./std(Rvect);
%nevectZScaled=(nevect-mean(nevect))./std(nevect);
%DataZScaled=[RvectZScaled,nevectZScaled];

%%%Note!! If you want the data to be Z-scaled, you should delete or comment
%%%out the next three lines and uncomment out the 3 lines above this!
RvectZScaled=Rvect;
nevectZScaled=nevect;
DataZScaled=[RvectZScaled,nevectZScaled];
N=length(Data(:,1));%number of data points


figure(1)
scatter(Rvect,nevect)
xlabel('Radius (um)');
ylabel('n_e (RIU)');
%We want to map out the nearest neighbors for each point in the z-scaled
%space. To start, we calculate the Euclidian distance of each point in this
%space.
EucDistance=zeros(N);%Preinitializing for speed.
for i=1:N
    for j=1:N
        EucDistance(i,j)=norm(DataZScaled(i,:)-DataZScaled(j,:));
        %EucDistance(i,j)=sqrt((1/100)*(Data(i,1)-Data(j,1))^2+(Data(i,2)-Data(j,2))^2);
    end
end

%Next, we find the K smallest elements in each collumn of EucDistance, i.e.
%the K nearest neighbors.

K=5%We find the K+1 nearest neighbors first because each point will always be
%closest to itself so we must eliminate this trivial occurrence.
[KNNs, KNNidx]=mink(EucDistance,(K+1),1);

KNNs(1,:)=[];%Deletes the top row wich contains the trivial value of the
%Distance between each point and itself (zero).
KNNidx(1,:)=[];%Removes the coordinate of the describing the distance
%between each point and itself.

%Next, we form the adjacency matrix(A)which consists of zeros or ones where
%zeros denote no connection between points of the row and collumn indicies
%(there will always be zeros allong the diagonal since there can never be a
%point that is connected to itself) and ones represent a connection between
%points which is determined by finding the K-nearest neighbors (KNNs).

A=zeros(N);
for i=1:N
A(KNNidx(:,i),i)=1;%Adjacency matrix
%Adjacency matrix is being filled with ones according to the locations of
%the min distances found above.
end

%Now we form the degree matrix which describes how many points each point
%is connected to. Each number of connections is listed allong the diagonal,
%with zeros everywhere else.
D=0*(1:N);%Preinitializing for speed
for i=1:N
%D=diag(K*ones(N,1));
D(i)=sum(A(i,:));
end
D=diag(D);

%Next we calculate the graph Laplacian by taking the difference between D
%and A: (L = D - A). The Laplacian's diagonal consists of the degree of
%each node and the offdiagonal is the negative edge weights.

L=D-A;% The Graph Laplacian

%Now, lets generate the symmetric normalized Laplacian

L=D^(-1/2)*L*D^(-1/2);

[REigVects, EigVals]=eig(L);
EigVals=diag(EigVals);
EigVals=real(EigVals);
REigVects=real(REigVects);%Just keeps the real component. the imaginary is
%zero anyway
[EigVals, ind] = sort(EigVals);%Rearranging in increasing order
REigVects=REigVects(:,ind);%Reordering



figure(2)
clf
gplot(A,Data,'g-*')
hold on
scatter(Rvect,nevect)
xlabel('Radius (um)');
ylabel('n_e (RIU)');
hold off

figure(3)
scatter(1:N,EigVals)
xlabel('Eigen Value Index')
ylabel('Eigen Value')
xlim([0 10])

Clusts=input('How many clusters do you want to search for? (#)');
Y=REigVects(:,1:Clusts);%Keeps the number of right eigenvectors
%specified in the user input above. There will be N rows of Y (which is of
%Reals^(number of clusters) space). %We will next use the coordinates in
%each row of Y to form clusters C1, ..., Ck with the k-means algorithm.

%Next, we normalize each row of Y
for i=1:N
Y(i,:)=Y(i,:)./norm(Y(i,:));
end
rng('shuffle')%Initiallizes random number generator
Mean_EigSpace=zeros(Clusts);
Ycopy=Y;
for i=1:Clusts
    RNG=randi([1,N-(i-1)]);%Stores random integer index for selecting
    %an initial cluster mean.
    Mean_EigSpace(i,:)=Ycopy(RNG,:);
    Ycopy(RNG,:)=[];%Removes data point just used as a mean so it is not
    %selected again.
end

%Now that the initial means are selected we need to iteratively calculate
%the distances between these means and each point then assign each point to
%a cluster and recalculate the mean and repeat.
Distance=zeros(size(Rvect));%Preinitiallizing for speed.
Mean_originalSpace=zeros(Clusts, 2);

for t=1:20
    ClusterAssignmentStorage=zeros(size(Y));
for i=1:N
    for k=1:Clusts
    Distance(i,k)=norm(Y(i,:)-Mean_EigSpace(k,:));
    end
    [MinVal, minidx]=min(Distance(i,:));
    ClusterAssignmentStorage(i,minidx)=1;
end
for j=1:Clusts
    Mean_EigSpace(j, :)=ClusterAssignmentStorage(:,j)'*Y./sum(ClusterAssignmentStorage(:,j));
    Mean_originalSpace(j,:)=ClusterAssignmentStorage(:,j)'*Data./sum(ClusterAssignmentStorage(:,j));
end
figure(4)
clf
hold on
for k=1:Clusts
scatter(Rvect(find(ClusterAssignmentStorage(:,k) ~=0)),nevect(find(ClusterAssignmentStorage(:,k) ~=0)))
scatter(Mean_originalSpace(:,1),Mean_originalSpace(:,2),'x')

end
hold off
end
Mean_originalSpace