% Analysis script to accompany: Santosh?Ansumali,?Aloke?Kumar,?Samarth?Agrawal,?H J?Shashank,?Meher K?Prakash, A steady trickle-down from metro districts and improving epidemic-parameters characterize the increasing COVID-19 cases in India, medRxiv (2020)
%doi:?https://doi.org/10.1101/2020.09.28.20202978

clear
% Create labels for the different states being analyzed
States={'Others','Kerala','Karnataka', 'Maharashtra','Gujarat','Rajasthan','Punjab','Chandigarh','Haryana','Delhi','Uttarakhand','Uttar Pradesh','Bihar','Jharkhand','West Bengal','Assam','Tripura','Odisha','Andhra Pradesh','Tamil Nadu','Puducherry','Others','Chattisgarh','Madhya Pradesh'};
States={'OT','KL','KA', 'MH','GJ','RJ','PB','CH','HR','DE','UK','UP','BR','JH','WB','AS','TR','OR','AP','TN','PY','OT','CT','MP'};

% Read data on traffic of trains, domestic and international flights to
% each district

Temp=csvread('Train_Flight_Traffic_Districts.csv',1,3);

state=zeros(length(Temp),1);
access=zeros(length(Temp),1);
group=zeros(length(Temp),1);
transitions=zeros(length(Temp),1);
traffic=zeros(length(Temp),1);

% Define a hierarchical level of access based on the train/flight traffic.
% Districts with international flight traffic belong to Group 1, with subgroups 1 to 3;
% Districts with domestic flight access belong to Group 2, with subgroups 4 to 6; and so on

for i=1:length(Temp)
  state(i) = 25-Temp(i,5);
  if Temp(i,3)>25000
      col = 1;
      group(i)=1;
  elseif Temp(i,3)>5000
      col = 2;      
      group(i)=1;
  elseif Temp(i,3)>0
      col = 3;
      group(i)=1;
  elseif Temp(i,4)>5000 
      col = 4;
      group(i)=2;
  elseif Temp(i,4)>0
      col = 5;
      group(i)=2;
  elseif Temp(i,2)>50
      col = 6;
      group(i)=3;
  elseif Temp(i,2)>20
      col = 7;
      group(i)=3;
  else
      col = 8;
      group(i)=4;
  end
  access(i)=col;
  traffic(i)=Temp(i,4)+Temp(i,2);
end

% Note how many districts belong to each group (international access,
% domestic flight access, train access
groupOccupancy = [length(group(group==1)) length(group(group==2)) length(group(group==3)) length(group(group==4))];

% Read data of population, infection, recovery rates

AnalysisFolder = 'Rates_Districts_MayToAugust/';
cd (AnalysisFolder)
% sanity check: do we have an output file for every settings file we created?
MasterDirectory=dir;
for f=1:length(MasterDirectory)
    AllFiles{f} = MasterDirectory(f).name;
end

i=1;
for f=1:length(MasterDirectory)
            if strfind(AllFiles{f},'.csv')
              %fprintf('%s\n',AllFiles{f})
              ValidFiles{i} = AllFiles{f};
              i=i+1;
        end
end

avgR=zeros(4,length(ValidFiles));

for i=1:4
    for j=1:8
        avgParams{i,j}={};
    end
end

clear i j
currentBurden=zeros(length(ValidFiles),1);

%%
% Perform the analysis for every fortnight from the files in the folder Rates_Districts_MayToAugust/
for Fortnight=1:length(ValidFiles)
num=0;den=0;newTemp=0;
data=zeros(24,8);
clear DistTemp;
currentBurdenMatrix=[];
DistTemp=csvread(ValidFiles{Fortnight},1,2);
Dist(Fortnight,:,:)=DistTemp;
for i=1:length(Temp)
%Inter = Temp(any(Temp(:,3) ~= 0,2),:);
  num = num + access(i)* Dist(Fortnight,i,2);
  den = den + Dist(Fortnight,i,2);
  currentBurdenMatrix(length(currentBurdenMatrix)+1,1:2)=[access(i) Dist(Fortnight,i,2)];

  avgR(group(i),Fortnight) = avgR(group(i),Fortnight)+Dist(Fortnight,i,8);
  avgParams{group(i),Fortnight}= cat(1,avgParams{group(i),Fortnight},Dist(Fortnight,i,5));
  checkCond= Dist(Fortnight,i,2); 
  %checkCond = Dist(Fortnight,i,2)/Dist(Fortnight,i,1)*1e4;
  
      %if and(checkCond> 0.1,transitions(i)==0)
      if and(checkCond > 8000,transitions(i)==0)
          transitions(i)=Fortnight+0.2*rand;
      end
      
      data(state(i),access(i))=data(state(i),access(i)) + round(Dist(Fortnight,i,2));

end

size(currentBurdenMatrix)
meanCurrentBurden(Fortnight) = currentBurdenMatrix(:,1)'*currentBurdenMatrix(:,2)/sum(currentBurdenMatrix(:,2));
for st=1:length(currentBurdenMatrix)
newTemp=newTemp+currentBurdenMatrix(st,1)^2*currentBurdenMatrix(st,2)/sum(currentBurdenMatrix(:,2)) ;
end
stdCurrentBurden(Fortnight) = (newTemp- meanCurrentBurden(Fortnight)^2)^0.5 ;
currentBurden(Fortnight) = num/den;


% Generate bullseye plots
bullseye(data,States); colorbar; caxis([0 20000]);
print(['Fig' num2str(Fortnight)],'-depsc')

for i=1:4
avgR(i,Fortnight)=avgR(i,Fortnight)/groupOccupancy(i); 
end
%pause
fprintf('Fortnight %d length %d\n', Fortnight, length(avgParams{1,Fortnight}))
end

%plot for currentBurden Distribution
errorbar(currentBurden,stdCurrentBurden/sqrt(188),'ro'); hold on; errorbar(currentBurden,stdCurrentBurden/sqrt(188)); hold off
axis([0 9 1 8])
%%
%plot for the access vs. switch over the high transitions
transNew=transitions(transitions>0);
accessNew=access(transitions>0);
for it=1:8
    nTrans(it)=length(transNew(accessNew==it));
    meanTrans(it)=mean(transNew(accessNew==it));
    stdTrans(it)=std(transNew(accessNew==it))/nTrans(it);
end
errorbar(meanTrans(:),stdTrans(:));axis([0 9 2 9])

%%
% Plots and trend analysis for the different rates, over the 4 months
    clear xtemp ytemp;
    xtemp1=[];xtemp2=[];xtemp3=[];xtemp4=[];
    ytemp1=[];ytemp2=[];ytemp3=[];ytemp4=[];

for i=1:8
    length(xtemp1)
for j=1:4
        for k=1:length(avgParams{j,i})
            arr(k)=avgParams{j,i}{k};
        end
            %arr(arr>5)=0;
            meanParam(j,i)=mean(arr);
            stdParam(j,i) =std(arr);

            
 if j==1
     fprintf('j %d arr %d',j,length(arr))
            xtemp1((length(xtemp1)+1):(length(xtemp1)+length(arr)),1:2) = [i*ones(length(arr),1),arr'];
 elseif j==2
     fprintf('j %d arr %d',j,length(arr))
            xtemp2((length(xtemp2)+1):(length(xtemp2)+length(arr)),1:2) = [i*ones(length(arr),1),arr'];
 elseif j==3
     fprintf('j %d arr %d',j,length(arr))
            xtemp3((length(xtemp3)+1):(length(xtemp3)+length(arr)),1:2) = [i*ones(length(arr),1),arr'];
 else
     
     fprintf('j %d arr %d',j,length(arr))
             xtemp4((length(xtemp4)+1):(length(xtemp4)+length(arr)),1:2) = [i*ones(length(arr),1),arr'];
 end
%             for xl=1:length(xtemp2)
%                 if (xtemp2(xl,2) >5) 
%                     xtemp2(xl,2)=0;
%                 end
%             end
        clear arr
    end
end

errorbar(meanParam(1,:),1.96*stdParam(1,:)/sqrt(groupOccupancy(1)),'r');hold on;
errorbar(meanParam(2,:),1.96*stdParam(2,:)/sqrt(groupOccupancy(2)),'b');hold on;
errorbar(meanParam(3,:),1.96*stdParam(3,:)/sqrt(groupOccupancy(3)),'g');hold on;
errorbar(meanParam(4,:),1.96*stdParam(4,:)/sqrt(groupOccupancy(4)),'k');hold off;
%axis([0 8.5 0.5 2.5])

%%
% trend analysis results
stats=regstats(xtemp1(:,1),xtemp1(:,2),'linear');
stats.tstat.pval
plot(xtemp1(:,1),xtemp1(:,2),'*');
pause

stats=regstats(xtemp2(:,1),xtemp2(:,2),'linear');
stats.tstat.pval
plot(xtemp2(:,1),xtemp2(:,2),'*');
pause

stats=regstats(xtemp3(:,1),xtemp3(:,2),'linear');
stats.tstat.pval
plot(xtemp3(:,1),xtemp3(:,2),'*');
pause

stats=regstats(xtemp4(:,1),xtemp4(:,2),'linear');
stats.tstat.pval
plot(xtemp4(:,1),xtemp4(:,2),'*');
pause
%%
% Centre of gravity of current infections
plot(currentBurden)

%compositeTransitions=[access transitions];
cTx = access(transitions>0);
cTy = transitions(transitions>0);
for i=1:Fortnight
mstdx(i)=i;
mstd(i)=mean(cTy(round(cTx)==i));
end
plot(cTx,cTy,'*',mstdx,mstd,'rs')