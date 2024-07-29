%Sparse peptide MSMS

%Input parent peptide parameters
sequence = 'TISFKDDGTYK';
pepLength = length(sequence);

%Labeled AA
lab = {'V'; 'I'; 'L'};

%LAbeled AA rules
pm = [1, 5; 
      1, 0;
      1, 0];

precursorCharge = 2;

ppm = 2; %ppm error tolerance for search

%mode for final plotting. Centroid plots MS data as stick spectra, profile
%plots MS data as full profile spectra
%mode = 'centroid';
mode = 'profile';

%% IMPORT MS DATA

%If there is already data uploaded, open a box and ask if the user would
%really like to upload new data, or proceed with the stuff that is already
%open

if exist('Data', 'var')

    quest = 'There is already data uploaded, would you like to upload new data, or proceed with the existing data?';
    answer = questdlg(quest, 'Data upload warning', 'Contintue with existing data', 'Upload new data', 'Contintue with existing data');

end

if ~exist('Data', 'var') || strcmp(answer, 'Upload new data')

    %input mass spec data
    [filename, path1 ] = uigetfile('.asc', 'import asc format mass spectrum', 'MultiSelect', 'off');
    filename2 = strcat(path1, filename);
        
    if iscell(filename2) == 1
        filename = filename';
        filename2 = filename2';
    end
    
    %split the filename from the file type to make the name in the data
    %structure
    file = cell(size(filename,1),1);
    
    for i = 1:size(filename,1)
    
        [path1, name, ext] = fileparts(filename);
        file(i,1) = {name};
    
    end
    
    %import each file into the structure
    for i = 1:length(file)
        
        Data.(file{i,1}) = importdata(filename2);
    
        if size(Data.(file{i,1}),2) == 3
            %remove the third column cause we don't need it
            Data.(file{i,1})(:,3) = [];
        end
    end
    
    % if "fullProfile" is selected, make the user input an .xy spectrum to make
    % prettier plots to compare the data to the simulations
    fullProfile = 1;
    if fullProfile == 1
    
        [filename, path1 ] = uigetfile('.xy', 'import full profile .xy spectrum', 'MultiSelect', 'off');
        filename2 = strcat(path1, filename);
        
        if iscell(filename2) == 1
            filename = filename';
            filename2 = filename2';
        end
        
            for i = 1:length(file)
        
                Data.profileData = importdata(filename2);
        
            end
    end
end

 clear ext quest answer

%% Generate all N terminal and C terminal sequences
Nterm = cell(length(sequence) - 1,1);
Cterm = cell(length(sequence) - 1,1);

for i = 1:pepLength-1

    Nterm(i,1) = {sequence(1:i)};
    Cterm(i,1) = {sequence(i+1:end)};

end

B = struct();
Y = struct();

%calculate formulas for all B and Y ions (singly charged and protonated)
for i = 1:length(Nterm)

    Bcleavage = num2str(i);
    Ycleavage = num2str(pepLength-i);

    Bion = strcat('B', Bcleavage);
    Yion = strcat('Y', Ycleavage);

    B.(Bion) = atomiccomp(Nterm{i,1});
    B.(Bion).O = B.(Bion).O - 1;
    B.(Bion).H = B.(Bion).H - 1;

    Y.(Yion) = atomiccomp(Cterm{i,1});
    Y.(Yion).H = Y.(Yion).H + 1; 

end

%Generate a list of all B ion monoisotopic masses and all charge states
[Blist, Bheaders] = getBmasses(Nterm, B, precursorCharge);
[Ylist, Yheaders] = getYmasses(Cterm, Y, precursorCharge);


%% COMPARE TO DATA

%Copy the B and Y ion lists and add some extra columns for ppm error and
%obs mass

lablHeaders = {'obs mass', 'ppm error', lab{:,1}, '15N', '13C', 'lab m/z', 'obs m/z', 'ppm error'};
Bheaders = [Bheaders, lablHeaders];
Yheaders = [Yheaders, lablHeaders];

zpad = zeros(size(Blist,1), size(lablHeaders, 2));
Blist2 = [Blist, num2cell(zpad)];
Ylist2 = [Ylist, num2cell(zpad)];
clear zpad

%Count labels for B ions
for i = 1:size(Blist2,1)

    AAcount = aacount(Blist2{i,1});
    
    %count each labeled AA in lab
    for j = 1:length(lab)
        %count number of AA in current fragment
        Blist2(i,j+6) = num2cell(AAcount.(lab{j,1}));
        %count 15N
        Blist2(i, length(lab)+7) = num2cell(Blist2{i, length(lab)+7} + AAcount.(lab{j,1})*pm(j,1));
        %count 13C
        Blist2(i, length(lab)+8) = num2cell(Blist2{i, length(lab)+8} + AAcount.(lab{j,1})*pm(j,2));

    end


    if (Blist2{i,10} ~=0) || (Blist2{i,11} ~=0) 
        
        f = B.(Blist2{i,2}); %get formula
        z = Blist2{i,4}; %get charge
        c = Blist2{i,11}; %get 13C labels
        n = Blist2{i,10}; %get 15N labels

        f.H = f.H + z-1; %add proton proportional to charge-1 (alreadt M+H)
        [MD, Info] = isotopicdist2(f, 'Csparse', c, 'Nsparse', n);
        Blist2{i,12} = Info.MonoisotopicMass/z; %scale to charge
    else
        Blist2{i,12} = Blist{i,3};
    end

end
clear MD Info f z c n 

%Count labels for Y ions
%Y ions
for i = 1:size(Ylist2,1)

    AAcount = aacount(Ylist2{i,1});
    
    %for each labeled AA in lab
    for j = 1:length(lab)
        %count number of AA in current fragment
        Ylist2(i,j+6) = num2cell(AAcount.(lab{j,1}));
        %count 15N
        Ylist2(i, length(lab)+7) = num2cell(Ylist2{i, length(lab)+7} + AAcount.(lab{j,1})*pm(j,1));
        %count 13C
        Ylist2(i, length(lab)+8) = num2cell(Ylist2{i, length(lab)+8} + AAcount.(lab{j,1})*pm(j,2));

    end

        if (Ylist2{i,10} ~=0) || (Ylist2{i,11} ~=0)  %iff there are labels, calc the labeled mass
        
            f = Y.(Ylist2{i,2});
            z = Ylist2{i,4};
            c = Ylist2{i,11};
            n = Ylist2{i,10};

            f.H = f.H + z - 1;
            [MD, Info] = isotopicdist2(f, 'Csparse', c, 'Nsparse', n );
            Ylist2{i,12} = Info.MonoisotopicMass/z;

        else %if there are no labels, just scoot the unlabeled mass over
            
            Ylist2{i,12} = Ylist2{i,3};

        end
end

%search for B ions
for i = 1:size(Blist,1)

    %get the current labeled m/z
    mass = Blist2{i,3}; %unlabeled mass
    mass2 = Blist2{i,12}; %labeled mass

    %truncate the m/z list to +-5 Da of mass
    mz = Data.(file{1,1})(:,1);
    mz(mz < mass - 10) = [];
    mz(mz > mass + 10) = [];

    %calc the element wise ppm difference between mass and m/z
    if ~isempty(mz)

        %unlabeled
        p = 1000000*(mass - mz)./mass;
        pabs = abs(p);
        
        %find nearest neighbor mass
        Blist2(i,5) = num2cell(mz(pabs == min(pabs),1));

        %put ppm error in table
        Blist2(i,6) = num2cell(p(pabs == min(pabs),1));


        %labeled
        p = 1000000*(mass2 - mz)./mass2;
        pabs = abs(p);
        
        %find nearest neighbor mass
        Blist2(i,13) = num2cell(mz(pabs == min(pabs),1));

        %put ppm error in table
        Blist2(i,14) = num2cell(p(pabs == min(pabs),1));

    elseif isempty(mz)
        
        %put an arbitrary large ppm in the ppm column
        Blist2(i,6) = num2cell(1000000);
        Blist2(i,14) = num2cell(1000000);

    end

end

%search for Y ions
for i = 1:size(Ylist,1)

    %get the current labeled m/z
    mass = Ylist2{i,3};
    mass2 = Ylist2{i,12}; %labeled mass

    %truncate the m/z list to +-5 Da of mass
    mz = Data.(file{1,1})(:,1);
    mz(mz < mass2 - 10) = [];
    mz(mz > mass2 + 10) = [];

    %calc the element wise ppm difference between mass and m/z
    if ~isempty(mz)

        %unlabeled
        p = 1000000*(mass - mz)./mass;
        pabs = abs(p);
        
        %find nearest neighbor mass
        Ylist2(i,5) = num2cell(mz(pabs == min(pabs),1));

        %put ppm error in table
        Ylist2(i,6) = num2cell(p(pabs == min(pabs),1));

        %labeled
        p = 1000000*(mass2 - mz)./mass2;
        pabs = abs(p);
        
        %find nearest neighbor mass
        Ylist2(i,13) = num2cell(mz(pabs == min(pabs),1));

        %put ppm error in table
        Ylist2(i,14) = num2cell(p(pabs == min(pabs),1));

    elseif isempty(mz)
        
        %put an arbitrary large ppm in the ppm column
        Ylist2(i,6) = num2cell(1000000);
        Ylist2(i,14) = num2cell(1000000);

    end

end

%B IONS: mark rows that did not have a match for deletion
for i = 1:size(Blist2,1)

    if (Blist2{i,6} == 1000000) && (Blist2{i,end} == 1000000) %if neither unlabeled nor labeled peaks were found, delete that row

        Blist2{i,end} = 6969; %rows marked with 6969 will be deleted

    elseif (Blist2{i,6} == 1000000) && (abs(Blist2{i,end}) > ppm)

        Blist2{i,end} = 6969; %rows marked with 6969 will be deleted

    elseif (abs(Blist2{i,6}) > ppm) && (Blist2{i,end} == 1000000)

        Blist2{i,end} = 6969; %rows marked with 6969 will be deleted

    elseif (abs(Blist2{i,6}) > ppm) && (abs(Blist2{i,end}) > ppm)
        
         Blist2{i,end} = 6969; %rows marked with 6969 will be deleted

    end

end

%Y IONS: mark rows with no matches for deletion.
for i = 1:size(Ylist2,1)

    if (Ylist2{i,6} == 1000000) && (Ylist2{i,end} == 1000000) %if neither unlabeled nor labeled peaks were found, delete that row

        Ylist2{i,end} = 6969; %rows marked with 6969 will be deleted

    elseif (Ylist2{i,6} == 1000000) && (abs(Ylist2{i,end}) > ppm)

        Ylist2{i,end} = 6969; %rows marked with 6969 will be deleted

    elseif (abs(Ylist2{i,6}) > ppm) && (Ylist2{i,end} == 1000000)

        Ylist2{i,end} = 6969; %rows marked with 6969 will be deleted

    elseif (abs(Ylist2{i,6}) > ppm) && (abs(Ylist2{i,end}) > ppm)
        
         Ylist2{i,end} = 6969; %rows marked with 6969 will be deleted

    end

end


%remove B ions that do not meet ppm criteria
Blist2(cell2mat(Blist2(:,end)) == 6969, :) = [];

%remove Y ions that do not meet ppm criteria
Ylist2(cell2mat(Ylist2(:,end)) == 6969, :) = [];



% lablHeaders = {lab{:,1}, '15N', '13C', 'lab m/z', 'ppm error'};
% zpad = zeros(size(Blist2,1), length(lablHeaders));
% Blist3 = [Blist2 num2cell(zpad)];
% Bheaders = [Bheaders, lablHeaders];
% 
% zpad = zeros(size(Ylist2,1), length(lablHeaders));
% Ylist2 = [Ylist2 num2cell(zpad)];
% Yheaders = [Yheaders, lablHeaders];


clear MD Info f z c n 

clear AAcount i j lablHeaders

%% PLOT


figCount = 1;
names = fieldnames(Data);



%Generate figures for all B ions
for i = 1:size(Blist2,1)

mass = Blist2{i,3};
centroidData = Data.(names{1,1});
profileData = Data.profileData;

%normalize Y axis to 100
centroidData(:,2) = 100*centroidData(:,2)./max(centroidData(:,2));
profileData(:,2) = 100*profileData(:,2)./max(profileData(:,2));

%truncate data to within +/- 10 m/z of current mass
centroidData(centroidData(:,1) > mass + 10, :) = [];
centroidData(centroidData(:,1) < mass - 10, :) = [];
profileData(profileData(:,1) > mass + 10, :) = [];
profileData(profileData(:,1) < mass - 10, :) = [];

%rescale profile data to centroidData
scale = max(centroidData(:,2))/max(profileData(:,2));
profileData(:,2) = scale*profileData(:,2);

%Get the formula for the current ion
formula = B.(Blist2{i,2});
N = Blist2{i, end-4}; %get number of 15N labels
C = Blist2{i, end-3}; %get number of 13C labels
charge = Blist2{i,4};

%add protons proportional to charge (formula should already be singly
%charged)
formula.H = formula.H + charge - 1;

%generate unlabeled distribution
MD = isotopicdist(formula);

%generate labeled distribution (if there are labels)
if N ~= 0 && C == 0
    MD2 = isotopicdist2(formula, 'Nsparse', N);
elseif N == 0 && C ~= 0
    MD2 = isotopicdist2(formula, 'Csparse', C);
elseif N ~= 0 && C ~= 0
    MD2 = isotopicdist2(formula, 'Csparse', C, 'Nsparse', N);
end

%scale both distributions to their respective charges
MD(:,1) = MD(:,1)./charge;

if exist('MD2', 'var')
    MD2(:,1) = MD2(:,1)./charge;
end

%scale MD and MD2 to their nearest neighbor peak
m = MD(1,1); %pull out the monoisotopic peak
p = abs(1000000*(profileData(:,1)-m)./m); %ppm error calc
int = profileData(p == min(p),2); %extract the intensity for the corresponding peak
scale = int/MD(1,2);
MD(:,2) = MD(:,2)*scale;

if exist('MD2', 'var')
    m = MD2(1,1); %pull out the monoisotopic peak
    %p = abs(1000000*(profileData(:,1)-m)./m); %ppm error calc
    p = abs(1000000*(centroidData(:,1)-m)./m); %ppm error calc
    %int = profileData(p == min(p),2); %extract the intensity for the corresponding peak
    int = centroidData(p == min(p),2); %extract the intensity for the corresponding peak
    scale = int/MD2(1,2);
    MD2(:,2) = MD2(:,2)*scale;
end

%generate a title for the plot with the peptide, ion sequence, ion ID and
%charge
formatSpec = 'Peptide: %s\n Fragment ion: %s, ion sequence: %s\n charge: %d';
t = sprintf(formatSpec, sequence, Blist2{i,2}, Blist2{i,1}, Blist2{i,4});

%PLOT
figure(figCount)
if strcmp(mode, 'profile')
    plot(profileData(:,1), profileData(:,2))
elseif strcmp(mode, 'centroid')
    stem(centroidData(:,1), centroidData(:,2))
end
hold on
scatter(MD(:,1), MD(:,2), 'filled')
if exist('MD2', 'var')
    scatter(MD2(:,1), MD2(:,2), 'filled')
end
title(t)
xlabel('m/z')
ylabel('Intensity')
if exist('MD2', 'var')
    legend('observed', 'sim. unlabeled', 'sim. labeled')
else
    legend('observed', 'sim. unlabeled')
end
figCount = figCount + 1;

clear MD MD2 scale p m int
end


%Repeat for Y fragments
for i = 1:size(Ylist2,1);

mass = Ylist2{i,3};
centroidData = Data.(names{1,1});
profileData = Data.profileData;

%normalize Y axis to 100
centroidData(:,2) = 100*centroidData(:,2)./max(centroidData(:,2));
profileData(:,2) = 100*profileData(:,2)./max(profileData(:,2));

%truncate data to within +/- 10 m/z of current mass
centroidData(centroidData(:,1) > mass + 10, :) = [];
centroidData(centroidData(:,1) < mass - 10, :) = [];
profileData(profileData(:,1) > mass + 10, :) = [];
profileData(profileData(:,1) < mass - 10, :) = [];

%Get the formula for the current ion
formula = Y.(Ylist2{i,2});
N = Ylist2{i, end-4}; %get number of 15N labels
C = Ylist2{i, end-3}; %get number of 13C labels
charge = Ylist2{i,4};

%add protons proportional to charge (formula should already be singly
%charged)
formula.H = formula.H + charge - 1;

%generate unlabeled distribution
MD = isotopicdist(formula);

%generate labeled distribution (if there are labels)
if N ~= 0 && C == 0
    MD2 = isotopicdist2(formula, 'Nsparse', N);
elseif N == 0 && C ~= 0
    MD2 = isotopicdist2(formula, 'Csparse', C);
elseif N ~= 0 && C ~= 0
    MD2 = isotopicdist2(formula, 'Csparse', C, 'Nsparse', N);
end

%scale both distributions to their respective charges
MD(:,1) = MD(:,1)./charge;

if exist('MD2', 'var')
    MD2(:,1) = MD2(:,1)./charge;
end

%scale MD and MD2 to their nearest neighbor peak
m = MD(1,1); %pull out the monoisotopic peak
p = abs(1000000*(profileData(:,1)-m)./m); %ppm error calc
int = profileData(p == min(p),2); %extract the intensity for the corresponding peak
scale = int/MD(1,2);
MD(:,2) = MD(:,2)*scale;

if exist('MD2', 'var')
    m = MD2(1,1); %pull out the monoisotopic peak
    p = abs(1000000*(profileData(:,1)-m)./m); %ppm error calc
    int = profileData(p == min(p),2); %extract the intensity for the corresponding peak
    scale = int/MD2(1,2);
    MD2(:,2) = MD2(:,2)*scale;
end

%generate a title for the plot with the peptide, ion sequence, ion ID and
%charge
formatSpec = 'Peptide: %s\n Fragment ion: %s, ion sequence: %s\n charge: %d';
t = sprintf(formatSpec, sequence, Ylist2{i,2}, Ylist2{i,1}, Ylist2{i,4});

%PLOT
figure(figCount)
if strcmp(mode, 'profile')
    plot(profileData(:,1), profileData(:,2))
elseif strcmp(mode, 'centroid')
    stem(centroidData(:,1), centroidData(:,2))
end
hold on
scatter(MD(:,1), MD(:,2), 'filled')
if exist('MD2', 'var')
    scatter(MD2(:,1), MD2(:,2), 'filled')
end
title(t)
xlabel('m/z')
ylabel('Intensity')
if exist('MD2', 'var')
    legend('observed', 'sim. unlabeled', 'sim. labeled')
else
    legend('observed', 'sim. unlabeled')
end
figCount = figCount + 1;

clear MD MD2 scale p m int
end


clear N C mz mass i pabs names t Ycleavage Bcleavage

%Generate output table 


%% FUNCTIONS
%FUNCTIONS
%Generate B fragment masses
%Input C terminal sequences and structure with B fragment formulas

function [Blist, Bheaders] = getBmasses(Nterm, B, maxCharge)
    f = fieldnames(B);
    fc = length(f)*maxCharge;
    
    Bheaders = {'Sequence', 'Ion', 'theo m/z', 'charge'};
    if maxCharge == 1
        Blist = cell(fc,4);
        Blist(:,1) = Nterm(:,1);
        Blist(:,2) = f(:,1);
    elseif maxCharge > 1
        Blist = cell(fc,4);
        j = 1;
        charge = 1;
        for i = 1:length(Blist)
            
            if (charge > maxCharge)
                charge = 1;
                j = j + 1;
            end
            
            Blist(i,1) = Nterm(j,1);
            Blist(i,2) = f(j,1);
            Blist(i,4) = {charge};
            charge = charge + 1;
    
        end
    end
    
    charge = 1;
    j = 1;
    
    for i = 1:fc
    
        if charge > maxCharge
            charge = 1;
            j = j+1;
        end
    
        h = B.(f{j}).H; %initial proton count for this fragment (should be singly charged)
    
        if charge > 1
            B.(f{j}).H = h + charge - 1; 
        end
        
        %simulate isotope distribution for monoisotopic mass
        MD = isotopicdist(B.(f{j}));
    
        %divide by charge
        MD(:,1) = MD(:,1)./charge;
    
        Blist(i,3) = {MD(1,1)};
        
        charge = charge + 1;

        %reset protons to original value
        B.(f{j}).H = h;
    end



    %calculate all labeled masses


end

%Generate Y fragment masses
function [Ylist, Yheaders] = getYmasses(Cterm, Y, maxCharge)
    f = fieldnames(Y);
    fc = length(f)*maxCharge;
    
    Yheaders = {'Sequence', 'Ion', 'theo m/z', 'charge'};
    if maxCharge == 1
        Ylist = cell(fc,4);
        Ylist(:,1) = Cterm(:,1);
        Ylist(:,2) = f(:,1);
    elseif maxCharge > 1
        Ylist = cell(fc,4);
        j = 1;
        charge = 1;
        for i = 1:length(Ylist)
            
            if (charge > maxCharge)
                charge = 1;
                j = j + 1;
            end
            
            Ylist(i,1) = Cterm(j,1);
            Ylist(i,2) = f(j,1);
            Ylist(i,4) = {charge};
            charge = charge + 1;
    
        end
    end
    
    charge = 1;
    j = 1;
    
    for i = 1:fc
    
        if charge > maxCharge
            charge = 1;
            j = j+1;
        end
    
        h = Y.(f{j}).H; %initial proton count for this fragment (should be singly charged)
    
        if charge > 1
            Y.(f{j}).H = h + charge - 1; 
        end
        
        %simulate isotope distribution for monoisotopic mass
        MD = isotopicdist(Y.(f{j}));
    
        %divide by charge
        MD(:,1) = MD(:,1)./charge;
    
        Ylist(i,3) = {MD(1,1)};
        
        charge = charge + 1;

        %reset protons to original count
        Y.(f{j}).H = h;
    end

end

%Optional generate C fragments

%optional generate Z fragments