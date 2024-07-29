%Calculate the mass of sparse labeled peptides

%This code is built to analyze a MALDI peptide spectrum where the protein
%has been labeled with sparse labeles. In the first section, input a list
%of peptides to be assesed, and which amino acids are labeled in the
%variable "lab" 

dN = 0.997034898; %15N - 14N
dC = 1.003354835; %13C - 12C

%Input peptides
sequence = {'FEGDTLVNR'};
%sequence = {'SAMPEGYVQER'};
%sequence = {'LEYNFNSHNVYITADK'};
%sequence = {'AVVFLEPQWYR'};
%sequence = {'FEGDTLVNR'};

% sequence = {'RHDFFK';
% 'NGIKANFK';
% 'YPDHMKR';
% 'DDGTYKTR';
% 'EDGNILGHK';
% 'GEGEGDATNGK';
% 'FEGDTLVNR';
% 'IEWHEMSK';
% 'IELKGIDFK';
% 'LTLKFICTTGK';
% 'SAMPEGYVQER';
% 'TISFKDDGTYK';
% 'AEVKFEGDTLVNR';
% 'GEGEGDATNGKLTLK';
% 'FSVRGEGEGDATNGK';
% 'FEGDTLVNRIELK';
% 'GIDFKEDGNILGHK';
% 'SAMPEGYVQERTISFK';
% 'LEYNFNSHNVYITADK';
% 'HDFFKSAMPEGYVQER';
% 'LEYNFNSHNVYITADKQK';
% 'DHMVLLEFVTAAGITHGEFR';
% 'RDHMVLLEFVTAAGITHGEFR';
% 'GEELFTGVVPILVELDGDVNGHK';
% 'MHHHHHHHHMSGLNDIFEAQK';
% 'EDGNILGHKLEYNFNSHNVYITADK';
% 'GEELFTGVVPILVELDGDVNGHKFSVR';
% 'IEWHEMSKGEELFTGVVPILVELDGDVNGHK';
% 'MSGLNDIFEAQKIEWHEMSK';
% 'DWKDHK';
% 'SSTRDWK';
% 'VTYLQNGK';
% 'DSGSYFCR';
% 'TNIRSSTR';
% 'EEDPIHLR';
% 'DHKFKWR';
% 'GMRTEDLPK';
% 'WRKDPQDK';
% 'DWKDHKFK';
% 'VLEKDSVTLK';
% 'VTYLQNGKGR';
% 'SSTRDWKDHK';
% 'VTYLQNGKGRK';
% 'CHSWKNTALHK';
% 'ATLKDSGSYFCR';
% 'TNIRSSTRDWK';
% 'AVVFLEPQWYR';
% 'DSGSYFCRGLFGSK';
% 'YFHHNSDFYIPK';
% 'WVFKEEDPIHLR';
% 'NTALHKVTYLQNGK';
% 'EEDPIHLRCHSWK';
% 'KYFHHNSDFYIPK';
% 'NTALHKVTYLQNGKGR';
% 'AVVFLEPQWYRVLEK';
% 'GRKYFHHNSDFYIPK';
% 'ATLKDSGSYFCRGLFGSK';
% 'YFHHNSDFYIPKATLK';
% 'TEDLPKAVVFLEPQWYR';
% 'KYFHHNSDFYIPKATLK';
% 'WVFKEEDPIHLRCHSWK';
% 'CHSWKNTALHKVTYLQNGK';
% 'EEDPIHLRCHSWKNTALHK';
% 'GMRTEDLPKAVVFLEPQWYR';
% 'AVVFLEPQWYRVLEKDSVTLK';
% 'TEDLPKAVVFLEPQWYRVLEK';
% 'YFHHNSDFYIPKATLKDSGSYFCR';
% 'CQTNLSTLSDPVQLEVHIGWLLLQAPR';
% 'CQTNLSTLSDPVQLEVHIGWLLLQAPRWVFK'};


%turn the incorporation optimization on or off. Recommend starting with
%this off and then pare down to good looking peptides then turn back on.
optimize = 0;
scrmbl = 1; %set to 1 to assess possible scrambling, set to 0 to turn off

%Input labeled AA
%lab = {'V'; 'I'; 'L'};
lab = {'K'; 'G'; 'S'};

%Input label parameters
%Column 1 = number of 15N labels
%Column 2 = number of 13C labels
%VIL parameters
% pm = [1, 5; 
%       1, 0;
%       1, 0];

%KGS parameters
pm = [1, 0; 
      1, 2;
      1, 0];

%incriments that N and C can be removed via scrambling
%inc = [1, 5]; %VIL
inc = [1, 2]; %KGS

%ppm tolerance to generate a figure later
ppm = 5;

maxCharge = 1;

Icutoff = 0.08; %reject any matches below this relative intensity 




%PLOT MODE, choose "profile" or "centroid"
mode = "profile";



%% IMPORT MASS SPEC DATA (UI IMPORT)

if exist('Data', 'var')

    quest = 'There is already data uploaded, would you like to upload new data, or proceed with the existing data?';
    answer = questdlg(quest, 'Data upload warning', 'Contintue with existing data', 'Upload new data', 'Contintue with existing data');

end

if ~exist('Data', 'var') || strcmp(answer, 'Upload new data')

    quest = 'Do you have a .asc file AND an .xy file? or just an .xy file';
    answer = questdlg(quest, 'Choose file types', '.asc and .xy', '.xy only', '.asc and .xy');

    switch answer

        case '.asc and .xy'
            clear Data

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

        case '.xy only'

            clear Data

            [filename, path1 ] = uigetfile('.xy', 'import full profile .xy spectrum', 'MultiSelect', 'off');
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


            for i = 1:length(file)

                Data.profileData = importdata(filename2);

            end

            %Calculate peak centroids for the full profile data
            for i = 1:length(file)

                 Data.(file{i,1}) = mspeaks(Data.profileData(:,1), Data.profileData(:,2));
                 Itemp = 100*Data.(file{i,1})(:,2)/max(Data.(file{i,1})(:,2)); %normalize ints
                 Data.(file{i,1})(Itemp < Icutoff, :) = []; %remove entries below intensity cutoff

                 I1 = max(Data.(file{i,1})(:,2)); %get the max int for the centroid data
                 I2 = max(Data.profileData(:,2));

                 if I1<I2
                    Data.(file{i,1})(:,2) = Data.(file{i,1})(:,2) + (I2-I1); 
                 end

            end


    end

end


%% Theo calculations

%Initialize report table
headers = {'Peptide', 'Unlabeled m/z', 'charge', lab{:,1}, '#15N', '#13C', 'Theo Labeled m/z', 'obs m/z', 'ppm error'};
resultTable = cell(size(sequence,1)*maxCharge, (length(headers)));
resultTable(:,2:end) = num2cell(0);

%Calculate the mass of each peptide (M + H)
for i = 1:size(sequence,1)
    i2 = maxCharge*(i-1);
    
    for k = 1:maxCharge %check each allowable charge state

        %put the current sequence in the result table
        resultTable(i2 +k, 1) = sequence(i,1);

        %calculate elemental formula for each charge state
        f = atomiccomp(sequence{i,1});

        %add protons equal to the charge (k)
        f.H = f.H + k;
        
        %Calc the monoisotopic mass of each peptide and put it in the table
        [MD, info] = isotopicdist(f); 

        %Dividee m/z in MD by the charge
        MD(:,1) = MD(:,1)./k;

        %get the monoisotopic m/z and put it in the first row of results
        %table
        resultTable(i2+k, 2) = num2cell(MD(1,1));

        %add the charge to the table
        resultTable(i2+k, 3) = num2cell(k);
        
        %count the number of AA in each peptide
        AAcount = aacount(sequence{i,1});
    
        %Extract the count for each label AA, calc #15N and #13C
        for j = 1:length(lab)
            %count the number of AA
            resultTable(i2+k, j + 3) = num2cell(AAcount.(lab{j,1}));
    
            %15N
            resultTable(i2+k, length(lab)+4) = num2cell(resultTable{i2+k, length(lab)+4} + resultTable{i2+k, j+3}*pm(j,1));
    
            %13C
            resultTable(i2+k, length(lab)+5) = num2cell(resultTable{i2+k, length(lab)+5} + resultTable{i2+k, j+3}*pm(j,2));
        end
    
        N = resultTable{i2+k, 7};
        C = resultTable{i2+k,8};

        %The old way for calculating the labeled m/z still works 11/21/23
        %calc labeled m/z: 
        %resultTable(i+k-1, length(lab)+5) = resultTable(i+k-1,1) + (dN*resultTable(i+k-1, length(lab)+3) + dC*resultTable(i+k-1, length(lab)+4))/k;
        
        %calculate labeled peptide monoisotopid mass
        if N ~= 0 && C == 0
            MD2 = isotopicdist2(f, 'Nsparse', N);
        elseif N == 0 && C ~=0
            MD2 = isotopicdist2(f,'Csparse', C);
        else    
            MD2 = isotopicdist2(f, 'Nsparse', N, 'Csparse', C);
        end

        %scale to charge
        MD2(:,1) = MD2(:,1)./k;
        
        %put the labeled "monoisotopic" m/z in the table
        resultTable(i2+k, length(lab)+6) = num2cell(MD2(1,1));



    end
end

clear MD MD2 f info N C i2 k 

%% Compare to data

for i = 1:size(resultTable,1)
    %get the current labeled theo m/z
    mass = resultTable{i, length(lab)+6};

    %truncate the m/z list to +-5 Da of mass
    %mz = Data.(file{1,1})(:,1);
    mz = Data.(file{1,1});

    %Normalize intensity to base peak
    mz(:,2) = (100*mz(:,2))/max(mz(:,2));

    mz(mz(:,1) < mass - 5,:) = [];
    mz(mz(:,1) > mass + 5,:) = [];

    %calc the element wise ppm difference between mass and m/z
    if ~isempty(mz)
        p = 1000000*(mass - mz(:,1))./mass;
        pabs = abs(p);
        
        %find nearest neighbor mass
        resultTable(i, length(lab) + 7) = num2cell(mz(pabs == min(pabs),1));

        %if the intensity falls below the cutoff, put 1000000 in the ppm
        %box
        if mz(pabs == min(pabs),2) < Icutoff

            %put an arbitrary large ppm in the ppm column
            resultTable(i, length(lab) + 8) = num2cell(1000000);

        elseif mz(pabs == min(pabs),2) >= Icutoff

            %put ppm error in table
            resultTable(i, length(lab) + 8) = num2cell(p(pabs == min(pabs),1));

        end

    elseif isempty(mz)
        
        %put an arbitrary large ppm in the ppm column
        resultTable(i, length(lab) + 8) = num2cell(1000000);

    end
end

 %updat the sequence vector to have multiples of each peptide for each
 %charge
if maxCharge > 1
    sequence2 = cell(size(sequence,1)*maxCharge,1);
    j = 1;
    count = 1;

    for i = 1:size(sequence2,1)

        if count > maxCharge
            j = j+1;
            count = 1;
        end
        
        sequence2(i,1) = sequence(j);

        count = count + 1;

        

    end

clear j count
elseif maxCharge == 1
    
    sequence2 = sequence;

end

%remove any rows from report table and result table that have 1000000 in
%the ppm column
%sequence2(cell2mat(resultTable(:,end)) == 1000000,:) = [];
resultTable(cell2mat(resultTable(:,end)) == 1000000,:) = [];

%remove any rows where the ppm error is greater than user ppm cutoff
%sequence2(cell2mat(resultTable(:,end)) > ppm, :) = [];
resultTable(cell2mat(resultTable(:,end)) > ppm, :) = [];

%sequence2(cell2mat(resultTable(:,end)) < -ppm, :) = [];
resultTable(cell2mat(resultTable(:,end)) < -ppm, :) = [];

%Generate a table for the final report with headers
finalReport = cell2table(resultTable, 'VariableNames', headers(1,:));



% %Generate Report Table
% rt = num2cell(resultTable);
% if maxCharge == 1
%     rt = [sequence2, rt];
% else
%     rt = [sequence2, rt];
% end
% reportTable = cell2table(rt);
% reportTable.Properties.VariableNames = headers(1,:);



clear mz p pabs

%% PLOT
figCount = 1;
%Generate a figure for each peptide
for i = 1:size(resultTable,1)

 
    %get the current mass
    mass = resultTable{i, length(lab)+6};
    

    %get the current sequence
    s = resultTable{i,1};

    N = resultTable{i, length(lab)+4}; %get the number of 15N
    C = resultTable{i, length(lab)+5}; %get the number of 13C
    p = num2str(resultTable{i,end},3); %get the ppm error 

    %extract the data +- 15 m/z around the current peptide
    d = Data.profileData; %full profile
    dc = Data.(file{1,1}); %centroided

    %normalize the Y axis to maximum peak
    d(:,2) = 100*d(:,2)./max(d(:,2));
    dc(:,2) = 100*dc(:,2)./max(dc(:,2));

    d(d(:,1)< mass - 15,:) = [];
    d(d(:,1)> mass + 15,:) = [];
 
    dc(dc(:,1)< mass - 15,:) = [];
    dc(dc(:,1)> mass + 15,:) = [];

    %Calculate the formula
    f = atomiccomp(s);

    %get the charge
    charge = resultTable{i,3};

    %add protons proportional to charge
    f.H = f.H + charge;

    %generate unlabeled distribution
    [MD, Info] = isotopicdist(f, 'NoiseThreshold', 100);

    if N ~= 0 || C ~= 0 %generate a labeled distribution if there are labeles
        [MD2, Info2] = isotopicdist2(f, 'Nsparse', N, 'Csparse', C, 'NoiseThreshold', 100);
    end

    %Scale each MD by the charge
    MD(:,1) = MD(:,1)/charge;

    if exist('MD2', 'var')
        MD2(:,1) = MD2(:,1)/charge;   
    end

    %normalize the monoisotopic peak of MD to the nearest neighbor in d
    int = dc(dc(:,1) == resultTable{i, length(lab)+6},2);
    scale = int/MD(1,2);
    %MD(:,2) = MD(:,2).*scale;

    %Normalize unlabeled distribution to nearest neighbor
    m = resultTable{i,2}; %get theo unlabeled mass
    p2 = abs(1000000*(dc(:,1)-m)./m);%ppm error calculation
    int = dc(p2 == min(p2),2); %extract nearest intensity value for mass within tolerance
    scale = int/MD(1,2);
    MD(:,2) = MD(:,2)*scale;

    if exist('MD2','var')
        
        m = MD2(1,1);%pull out the monoisotopic peak
        p2 = abs(1000000*(dc(:,1)-m)./m); %ppm error calc
        int = dc(p2 == min(p2),2); %get the intensity for the peak closet to the theo peak
        scale = int/MD2(1,2);
        MD2(:,2) = MD2(:,2)*scale;

    end


%         if min(p2) < ppm
%             scale2 = int/MD2(1,2);
%             MD2(:,2) = MD2(:,2).*scale2;
%         else
%             MD2(:,2) = MD2(:,2).*scale2; %if there is no unlabeled peak, scale to labeled peak
%         end
%         
%         MD(:,2) = MD(:,2).*scale2;

    figure(figCount)
    if strcmp(mode, "profile")
        plot(d(:,1), d(:,2))
    elseif strcmp(mode,"centroid")
        stem(dc(:,1), dc(:,2), 'Marker', '|')
    end
    hold on
    scatter(MD(:,1), MD(:,2), 'filled')
    if N ~= 0 || C ~=0 %only add the labeled dist if there are labels.
        scatter(MD2(:,1), MD2(:,2), 'filled')
    end
    xlabel('m/z')
    ylabel('Intensity')
    legend('observed', 'simulated, unlabeled', 'simulated, labeled')
    title(sprintf(s))
    formatSpec = 'Labeled AA: %s, %s, %s, #15N: %d, #13C: %d \nlabeled m/z: %s, charge: %d, ppm error: %s';
    subtitle(sprintf(formatSpec, lab{1,:}, lab{2,:}, lab{3,:}, N, C, num2str(mass), charge, p))
    drawnow

    figCount = figCount + 1;

end
    
    if optimize == 1 && ((N + C) ~=0)
        %[incorp, figCount] = sIsotopeRMSE(sequence(i), d, dc, MD, MD2, lab, resultTable, figCount);
    end

%perform the fitting process for each peptide in the results table
if scrmbl == 1
    dnd = 0; %switch to disable pop-up prompt

    for i = 1:size(resultTable,1)

        %get the current mass
        mass = resultTable{i, length(lab)+6};


        %get the current sequence
        s = resultTable{i,1};

        %extract the data +- 15 m/z around the current peptide
        d = Data.profileData; %full profile
        dc = Data.(file{1,1}); %centroided

        %normalize the Y axis to maximum peak
        d(:,2) = 100*d(:,2)./max(d(:,2));
        dc(:,2) = 100*dc(:,2)./max(dc(:,2));

        d(d(:,1)< mass - 15,:) = [];
        d(d(:,1)> mass + 15,:) = [];

        dc(dc(:,1)< mass - 15,:) = [];
        dc(dc(:,1)> mass + 15,:) = [];

        %Calculate the formula
        f = atomiccomp(s);

        %get the charge
        charge = resultTable{i,3};

        %add protons proportional to charge
        f.H = f.H + charge;


        %run the fitting function
        [incorp, figCount, dnd] = sIsotopeRMSEv3(s, d, dc, inc, resultTable(i,:), figCount, ppm, dnd);

        %Convert incorporations to percent
        incorp(:,end) =  incorp(:,end)*100;

    end


end
    
    



%clear N C d dc f m mass int Info Info2 i p p2 MD 
clear C dc f formatSpec headers Info Info2 int m mass MD MD2 scale
clear AAcount answer b c d d2 DF DFtemp dI dists e fields fp fspec fspec2 i l l3 Icomb idx1 idx2 Itemp
clear n N numFields o p P p2 P2 pkInfo poop ppmPoop ppmTemp quest step sz w windowCenter xHighLim xLowLim xmax ymax

%% save report table
quest = 'Do you want to save the output table? A copy will always be available in the workspace in the variable "finalReport"';
answer = questdlg(quest);

if strcmp (answer, 'Yes')
    [f, p] = uiputfile('.csv', 'Specify directory and filename for results table');
    fp = strcat(p, f);
    writetable(finalReport, fp);
end

clear quest answer



%% Calculate combination distribtion from linear combinations of the labeled and unlabeled distributions
%This function determines the percent incorporation for a heavy isotope
%label when there is only one label present 

function [incorp, figCount] = sIsotopeRMSE(sequence, d, dc, MD, MD2, lab, resultTable, figCount)

    
    A = [1,0]; %weights for linear combinations. Each element must sum to 1
    numSteps = 100; %how many intervals between 0 and 1 to test
    steps = linspace(0,1,numSteps);

    %normalize both MD and MD2 so they start from equal intensity 
    MD(:,2) = 100*MD(:,2)/max(MD(:,2));
    MD2(:,2) = 100*MD2(:,2)/max(MD2(:,2));
    
    dc2 = dc; %make a copy of the centroided data
    
    %generate a new matrix that contains MD and MD2, align each peak by nominal
    %mass
    %calculate the number of unique nominal masses in MD and MD2
    nom = unique(round([MD(:,1); MD2(:,1)]));
    %calculate how many Da labeled shift to index it into dists
    n = resultTable(1, length(lab) + 3) + resultTable(1, length(lab) + 4); %add values from the %15N and #13C columns
    
    dists = zeros(length(nom),4);
    dists(1:size(MD,1), 1:2) = MD; %enter unlabeled distribution
    dists((1+n):(size(MD2,1)+n), 3:4) = MD2; %enter labeled distribution
    
    %Initialize a linear combination distribution matrix 
    MDlc = zeros(size(dists,1),2);
    
    %get the peaks in dc that line up with the theo peaks
    
    for j = 1:size(MDlc,1) %populate each m/z value
    
        if dists(j,1) ~= 0 && dists(j,3) ~= 0 %average masses if both distributions contribute
            MDlc(j,1) = mean([dists(j,1), dists(j,3)]);
        elseif dists(j,1) == 0 %if only one distribution has a mass at this value, just move it over
            MDlc(j,1) = dists(j,3);
        elseif dists(j,3) == 0
            MDlc(j,1) = dists(j,1);
        end
    
    end
    
    if any(~any(MDlc,2)) %if any rows contain all zeros, remove the row
        MDlc(~any(MDlc,2),:) = [];
        dists(~any(dists,2),:) = [];
    end

    %remove values from dc2 greater than and less than the theo dists
    dc2(dc2(:,1) < (min(MDlc(:,1))-0.1), :) = [];
    dc2(dc2(:,1) > (max(MDlc(:,1))+ 0.1), :) = [];
    
    ppmMat = abs(1000000*((MDlc(:,1) - dc2(:,1)')./MDlc(:,1))); %find the nearest neighbor isotope peak for each theo peak
    [M,I] = min(ppmMat,[],2 ); %get the index of the smallest ppm in each row
    
    
    dc2 = dc2(I,:); %remove any values that are not part of the isotope distribution
    
    RMSE = zeros(numSteps,1);
    
    %LOOP
    for j = 1:numSteps

        %Set weighting coefficients
        A(1,2) = steps(j); %weight for labeled dist
        A(1,1) = 1-steps(j); %weight for unlabeled dist
    
        %apply weights to intensities
        MDlc(:,2) = A(1,1)*dists(:,2) + A(1,2)*dists(:,4);

        %Scale to the unlabeled distribution
        I = MDlc(1,2);
        I2 = dc2(1,2);
        scale3 = I2/I;
        MDlc(:,2) = scale3*MDlc(:,2);

        %If any intensity exceeds 100%, scale the entire distribution back
        %down so the max is 100%
        if max(MDlc(:,2)) > 100
            I = max(MDlc(:,2));
            scale = 100/I;
            MDlc(:,2) = MDlc(:,2)*scale;

        end
    
        %Find the largest scale coefficient and scale the new distribution
        %(CHANGE THIS TO NEAREST NEIGHBOR TO THE MAX PEAK)
%         maxI = max(MDlc(:,2));
%         maxI2 = max(dc2(:,2));
%         scale3 = maxI2/maxI;
%         MDlc(:,2) = scale3*MDlc(:,2);
    
        %RMSE calc
        RMSE(j,1) = rmse(MDlc(:,2), dc2(:,2));
    
        f = figure(figCount);
        formatspec = 'Theoretical vs. observed mass spectrum for %s';
        t = sprintf(formatspec, sequence{:});
        
        subplot(1,2,1)
        plot(d(:,1), d(:,2))
        title(t)
        hold on
        scatter(MDlc(:,1), MDlc(:,2), 'filled')
        xlabel('m/z')
        ylabel('Intensity')
        hold off

        subplot(1,2,2)
        title('RMSE')
        plot(100*steps,RMSE(:,1))
        xlabel('Percent Incorpoation')
        ylabel('RMSE')
        hold off
        
        
    end
    %close(f)
    [pks, locs] = findpeaks(-RMSE);

    %recalculate distribution at best fitting incorporation
    %Set weighting coefficients
    A(1,2) = steps(locs); %weight for labeled dist
    A(1,1) = 1-steps(locs); %weight for unlabeled dist

    %apply weights to intensities
    MDlc(:,2) = A(1,1)*dists(:,2) + A(1,2)*dists(:,4);
    
    %Find the largest scale coefficient and scale the new distribution
    maxI = max(MDlc(:,2));
    maxI2 = max(dc2(:,2));
    scale3 = maxI2/maxI;
    MDlc(:,2) = scale3*MDlc(:,2);

    figure(figCount)
    subplot(1,2,1)
    plot(d(:,1), d(:,2))
    title(t)
    hold on
    scatter(MDlc(:,1), MDlc(:,2), 'filled')
    xlabel('m/z')
    ylabel('Intensity')
    hold off
    ylim([0,105])


    subplot(1,2,2)
    title('RMSE')
    plot(100*steps,RMSE(:,1))
    hold on
    xlabel('Percent Incorpoation')
    ylabel('RMSE')
    labl = num2str(100*steps(locs), 2);
    xline(locs, '--b', labl);
    hold off

    figCount = figCount + 1; 
    incorp = labl;


 
end

function [incorp, figCount] = sIsotopeRMSE3Dist(sequence, d, dc, MD, MD2, MD3, lab, resultTable, figCount)



    
    A = [1,0]; %weights for linear combinations. Each element must sum to 1
    numSteps = 100; %how many intervals between 0 and 1 to test
    steps = linspace(0,1,numSteps);

    %normalize both MD and MD2 so they start from equal intensity 
    MD(:,2) = 100*MD(:,2)/max(MD(:,2)); %unlabeled
    MD2(:,2) = 100*MD2(:,2)/max(MD2(:,2)); %scrambled
    MD3(:,2) = 100*MD3(:,2)/max(MD3(:,2)); %Fully labeled
    
    dc2 = dc; %make a copy of the centroided data

    %remove values from dc2 greater than and less than the theo dists
    dc2(dc2(:,1) < (min(MD(:,1))-0.1), :) = [];
    dc2(dc2(:,1) > (max(MD3(:,1))+ 0.1), :) = [];

    %FIT THE UNLABELED AND SCRAMBLED DIST
    %I am going to just weight the monoisotopic peaks for now, SAMPEGVER
    %has no overlap between the two
    m1 = MD(1,:);
    m2 = MD2(1,:);

    %Scale to the unlabeled distribution
        I = m1(1,2);
        I2 = dc2(1,2);
        scale = I2/I;
        MD(:,2) = scale*MD(:,2);
    %Scale the scramble distribution
        I = m2(1,2);
        p = abs((1000000*(dc2(:,1)-m2(1,1)))/m2(1,1));
        I2 = dc2(p == min(p),2);
        scale = I2/I;
        MD2(:,2) = MD2(:,2)*scale;

        figure(100)
        stem(dc2(:,1), dc2(:,2))
        hold on
        scatter(MD(:,1), MD(:,2))
        scatter(MD2(:,1), MD2(:,2))

    %Combine unlabeled and scrambled distribution and re-normalize to 100
    %YOU MIGHT NEED TO CHANGE THIS LATER CAUSE IT ONLY WORKS IF THE
    %DISTRIBUTIONS DON'T OVERLAP
    MDcomb = [MD; MD2];
    MDcomb(:,2) = 100*MDcomb(:,2)/max(MDcomb(:,2));
    
    %generate a new matrix that contains MDcomb and MD3, align each peak by nominal
    %mass
    %calculate the number of unique nominal masses in MD and MD2
    nom = unique(round([MDcomb(:,1); MD3(:,1)]));
    %calculate how many Da labeled shift to index it into dists
    n = resultTable(1, length(lab) + 3) + resultTable(1, length(lab) + 4); %add values from the %15N and #13C columns
    
    dists = zeros(length(nom),4);
    dists(1:size(MDcomb,1), 1:2) = MDcomb; %enter unlabeled distribution
    dists((1+n):(size(MD3,1)+n), 3:4) = MD3; %enter labeled distribution


    
    %Initialize a linear combination distribution matrix 
    MDlc = zeros(size(dists,1),2);
    

    
    for j = 1:size(MDlc,1) %populate each m/z value
    
        if dists(j,1) ~= 0 && dists(j,3) ~= 0 %average masses if both distributions contribute
            MDlc(j,1) = mean([dists(j,1), dists(j,3)]);
        elseif dists(j,1) == 0 %if only one distribution has a mass at this value, just move it over
            MDlc(j,1) = dists(j,3);
        elseif dists(j,3) == 0
            MDlc(j,1) = dists(j,1);
        end
    
    end
    
    if any(~any(MDlc,2)) %if any rows contain all zeros, remove the row
        MDlc(~any(MDlc,2),:) = [];
        dists(~any(dists,2),:) = [];
    end

    
    ppmMat = abs(1000000*((MDlc(:,1) - dc2(:,1)')./MDlc(:,1))); %find the nearest neighbor isotope peak for each theo peak
    [M,I] = min(ppmMat,[],2 ); %get the index of the smallest ppm in each row
    
    
    dc2 = dc2(I,:); %remove any values that are not part of the isotope distribution
    
    RMSE = zeros(numSteps,1);
    
    %LOOP
    for j = 1:numSteps

        %Set weighting coefficients
        A(1,2) = steps(j); %weight for labeled dist
        A(1,1) = 1-steps(j); %weight for unlabeled dist
    
        %apply weights to intensities
        MDlc(:,2) = A(1,1)*dists(:,2) + A(1,2)*dists(:,4);

        %Scale to the unlabeled distribution
        I = MDlc(1,2);
        I2 = dc2(1,2);
        scale3 = I2/I;
        MDlc(:,2) = scale3*MDlc(:,2);

        %If any intensity exceeds 100%, scale the entire distribution back
        %down so the max is 100%
        if max(MDlc(:,2)) > 100
            I = max(MDlc(:,2));
            scale = 100/I;
            MDlc(:,2) = MDlc(:,2)*scale;

        end
    
        %Find the largest scale coefficient and scale the new distribution
        %(CHANGE THIS TO NEAREST NEIGHBOR TO THE MAX PEAK)
%         maxI = max(MDlc(:,2));
%         maxI2 = max(dc2(:,2));
%         scale3 = maxI2/maxI;
%         MDlc(:,2) = scale3*MDlc(:,2);
    
        %RMSE calc
        RMSE(j,1) = rmse(MDlc(:,2), dc2(:,2));
    
        f = figure(figCount);
        formatspec = 'Theoretical vs. observed mass spectrum for %s';
        t = sprintf(formatspec, sequence{:});
        
        subplot(1,2,1)
        plot(d(:,1), d(:,2))
        title(t)
        hold on
        scatter(MDlc(:,1), MDlc(:,2), 'filled')
        xlabel('m/z')
        ylabel('Intensity')
        hold off

        subplot(1,2,2)
        title('RMSE')
        plot(100*steps,RMSE(:,1))
        xlabel('Percent Incorpoation')
        ylabel('RMSE')
        hold off
        
        
    end
    %close(f)
    [pks, locs] = findpeaks(-RMSE);

    %recalculate distribution at best fitting incorporation
    %Set weighting coefficients
    A(1,2) = steps(locs); %weight for labeled dist
    A(1,1) = 1-steps(locs); %weight for unlabeled dist

    %apply weights to intensities
    MDlc(:,2) = A(1,1)*dists(:,2) + A(1,2)*dists(:,4);
    
    %Find the largest scale coefficient and scale the new distribution
    maxI = max(MDlc(:,2));
    maxI2 = max(dc2(:,2));
    scale3 = maxI2/maxI;
    MDlc(:,2) = scale3*MDlc(:,2);

    figure(figCount)
    subplot(1,2,1)
    plot(d(:,1), d(:,2))
    title(t)
    hold on
    scatter(MDlc(:,1), MDlc(:,2), 'filled')
    xlabel('m/z')
    ylabel('Intensity')
    hold off
    ylim([0,105])


    subplot(1,2,2)
    title('RMSE')
    plot(100*steps,RMSE(:,1))
    hold on
    xlabel('Percent Incorpoation')
    ylabel('RMSE')
    labl = num2str(100*steps(locs), 2);
    xline(locs, '--b', labl);
    hold off

    figCount = figCount + 1; 
    incorp = labl;
 
end

function [incorp, figCount, dnd] = sIsotopeRMSEv3(sequence, d, dc, inc, resultTable, figCount, ppm, dnd)
%% sequence is a char array protein sequence
sequence = resultTable{1,1};


%d is full profile MS data (1st column is m/z and 2nd column is intensity)

%dc is centroided MS data 

%inc is a 1x2 matrix that contains the incriments in which labels are
%stripped, must be integer values. ex. [1,5] First value is for 15N and
%second is for 13C

%resultTable is a single row from the table of matches from SparseDigest

%FigCount is the current figure number

%plotVersion = 'v1';
plotVersion = 'v2';

    a = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
    

    %make a list of all the peaks in the observed data and their widths in
    %m/z
    [pks, locs, w, pr] = findpeaks(d(:,2), d(:,1), 'NPeaks', 20, 'SortStr','descend'); 
    pkInfo = [locs, pks, w, pr];
    clear locs pks w pr

%     figure(1)
%     plot(d(:,1),d(:,2))
%     hold on
%     %xline(locs)

%     figure(1)
%     hold on
%     plot(d(:,1), d(:,2))
%     hold off


    %find the peak width of the monoisotopic peak
    m = resultTable{1,end-1};
    ppmTemp = abs(1000000*(m- pkInfo(:,1))./m);
    w = pkInfo(ppmTemp == min(ppmTemp), 3); %Peak width to use for simulations

    f = atomiccomp(sequence);

    charge = resultTable{1,3};
    f.H = f.H + charge;

    NC = [resultTable{1,7}, resultTable{1,8}];

    n = (1+NC(1,1)/inc(1,1))*(1+NC(1,2)/inc(1,2)); %calculate number of sparse label combintations
    NCscram = zeros(n,2);

    %Generate a list of all possible N and C incorporation combinations
    idx = 1;
    for i = 1:(1+NC(1,2)/inc(1,2))
        for j = 1:(1+NC(1,1)/inc(1,1))
            %i is carbon j is nitrogen

            NCscram(idx,1) = (j-1)*inc(1,1); %
            NCscram(idx,2) = (i-1)*inc(1,2);
            idx = idx + 1;
            
        end
    end
    clear idx

    if dnd == 0
        q = 'The fitting process is about to begin, would you like to review the possible labeling combinations?';
        ttl = 'Determining label incorporation';
        answer = questdlg(q, ttl, 'Yes', 'No', 'Don''t ask again', 'No');
    
        %Handle response
        switch answer
            case 'Yes'
    
              %fig = uifigure('CloseRequestFcn', @(varargin)uiresume(gcbf));
              fig = uifigure('CloseRequestFcn', @(src,event) uiresume(src));
              vars = {'15N','13C'};
              tbl = table(NCscram(:,1), NCscram(:,2), 'VariableNames', vars);
              UItbl = uitable(fig, "Data", tbl, "CellSelectionCallback", @(h,e) set(h, 'UserData', e));
              UItbl.ColumnEditable = true;
              b = uibutton(fig, "Text", "Add Row", ...
                            "ButtonPushedFcn", @(src,event) buttonPushed(UItbl, fig), ...
                            "Position", [400, 250, 100, 22]);
              b2 = uibutton(fig, "Text", "Delete Row", ...
                            "ButtonPushedFcn", @(src,event) buttonPushed2(UItbl, fig), ...
                            "Position", [400, 200, 100, 22]);
              b3 = uibutton(fig, "Text", "Close", ...
                          "ButtonPushedFcn", @(src,event) buttonPushed3(UItbl, fig), ...
                          "Position", [400, 150, 100, 22]);
    
    
              uiwait(fig);
    
              test = UItbl.Data;
              close(fig);
              delete(fig);
    
              NCscram = table2array(test);
              n = length(NCscram);
              clear vars test
    
            case 'no'
    
                %close fig;
    
            case 'Don''t ask again'
    
                %close fig;
                dnd = 1;
        end
    end

    %get the peak width (in Da) of the monoisotopic peak (if its there)
    %Else, get the width of the fully labeled peak 
    
    dists = struct();
    %Generate all possible isotope distributions
    for i = 1:n
        %fd = num2str(i);

        if w < 0.001 %increase point density of isotope simulation for fine structure spectra
            [MD, info, DF] = isotopicdist2(f, 'Nsparse', NCscram(i,1), 'Csparse', NCscram(i,2),'FFTResolution', 10000, 'Resolution', w, 'FFTrange', 20);
        else
            [MD, info, DF] = isotopicdist2(f, 'Nsparse', NCscram(i,1), 'Csparse', NCscram(i,2), 'Resolution', w/2, 'FFTrange', 20);
        end

        clear info
       
        %Normalize values in DF to 1
        DF(:,2) = DF(:,2)/max(DF(:,2));
        MD(:,2) = MD(:,2)/max(MD(:,2));

        %scale distribution to charge
        DF(:,1) = DF(:,1)./charge;
        MD(:,1) = MD(:,1)./charge;
        
        dists.(a(i)).DF = DF;
        dists.(a(i)).MD = MD;

        if i == n
            xmax = max(MD(:,1) + 1);
        end
      
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Applying single point calibration to d & dc
    mono2 = dists.(a(1)).MD(1,1);
    p2 = abs(1000000*(dc(:,1)-mono2)./mono2);%ppm error calculation
    dcTemp = dc(p2(:,1)<5,:);
    CalPoint = dcTemp(dcTemp(:,2)==max(dcTemp(:,2)));
    
    if ~isempty(CalPoint)
        CalFactor = mono2 - CalPoint;
        dc(:,1) = dc(:,1) + CalFactor;
        d(:,1) = d(:,1) + CalFactor;        
        %error('No monoisotopic peak found. Please check your input sequence, charge, and dataset\n If none of those fix it, there might be no monoisotopic peak in the data.')
    end
    clear CalFactor CalPoint mono2 p2 dcTemp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mz = dists.(a(1)).DF(:,1); %Get the vector of m/z 
    fields = fieldnames(dists);


    if w < 0.001 %data trimming for fine structure stuff

        fl = length(fields);
    
        %Find the largest m/z across all dists and remove it to save memory
        pks = mspeaks(dists.(a(fl)).DF(:,1), dists.(a(fl)).DF(:,2), 'HeightFilter', 0.01);
        pks = max(pks(:,1))+1;

        %Trim the unified mz vector
         mz(mz>pks) = [];

         %get its length
         l = length(mz);

        for i = 1:length(fields) %trim each simulation so they are the same length
            
            l1 = length(dists.(a(i)).DF); %calc initial length
            dists.(a(i)).DF(dists.(a(i)).DF > max(mz),:) = []; %remove values above maximum value
            l2 = length(dists.(a(i)).DF); %length after trimming
            %tempmz = zeros(l1 - l2,2); %zero fill vector
            %tempmz(:,1) = mz(1:(l1-l2)); %add m/z to zero fill vector
            %dists.(a(i)) = [tempmz;  dists.(a(i))]; %concatenate so all distributions are same length
    
            l3 = length(dists.(a(i)).DF);
            l = max(l, l3);
    
            if i > 1 %Add pad to left size with 0 intensity
                tempmz = zeros(l - l3,2);
                tempmz(:,1) = mz(1:(l-l3));
                dists.(a(i)).DF = [tempmz;  dists.(a(i)).DF]; %concatenate zero pad to left side so all distributions are same length
            end
    
        end
    
    else %data trimming for standard isotopic resolution

            l = length(dists.(a(1)).DF);

            
            
            for i = 1:length(fields) %trim each simulation so they are the same length
            
                l1 = length(dists.(a(i)).DF); %calc initial length
                dists.(a(i)).DF(dists.(a(i)).DF > max(mz),:) = []; %remove values above maximum value
                l2 = length(dists.(a(i)).DF); %length after trimming
                tempmz = zeros(l1 - l2,2); %zero fill vector
                tempmz(:,1) = mz(1:(l1-l2)); %add m/z to zero fill vector
                dists.(a(i)).DF = [tempmz;  dists.(a(i)).DF]; %concatenate so all distributions are same length
        
                l3 = length(dists.(a(i)).DF);
                l = max(l, l3);
        
                if l > l3 %Add pad to left size with 0 intensity
                    tempmz = zeros(l - l3,2);
                    tempmz(:,1) = mz(1:(l-l3));
                    dists.(a(i)).DF = [tempmz;  dists.(a(i)).DF]; %concatenate zero pad to left side so all distributions are same length
                end
            end
    
    end

    clear l l1 l2 l3 tempmz

    n = 0:0.01:1;
    A = zeros(length(fields),1); %weighting coefficients for each distribution
    A(1,1) = 1; %All ratios will be calculated relative to unlabeled. 

    mono = resultTable{1,2};
    DFcomb = zeros(length(mz),2);
    DFcomb(:,1) = mz; %simulation of the combined distributions
    %DFcomb = [mz, zeros(length(mz))]; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START THE FITTING PROCESS

    for i = 1:length(fields) %loop through each of the labeled distributions
%         if i == 4
%             sprintf('testing');
%         end
        
        
        %SCALE THE FIRST DISTRIBUTION TO THE MONOISOTOPIC PEAK IF IT IS
        %THERE
        if i == 1
            %scale the height of all distributions to the unlabeled, if no unlabeled present, scale to fully labeled
            ppmTemp = abs(1000000*(dc(:,1)-mono)./mono);

            if any(ppmTemp < ppm) %if any peaks meet ppm criteria, scale to mono

                I = dc(ppmTemp == min(ppmTemp), 2);
                scale = I/max(dists.(a(1)).DF(:,2));
                A(1,1) = scale;
                DFcomb(:,2) = A(1,1)*dists.(a(1)).DF(:,2); %"A" contains scaling coefficients for all distributions

            elseif ~any(ppmTemp(ppmTemp < ppm)) %if there is no mono peak, set scale to 0

                A(1,1) = 0;

            end

            %generate a peak list for the combined distribution
            %P1 = mspeaks(DFcomb(:,1), DFcomb(:,2), 'HeightFilter', 0.1);

        elseif i > 1 %PRECALCULATIONS BEFORE FITTING

            I = 0; %Reset the Intensity threshold to 0

            %Get the peak values and heights from the combined distribution
            P = mspeaks(DFcomb(:,1), DFcomb(:,2), 'HeightFilter', 0.01); %Get the peak list for the combined distribution
            MD = dists.(a(i)).MD; %get the peak list for the current isotope distribution
            %poop = mspeaks(dists.(a(i)).DF(:,1), dists.(a(i)).DF(:,2), 'HeightFilter', 0.1); %Get the peak list for the nth isotope distribution
            n = sum(NCscram(i,:)) + 1; %calculate number of neutrons added +1 to get nth peak in distribution

            %calculate the index of the first peak that would change
            Atemp = A;
            Atemp(i,1) = Atemp(i,1) + 10; %arbitrarily increase A to find the index of the first peak that changes

            Itemp = zeros(length(DFcomb),1);

            for k = 1:length(Atemp)
                Itemp = Itemp + Atemp(k,1)*dists.(a(k)).DF(:,2);
            end

            %temporary peak list
            P2 = mspeaks(DFcomb(:,1), Itemp, 'HeightFilter', 0.01);

            %Calculate the difference in intensities between P and P2
            l1 = length(P);
            l2 = length(P2);
            
            if l2-l1 == 0 %if the vectors are the same lenght, subtract
                dI = P2(:,2) - P(:,2);
            elseif l2-l1 > 0 % add zero pad if the vectors are different lenghts
                dI = P2(:,2) - [P(:,2); zeros(l2-l1,1)]; 
            end

            %replace negative values in dI with zero
            dI(dI < 0) = 0;

            %replace small intensities (I<0.5) with zero
            dI(dI < 0.5) = 0;

            %get the index of the first positive value in dI adn set it as
            %n
            idx = find(dI, 1, 'first');

            n = idx;

            clear l1 l2 idx Atemp Itemp 

            %if the index n is out of range, recalculate the index for the
            %nth isotope for the current isotope pattern
            if (~isempty(P)) && (n > length(P)) || exist("N", "var")

                l = length(P);
                nl = n - l;
                P = [P; zeros(nl,2)];
                P(n,1) = P2(n,1);

            end

            ppmMD = abs(1000000*(dc(:,1)-MD(1,1))./MD(1,1));

            if ~isempty(P)
                ppmTemp = abs(1000000*(dc(:,1)-P(n,1))./P(n,1));
            elseif isempty(P)
                ppmTemp = zeros(size(ppmMD)) + 100;
            end

            %Get the indices of the match between the combined dist and the
            %data as well as the nth distribution and the data, if any
            idx1 = ppmTemp < ppm;
            idx2 = ppmMD < ppm;

            if any(ppmTemp < 2*ppm) && (~any(idx1) && ~any(idx2))
                fspec = 'Error in the fitting process: no peaks mathcing within %d ppm, but matches were found at %d ppm.';
                error(fspec, ppm, 2*ppm)
            end

            if any(idx1) && ~any(idx2) %Get the peak closest to the combination dist

                I = dc(ppmTemp == min(ppmTemp), 2);

            elseif ~any(idx1) && any(idx2) %Get the peak closest to the individual dist mono peak

                I = dc(ppmMD == min(ppmMD), 2);
 

            %if the both distributions have a match but they line up with
            %different peaks, get the intensity that lines up with the nth
            %distribution monoisotopic peak
            elseif any(idx1) && any(idx2) && ~all(idx1 == idx2)

                I = dc(ppmMD == min(ppmMD), 2);

            elseif any(idx1) && any(idx2) && all(idx1 == idx2) %if both have hits and they match the same peak, use P
                
                I = dc(ppmTemp == min(ppmTemp), 2);

            elseif ~any(idx1) && ~any(idx2) %if there are no matches, set I to 0;
                I = 0;
            end

            DFtemp = DFcomb; %temporary copy of combined dist

            %If P hasn't been established and there is a match for the
            %current dist
            if isempty(P) && any(idx2)
                DFcomb(:,2) = I*dists.(a(i)).DF(:,2);
                P = mspeaks(DFcomb(:,1), DFcomb(:,2), 'HeightFilter', 0.01);
                %if (~isempty(P)) && (n > length(P)) 
                    N = round((poop(1,1) - mono)*charge);
                    n = n - N;
               % end
            end

            %increase the weighting coefficient in steps until the nth peak matches the intensity in the data
            if ~isempty(P)
                while P(n,2) < I
                    
%                     if(P(n,2) > 80)
%                         sprintf('testing')
%                     end

                    e = 100*P(n,2)/I;
                    
                    if e < 90
                        step = I/20;
                    else
                        step = I/50;
                    end
                    
                    %Increase the current weighting coefficient
                    A(i,1) = A(i,1) + step;

                    %sum the intensities of all distributions
                    Icomb = zeros(length(DFcomb),1); %intensity vector for combined distribution

                    for k = 1:length(A)
                        Icomb = Icomb + A(k,1)*dists.(a(k)).DF(:,2);
                    end
    
                    P = mspeaks(DFcomb(:,1), Icomb(:,1), 'HeightFilter', 0.01);
                    

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % LIVE PLOT

                    if strcmp(plotVersion, 'v1') %PLOT VERSION 1

                        %Plot the individual distributions
                        lf = figure(400);
                        tiledlayout(1,2,'TileSpacing','tight')
                        set(lf, 'Position', [100, 100, 1120, 420])
                        set(lf, 'Color', [1,1,1])
                        nexttile
                        plot(dists.(a(1)).DF(:,1),A(1,1)*dists.(a(1)).DF(:,2))
                        I1 = max(A(1,1)*dists.(a(1)).DF(:,2));
                        I2 = 0;
                        hold on
    
                        for k = 2:length(A)
                            plot(dists.(a(k)).DF(:,1), A(k,1)*dists.(a(k)).DF(:,2))
                            In = max(A(k,1)*dists.(a(k)).DF(:,2));
                            I2 = max(I2, In);
                        end
    
                        Ilim = max(I1, I2) + 5;
                        hold off
                        title('Individual Simulations')
                        ylim([0,Ilim])
                        xlim([mono-1, xmax ])
    
                        %plot the combined distributions overlaid with the data
                        nexttile
                        plot(d(:,1),d(:,2))
                        hold on
                        plot(DFcomb(:,1), Icomb(:,1))
                        I3 = max(d(:,2))+5;
                        if ~isempty(P)
                            stem(P(:,1), P(:,2))
                        end
                        ylim([0, I3])
                        xlim([mono-1, xmax ])
                        hold off
                        title('Summed simulations')
                        drawnow

                    elseif strcmp(plotVersion, 'v2') %PLOTTING VERSION 2 WITH ZOOM IN ON EACH ISOTOPE

                        lf = figure(figCount);
                        
                        %, 'units', 'normalized', 'outerposition', [0 0 1 1]
                        titles = {'Peptide Spectrum'; "Monoisotopic" ;'1st Isotope'; '2nd Isotope'; '3rd Isotope' ;'4th Isotope'; '5th Isotope' 
                                    '6th Isotope' ;'7th Isotope'; '8th Isotope' ;'9th Isotope' ;'10th Isotope'};

                        numFields = length(fieldnames(dists));

                        c = size(NCscram,1) + 1;

                        t = tiledlayout(1,c, "TileSpacing", "Tight");
                        h = gobjects(1,c);
                        lf.WindowState = "maximized";

%                         if numFields<= 4
%                             t = tiledlayout(1,c);
%                             h = gobjects(1,c);
%                         elseif numFields > 4 && numFields <= 8
%                             t = tiledlayout(2,5);
%                             h = gobjects(1,10);
%                         elseif numFields > 8
%                             t = tiledlayout(3,4);
%                             h = gobjects(1,12);
%                         end

                        for k = 1:length(h)

                            if k == 1 %plot the first pane, full dist

                                nexttile
                                plot(d(:,1),d(:,2))
                                hold on
                                plot(DFcomb(:,1), Icomb(:,1))
                                I3 = max(d(:,2))+5;
                                if ~isempty(P)
                                    stem(P(:,1), P(:,2), "Color", "#D95319")
                                end
                                ylim([0, I3])
                                xlim([mono-1, xmax ])
                                hold off
                                title(sequence)

                            else %plot the nth pane

                                d2 = d;

                                %windowCenter = mono + (k-2)/charge;
                                if k == 2
                                    windowCenter = dists.(a(k-1)).MD(1,1);
                                elseif k > 2
                                    windowCenter = dists.(a(k-1)).MD(1,1); % + 0.006/charge;
                                end

                                if w < 0.001
                                    sz = 15*w; %peak width x 2
                                else
                                    sz = 4*w;
                                end
                                %sz = 1/(charge*25);
                                %sz = 1/(charge*50);
                                xLowLim = windowCenter - sz;
                                xHighLim = windowCenter + sz;
                                
                                d2(d2(:,1) < xLowLim,:) = [];
                                d2(d2(:,1) > xHighLim ,:) = [];

                                ymax = max(d2(:,2))+0.05*max(d2(:,2));

                                ax = nexttile(t);
                                h(k) = plot(d2(:,1), d2(:,2), 'LineWidth',2.0); %plot the data
                                hold on

%                                 for o = 1:numFields
%                                     stem(dists.(a(o)).MD(:,1), A(o,1)*dists.(a(o)).MD(:,2), 'Marker', 'none', 'LineWidth', 2.0);
%                                     %plot(dists.(a(o)).DF(:,1), A(o,1)*dists.(a(o)).DF(:,2), 'LineWidth', 2.0, 'LineStyle', ":")
%                                 end

                                plot(DFcomb(:,1), Icomb(:,1), 'LineWidth', 2.0, "Color", "#D95319") %plot the combined distribution
                                
                                if k <= 2
                                    title(titles(k,1));

                                    if k == 2
                                        subtitle('Unlabeled');
                                    end

                                elseif k > 2
                                    %calc the number of isotopes added
                                    b = sum(NCscram(k-1,:));
                                    if b == 1
                                        title('1st isotope')
                                    elseif b == 2
                                        title('2nd isotope')
                                    elseif b == 3
                                        title('3rd isotope')
                                    elseif b >= 4
                                        fspec = '%dth isotope';
                                        title(sprintf(fspec, b));
                                    end


                                    if k > 2
                                        fspec2 = 'Labels: %d 15N, %d 13C, %sx';
                                        poop = num2str(A(k-1), 2);
                                        subtitle(sprintf(fspec2, NCscram(k-1,1), NCscram(k-1,2), poop))
                                    end
                                    
                                end
                                
                                xlabel('m/z');
                                ylabel('Intensity');

                      
                                for o = 1:numFields
                                    %stem(dists.(a(o)).MD(:,1), A(o,1)*dists.(a(o)).MD(:,2), 'Marker', 'none', 'LineWidth', 2.0);
                                    plot(dists.(a(o)).DF(:,1), A(o,1)*dists.(a(o)).DF(:,2), 'LineWidth', 2.0, 'LineStyle', "--")
                                end

                                
                        
                                xlim([xLowLim, xHighLim])
                                ylim([0, ymax])
                                %ylim([0,7])

                                hold off
                                
                            end
                              
                        end

                        drawnow

                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % END PLOTTING SECTION

                    if isempty(P)
                        for k = 1:length(A)
                            Icomb = Icomb + A(k,1)*dists.(a(k)).DF(:,2);
                        end
                        
                        %update the peak list 
                        P = mspeaks(DFcomb(:,1), Icomb(:,1), 'HeightFilter', 0.01);
                    end
                end

                %Update intensities in DFcomb if it has any values
                if any(Icomb)
                    DFcomb(:,2) = Icomb;
                    Icomb = zeros(length(DFcomb),1);
                end


            end
                

        end
            

    end

    %close(lf)

    %Normalize weighting coefficients so they add up to 1
    Anorm = A(:,1)./sum(A(:,1));

    incorp = [NCscram Anorm];

    %Output and display percentages

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %OUTPUT FINAL PLOTS

    %Plot 1: Individual distributions scaled to their appropriate weights
    if strcmp(plotVersion, 'v1')
        figure(figCount)
        tiledlayout(1,2,'TileSpacing','tight')
        set(gcf, 'Position', [100, 100, 1120, 420])
        nexttile
        plot(dists.(a(1)).DF(:,1),A(1,1)*dists.(a(1)).DF(:,2))
        I1 = max(A(1,1)*dists.(a(1)).DF(:,2));
        I2 = 0;
        hold on

        for k = 2:length(A)
            plot(dists.(a(k)).DF(:,1), A(k,1)*dists.(a(k)).DF(:,2))
            In = max(A(k,1)*dists.(a(k)).DF(:,2));
            I2 = max(I2, In);
        end

        Ilim = max(I1, I2) + 5;
        hold off
        title('Individual Simulations')
        ylim([0,Ilim])
        xlim([mono-1, xmax ])

        %plot the combined distributions overlaid with the data
        nexttile
        plot(d(:,1),d(:,2))
        hold on
        %plot(DFcomb(:,1), Icomb(:,1))
        I3 = max(d(:,2))+5;
        if ~isempty(P)
            stem(P(:,1), P(:,2))
        end
        ylim([0, I3])
        xlim([mono-1, xmax ])
        hold off
        title('Summed simulations')

        figCount = figCount + 1;

    elseif strcmp(plotVersion, 'v2')

        %Plot 2: Combined distribution overlayed with the observed data.
        figure(figCount);
        set(gcf, 'WindowState', 'maximized');
        titles = {'Peptide Spectrum'; "Monoisotopic" ;'1st Isotope'; '2nd Isotope'; '3rd Isotope' ;'4th Isotope'; '5th Isotope'
            '6th Isotope' ;'7th Isotope'; '8th Isotope' ;'9th Isotope' ;'10th Isotope'};

        numFields = length(fieldnames(dists));

        c = size(NCscram,1) + 1;

        t = tiledlayout(1,c);
        h = gobjects(1,c);
        %
        % if numFields<= 4
        %     t = tiledlayout(1,c);
        %     h = gobjects(1,c);
        % elseif numFields > 4 && numFields <= 8
        %     t = tiledlayout(2,5);
        %     h = gobjects(1,10);
        % elseif numFields > 8
        %     t = tiledlayout(3,4);
        %     h = gobjects(1,12);
        % end

        for k = 1:length(h)

            if k == 1 %plot the first pane, full dist

                nexttile
                plot(d(:,1),d(:,2))
                hold on
                plot(DFcomb(:,1), DFcomb(:,2))
                I3 = max(d(:,2))+5;
                if ~isempty(P)
                    stem(P(:,1), P(:,2), "Color", "#D95319")
                end
                ylim([0, I3])
                xlim([mono-1, xmax ])
                hold off
                title(sequence)

            else %plot the nth pane

                d2 = d;

                if k == 2
                    %windowCenter = mono + (k-2)/charge;
                    windowCenter = dists.(a(k-1)).MD(1,1);
                elseif k > 2
                    windowCenter = dists.(a(k-1)).MD(1,1)  + 0.006/charge;
                end

                if w < 0.001
                    sz = 12*w; %peak width x 2
                else
                    sz = 4*w;
                end
                %sz = 1/(charge*25);
                %sz = 1/(charge*50);
                xLowLim = windowCenter - sz;
                xHighLim = windowCenter + sz;

                d2(d2(:,1) < xLowLim,:) = [];
                d2(d2(:,1) > xHighLim ,:) = [];

                ymax = max(d2(:,2))+0.05*max(d2(:,2));

                ax = nexttile(t);
                h(k) = plot(d2(:,1), d2(:,2), 'LineWidth',2.0); %plot the data
                hold on

                %         for o = 1:numFields
                %             stem(dists.(a(o)).MD(:,1), A(o,1)*dists.(a(o)).MD(:,2), 'Marker', 'none', 'LineWidth', 2.0);
                %             %plot(dists.(a(o)).DF(:,1), A(o,1)*dists.(a(o)).DF(:,2), 'LineWidth', 2.0, 'LineStyle', ":")
                %         end

                plot(DFcomb(:,1), DFcomb(:,2), 'LineWidth', 2.0, "Color", "#D95319") %plot the combined distribution

                if k <= 2
                    title(titles(k,1));

                    if k == 2
                        fspec2 = 'Unlabeled: %s%%';
                        poop = num2str(100*Anorm(1,1), 2);
                        subtitle(sprintf(fspec2, poop));
                    end

                elseif k > 2
                    %calc the number of isotopes added
                    b = sum(NCscram(k-1,:));
                    if b == 1
                        title('1st isotope')
                    elseif b == 2
                        title('2nd isotope')
                    elseif b == 3
                        title('3rd isotope')
                    elseif b >= 4
                        fspec = '%dth isotope';
                        title(sprintf(fspec, b));
                    end


                    if k > 2
                        fspec2 = 'Labels: %d 15N, %d 13C, %s%%';
                        poop = num2str(100*Anorm(k-1), 2);
                        subtitle(sprintf(fspec2, NCscram(k-1,1), NCscram(k-1,2), poop))
                    end

                end

                xlabel('m/z');
                ylabel('Intensity');


                for o = 1:numFields
                    %stem(dists.(a(o)).MD(:,1), A(o,1)*dists.(a(o)).MD(:,2), 'Marker', 'none', 'LineWidth', 2.0);
                    plot(dists.(a(o)).DF(:,1), A(o,1)*dists.(a(o)).DF(:,2), 'LineWidth', 2.0, 'LineStyle', "--")
                end



                xlim([xLowLim, xHighLim])
                ylim([0, ymax])
                hold off

                figCount = figCount + 1;

            end

        end
    end
end

function [UItbl] = buttonPushed(UItbl, fig)

  data = table2array(UItbl.Data);
  vNames = UItbl.Data.Properties.VariableNames;
  sz = size(data);
  data = [data; zeros(1,sz(1,2))];
  data = table(data(:,1), data(:,2), 'VariableNames', vNames);

  UItbl.Data = data;


end

function [UItbl] = buttonPushed2(UItbl, fig)

Index = get(UItbl, 'UserData');
idx = unique(Index.Indices(:,1));
D = table2array(UItbl.Data);
vars = UItbl.Data.Properties.VariableNames;
if ~isempty(Index)
    D(idx, :) = [];
end
UItbl.Data = table(D(:,1), D(:,2), 'VariableNames', vars);

end

function buttonPushed3(UItbl, fig)

    %close 
    uiresume(fig)
    
end
