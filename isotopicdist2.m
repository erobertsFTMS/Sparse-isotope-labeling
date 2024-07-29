function [md, info, df]=isotopicdist2(p,varargin)
%  ISOTOPICDIST high resolution isotope mass distribution and density function.
%
%  [MD,INFO,DF] = ISOTOPICDIST(P) returns the expected mass distribution in
%  MD, the expected density function in DF, and the monoisotopic mass,
%  average mass, most abundant mass, nominal mass, and empirical formula in
%  the structure INFO for a peptide sequence given by the input character 
%  vector or string P. P may also be a string vector or cell array of 
%  character vectors with multiple peptides.
%
%  ISOTOPICDIST2(B) accepts a numeric vector of the form [C H N O S] where C,
%  H, N, O, and S are non-negative values that represent the number of
%  atoms of Carbon, Hydrogen, Nitrogen, Oxygen and Sulfur respectively in
%  the compound. B may be a M by 5 numeric matrix with multiple compounds
%  along the rows.
%
%  ISOTOPICDIST2(F) accepts an empirical chemical formula represented by the
%  structure F. Fieldnames in F are valid element symbols (case sensitive)
%  and the respective values are the number of atoms for each element. F
%  may be a structure array with multiple chemical formulas.
%
%  ISOTOPICDIST2(...,'NTERMINAL',NTER) modifies the N-terminal of the peptide.
%  Valid options are: 'None', 'Amine' (default), 'Formyl', and 'Acetyl'.
%  NTER may also be a structure with the same format as F containing an
%  empirical formula of a custom modification.
%
%  ISOTOPICDIST2(...,'CTERMINAL',CTER) modifies the C-terminal of the peptide.
%  Valid options are: 'None', 'Free Acid' (default), and 'Amide'. CTER may
%  also be a structure with the same format as F containing an empirical
%  formula of a custom modification. 
%
%  ISOTOPICDIST2(...,'RESOLUTION',R) sets the approximate resolution of the
%  instrument, resolution is given as the Gaussian width at FWHH. R
%  defaults to 1/8 Daltons.
%
%  ISOTOPICDIST2(...,'FFTRESOLUTION',N) sets the number of points used to
%  compute the FFT. N is given in number of points per Dalton. N defaults to
%  1000. 
%
%  ISOTOPICDIST2(...,'FFTRANGE',Q) sets the absolute range (in Daltons) for
%  the FFT and the output density function. By default, Q is automatically
%  estimated based on the weight of the molecule. 
%
%  ISOTOPICDIST2(...,'C13abund',Q) sets relative contribution of 13C to the
%  value specified by the user so simulate uniformly labeled molecules. The 
%  value of 12C is automatically adjusted
%  so that %12C + %13C = 100. Acceptable values for Q are 0 - 100.
%
%  ISOTOPICDIST2(...,'N15abund',Q) sets relative contribution of 15N to the
%  value specified by the user. The value of 15N is automatically adjusted
%  so that %14N + %14N = 100. Acceptable values for Q are 0 - 100.
%
%  ISOTOPICDIST2(...,'Nsparse',Q) Adds Q 15N to the simulation. Q can
%  range from 1 to the number of nitrogen in the current input formula. The
%  number of natural abundance nitrogen is automatically subtracted by Q.
%  For the sparce labeled 15N, the A entries should be [15.0001088984 0 0 0
%  0] and the B entries should be [100 0 0 0 0]. 
%
%  ISOTOPICDIST2(...,'Csparse',Q) Adds Q 13C to the simulation. Q can
%  range from 1 to the number of nitrogen in the current input formula. The
%  number of natural abundance nitrogen is automatically subtracted by Q.
%  For the sparce labeled 13C, the A entries should be [13.0033548378 0 0 0
%  0] and the B entries should be [100 0 0 0 0]. 
%
%
%  Notes: 
%  1. The actual FFT range used internally by ISOTOPICDIST is further
%  increased such that Q*N is a power of two. 
%  2. Increase N and reduce R to achieve ultrahigh resolution but ensure
%  that Q*N is within the available memory. Ultrahigh resolution allows you
%  to resolve micropeaks with the same nominal mass but slightly different
%  exact masses. 
%
%  ISOTOPICDIST2(...,'FFTLOCATION',L) specifies the location of the FFT
%  range (window) defined by Q. L defaults to 1/16. It does this by setting
%  the location of the lower limit of the FFT range, relative to the
%  location of the monoisotopic peak, which is computed by ISOTOPICDIST.
%  The location of the lower limit of the FFT range is set to 
%          FFT_range_low = mass of the monoisotopic peak - (L * Q)
%
%  ISOTOPICDIST2(...,'NOISETHRESHOLD',T) removes any point in the mass
%  distribution that is smaller than 1/T times the most abundant mass. T
%  defaults to 1e6. 
%
%  ISOTOPICDIST2(...,'SHOWPLOT',SP) plots the isotopic density function and
%  the scaled mass distribution points that overlaps the density function
%  when SP is TRUE. SP defaults to FALSE when ISOTOPICDIST is called with
%  output arguments, otherwise, SP defaults to true. In the case of
%  multiple elements in one of the inputs P, F, or B, SP can also contain
%  an index to indicate which element ISOTOPICDIST uses for the plot. 
%
%
%
%   Examples:
%
%        % Calculate and display the isotopic distribution of Glutamine:
%        MD = isotopicdist([5 10 2 3 0],'SHOWPLOT',true);
%
%        % Calculate and display the isotopic distribution of the peptide
%        % 'MATLAP' with Acetyl N-terminal and Amide C-terminal:
%        MD = isotopicdist('MATLAP','NTERM','Acetyl','CTERM','Amide','SHOWPLOT',true);
%
%        % Display the isotopic distribution of the "averagine" model 
%        % with a given element composition:
%        isotopicdist([4.9384 7.7583 1.3577 1.4773 0.0417])

%   Copyright 2009-2012 The MathWorks, Inc.

% References:
% [1] Alan L. Rockwood, Steven L. Van Orden, and Richard D. Smith, "Rapid
%     Calculation of Isotope Distributions" Anal. Chem. 67:15 pp 26992704
%     (1995).
% [2] Alan L. Rockwood, Steven L. Van Orden, Richard D. Smith "Ultrahigh
%     Resolution Isotope Distribution Calculations" Rapid Commun. Mass
%     Spectrum 10, pp 5459 (1996).
% [3] Senko, M.W., et al., "Automated assignment of charge states from
%     resolved isotopic peaks for multiply charged ions" J. Am. Soc. Mass
%     Spectrom. 6 pp 5256 (1995).
% [4] Senko, M.W., et al., "Determination of monoisotopic masses and ion
%     populations for large biomolecules from resolved isotopic
%     distributions" J. Am. Soc. Mass Spectrom. 6, pp 229233 (1995).

% check input type (structure, string or array of 5 elements)
if nargin > 0
    p = convertStringsToChars(p);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if ~isstruct(p) && ~ischar(p) && ~iscellstr(p) && (~isnumeric(p) || size(p,2)~=5)
    error(message('bioinfo:isotopicdist:IncorrectInputArguments'));
end

bioinfochecknargin(nargin,1,mfilename);

ele = {'C','H','N','O','S','He','Li','Be','B','F','Ne','Na','Mg','Al','Si',...
       'P','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',...
       'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',...
       'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',...
       'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',...
       'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At',...
       'Rn','Fr','Ra','Ac','Th','Pa','U'};

[w, massRange, L, K, R, nt, plotId,Mnter,Anter,Bnter,symnterm,Mcter,Acter,Bcter,symcterm, Cchange, Cabund, Nchange, Nabund, Csparse, Nsparse] = parse_inputs(ele,nargout,varargin{:});
%INSERT [13C] and/or [15N] into the ele variable if the user specified
%sparse labeling

% When one sequence peptide or many sequence peptides are given convert that
% to the CHNOS matrix with element composition for each peptide:
if ischar(p) || iscellstr(p)
    M = peptide2chnos(p);
    syminM = (1:5)';
    [A,B] = chnosIsotopicComposition(Cchange, Cabund, Nchange, Nabund, Csparse, Nsparse);
elseif isnumeric(p) && size(p,2) == 5;
    M = p;
    syminM = (1:5)';
    [A,B] = chnosIsotopicComposition(Cchange, Cabund, Nchange, Nabund, Csparse, Nsparse);
elseif isstruct(p)
    try
        [M,A,B,syminM] = struct2formula(p(:), ele, Cchange, Cabund, Csparse, Nchange, Nabund, Nsparse); %%%%%%%%%%%%%
    catch ME
        if strcmp(ME.identifier,'MATLAB:cell2mat:MixedDataTypes') || strcmp(ME.identifier,'MATLAB:badsubscript')
            error(message('bioinfo:isotopicdist:invalidFormula'))
        else
            rethrow(ME)
        end
    end
end

% To update the M, A and B array with modification of the n-term
if any(Mnter)
    [syminM,~,h] = unique([syminM;symnterm]);
    A(h,:) = [A;Anter];
    B(h,:) = [B;Bnter];
    M = accumarray(h,[M Mnter])';
end

% To update the M, A and B array with modification of the c-term
if any(Mcter)
    [syminM,~,h] = unique([syminM;symcterm]);
    A(h,:) = [A;Acter];
    B(h,:) = [B;Bcter];
    M = accumarray(h,[M Mcter])';
end

%IDIOT CHECK: The user cannot add more sparse labeles than there are in the
%curernt elemental formula
if Csparse >= M(1)
    Csparse
    M
    %error('The number of labeled 13C must be less than the number of carbon in the submitted formula')
end

if Nsparse >= M(3)
    Nsparse
    M
    %error('The number of labeled 15N must be less than the number of carbon in the submitted formula')
end

if (size(M,2) > 5) && (Csparse > 0 || Nsparse > 0)
    error('Hey buddy, it looks like you are doing sparse labeling and using elements other than CHNOS, sorry but I didnt have time to figure out how to code that, sorry! ~Eli')
end

%Update the M array with any sparse labels
if Csparse > 0 && Nsparse == 0 % Add 13C labels
    M = [M(1) - Csparse, Csparse, M(2:end)]; %Insert 13C into formula and subtract the number of labeled carbons from total carbons
end

if Nsparse > 0 && Csparse == 0 % Add 15N labels
    M = [M(1:2), M(3) - Nsparse, Nsparse, M(4:end)]; %Insert 13C into formula and subtract the number of labeled carbons from total carbons
end

if Nsparse > 0 && Csparse > 0 % Add BOTH 13C and 15N labels 
    M = [M(1) - Csparse, Csparse, M(2), M(3) - Nsparse, Nsparse, M(4:end)];%Do both 13C and 15N
end

% Removing columns with zeros
h = any([A;B],1);
A = A(:,h);
B = B(:,h);

% m is the number of elements
% n is the maximum number of isotopes for the current elements
% q = number of compounds to process
[m n] = size(A);
q = size(M,1);

% Make sure that the isotopic abundance of every element adds up to 1
for i = 1:m
    B(i,:) = B(i,:) ./ sum(B(i,:));
end


% Mass range required for low resolution Gaussian pulse:
mrgp = sqrt(-log(.001)./log(2)).*w./L;

% Mass range required for very large molecules (assuming a symmetric
% distribution around the average mass)


    %mrlm = max(M*(sum(A.*B,2)-A(:,1))); %MODIFIED MASS RANGE FOR 99%
    %enriched stuff
    mrlm = max(1.2*M*(sum(A.*B,2)-A(:,1)));

    %mrlm = max(2*M*(sum(A.*B,2)-A(:,1))); %avg mass - mono mass ORIGINAL



% Estimate the best mass range for the FFT by following [1]
if massRange == 0 % if different from 0, user already set the value by a PVP
    massRange = zeros(q,1);
    % Find out what would be the mass range required for this formula
    % Distribution mean for each element:
    dm = sum(B.*A,2);
    % Distribution second moment for every element:
    dsm = sum((bsxfun(@plus,dm,-A).^2).*B,2);
    % Distribution standard deviation for every element:
    % dsd = sqrt(dsm); % Not used, but left here for documentation purposes
    for i = 1:q
        % Standard deviation
        sigma = sqrt(M(i,:)*dsm);
        % Eq 9. in [1]
        massRange(i) = 1.2*K*sqrt((1+sigma^2));
    end
    % The worst case scenario, since the FFT for every element is to be
    % computed only once:
    massRange = max(massRange);
    % Safe lower-bound for very large molecules:
    massRange = max(massRange,mrlm);

    if mrlm > massRange
        fprintf ('%d', Nchange)
    end
    % Automatically correct for low resolution instruments:
    massRange = max(massRange,mrgp);

else % if the user set the range manually we hint an improper selection
    % Low resolution Gaussian pulse needs a larger range:
    if massRange < mrgp(1)
        warning(message('bioinfo:isotopicdist:InadequateMassRangeForLowRes', sprintf( '%g', mrgp )))
    end
    % Very large molecules needs a larger range:
    if massRange < mrlm
        warning(message('bioinfo:isotopicdist:InadequateMassRangeForLargeMolecules', sprintf( '%g', mrlm )))
    end
end

% Find what would be the size of the vector for the FFT transform:
 N = pow2(nextpow2(R*massRange));

% Find the left most istotpe for very element 

    Alm = min((A==0)./eps+A,[],2); 
    
    %Find the Average mass for each element
    if Cchange == 1 || Nchange == 1
        Avg = sum(A.*B,2);
    end
    % Eq 1. in [1]. Note that in this implementation masses are shifted to the
    % first isotope of every element (heterodyning), this permits using a
    % smaller vector for the FFT, the resulting final isotopic mass
    % distribution will be shifted to recover the correct masses:
    
    v = zeros(m,N);
    for i=1:m
        for j=1:n
            if A(i,j) > 0
                idx = round((A(i,j)-Alm(i))*R)+1;
                v(i,idx) = B(i,j); %put the intensities in B in a matrix B, spaced out enough to FFT later
            end
        end
    end

%     %PUT AN AVERAGE MASS HETERODYNE HERE?
%     if Cchange == 1 || Nchange == 1
%         offset = zeros(m,1);
% 
%         for i = 1:m
%             offset(i,1) = round((A(i,1)-Avg(i))*R);
%         end
%         offset = max(abs(offset));
% 
%         v = zeros(m,N+offset);
%         for i=1:m
%             for j=1:n
%                 if A(i,j) > 0
%                     idx = round((A(i,j)-Avg(i))*R)+1;
%                     v(i,idx+offset) = B(i,j); %put the intensities in B in a matrix B, spaced out enough to FFT later
%                     %idx = round((Avg(i)-A(i,j))*R)+1 - offset;
%                     
%                 end
%             end
%         end
% 
%     end



% Eq 2. in [1]
V = fft(v,N,2); %Calculate f(u) for each individual element

% Set the Gaussian pulse in the mass domain shifted by N*L (this factor is
% corrected later with the final shift)
sig = w.*R./sqrt(log(2));
s = exp(-((-N*L:(1/L-1)*N*L-1)./sig).^2);
S = fft(s); %fourier transform the peak shape function into u domain

% Allocate output variables
info = repmat(struct,q,1);
md = cell(q,1);
df = cell(q,1);

for i=1:q
    
    % Eq 5. in [1]
    F = ones(1,N);
    for j = 1:m
        F = F.*(V(j,:).^M(i,j)); %Generate f(u)chnos = f(u)c^a * f(u)h^b * f(u)n^c * f(u)o^d * f(u)s^e
    end
    
    % find the left most mass term ELI THIS IS THE HETERODYNE CORRECTION
    %Alm(1,1) = 13.00335;
    mzlm = M(i,:)*Alm;

%     if Cchange == 1 || Nchange == 1 %shift to the average mass if enriched
%         mzlm = M(i,:)*sum(A.*B,2);
%     end
    
    f = real(ifft(F.*S)); %Multiply f(u) times s(u) and ifft to convolute
    mz = mzlm +(0:N-1)/R-N*L/R; %Generate m/z values, add monoisotopic mass to correct for heterodyne
    
    P = mspeaks(mz,f,'denoising','false');
    
    % Remove any peak that appears left to the leftmost mass, these are
    % peaks that wrapped around due to the circular convolution of the fft
    P(P(:,1)<mzlm-w,:)=[];


    
    % Remove any peak below the noise threshold
    h = P(:,2)>max(P(:,2))/nt;
    if any(h)
       P = P(h,:);
    else
      P = P(P(:,2)==max(P(:,2)),:);
    end
    
    % Improve mass accuracy by interpolating peak maxima with a spline:
    for j = 1:size(P,1)
        [~,h] = min(abs(mz-P(j,1)));
        p = polyfit(mz(h+[-1 0 1])-mz(h),f(h+[-1 0 1])./f(h),2);
        xchat = roots(p*[2 0;0 1;0 0]);
        newInt = polyval(p,xchat)*f(h);
        P(j,:) = [xchat+mz(h),newInt];
    end
    
    % Normalize the area of the isotopic mass distribution
    P(:,2) = P(:,2)./sum(P(:,2));
    
    % Normalize the area of the isotopic density distribution
    f = f*(R/sum(f));
    
    if Csparse > 0 && Nsparse == 0
        syminM = (1:6)';
        ele = [ele(1), {'C13'} ele(2:end)]; %insert [13c] into the ele variable. Had to use C13 because of fieldname rules
    end

    if Nsparse > 0 && Csparse == 0
        syminM = (1:6)';
        ele = [ele(1:3), {'N15'} ele(4:end)]; %insert [15N] into the ele variable. Had to use N15 because of fieldname rules
    end

    if Csparse > 0 && Nsparse > 0
        syminM = (1:7)';
        ele = [ele(1), {'C13'}, ele(2:3), {'N15'}, ele(4:end)]; %insert C13 AND 15N into the ele variable
    end
    
    if (nargout>1) || (i==plotId)
        info(i).NominalMass = M(i,:)*round(A(:,1));
        info(i).MonoisotopicMass = M(i,:)*A(:,1);
        info(i).ObservedAverageMass = sum(mz.*f)./sum(f);
        info(i).CalculatedAverageMass = M(i,:)*sum(A.*B,2);
        info(i).MostAbundantMass = P(max(P(:,2))==P(:,2),1);
        info(i).Formula = cell2struct(num2cell(M(i,:)),ele(syminM),2);
    end
    
    if i==plotId
        figure('Tag','isotopicdist')

        h7 = plot(mz,f./(max(f)./max(P(:,2))));
        hold on
        h6 = plot(P(:,1),P(:,2),'rx');
        h1 = plot(info(i).NominalMass,0,'^g','markerfacecolor',[.1 1 .9]);
        h2 = plot(info(i).MonoisotopicMass,0,'^g','markerfacecolor',[1 .4 0]);
        h3 = plot(info(i).MostAbundantMass,max(P(:,2)),'or');
        h4 = plot(info(i).ObservedAverageMass,max(P(:,2)),'--w','visible','off');
        h5 = plot(info(i).CalculatedAverageMass,max(P(:,2)),'--w','visible','off');
        
        legend([h1 h2 h3 h4 h5 h6 h7], ...
            {sprintf('Nominal mass = %6.10g',info(i).NominalMass),...
            sprintf('Monoisotopic mass = %6.10g',info(i).MonoisotopicMass),...
            sprintf('Most Abundant mass = %6.10g',info(i).MostAbundantMass),...
            sprintf('Observed avg. mass = %6.10g',info(i).ObservedAverageMass),...
            sprintf('Calculated avg. mass = %6.10g',info(i).CalculatedAverageMass),...
            'Mass Distribution Function','Prob. Density Function (scaled)'},...
            'Location','NorthOutside');
        ts = sprintf(sprintf('%s_{%%g}',ele{syminM}),M(i,:));
        ts = regexprep(ts,'_\{1\}','');
        ts = regexprep(ts,'[A-Z][a-z]?_\{0\}','');
        title(ts)
    end
    
    if q==1
        md = P;
        if nargout>2
            df = [mz',f'];
        end
    else
        md{i} = P;
        if nargout>2
            df{i} = [mz',f'];
        end
    end
end

if nargout==0
    clear md
end

end

function M = peptide2chnos(p)
% Convert a peptide string (or a cellstr with peptides) to its empirical
% formula:

% Aminoacid element composition:
% Columns ordered by "C H N O S" and rows by int2aa
CHNOS = [3  7  1 2 0
    6  14 4 2 0
    4  8  2 3 0
    4  7  1 4 0
    3  7  1 2 1
    5  10 2 3 0
    5  9  1 4 0
    2  5  1 2 0
    6  9  3 2 0
    6  13 1 2 0
    6  13 1 2 0
    6  14 2 2 0
    5  11 1 2 1
    9  11 1 2 0
    5  9  1 2 0
    3  7  1 3 0
    4  9  1 3 0
    11 12 2 2 0
    9  11 1 3 0
    5  11 1 2 0 ];

try
    if ischar(p)
        M = accumarray(aa2int(p)',1,[20,1])'*CHNOS;
        M =  M - [0 2 0 1 0]*(numel(p)-1);
    else
        q = numel(p);
        M = zeros(q,5);
        for i = 1:q
            M(i,:) =  accumarray(aa2int(p{i})',1,[20,1])'*CHNOS;
            M(i,:) =  M(i,:) - [0 2 0 1 0]*(numel(p{i})-1);
        end
    end
catch ME
    if strcmp(ME.identifier,'MATLAB:accumarray:nonPosIntIndValues') || ...
            strcmp(ME.identifier,'MATLAB:accumarray:subsValuesTooLarge')
        error(message('bioinfo:isotopicdist:invalidAminoacidSymbols'))
    else
        rethrow(ME)
    end
end
end


function [A,B] = chnosIsotopicComposition(Cchange, Cabund, Nchange, Nabund, Csparse, Nsparse)

%check to see if the user asked for sparse labeling AND uniform labeling.
%If they did, politely tell them to fuck off
if (Cchange == 1 && (Csparse > 0 || Nsparse > 0)) || (Nchange == 1 && (Csparse > 0 || Nsparse > 0))
    error('HEY! THATS ILLEGAL! You can''t do uniform labeling AND sparse labeling at the same time yet! Pick one bucko ')
end


% Isotopic composition for the five basic elements in the amonoacids
A = [12           13.0033548378  14.003241988  0 0 0 0 0 0 0
    1.0078250321  2.014101778    3.0160492675  0 0 0 0 0 0 0
    14.0030740052 15.0001088984  0             0 0 0 0 0 0 0
    15.9949146221 16.9991315     17.9991604    0 0 0 0 0 0 0
    31.97207069   32.9714585     33.96786683   35.96708088 0 0 0 0 0 0];

B = [98.93   1.07    0     0     0 0 0 0 0 0
    99.9885  0.0115  0     0     0 0 0 0 0 0
    99.632   0.368   0     0     0 0 0 0 0 0
    99.757   0.038   0.205 0     0 0 0 0 0 0
    94.93    0.76    4.29  0.02  0 0 0 0 0 0];

if Cchange == 1
    %replace Carbon isotope abundances with user input abundances
    B(1, 1:2) = Cabund; 

end

if Nchange == 1
    %replace Nitrogen isotope abundances with user input abundances
    B(3, 1:2) = Nabund;

end

if Csparse > 0 && Nsparse == 0 %Add sparse labeled 13C
    A = [A(1,:); [13.0033548378 0 0 0 0 0 0 0 0 0]; A(2:end, :)]; %insert a row with the mass of 13C
    B = [B(1,:); [100 0 0 0 0 0 0 0 0 0]; B(2:end, :)]; %insert a row with 100% abundance for 13C

end

if Nsparse > 0 && Csparse == 0 %Add sparse labeled 15N
    A = [A(1:3,:); [15.0001088984 0 0 0 0 0 0 0 0 0]; A(4:end, :)];  %insert a row with the mass of 15N
    B = [B(1:3,:); [100 0 0 0 0 0 0 0 0 0]; B(4:end, :)]; %insert a row with 100% abundance for 15N
end

if Csparse > 0 && Nsparse > 0 %Add sparse labeled 13C AND 15N
    A = [A(1,:); [13.0033548378 0 0 0 0 0 0 0 0 0]; A(2:3,:); [15.0001088984 0 0 0 0 0 0 0 0 0]; A(4:end, :)]; %insert a row with the mass of 13C AND a row with the mass of 15N
    B = [B(1,:); [100 0 0 0 0 0 0 0 0 0]; B(2:3,:); [100 0 0 0 0 0 0 0 0 0]; B(4:end, :)]; 
end

end

function [M,A,B,syminM] = struct2formula(st,ele, Cchange, Cabund, Csparse, Nchange, Nabund, Nsparse)

c = struct2cell(st);
c(cellfun(@isempty,c)) = {0};
M = cell2mat(c');
ws = warning('off','bioinfo:seqmatch:StringNotFound');
syminM = seqmatch(fieldnames(st),ele,'exact',true);
warning(ws);

% element's isotopes and their natural abundance (Ref:
% http://physics.nist.gov/PhysRefData/Compositions/index.html)

elements = [...
    12         13.0033548378000         14.0032419880000         0      0      0      0      0      0      0            98.9300000000000         1.07000000000000         0      0      0      0      0      0      0      0
    1.00782503210000         2.01410177800000         3.01604926750000         0      0      0      0      0            0      0      99.9885000000000         0.0115000000000000     0      0      0      0      0      0      0            0
    14.0030740052000         15.0001088984000         0      0      0      0      0      0      0      0            99.6320000000000         0.368000000000000       0      0      0      0      0      0      0      0
    15.9949146221000         16.9991315000000         17.9991604000000         0      0      0      0      0      0            0      99.7570000000000         0.0380000000000000     0.205000000000000       0      0      0      0      0            0      0
    31.9720706900000         32.9714585000000         33.9678668300000         35.9670808800000         0      0      0            0      0      0      94.9300000000000         0.760000000000000       4.29000000000000         0.0200000000000000            0      0      0      0      0      0
    4.00260324970000         3.01602930970000         0      0      0      0      0      0      0      0            99.9998630000000         0.000137000000000000  0      0      0      0      0      0      0      0
    7.01600400000000         6.01512230000000         0      0      0      0      0      0      0      0            92.4100000000000         7.59000000000000         0      0      0      0      0      0      0      0
    9.01218210000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    11.0093055000000         10.0129370000000         0      0      0      0      0      0      0      0            80.1000000000000         19.9000000000000         0      0      0      0      0      0      0      0
    18.9984032000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    19.9924401759000         20.9938467400000         21.9913855100000         0      0      0      0      0      0            0      90.4800000000000         0.270000000000000       9.25000000000000         0      0      0      0      0            0      0
    22.9897696700000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    23.9850419000000         24.9858370200000         25.9825930400000         0      0      0      0      0      0            0      78.9900000000000         10         11.0100000000000         0      0      0      0      0      0      0
    26.9815384400000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    27.9769265327000         28.9764947200000         29.9737702200000         0      0      0      0      0      0            0      92.2297000000000         4.68320000000000         3.08720000000000         0      0      0      0      0            0      0
    30.9737615100000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    34.9688527100000         36.9659026000000         0      0      0      0      0      0      0      0            75.7800000000000         24.2200000000000         0      0      0      0      0      0      0      0
    39.9623831230000         35.9675462800000         37.9627322000000         0      0      0      0      0      0            0      99.6003000000000         0.336500000000000       0.0632000000000000     0      0      0      0      0            0      0
    38.9637069000000         39.9639986700000         40.9618259700000         0      0      0      0      0      0            0      93.2581000000000         0.0117000000000000     6.73020000000000         0      0      0      0      0            0      0
    39.9625912000000         41.9586183000000         42.9587668000000         43.9554811000000         45.9536928000000            47.9525340000000         0      0      0      0      96.9410000000000         0.647000000000000            0.135000000000000       2.08600000000000         0.00400000000000000    0.187000000000000       0      0      0            0
    44.9559102000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    47.9479471000000         45.9526295000000         46.9517638000000         48.9478708000000         49.9447921000000            0      0      0      0      0      73.7200000000000         8.25000000000000         7.44000000000000            5.41000000000000         5.18000000000000         0      0      0      0      0
    50.9439637000000         49.9471628000000         0      0      0      0      0      0      0      0            99.7500000000000         0.250000000000000       0      0      0      0      0      0      0      0
    51.9405119000000         49.9460496000000         52.9406538000000         53.9388849000000         0      0      0            0      0      0      83.7890000000000         4.34500000000000         9.50100000000000         2.36500000000000            0      0      0      0      0      0
    54.9380496000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    55.9349421000000         53.9396148000000         56.9353987000000         57.9332805000000         0      0      0            0      0      0      91.7540000000000         5.84500000000000         2.11900000000000         0.282000000000000            0      0      0      0      0      0
    58.9332002000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    57.9353479000000         59.9307906000000         60.9310604000000         61.9283488000000         63.9279696000000            0      0      0      0      0      68.0769000000000         26.2231000000000         1.13990000000000            3.63450000000000         0.925600000000000       0      0      0      0      0
    62.9296011000000         64.9277937000000         0      0      0      0      0      0      0      0            69.1700000000000         30.8300000000000         0      0      0      0      0      0      0      0
    63.9291466000000         65.9260368000000         66.9271309000000         67.9248476000000         69.9253250000000            0      0      0      0      0      48.6300000000000         27.9000000000000         4.10000000000000            18.7500000000000         0.620000000000000       0      0      0      0      0
    68.9255810000000         70.9247050000000         0      0      0      0      0      0      0      0            60.1080000000000         39.8920000000000         0      0      0      0      0      0      0      0
    73.9211782000000         69.9242504000000         71.9220762000000         72.9234594000000         75.9214027000000            0      0      0      0      0      36.2800000000000         20.8400000000000         27.5400000000000            7.73000000000000         7.61000000000000         0      0      0      0      0
    74.9215964000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    79.9165218000000         73.9224766000000         75.9192141000000         76.9199146000000         77.9173095000000            81.9167000000000         0      0      0      0      49.6100000000000         0.890000000000000            9.37000000000000         7.63000000000000         23.7700000000000         8.73000000000000         0      0      0            0
    78.9183376000000         80.9162910000000         0      0      0      0      0      0      0      0            50.6900000000000         49.3100000000000         0      0      0      0      0      0      0      0
    83.9115070000000         77.9203860000000         79.9163780000000         81.9134846000000         82.9141360000000            85.9106103000000         0      0      0      0      57         0.350000000000000       2.28000000000000            11.5800000000000         11.4900000000000         17.3000000000000         0      0      0      0
    84.9117893000000         86.9091835000000         0      0      0      0      0      0      0      0            72.1700000000000         27.8300000000000         0      0      0      0      0      0      0      0
    87.9056143000000         83.9134250000000         85.9092624000000         86.9088793000000         0      0      0            0      0      0      82.5800000000000         0.560000000000000       9.86000000000000         7          0      0            0      0      0      0
    88.9058479000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    89.9047037000000         90.9056450000000         91.9050401000000         93.9063158000000         95.9082760000000            0      0      0      0      0      51.4500000000000         11.2200000000000         17.1500000000000            17.3800000000000         2.80000000000000         0      0      0      0      0
    92.9063775000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    97.9054078000000         91.9068100000000         93.9050876000000         94.9058415000000         95.9046789000000            96.9060210000000         99.9074770000000         0      0      0      24.1300000000000         14.8400000000000            9.25000000000000         15.9200000000000         16.6800000000000         9.55000000000000         9.63000000000000            0      0      0
    96.9063650000000         97.9072160000000         98.9062546000000         0      0      0      0      0      0            0      0      0      0      0      0      0      0      0      0      0
    101.904349500000         95.9075980000000         97.9052870000000         98.9059393000000         99.9042197000000            100.905582200000         103.905430000000         0      0      0      31.5500000000000         5.54000000000000            1.87000000000000         12.7600000000000         12.6000000000000         17.0600000000000         18.6200000000000            0      0      0
    102.905504000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    105.903483000000         101.905608000000         103.904035000000         104.905084000000         107.903894000000            109.905152000000         0      0      0      0      27.3300000000000         1.02000000000000            11.1400000000000         22.3300000000000         26.4600000000000         11.7200000000000         0      0      0            0
    106.905093000000         108.904756000000         0      0      0      0      0      0      0      0            51.8390000000000         48.1610000000000         0      0      0      0      0      0      0      0
    113.903358100000         105.906458000000         107.904183000000         109.903006000000         110.904182000000            111.902757200000         112.904400900000         115.904755000000         0      0      28.7300000000000            1.25000000000000         0.890000000000000       12.4900000000000         12.8000000000000         24.1300000000000            12.2200000000000         7.49000000000000         0      0
    114.903878000000         112.904061000000         0      0      0      0      0      0      0      0            95.7100000000000         4.29000000000000         0      0      0      0      0      0      0      0
    119.902196600000         111.904821000000         113.902782000000         114.903346000000         115.901744000000            116.902954000000         117.901606000000         118.903309000000         121.903440100000         123.905274600000            32.5800000000000         0.970000000000000       0.660000000000000       0.340000000000000       14.5400000000000            7.68000000000000         24.2200000000000         8.59000000000000         4.63000000000000         5.79000000000000
    120.903818000000         122.904215700000         0      0      0      0      0      0      0      0            57.2100000000000         42.7900000000000         0      0      0      0      0      0      0      0
    129.906222800000         119.904020000000         121.903047100000         122.904273000000         123.902819500000            124.904424700000         125.903305500000         127.904461400000         0      0      34.0800000000000            0.0900000000000000     2.55000000000000         0.890000000000000       4.74000000000000         7.07000000000000            18.8400000000000         31.7400000000000         0      0
    126.904468000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    131.904154500000         123.905895800000         125.904269000000         127.903530400000         128.904779500000            129.903507900000         130.905081900000         133.905394500000         135.907220000000         0            26.8900000000000         0.0900000000000000     0.0900000000000000     1.92000000000000         26.4400000000000            4.08000000000000         21.1800000000000         10.4400000000000         8.87000000000000         0
    132.905447000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    137.905241000000         129.906310000000         131.905056000000         133.904503000000         134.905683000000            135.904570000000         136.905821000000         0      0      0      71.6980000000000         0.106000000000000            0.101000000000000       2.41700000000000         6.59200000000000         7.85400000000000         11.2320000000000            0      0      0
    138.906348000000         137.907107000000         0      0      0      0      0      0      0      0            99.9100000000000         0.0900000000000000     0      0      0      0      0      0      0      0
    139.905434000000         135.907140000000         137.905986000000         141.909240000000         0      0      0            0      0      0      88.4500000000000         0.185000000000000       0.251000000000000       11.1140000000000            0      0      0      0      0      0
    140.907648000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    141.907719000000         142.909810000000         143.910083000000         144.912569000000         145.913112000000            147.916889000000         149.920887000000         0      0      0      27.2000000000000         12.2000000000000            23.8000000000000         8.30000000000000         17.2000000000000         5.70000000000000         5.60000000000000            0      0      0
    144.912744000000         146.915134000000         0      0      0      0      0      0      0      0      0            0      0      0      0      0      0      0      0      0
    151.919728000000         143.911995000000         146.914893000000         147.914818000000         148.917180000000            149.917271000000         153.922205000000         0      0      0      26.7500000000000         3.07000000000000            14.9900000000000         11.2400000000000         13.8200000000000         7.38000000000000         22.7500000000000            0      0      0
    152.921226000000         150.919846000000         0      0      0      0      0      0      0      0            52.1900000000000         47.8100000000000         0      0      0      0      0      0      0      0
    157.924101000000         151.919788000000         153.920862000000         154.922619000000         155.922120000000            156.923957000000         159.927051000000         0      0      0      24.8400000000000         0.200000000000000            2.18000000000000         14.8000000000000         20.4700000000000         15.6500000000000         21.8600000000000            0      0      0
    158.925343000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    163.929171000000         155.924278000000         157.924405000000         159.925194000000         160.926930000000            161.926795000000         162.928728000000         0      0      0      28.1800000000000         0.0600000000000000            0.100000000000000       2.34000000000000         18.9100000000000         25.5100000000000         24.9000000000000            0      0      0
    164.930319000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    165.930290000000         161.928775000000         163.929197000000         167.932368000000         166.932045000000            169.935460000000         0      0      0      0      33.6100000000000         0.140000000000000            1.61000000000000         26.7800000000000         22.9300000000000         14.9300000000000         0      0      0            0
    168.934211000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    173.938858100000         167.933894000000         169.934759000000         170.936322000000         171.936377700000            172.938206800000         175.942568000000         0      0      0      31.8300000000000         0.130000000000000            3.04000000000000         14.2800000000000         21.8300000000000         16.1300000000000         12.7600000000000            0      0      0
    174.940767900000         175.942682400000         0      0      0      0      0      0      0      0            97.4100000000000         2.59000000000000         0      0      0      0      0      0      0      0
    179.946548800000         173.940040000000         175.941401800000         176.943220000000         177.943697700000            178.945815100000         0      0      0      0      35.0800000000000         0.160000000000000            5.26000000000000         18.6000000000000         27.2800000000000         13.6200000000000         0      0      0            0
    180.947996000000         179.947466000000         0      0      0      0      0      0      0      0            99.9880000000000         0.0120000000000000     0      0      0      0      0      0      0      0
    183.950932600000         179.946706000000         181.948206000000         182.950224500000         185.954362000000            0      0      0      0      0      30.6400000000000         0.120000000000000       26.5000000000000            14.3100000000000         28.4300000000000         0      0      0      0      0
    186.955750800000         184.952955700000         0      0      0      0      0      0      0      0            62.6000000000000         37.4000000000000         0      0      0      0      0      0      0      0
    191.961479000000         183.952491000000         185.953838000000         186.955747900000         187.955836000000            188.958144900000         189.958445000000         0      0      0      40.7800000000000         0.0200000000000000            1.59000000000000         1.96000000000000         13.2400000000000         16.1500000000000         26.2600000000000            0      0      0
    192.962924000000         190.960591000000         0      0      0      0      0      0      0      0            62.7000000000000         37.3000000000000         0      0      0      0      0      0      0      0
    194.964774000000         189.959930000000         191.961035000000         193.962664000000         195.964935000000            197.967876000000         0      0      0      0      33.8320000000000         0.0140000000000000            0.782000000000000       32.9670000000000         25.2420000000000         7.16300000000000         0      0      0            0
    196.966552000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    201.970626000000         195.965815000000         197.966752000000         198.968262000000         199.968309000000            200.970285000000         203.973476000000         0      0      0      29.8600000000000         0.150000000000000            9.97000000000000         16.8700000000000         23.1000000000000         13.1800000000000         6.87000000000000            0      0      0
    204.974412000000         202.972329000000         0      0      0      0      0      0      0      0            70.4760000000000         29.5240000000000         0      0      0      0      0      0      0      0
    207.976636000000         203.973029000000         205.974449000000         206.975881000000         0      0      0            0      0      0      52.4000000000000         1.40000000000000         24.1000000000000         22.1000000000000            0      0      0      0      0      0
    208.980383000000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    208.982416000000         209.982857000000         0      0      0      0      0      0      0      0      0            0      0      0      0      0      0      0      0      0
    209.987131000000         210.987481000000         0      0      0      0      0      0      0      0      0            0      0      0      0      0      0      0      0      0
    210.990585000000         220.011384100000         222.017570500000         0      0      0      0      0      0            0      0      0      0      0      0      0      0      0      0      0
    223.019730700000         0      0      0      0      0      0      0      0      0      0      0      0            0      0      0      0      0      0      0
    223.018497000000         224.020202000000         226.025402600000         228.031064100000         0      0      0            0      0      0      0      0      0      0      0      0      0      0      0      0
    227.027747000000         0      0      0      0      0      0      0      0      0      0      0      0            0      0      0      0      0      0      0
    230.033126600000         232.038050400000         0      0      0      0      0      0      0      0      100            0      0      0      0      0      0      0      0      0
    231.035878900000         0      0      0      0      0      0      0      0      0      100       0      0            0      0      0      0      0      0      0
    238.050782600000         233.039628000000         234.040945600000         235.043923100000         236.045561900000            0      0      0      0      0      99.2745000000000         0      0.00550000000000000    0.720000000000000            0      0      0      0      0      0];

%PUT THE C CHANGE AND N CHANGE STUFF HERE
if Cchange == 1
    %replace Carbon isotope abundances with user input abundances
    elements(1, 11:12) = Cabund; 

end

if Nchange == 1
    %replace Nitrogen isotope abundances with user input abundances
    elements(3, 11:12) = Nabund;

end

A = elements(syminM,1:10);
B = elements(syminM,11:20);

if Csparse > 0 && Nsparse == 0 %Add sparse labeled 13C
    A = [A(1,:); [13.0033548378 0 0 0 0 0 0 0 0 0]; A(2:end, :)]; %insert a row with the mass of 13C
    B = [B(1,:); [100 0 0 0 0 0 0 0 0 0]; B(2:end, :)]; %insert a row with 100% abundance for 13C

end

if Nsparse > 0 && Csparse == 0 %Add sparse labeled 15N
    A = [A(1:3,:); [15.0001088984 0 0 0 0 0 0 0 0 0]; A(4:end, :)];  %insert a row with the mass of 15N
    B = [B(1:3,:); [100 0 0 0 0 0 0 0 0 0]; B(4:end, :)]; %insert a row with 100% abundance for 15N
end

if Csparse > 0 && Nsparse > 0 %Add sparse labeled 13C AND 15N
    A = [A(1,:); [13.0033548378 0 0 0 0 0 0 0 0 0]; A(2:3,:); [15.0001088984 0 0 0 0 0 0 0 0 0]; A(4:end, :)]; %insert a row with the mass of 13C AND a row with the mass of 15N
    B = [B(1,:); [100 0 0 0 0 0 0 0 0 0]; B(2:3,:); [100 0 0 0 0 0 0 0 0 0]; B(4:end, :)]; 
end

end

function [w, massRange, L, K, R, nt, plotId,Mnter,Anter,Bnter,symnterm,Mcter,Acter,Bcter,symcterm, Cchange, Cabund, Nchange, Nabund, Csparse, Nsparse] = parse_inputs(ele,nout,varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(numel(varargin),2) == 1
    error(message('bioinfo:isotopicdist:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'resolution','fftrange','fftlocation','fftresolution','noisethreshold','showplot','nterminal','cterminal','charge', 'C13abund', 'N15abund', 'Csparse', 'Nsparse'};

% Setting up defaults
R = 1000; % Isotopic resolution: number of samples per Dalton
K = 10;   % Small-large molecules smooth factor used in Eq 9. Roockwood 1995
L = 1/16; % N*L is the fraction that the Gaussian pulse is shifted within the
% mass range, this allows to safely reduce the resolution and
% avoid that the left tail appears on the right of the spectra
w = 1/8;  % Width of the Gaussian pulse (in Daltons), usually this value should
% match the resolution of the instrument.
nt = 1e6; % noise removal threshold
massRange = 0; % to initialize (default value assigned in main function)
Cchange = 0; %binary switch to indicate if uniform 13C labeling will happen
Cabund = []; %Abundance values for 12C and 13C set by user
Nchange = 0;%binary switch to indicate if uniform 15N labeling will happen
Nabund = []; %abundance values for 14N and 15N set by user
Nsparse = 0; %number of sparse 15N to incorporate
Csparse = 0; %number of sparse 13C to incorporate

if nout == 0
    plotId = 1; 
else
    plotId = 0;
end
Mnter = [0 0 0 0 0]; % n-tern defaults to Amine
Mcter = [0 0 0 0 0]; % c-tern defaults to Free acid
[Anter,Bnter] = chnosIsotopicComposition(Cchange, Cabund, Nchange, Nabund, Csparse, Nsparse);
[Acter,Bcter] = chnosIsotopicComposition(Cchange, Cabund, Nchange, Nabund, Csparse, Nsparse);
symnterm = (1:5)';
symcterm = (1:5)';

% Loop over the values
for j=1:2:numel(varargin)-1
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % Resolution
            if (~isscalar(pval) || pval<=0)
                error(message('bioinfo:isotopicdist:IncorrectResolution'))
            end
            w = pval;
        case 2  % FFT Range
            if (~isscalar(pval) || pval<=0)
                error(message('bioinfo:isotopicdist:IncorrectFFTRange'))
            end
            massRange = pval;
        case 3  % FFTLocation
            if (~isscalar(pval) || pval<0 || pval>1)
                error(message('bioinfo:isotopicdist:IncorrectFFTLocation'))
            end
            L = pval;
        case 4  % FFT Resolution
            if (~isscalar(pval) || rem(pval,1) || pval<=0)
                error(message('bioinfo:isotopicdist:IncorrectFFTResolution'))
            end
            R = pval;
        case 5  % Noise Threshold
            if (~isscalar(pval) || pval<0)
                error(message('bioinfo:isotopicdist:IncorrectnoiseThreshold'))
            end
            nt = pval;
        case 6  % Show Plot
            if bioinfoprivate.opttf(pval,okargs{k},mfilename);
                if isnumeric(pval)
                    if isscalar(pval)
                        plotId = double(pval);
                    else
                        plotId = 1;
                        warning(message('bioinfo:isotopicdist:SPNoScalar'))
                    end
                else
                    plotId = 1;
                end
            else
                plotId = 0;
            end
        case 7  % N Terminal Modification
            if isstruct(pval)
                try
                    s = pval(1);
                    if isfield(s,'H') 
                        s.H = s.H - 1;
                    else
                        s.H = -1;
                    end
                    [Mnter,Anter,Bnter,symnterm] = struct2formula(s,ele);
                catch ME
                    if strcmp(ME.identifier,'MATLAB:cell2mat:MixedDataTypes')
                        error(message('bioinfo:isotopicdist:invalidFormulaNTerminus'))
                    else
                        rethrow(ME)
                    end
                end
            elseif ischar(pval)
                ValidNInputs = {'None', 'Amine', 'Formyl', 'Acetyl'}; 
                Nvalues = [0 -1 0 0 0; 0 0 0 0 0; 1 0 0 1 0; 2 2 0 1 0];
                h = find(strncmpi(pval,ValidNInputs,numel(pval)));
                if isempty(h)
                    error(message('bioinfo:isotopicdist:IncorrectNTerminus'))
                elseif numel(h)>1
                    error(message('bioinfo:isotopicdist:AmbiguousNTerminus'))
                end
                Mnter = Nvalues(h,:);
            end
        case 8  % C Terminal Modification
            if isstruct(pval)
                try                    
                    s = pval(1);
                    if isfield(s,'H') 
                        s.H = s.H - 1;
                    else
                        s.H = -1;
                    end
                    if isfield(s,'O') 
                        s.O = s.O - 1;
                    else
                        s.O = -1;
                    end
                    [Mcter,Acter,Bcter,symcterm] = struct2formula(s,ele);
                catch ME
                    if strcmp(ME.identifier,'MATLAB:cell2mat:MixedDataTypes')
                        error(message('bioinfo:isotopicdist:invalidFormulaCTerminus'))
                    else
                        rethrow(ME)
                    end
                end
            elseif ischar(pval)
                ValidCInputs = {'None', 'Free Acid', 'Amide'};
                Cvalues = [0 -1 0 -1 0; 0 0 0 0 0; 0 1 1 -1 0];
                h = find(strncmpi(pval,ValidCInputs,numel(pval)));
                if isempty(h)
                    error(message('bioinfo:isotopicdist:IncorrectCTerminus'))
%               elseif numel(h)>1
%                   error('bioinfo:isotopicdist:AmbiguousCTerminus', ...
%                       'Ambiguous C modification type. Options are ''None'', ''Free Acid'', or ''Amide''. Or use a structure to set a different formula modification.')
                end
                Mcter = Cvalues(h,:);
            end
        case 10 %modify C isotopic abundances
            
            if size(pval) ~= [1,1]
                error('Arrays are not accepted as input. Please input a double representing the %13C abundance')
            end
            Cchange = 1;
            Cabund = [100-pval, pval];
            

        case 11 %modify N isotopic abundance
            
            if size(pval) ~= [1,1]
                error('Arrays are not accepted as input. Please input a double representing the %15N abundance')
            end

            Nchange = 1;
            Nabund = [100-pval, pval];

        case 12 %Add sparse labeled 13C
            if size(pval) ~= [1,1]
                error('Arrays are not accepted as input. Please input an integer representing the %15N abundance')
            end
            if mod(pval,1) ~=0
                error('Only Integer values are accepted. Please enter an integer.')
            end

            Csparse = pval;

        case 13 %Add sparse labeled 15N
            if size(pval) ~= [1,1]
                error('Arrays are not accepted as input. Please input an integer representing the %15N abundance')
            end
            if mod(pval,1) ~=0
                error('Only Integer values are accepted. Please enter an integer.')
            end

            Nsparse = pval;

            
    end
end
end
