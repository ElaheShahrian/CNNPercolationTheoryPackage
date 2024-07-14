clear all
close all
clc
p = 3 ;
for iii = 1:5

    %General Percolate matrix Permeability
    exportmatrix = percolation_matrix_generator(64,0, 1, 80,p);
    
    % Preparation data for inputing to the eclipse, 1)trasnform
    transformed_exportmatrix = transformed_inc_file(exportmatrix,'PERMX','PERMX.INC');
    
    %General Percolate matrix Porosity
    exportporosity = percolation_porosity_generator(exportmatrix ,1, 0.0001, 1);
    
    % Preparation data for inputing to the eclipse, 1)trasnform
    transformed_exportporosity = transformed_inc_file(exportporosity,'PORO','PORO.INC');
    
    % sitation percolation
    % pOSites = sitation(64,exportmatrix,p,iii);
    
    % graphic, binary image(gray)
    figure_percolation(exportmatrix,iii,p)
    
    % run(bat) eclipse
    % run_eclipse()
    
    
    % load RSM
    DataS = load_RSM();
    
    % prepare DataS
    Inj_BPR = get_property_RSM(DataS,1,64,1,260,'BPR') ;
    Pro_BPR = get_property_RSM(DataS,65,128,1,260,'BPR') ;
    
    % darcy calculation
    dx=1000;
    dy=1000;
    dz=1000;
    visco = 0.633556;
    kx(iii) = darcy_equation(mean(Pro_BPR)-mean(Inj_BPR),DataS.FOPR(end),visco,64,dx,dy,dz) ;
end

xlswrite(['kx' num2str(p) '0.xlsx'],kx');


%%
function RN = percolation_matrix_generator(L,defaultvalue, othervalue, boundryvalue,p)
    %L= size of model, defaultvalue=1, othervalue=0, boundryvalue=80, p=probability percolation 
    RN =1+(10-1) * rand(L,L);
    a(1:L,1)=[boundryvalue];
    b(1:L,1)=[boundryvalue];
    RN(RN<p)=othervalue;
    RN(RN>=p)=defaultvalue;
    RN =[a RN b];
end


function transformed_RN = transformed_inc_file(RN,keyword,file_name) 
    transformed_RN = fliplr(transpose(RN));
    
    string = cell(1,1);
    string{1,1} = keyword;
    transformed_RN_reshape=reshape(transformed_RN,[numel(transformed_RN),1]); % reshape transformed_cc
    transformed_RN_str=num2str(transformed_RN_reshape);
    string(2:1+numel(transformed_RN_reshape),1) = cellstr(transformed_RN_str) ;
    string{numel(string)+1,1} = '/' ;
    writecell(string,file_name,'FileType','text') ;

end

function porosity = percolation_porosity_generator(RN,defaultvalue, othervalue, boundryvalue)
     %RN=input, defaultvalue=1, othervalue=0.0001, boundryvalue=1
    a = RN(1:end,1) ;
    a(:) = boundryvalue;
    b = RN(1:end,end);
    b(:) = boundryvalue;
    c = RN(1:end,1+1:end-1);
    u = unique(c);
    c(c==u(1)) = othervalue;
    c(c==u(2)) = defaultvalue;
    porosity = [a c b];

end


function pOSites = sitation(L,exportmatrix,p,iii)
    a = exportmatrix(:,1+1:end-1);
    nSites  = numel(a);
    nOSites = nnz(logical(a==1));
    pOSites = nOSites/nSites ;
    


    % Need to fix
    [cls, numC] = bwlabel(exportmatrix,4);
    nC = zeros(numC,1);
    for c = 1 : numC
        nC(c) = length(cls(cls == c));
    end
    % s   counter from 1 to max number of occupied sites in a cluster
    s = 1:max(nC);
    % nS  cluster size distribution
    nS = zeros(max(nC),1);
    for c = s
        nS(c) = length(nC(nC == c));
    end
    % Probability that an occupied site chosen at random is part of an s-site cluster
    probROS = nS ./ sum(nS);
    % mean clustert size = mean number of occupied sites in cluster
    meanClusterSize = mean(nC);
    S = sum(s.*probROS');
    % CHECKING FOR SPANNING CLUSTERS
    %    X spanning
    isSpanX = zeros(numC,1); isSpanY = isSpanX;
    %isSpanX = zeros(numC,1); isSpanY = isSpanX;
    for c = 1 : numC
        col1 = isempty(find(cls(cls(:,1)==c)));
        colN = isempty(find(cls(cls(:,end)==c)));
        col1N = col1 + colN;
        if col1N == 0
            isSpanX(c) = 1;
        end
        col1 = isempty(find(cls(cls(1,:)==c)));
        colN = isempty(find(cls(cls(end,:)==c)));
        col1N = col1 + colN;
        if col1N == 0
            isSpanY(c) = 1;
        end
    end
    disp('  ');
    fprintf('Linear dimension  L = %2.0f \n',L);
    fprintf('Total number of sites  nSties = %2.0f \n',nSites);
    fprintf('Probability of a site being occupied   p = %2.2f \n',p);
    disp('  ');
    fprintf('Number of occupied sites nOSites = %2.0f \n',nOSites);
    fprintf('Simulation: probability of a site being occupied pOSites = %2.2f \n',pOSites);
    fprintf('Number of clusters  numC = %2.0f \n',numC);
    disp('Cluster Size: number of occupied sites in each cluster  nC ');
    fprintf('   %2.0f ',nC');
    disp('  ')
    disp('Cluster Size Distribution  ');
    disp('    s      nS    probROS');
    for c = s
        fprintf('   %2.0f     %2.0f     %2.3f   \n',s(c), nS(c), probROS(c));
    end
    disp('  ');
    disp('Mean Cluster Size  S:  sites / clusters / prob   ')
    fprintf('S = mean(nC) = %2.3f \n',meanClusterSize);
    fprintf('S = SUM(s .* probROS) = %2.3f \n',S);
    disp('  ');
    disp('Spanning clusters = 1')
    disp('cluster #   rows   cols');

    for c = 1: numC
        fprintf('     %2.0f      %u2.0f       %2.0f     \n',c, isSpanY(c), isSpanX(c));
    end

    for c = 1 : numC;
        [rN, cN] = find(cls == c);
        %       text(cN(1),rN(1),num2str(c,'%2.0f'),'Color','w');
    end


    % max cluster
    maxC= max(nC);
    P_big(iii,1)= maxC/nOSites;
    xlswrite(['P_big' num2str(p) '0.xlsx'],P_big);
    
end


function figure_percolation(RN_transformed,iii,p)
    
    % Figure
    f=figure(iii);
    set(f,'units','normalized','position',[0.1 0.52 0.23 0.32],'color','w');
    
    % Axes
    p_handle = imagesc(~RN_transformed);
    axis equal
    axis tight
    colormap gray
    set(gca,'XColor','none','YColor','none')

    if exist(['Figure' num2str(p) '0']) ~= 7  % 7 = folder
        mkdir(['Figure' num2str(p) '0'])  % make folder
    end
    filename = fullfile([cd '\Figure' num2str(p) '0'], ['figure' num2str(iii)]);
    
    % Way 1
    Frame = getframe(f);
    imwrite(Frame.cdata,[filename '.png'])
    
    % Way 2
    %     print(filename ,'-dpng');
    close(f);
    pause(2)
end


function run_eclipse()
   
    if exist([cd '\A'])== 7 % 7 = folder   2 = file
        rmdir([cd '\A'],'s') % remove folder
    end
    mkdir(['A']); % make folder
    
    % a simple way
%     BF= [cd '\A.data']; % base file
%     NF = [cd '\A\Sample.data']; % new file
%     copyfile(fullfile(BF), NF);

    % copy files to A folder
    copyfile(fullfile([cd '\A.data']), [cd '\A\Sample.data']);
    copyfile(fullfile([cd '\PERMX.INC']),[cd '\A\PERMX.INC']);
    copyfile(fullfile([cd '\PORO.INC']),[cd '\A\PORO.INC']);

    %make .bat file
    ID = fopen([cd '\A\File_Run.bat'],'w+');
    fprintf(ID,['eclrun eclipse %s\\A\\Sample.DATA \n'],cd);
    fclose(ID);

    %run bat file with eclipse
    system([cd '\A\File_Run.bat']); 
    clc
end
  
function DataS = load_RSM()
    if exist([cd '\A\Sample.RSM']) == 2
        IDr = fopen([cd '\A\Sample.RSM'],'r');
        
        n=1;
        while true
            n=n+1;
            RSMFile{n,:}=fgetl(IDr);
            if RSMFile{n,:}==-1
                End_row_num=n-1;
                break
            end
        end
        fclose(IDr);

        DataS = To_Get_RSM_File(RSMFile);


    end
end

function S = To_Get_RSM_File(RSMFile)
    % Find Head and End of Table
    T_index_H = [] ;
    T_index_E = [] ;
    H_index_Name = [] ;
    T_index_A = [] ;
    for i = 1:size(RSMFile,1)-1
        if ~isempty(RSMFile{i,:})
            if ~isempty([regexp(RSMFile{i,:},'TIME')])
                if isempty(T_index_H)
                    T_index_H = [T_index_H , i+4] ;
                    H_index_Name = [H_index_Name , i] ;
                else
                    T_index_H = [T_index_H , i+4] ;
                    T_index_E = [T_index_E ,T_index_H(end)-7];
                    H_index_Name = [H_index_Name , i] ;
                end
            end
        end
    end
    T_index_E = [T_index_E , size(RSMFile,1)-1] ;
    T_index_O = zeros(size(T_index_E)) ;
    for i = 1:length(T_index_E)
        for j = H_index_Name(i):T_index_H(i)
            if ~isempty([regexp(RSMFile{j,:},'*10**')])
                T_index_O(i) = 1 ;
            end
        end
    end
    T_index_H = T_index_H + T_index_O ;
    
    % Processing
    % P1:
    Header_Name = {} ;
    for i = 1:length(H_index_Name)
        splited_headers = strsplit(RSMFile{H_index_Name(i),:}) ;
    
        for j = 1:length(splited_headers)
            if ~isempty(splited_headers{j})
                Header_Name{1,end+1} = splited_headers{1,j} ;
            end
        end
    
    end
    
    %P2:
    for i = 1:length(T_index_H)
        AAA(:,i) = RSMFile(T_index_H(i):T_index_E(i),1) ;
    end
    
    %P3:
    for i = 1:size(AAA,1)
        Experssion = '' ;
        for j = 1:length(T_index_H)
            if j ~= length(T_index_H)
                separator = '' ;
            else
                separator = ' ' ;
            end
            Experssion = [Experssion , AAA{i,j} , separator] ;
        end
        Table{i,1} = Experssion ;
    end
    
    %P4:
    for i = 1:length(Table)
        Eval(i,:) = str2num(Table{i,:}) ;
    end
    
    %P5:
    New_Header_Name = CountHeaders(Header_Name) ;
    
    %P6:
    S = struct ;
    for i = 1:size(New_Header_Name,2)
        eval(['S.' New_Header_Name{i} ' = Eval(:,' num2str(i) ') ;']);
    end
end


function newdata = CountHeaders(data)
    newdata = data;
    counts = containers.Map();
    for i = 1:numel(data)
        if isKey(counts, data{i})
            counts(data{i}) = counts(data{i}) + 1;
            newdata{i} = [data{i}, num2str(counts(data{i}))];
        else
            counts(data{i}) = 1;
        end
    end
end


function Result = get_property_RSM(DataS,Start_Grid,End_Grid,Start_Row,End_Row,StringVar)
    try
        % StringVar is a vector in DataS like BPR
        for i = Start_Grid:End_Grid
            if i==1
                eval(['BPRT(:,i) = DataS.' StringVar ';'])
            else
                eval(['BPRT(:,i-Start_Grid+1) = DataS.' StringVar  num2str(i) ';'])
            end
        end
        Result = BPRT(Start_Row:End_Row,:);
    catch
        Result = [] ;
    end
end

function kx = darcy_equation(dp,q,visco,Length,dx,dy,dz)
    L = Length+2*dx;
    A = (dy*L)*(dz*1) ;
    k = -((q*visco*L)./(1.127e-3*A*dp));
    kx = mean(mean(k)) ;
end


