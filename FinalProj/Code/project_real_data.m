%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Real Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
flag_save = 1;

addpath(genpath('TLEs'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          Echo 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 'SGP4';
%     model = 'SDP4';

% noradID = '00049';      % Echo 1A (LEO)
noradID = '02253';      % PAGEOS-A (Polar)
% noradID = '02324';      % PasComSat/OV1-8 (LEO)
filename = sprintf([ noradID '.txt']);
fid = fopen(filename);
res={};
while ~feof(fid)
    res{end+1,1} = fgetl(fid);
end
fclose(fid);
no_of_lines = numel(res);
kmax = no_of_lines/2;

range = 1:kmax;
fprintf('Total entries to be processed: %.0f\n',length(range));

% Retrieve UTC_time from TLE entries
[~,~,jd_cache,~,~] = LoadNORAD(filename,model,kmax);
jd = jd_cache(range);
clear jd_cache
for idx = 1:length(jd)
    if idx == 1
        UTC_time = zeros(length(jd),6);
    end
    UTC_time_scratch = JD2Date(jd(idx));
    UTC_time(idx,1:6) = UTC_time_scratch;
end
clear jd UTC_time_scratch
UTC_datetime = datetime(UTC_time);

% Convert TLE to [r,v] using SCT Toolbox
[r,v,jD,coe,x] = LoadNORAD(filename,model,kmax);
rv_eci = [r' v']*1e3; % [r,v] (meters)


%% Convert TLE to classical orbital elements using SCT Toolbox
[el, name] = NORADToEl( [], [], [], model, filename );
% el            (:,6)   Elements vector [a,i,W,w,e,M]

%% Save data to .mat file
dataname = sprintf([ noradID '_data.mat']);
realdata = [UTC_time el rv_eci];
save(fullfile(pwd,'Data',dataname),'realdata');

% %% Sort data in chronological order 
% % (Some data (not all) retrieved from spacetrack is out of chronological
% % order.)
% [UTC_datetime,I] = sort(UTC_datetime);
% UTC_time = UTC_time(I,:);
% el = el(I,:);

%% Plotting
f = figure();
f.Units = 'centimeters';
f.Position = [1 1 55 10];

subplot(1,3,1);
plot(UTC_datetime,el(:,1),':ok','MarkerSize',2,'MarkerFaceColor','k'); grid on; set(gca,'FontSize',12);
title('Semi-major Axis');
ylabel('a [km]');

subplot(1,3,2);
plot(UTC_datetime,el(:,5),':ok','MarkerSize',2,'MarkerFaceColor','k'); grid on; set(gca,'FontSize',12);
title('Eccentricity');
ylabel('e');

subplot(1,3,3);
plot(UTC_datetime,rad2deg(el(:,2)),':ok','MarkerSize',2,'MarkerFaceColor','k'); grid on; set(gca,'FontSize',12);
title('Inclination');
ylabel('i [^\circ]');

%% Save images
if flag_save
    q = get(0,'children');
    for f = 1:length(q)
        figure(f);
        set(gca,'FontSize',12);
        figname = sprintf([ 'Real Data ' noradID ' Figure ' num2str(length(q)+1-f)]);
        print(q(f),fullfile(pwd,'Figures',figname),'-dpng','-r300');
    end
end

%{
% Define filename
filename1 = 'Tiangong-1mod.txt';

% Load data from filename
fid1 = fopen(filename1);
res={};
while ~feof(fid1)
    res{end+1,1} = fgetl(fid1);
end

no_of_lines = numel(res);
kmax = no_of_lines;
range = 1:kmax;

k = 0;
for idx = 1:kmax
    str = res{idx,1};
    
    regexp(str,' ','split');
    
    
    pattern = 'EC=';
    keep_idx(idx,1) = contains(str,pattern);
    if keep_idx(idx,1) == 1
        k = k+1;
        datasinglecell{k,1} = res{idx,1};
    end
end
datacells = regexp(datasinglecell,' ','split');
datacells = vertcat(datacells{:});
kindex = 1:length(datacells);
datacells = regexp(datacells(:,2),',','split');

for idx = 1:length(datacells)
    dataline = datacells{idx,1};
    dataline(1) = [];
    no_of_meas = length(dataline)/10;
    range = 10:10:length(dataline);
    dataline(range) = [];
    A = reshape(dataline,9,[])';
    datacells{idx,1} = cellfun(@str2num, A);
    datacells{idx,1} = [datacells{idx,1} kindex(idx)*ones(length(A),1)];
end

datatable = cell2mat(datacells);
fclose(fid1);
save(filename1a,'datatable'); fprintf('\nmatfile saved.\n');


%}

