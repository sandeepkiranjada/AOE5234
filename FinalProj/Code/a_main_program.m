%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%            Main Program for Computing Errors between brdc and TLE Data             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                                                    %%%%%%%
%%%%%%%     Dependencies:                                                                  %%%%%%%
%%%%%%%                     Toolboxes:                                                     %%%%%%%
%%%%%%%                                 SCT Toolbox                                        %%%%%%%
%%%%%%%                                 constell Toolbox                                   %%%%%%%
%%%%%%%                     Functions:                                                     %%%%%%%
%%%%%%%                                main_brdc.m                                         %%%%%%%
%%%%%%%                                match_brdc2tle.m                                    %%%%%%%
%%%%%%%                                read_rinex_nav.m                                    %%%%%%%
%%%%%%%                                parsef.m                                            %%%%%%%
%%%%%%%                                                                                    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close('all'), tic

% for prn = [1:3 5:17 19:32]
for prn = 1
    
    clearvars -except prn
    model = 'SGP4';
%     model = 'SDP4';
    fprintf('\nProcessing PRN%02.0f\n',prn);
    
    fname = sprintf(['tle_prn' num2str(prn,'%02.0f') '.txt']);
    fid = fopen(fname);
    res={};
    while ~feof(fid)
        res{end+1,1} = fgetl(fid);
    end
    fclose(fid);
    no_of_lines = numel(res);
    kmax = no_of_lines/2;
    
    dataname = sprintf(['prn' num2str(prn,'%02.0f') '_data.mat']);
    seed = floor(kmax/3);
    range = seed:seed:kmax;
%     range = 1:5000;
%     range = 1:kmax;
    fprintf('Total entries to be processed: %.0f\n',length(range));
    
    %     %{
    
    % Retrieve UTC_time from TLE entries
    [~,~,jd_cache,~,~] = LoadNORAD(fname,model,kmax);
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
    
    % UTC to GPS time conversion
    [GPS_week,GPS_sec,~] = utc2gps(UTC_time);%,leap_sec); % uses lookup table of leap seconds
    GPS_time = [GPS_week,GPS_sec];%-leap_sec+1];
    
    %% Propagate broadcast data to TLE time
    yyyy = UTC_time(:,1);
    doy = day(datetime(UTC_time),'dayofyear'); % compute for fractional day of year
    t_start = GPS_time;
    popu = 1;   % 1 if brdc files need to be downloaded from CDDIS
                % 0 if brdc files are already populated
    brdc_ecef = main_brdc(prn,yyyy,doy,t_start,popu); % brdc [r,v] (meters)
    
    % Convert brdc data from ECEF to ECI using GPS_time
    [x_eci_cache,v_eci_cache] = ecef2eci(GPS_time,brdc_ecef(:,1:3),brdc_ecef(:,4:6));
    brdc_eci = [x_eci_cache,v_eci_cache];
    clear x_eci_cache v_eci_cache
    
    % Solve for brdc range in ECI and ECEF
    for idx = 1:length(range)
        if idx == 1
            brdc_eci_range = zeros(length(range),1);
            brdc_ecef_range = zeros(length(range),1);
        end
        brdc_eci_range(idx) = norm(brdc_eci(idx,1:3));
        brdc_ecef_range(idx) = norm(brdc_ecef(idx,1:3));
    end
    
    %% Convert TLE to [r,v] using SCT Toolbox
    [r,v,~,~,~] = LoadNORAD(fname,model,kmax);
    sct_eci_cache = [r' v']*1e3; % [r,v] (meters)
    sct_eci = sct_eci_cache(range,:);
    clear sct_eci_cache r v
    
    % Convert SCT data from ECI to ECEF using GPS_time
    [x_ecef,v_ecef] = eci2ecef(GPS_time,sct_eci(:,1:3),sct_eci(:,4:6));
    sct_ecef = [x_ecef,v_ecef]; % [r,v] (meters)
    clear x_ecef v_ecef
    
    % Solve for SCT range in ECI and ECEF
    for idx = 1:length(range)
        if idx == 1
            sct_eci_range = zeros(length(range),1);
            sct_ecef_range = zeros(length(range),1);
        end
        sct_eci_range(idx) = norm(sct_eci(idx,1:3));
        sct_ecef_range(idx) = norm(sct_ecef(idx,1:3));
    end
    
    % Get ECEF error vectors
    er_sct = sct_ecef-brdc_ecef;
    
    % Convert data from ECEF to Local Level (LL) coordinates
    [ll_vec_sct] = ecef2ll(er_sct(:,1:3), brdc_ecef(:,1:3), brdc_ecef(:,4:6));
    
    save(dataname);
    
    %}
    
    load(dataname,'UTC_time','range','prn','brdc_eci','brdc_eci_range','brdc_ecef','brdc_ecef_range', ...
        'sct_eci','sct_eci_range','sct_ecef','sct_ecef_range','ll_vec_sct');
    
    %% Plotting
    x0 = 5; y0 = 600;
    lx = 550; ly = 400;
    space1 = 560;
    
    %% Convert UTC time into MATLAB's serial date numbers
    t_sd = datenum(UTC_time);
    t_cd = datetime(t_sd,'ConvertFrom','datenum');
    
    %% Set colors
    blue = [0,0.447,0.741];
    red = [0.85,0.325,0.098];
    yellow = [0.929,0.694,0.125];
    violet = [0.494,0.184,0.556];
    green = [0.466,0.674,0.188];
    lblue = [0.301,0.745,0.933];
    dred = [0.635,0.078,0.184];
    
    %% Plot results from SCT's LoadNORAD
    PlotOpts.Color = yellow;
    PlotOpts.MarkerFaceColor = yellow;
    PlotOpts.MarkerSize = 2;
    PlotOpts4sct(1) = {PlotOpts};
    picsize = [x0+0*space1 y0 lx ly]; figure('rend','painters','pos',picsize)
    plot(t_sd,sct_eci(:,1)-brdc_eci(:,1),'o',PlotOpts4sct{1}); hold on; ylabel('[X-error]_{ECI} (m)');
    picsize = [x0+1*space1 y0 lx ly]; figure('rend','painters','pos',picsize)
    plot(t_sd,sct_eci(:,2)-brdc_eci(:,2),'o',PlotOpts4sct{1}); hold on; ylabel('[Y-error]_{ECI} (m)');
    picsize = [x0+2*space1 y0 lx ly]; figure('rend','painters','pos',picsize)
    plot(t_sd,sct_eci(:,3)-brdc_eci(:,3),'o',PlotOpts4sct{1}); hold on; ylabel('[Z-error]_{ECI} (m)');
    picsize = [x0+0*space1 y0-ly-150 lx ly]; figure('rend','painters','pos',picsize)
    figure(4), plot(t_sd,sct_eci_range-brdc_eci_range,'o',PlotOpts4sct{1}); hold on; ylabel('Range-error (m)')
    
    for idx = 1:4
        figure(idx)
        if t_sd(length(range)) < t_sd(1)
            tmax = datenum([2018 04 01]);
        else
            tmax = t_sd(length(range));
        end
        axis([t_sd(1) tmax -5000 5000]);
        legendname = sprintf(['Errors for PRN' num2str(prn,'%02.0f')]);
        legend(legendname);
    end
    
    %% for SCT
    picsize = [x0+1*space1 y0-ly-150 lx+200 ly+100]; figure('rend','painters','pos',picsize)
    plot(t_sd,ll_vec_sct(:,1),'o','MarkerSize',2,'MarkerEdgeColor',green,'MarkerFaceColor',green); hold on; grid on
    plot(t_sd,ll_vec_sct(:,2),'o','MarkerSize',2,'MarkerEdgeColor',blue,'MarkerFaceColor',blue); hold on; grid on
    plot(t_sd,ll_vec_sct(:,3),'o','MarkerSize',2,'MarkerEdgeColor',red,'MarkerFaceColor',red); hold on; grid on
    axis([t_sd(1) tmax -5000 5000]);
    titlename = sprintf(['Errors for PRN' num2str(prn,'%02.0f')]);
    title(titlename);
    legend('In-track','Cross-track','Radial');
    ylabel('Error (m)');
    xlabel('Elapsed time');
    
    %% Save images
    q = get(0,'children');
    fig_date_time = datestr(now,'dd-mmm-yyyy HH-MM');
    fig_date = fig_date_time(1:11);
    fig_time = fig_date_time(13:end);
    for f = 1:length(q)
        figure(f);
        grid on
        set(gca,'FontSize',12);
        datetick('x','ddmmmyyyy','keeplimits');
    end
    for f = 1:length(q)
        figure(f);
        %         print(q(f),['PRN' num2str(prn,'%02.0f') ' ' fig_date ' ' fig_time ' Figure ' num2str(length(q)+1-f)],'-dpng','-r300');
        %         print(q(f),['PRN' num2str(prn,'%02.0f') ' ' fig_date ' Figure ' num2str(length(q)+1-f)],'-dpng','-r300');
        %         print(q(f),['PRN' num2str(prn,'%02.0f') ' ' model '
        %         Figure ' num2str(length(q)+1-f)],'-dpng','-r300'); % With model type (SGP4 vs SDP4)
        print(q(f),['PRN' num2str(prn,'%02.0f') ' Figure ' num2str(length(q)+1-f)],'-dpng','-r300');
    end
    
    disp(' ');
    toc
    
%     close('all');
    fclose('all');
end