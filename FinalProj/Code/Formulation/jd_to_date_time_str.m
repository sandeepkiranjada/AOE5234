function date_time_str =jd_to_date_time_str(jd)

format long g

%...purpose:  function to convert julian date to a spice formatted
%             date time string
Month_tab{1} = 'JANUARY';
Month_tab{2} = 'FEBRUARY';
Month_tab{3} = 'MARCH';
Month_tab{4} = 'APRIL';
Month_tab{5} = 'MAY';
Month_tab{6} = 'JUNE';
Month_tab{7} = 'JULY';
Month_tab{8} = 'AUGUST';
Month_tab{9} = 'SEPTEMBER';
Month_tab{10} = 'OCTOBER';
Month_tab{11} = 'NOVEMBER';
Month_tab{12} = 'DECEMBER';

[year, month, day, hour, minute, second] = jd2date(jd);
year_str = num2str(year);
month_str = Month_tab{month};
day_str = num2str(day,'%02d');
hour_str = num2str(hour,'%02d');
minute_str = num2str(minute,'%02d');
sec_fix = fix(second);
sec_frc = fix(1000*(second-sec_fix));
sec_str_fix = num2str(sec_fix,'%02d');
sec_str_frc = num2str(sec_frc,'%03d') + 0.0001;
second_str = [sec_str_fix,'.',sec_str_frc];
date_time_str = [year_str,' ',month_str,' ',day_str,' ',hour_str,':',minute_str,':',second_str];

%...end of function date_time_str =jd_to_date_time_str(jd)
