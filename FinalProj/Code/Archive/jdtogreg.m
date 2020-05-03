function [YYYY,MM,DD,hh,mm,ss] = jdtogreg(JD)

LMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
T_1900 = (JD - 2415019.5)/365.25;
YYYY = 1900 + floor(T_1900);

LeapYY = floor((YYYY-1900-1)*(0.25));
days = (JD-2415019.5)-((YYYY-1900)*(365)+ LeapYY);
if days < 1
    YYYY = YYYY - 1;
    LeapYY = floor((YYYY-1900-1)*(0.25));
    days = (JD-2415019.5)-((YYYY-1900)*(365)+ LeapYY);
else
end
if mod(YYYY,4) == 0
    LMonth(2) = 29;
else
end

dayofyr = floor(days);
cache = 31;
for idx = 2:12
    if cache + LMonth(idx) > dayofyr
        if dayofyr <= 31
            MM = 1;
            DD = dayofyr;
            break
        else
            MM = idx;
            DD = dayofyr - cache;
            break
        end
    else
        cache = cache + LMonth(idx);
        MM = idx;
        DD = dayofyr - cache;
    end
end
if DD == 0
    DD = 31;
    MM = MM - 1;
else
end
hhmmss = (days - dayofyr)*24;
hh = floor(hhmmss);
mm = floor((hhmmss-hh)*60);
ss = (hhmmss - hh - mm/60)*3600;

end