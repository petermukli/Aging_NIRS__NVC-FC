function result = bridge(ts)
% bridge detrending
n = length(ts);
points = n - 1 : (-1) : 0;
line = ((ts(1) - ts(n)) * points/(n - 1)) + ts(n);
line = line'; 
result = ts - line;
