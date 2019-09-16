% Katy RIojas
% Trim to start and end points of insertion
% 9/8/19
function [p1_tvec,p1_vec,p2_tvec,p2_vec,p3_tvec,p3_vec,...
    p4_tvec,p4_vec,c1_tvec,c1_vec,c2_tvec,c2_vec,c3_tvec,c3_vec,releaseTimes] = ...
    trimManualTrials(data_pman1,data_pman2,data_pman3,data_pman4,...
    data_cman1,data_cman2,data_cman3)

%Find trim points for time vectors
LL = 5;
Fthresh1 = 50;
Fthresh2 = 125;
pman1_peak_idx = find(data_pman1.Fmag>Fthresh1,1);
tPP_pman1_idx = find(data_pman1.Fmag(pman1_peak_idx:end)<LL,1);
tPP_pman1 = data_pman1.time(tPP_pman1_idx+pman1_peak_idx);

% Video 1:
pman1_vid_tap = 5.93; %[s]
pman1_vid_start = 10.78; % [s]
pman1_vid_depth = 55.05; % [s]
pman1_vid_release = 73.78; % [s]

% Calculations
diff_pman1_start = pman1_vid_start - pman1_vid_tap;
diff_pman1_end = pman1_vid_depth - pman1_vid_start;
diff_pman1_release = pman1_vid_release - pman1_vid_start;
p1_tstart = tPP_pman1 + diff_pman1_start;
p1_tend = p1_tstart + diff_pman1_end;
p1_release = p1_tstart + diff_pman1_release;

p1_idx_start = find(data_pman1.time>p1_tstart,1);
p1_idx_end = find(data_pman1.time>p1_tend,1);
p1_idx_release = find(data_pman1.time>p1_release,1);

pman1_release = data_pman1.time(p1_idx_release);

p1_vec = p1_idx_start:p1_idx_end;
p1_tvec = data_pman1.time(p1_vec);

% Video 2:
% First we find the time after the peak in the force data
% pman2_peak_idx = find(data_pman2.Fmag>Fthresh2,1);
% tPP_pman2_idx = find(data_pman2.Fmag(pman2_peak_idx:end)<LL,1);
% tPP_pman2 = data_pman2.time(tPP_pman2_idx + pman2_peak_idx);
tPP_pman2 = 5.28; % time directly after the 3rd peak of second tap
% then we pull in the video data times
pman2_vid_tap = 7.80; %[s]
pman2_vid_start = 11.93; % [s]
pman2_vid_depth = 141.32; % [s]
pman2_vid_release = 144.3; % [s]

diff_pman2_start = pman2_vid_start - pman2_vid_tap;
diff_pman2_end = pman2_vid_depth - pman2_vid_start;
diff_pman2_release = pman2_vid_release - pman2_vid_start;
p2_tstart = tPP_pman2 + diff_pman2_start;
p2_tend = p2_tstart + diff_pman2_end;
p2_release = p2_tstart + diff_pman2_release;

% find indexes for exporting
p2_idx_start = find(data_pman2.time>p2_tstart,1);
p2_idx_end = find(data_pman2.time>p2_tend,1);
p2_idx_release = find(data_pman2.time>p2_release,1);
pman2_release = data_pman2.time(p2_idx_release);

p2_vec = p2_idx_start:p2_idx_end;
p2_tvec = data_pman2.time(p2_vec);

% Video 3
pman3_peak_idx = find(data_pman3.Fmag>Fthresh1,1);
tPP_pman3_idx = find(data_pman3.Fmag(pman3_peak_idx:end)<LL,1);
tPP_pman3 = data_pman3.time(tPP_pman3_idx+pman3_peak_idx);
% then we pull in the video data times
pman3_vid_tap = 8.00; %[s]
pman3_vid_start = 13.15; % [s]
pman3_vid_depth = 81.53; % [s]
pman3_vid_release = 87.35; % [s]

diff_pman3_start = pman3_vid_start-pman3_vid_tap;
diff_pman3_end = pman3_vid_depth - pman3_vid_start;
diff_pman3_release = pman3_vid_release - pman3_vid_start;

p3_tstart = tPP_pman3+diff_pman3_start;
p3_tend = p3_tstart + diff_pman3_end;
p3_release = p3_tstart + diff_pman3_release;

% find indexes for exporting
p3_idx_start = find(data_pman3.time>p3_tstart,1);
p3_idx_end = find(data_pman3.time>p3_tend,1);
p3_idx_release = find(data_pman3.time>p3_release,1);
pman3_release = data_pman3.time(p3_idx_release);

p3_vec = p3_idx_start:p3_idx_end;
p3_tvec = data_pman3.time(p3_vec);

% Video 4
pman4_peak_idx = find(data_pman4.Fmag>Fthresh1,1);
tPP_pman4_idx = find(data_pman4.Fmag(pman4_peak_idx:end)<LL,1);
tPP_pman4 = data_pman4.time(tPP_pman4_idx+pman4_peak_idx);
% then we pull in the video data times
pman4_vid_tap = 3.38; %[s]
pman4_vid_start = 9.48; % [s]
pman4_vid_depth = 74.47; % [s]
pman4_vid_release = 106.12; % [s]

diff_pman4_start = pman4_vid_start-pman4_vid_tap;
diff_pman4_end = pman4_vid_depth - pman4_vid_start;
diff_pman4_release = pman4_vid_release - pman4_vid_start;

p4_tstart = tPP_pman4 + diff_pman4_start;
p4_tend = p4_tstart + diff_pman4_end;
p4_release = p4_tstart + diff_pman4_release;

% find indexes for exporting
p4_idx_start = find(data_pman4.time>p4_tstart,1);
p4_idx_end = find(data_pman4.time>p4_tend,1);
p4_idx_release = find(data_pman4.time>p4_release,1);
pman4_release = data_pman4.time(p4_idx_release);

p4_vec = p4_idx_start:p4_idx_end;
p4_tvec = data_pman4.time(p4_vec);

% Cadaver 1:
% cman1_peak_idx = find(data_cman1.Fmag>50,1);
% tPP_cman1_idx = find(data_cman1.Fmag(cman1_peak_idx:end)<LL,1);
% tPP_cman1 = data_cman1.time(tPP_cman1_idx+cman1_peak_idx);
tPP_cman1 = 81.8599; %[s] manually segmented
% then we pull in the video data times
cman1_vid_tap = 69; %[s] % this is the end tap
cman1_vid_start = 8.7; % [s] %unclear on vid
cman1_vid_depth = 54; % [s]
cman1_vid_release = 106.12; % [s]

diff_cman1_start = cman1_vid_start - cman1_vid_tap;
diff_cman1_end = cman1_vid_depth - cman1_vid_start;
diff_cman1_release = cman1_vid_release - cman1_vid_start;

c1_tstart = tPP_cman1 + diff_cman1_start;
c1_tend = c1_tstart + diff_cman1_end;
c1_release = c1_tstart + diff_cman1_release;

% find indexes for exporting
c1_idx_start = find(data_cman1.time>c1_tstart,1);
c1_idx_end = find(data_cman1.Fmag(c1_idx_start:end)>Fthresh2,1);
c1_idx_end = c1_idx_end + c1_idx_start -1;

if isempty(c1_idx_end)
    c1_idx_end = find(data_cman1.time>c1_tend,1);
end

c1_idx_release = find(data_cman1.time>c1_release,1);
cman1_release = data_cman1.time(c1_idx_release);

c1_vec = c1_idx_start:c1_idx_end;
c1_tvec = data_cman1.time(c1_vec);

% Video 2:
cman2_peak_idx = find(data_cman2.Fmag>Fthresh1,1);
tPP_cman2_idx = find(data_cman2.Fmag(cman2_peak_idx:end)<LL,1);
tPP_cman2 = data_cman2.time(tPP_cman2_idx+cman2_peak_idx);

cman2_vid_tap = 5.17; %[s] % this is the end tap
cman2_vid_start = 12.2; % [s] %unclear on vid
cman2_vid_depth = 83; % [s]
% cman1_release = 106.12; % [s]
diff_cman2_start = cman2_vid_start - cman2_vid_tap;
diff_cman2_end = cman2_vid_depth - cman2_vid_start;
c2_tstart = tPP_cman2 + diff_cman2_start;
c2_tend = c2_tstart + diff_cman2_end;
% find indexes for exporting
c2_idx_start = find(data_cman2.time>c2_tstart,1);
c2_idx_end = find(data_cman2.Fmag(c2_idx_start:end)>Fthresh2,1);
c2_idx_end = c2_idx_end + c2_idx_start - 1;

if isempty(c2_idx_end)
    c2_idx_end = find(data_cman2.time>c2_tend,1);
end

c2_vec = c2_idx_start:c2_idx_end;
c2_tvec = data_cman2.time(c2_vec);

% Video 3:
cman3_peak_idx = find(data_cman3.Fmag>Fthresh1,1);
tPP_cman3_idx = find(data_cman3.Fmag(cman3_peak_idx:end)<LL,1);
tPP_cman3 = data_cman3.time(tPP_cman3_idx+cman3_peak_idx);

cman3_vid_tap = 3.68; %[s] % this is the end tap
cman3_vid_start = 11.68; % [s] %unclear on vid
cman3_vid_depth = 96.8; % [s]
% cman1_vid_release = 106.12; % [s]
diff_cman3_start = cman3_vid_start - cman3_vid_tap;
diff_cman3_end = cman3_vid_depth - cman3_vid_start;
c3_tstart = tPP_cman3 + diff_cman3_start;
c3_tend = c3_tstart + diff_cman3_end;
% find indexes for exporting
c3_idx_start = find(data_cman3.time>c3_tstart,1);
c3_idx_end = find(data_cman3.Fmag(c3_idx_start:end)>Fthresh2,1);
c3_idx_end = c3_idx_end + c3_idx_start - 1;

if isempty(c3_idx_end)
    c3_idx_end = find(data_cman3.time>c3_tend,1);
end

c3_vec = c3_idx_start:c3_idx_end;
c3_tvec = data_cman3.time(c3_vec);

releaseTimes = [pman1_release;pman2_release;...
    pman3_release;pman4_release;cman1_release];
end