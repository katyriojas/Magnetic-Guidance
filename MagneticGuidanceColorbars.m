%% Load data

% guided/unguided phantom data
% load('data\phantom\data_phantom.mat');

% manual phantom data
load('data\phantom\ManualPhantomData.mat');


%% Create base profile (circle/square)

base_type = 'circle';
width = 20;
spacing = 1.02*width;

if strcmp(base_type, 'square') 
    q = linspace(0,2*pi,5) + pi/4;
    base = (width/2)*[cos(q); sin(q)];   % Base curve is a square
elseif strcmp(base_type, 'circle')
    % Base curve is a circle
    cir_res = 40;   % circle resolution (# pts)
    q = linspace(0, 2*pi, cir_res);
    base = (width/2)*[cos(q); sin(q)]; 
else
    error('base_type must be ''circle'' or ''square''')
end
  
%% create paths

% trim 10 degrees back from final angular depth
% trim_deg = 10;
% mag1_end = mag1.interp_angdepth(end)-trim_deg;
% nomag1_end = nomag1.interp_angdepth(end)-trim_deg;
% mag1_end_ind = find(mag1.interp_angdepth >= mag1_end, 1);
% nomag1_end_ind = find(nomag1.interp_angdepth >= nomag1_end, 1);
% max_ang_index1 = find(st_theta > mag1_end, 1);
% max_ang_index2 = find(st_theta > nomag1_end, 1);


for ii=1:length(data)
    path(ii).mag   = [data(ii).mag_interp_angdepth,  (ii*spacing-width/3)*ones(size(data(ii).mag_interp_angdepth)), zeros(size(data(ii).mag_interp_angdepth))]';
    path(ii).nomag = [data(ii).nomag_mea_interp_angdepth,  (-ii*spacing+width/3)*ones(size(data(ii).nomag_mea_interp_angdepth)), zeros(size(data(ii).nomag_mea_interp_angdepth))]';
end

%% extrude path to create surface

for ii=1:length(data)

    [X(ii).mag, Y(ii).mag, Z(ii).mag] = extrude(base, path(ii).mag, 0);
    [X(ii).nomag, Y(ii).nomag, Z(ii).nomag] = extrude(base, path(ii).nomag, 0);
end

%% plot
figure(13); clf(13); 
axis equal; 
hold on; 
% grid off;

% plot extrusions
for ii=1:length(data)
    h(ii).mag = surf(X(ii).mag, Y(ii).mag, Z(ii).mag);
    h(ii).mag.CData = repmat(data(ii).mag.Fmag_smooth', [length(base), 1]);

    h(ii).nomag = surf(X(ii).nomag, Y(ii).nomag, Z(ii).nomag);
    h(ii).nomag.CData = repmat(data(ii).nomag_mea.Fmag_smooth', [length(base), 1]);
end



shading flat
view(0,90)

% rescale colormap
cMap = hot(256+70);
% cMap = cool(256);
cMap = cMap(1:end-70,:);
dataMin = 0;
dataMax = 100;
% dataMax = max(data_mag1.Fmag_smooth);
% dataMin = min(data_mag1.Fmag_smooth);
centerPoint = 5;
scalingIntensity = 5;

x = 1:length(cMap); 
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*511/max(x)+1; 
newMap = interp1(x, cMap, 1:512);
cmap_scaled = colormap(newMap);

xlabel('Angular Insertion Depth (deg)')
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])

cbar = colorbar;
cbar.Label.String = 'Force (mN)';
cbar.Label.FontSize = 14';
