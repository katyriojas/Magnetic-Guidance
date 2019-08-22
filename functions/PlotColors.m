colors.light_blue = [114, 147, 203];
colors.light_orange = [225, 151, 76];
colors.light_green = [158, 223, 109];
colors.light_red = [211, 94, 96];
colors.light_gray = [128, 133, 133];
colors.light_purple = [144, 103, 167];
colors.light_brown = [171, 104, 87];
% colors.blue = [57, 106, 177];
colors.blue = [45, 85, 177];
colors.orange = [210, 110, 39];
colors.green = [50, 160, 70];
colors.red = [204, 37, 41];
colors.gray = [83, 81, 84];
colors.purple = [107, 76, 154];
colors.brown = [146, 36, 40];

% rescale from [0->255] to [0->1]
F=fieldnames(colors);
for ii=1:numel(F)
   colors.(F{ii}) = colors.(F{ii}) / 255;
end

clear F G