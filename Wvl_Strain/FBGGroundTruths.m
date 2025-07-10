function [xyz, xyze] = FBGGroundTruths(directory, filetype, sample, xyze, Box2Needle_frame)

path = directory + filesep +  filetype;
list = dir(path);
names = {list(:).name};
names = cellfun(@convertCharsToStrings, names);
filepaths = directory + filesep + names;

% 3D coordinates
switch filetype
    
    case "*.csv"
        xyz = readtable(filepaths(sample));
        xyz = table2array(xyz);
    
    case "*.xlsx"
        datafile = filepaths(sample);
        shts = sheetnames(datafile);
        xyzCT = readtable(datafile,'Sheet',shts(1));
        xyzCT = table2array(xyzCT);

        % Transformation matrix
        Box2CT_frame = readtable(datafile,'Sheet',shts(3));
        Box2CT_frame = table2array(Box2CT_frame);
        CT2Needle_frame = Box2Needle_frame*inv(Box2CT_frame);

        % GT in needle frame
        R = CT2Needle_frame(1:3,1:3)';
        t = CT2Needle_frame(1:3,4);
        xyz = R*xyzCT' - R*t;
        xyz = xyz';

    case "*.dat"
        xyz = readtable(filepaths(sample));
        xyz = xyz{1:end,1:3};
end

% [xyze, xyz, Ls] = lengthcorr(xyze, xyz);

% Tip removal
xyz = tipremove(xyz');
xyz = flip(xyz,2);

% GT s array (initial)
sGT = scal(xyz);
[sGT, index] = unique(sGT);
xyz = xyz(:,index);

% Interpolation
s = scal(flip(xyze,1)');
s = s(s <= sGT(end));
xyz = interp1(sGT, xyz', s, 'linear', 'extrap');
sGT = scal(xyz');
error = mean(abs(s - sGT));

% Loop
while error > 1e-5
    xyz = interp1(sGT, xyz, s, 'linear', 'extrap');
    sGT = scal(xyz');
    error = mean(abs(s - sGT));
end
xyz = flip(xyz,1);

%%%%%%%%
function [newcoord] = tipremove(coord)

    subs = diff(coord, 1, 2);
    lx = vecnorm(subs, 2, 1);
    L = sum(lx);
    x = [0, cumsum(lx)];

    tip = 4;
    idx = find(L - x >= tip, 1, 'last');
    l = tip - (L - x(idx));
    p1 = coord(:, idx);
    p2 = coord(:, idx + 1);
    Lseg = lx(idx);
    t = l / Lseg;
    newPt = p2 - t * (p2 - p1);
    newcoord = [coord(:, 1:idx), newPt];

end % tipremove
%%%%%%%%

%%%%%%%%
function [s] = scal(coord)

    subs = diff(coord, 1, 2);
    ls = vecnorm(subs, 2, 1);
    s = [0, cumsum(ls)];   

end % function scal
%%%%%%%%

end % function FBGGroundTruths