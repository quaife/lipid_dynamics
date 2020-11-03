function [velx,vely,torq] = loadVelFile(file)
fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
N = val(1);
nb = val(2);
val = val(3:end);

ntime = numel(val)/(3*nb+1);
% 2 velocities, 1 torque, 1 time
if ntime ~= ceil(ntime);
  disp('PROBLEM WITH VALUES FOR N AND nb');
end

time = zeros(ntime,1);
velx = zeros(1,nb,ntime);
vely = zeros(1,nb,ntime);
torq = zeros(1,nb,ntime);

istart = 1;
for m = 1:ntime
  % load x component of velocity
  for k=1:nb
    iend = istart;
    velx(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end

  % load y component of velocity
  for k=1:nb
    iend = istart;
    vely(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end

  % load torque
  for k=1:nb
    iend = istart;
    torq(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end

  time(m) = val(istart);
  istart = istart + 1;
end



