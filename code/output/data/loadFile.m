function [posx,posy,xc,tau,time] = loadFile(file)
fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
N = val(1);
nb = val(2);
val = val(3:end);

ntime = numel(val)/(2*N*nb+3*nb+1);
% 2 positions, 1 time, 2 centers, 1 angle
if ntime ~= ceil(ntime);
  disp('PROBLEM WITH VALUES FOR N AND nb');
end

posx = zeros(N,nb,ntime);
posy = zeros(N,nb,ntime);
time = zeros(ntime,1);
xc = zeros(2,nb,ntime);
tau = zeros(1,nb,ntime);

istart = 1;
for m = 1:ntime
  % load x positions
  for k=1:nb
    iend = istart + N - 1;
    posx(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end

  % load y positions
  for k=1:nb
    iend = istart + N - 1;
    posy(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end

  % load centers
  for k=1:nb
    iend = istart + 1;
    xc(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end

  % load angles
  for k=1:nb
    iend = istart;
    tau(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end

  time(m) = val(istart);
  istart = istart + 1;
end



