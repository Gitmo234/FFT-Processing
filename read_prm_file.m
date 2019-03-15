% script to read a basic midas file, and save the each fft (or psf)

%% config variable ------------------------------------------------
%% values to change for data processing
skip_seconds = 1; % number of seconds to skip at beginning of data
fft_seconds  = 1; % number of seconds to use for input for fft 

%% values for the output
outputname_prefix = 'output'; % prefix for the output file name

%% values for the plot
showplot = 1; % if showplot, only show a single plot

% fft size
fft_pwr2 = 16;
nfft = 2^fft_pwr2;
interp_f = 1;

% window
win = hanning(nfft);
n_win = length(win);


%% main program ------------------------------------------------

if exist ('OCTAVE_VERSION', 'builtin')
  pkg load signal 
end

if ~exist('stdout')
  stdout = 1;
end

last_path = pwd();

[file,path] = uigetfile([last_path filesep '*.prm'], 'Select the .PRM file');

filename = fullfile(path,file);

% open file
fid = fopen (filename, 'r');

if fid == -1
  error('Couldnt open file "%s"', filename);
end

try

  % read header
  fseek(fid, 0, 'bof');
  hcb = {};
  hcb.version = fread(fid, 4, '*char')';
  hcb.head_rep = fread(fid, 4, '*char')';
  hcb.data_rep = fread(fid, 4, '*char')';
  hcb.detached = fread(fid, 1, 'int32');
  hcb.protected = fread(fid, 1, 'int32');
  hcb.pipe = fread(fid, 1, 'int32');
  hcb.ext_start = fread(fid, 1, 'int32') * 512;
  hcb.ext_size = fread(fid, 1, 'int32');
  hcb.data_start = fread(fid, 1, 'double');
  hcb.data_size = fread(fid, 1, 'double');
  hcb.type = fread(fid, 1, 'int32');  
  hcb.format = fread(fid, 2, '*char')';
  hcb.flagmask = fread(fid, 1, 'int16');
  hcb.timecode = fread(fid, 1, 'double');
  hcb.inlet = fread(fid, 1, 'int16');
  hcb.outlets = fread(fid, 1, 'int16');
  hcb.outmask = fread(fid, 1, 'int32');
  hcb.pipeloc = fread(fid, 1, 'int32');
  hcb.pipesize = fread(fid, 1, 'int32');
  hcb.in_byte = fread(fid, 1, 'double');
  hcb.out_byte = fread(fid, 1, 'double');
  hcb.outbytes = fread(fid, 8, 'double');
  hcb.keylength = fread(fid, 1, 'int32');
  hcb.keywords = fread(fid, 92, '*char')';
            
  % check looks like a midas file
  if ~strcmp(hcb.version,'BLUE')
    error ('Not a BLUE formatted midas file');
  end
  if ~strcmp(hcb.head_rep,'EEEI') || ~strcmp(hcb.data_rep,'EEEI')
    error ('Header or Data representation is not EEEI');
  end  
  if hcb.type ~= 1000
    error ('Can only process type 1000 files');
  end
  if hcb.format(1) ~= 'S'
    error ('Can only process scalar type files');
  end  
  
  % read the type 1000 adjunct
  hcb.adjunct = struct(...
   'xstart', fread(fid, 1, 'double'), ...
   'xdelta', fread(fid, 1, 'double'), ...
   'xunits', fread(fid, 1, 'int32') ...
  ); 
 
  % read ext header (currently not needed)
  % do nothing
 
  % get info we need from the header 
  switch hcb.format(2)
    case 'I'
      sampletype = 'int16';
      samplesize = 2;
    case 'L'
      sampletype = 'int32';
      samplesize = 4; 
    case 'F'
      sampletype = 'single';
      samplesize = 4; 
    otherwise
      error ('unsupported format type of "%c"',  hcb.format(2));  
  end
  
  Fs = 1.0/hcb.adjunct.xdelta
  sampletype
  datasize = hcb.data_size/samplesize  
  
  skip_samples = skip_seconds * Fs
  fft_samples = fft_seconds * Fs
  
  % go start of data
  fseek(fid, hcb.data_start, 'bof');
  
  % skip initial data that we wanted to skip
  fseek(fid, skip_samples*samplesize, 'cof');
  datasize = datasize - skip_samples;
  
  % loop reading the data in 
  cnt = 1;
  
  while datasize >= fft_samples
    x = single(fread(fid, [1, fft_samples], sampletype));

    datasize = datasize - length(x);

    % convert to  frequency domain data
    if exist ('OCTAVE_VERSION', 'builtin')
      overlap = 0.5
    else
      overlap = nfft/2
    end
    
    [psd, f] = pwelch(x-mean(x), win, overlap, interp_f*nfft, Fs);
  
    % on first time print out some stats
    if cnt == 1
      fprintf(stdout, 'psd size = %d %s samples \n', length(psd), class(psd)); 
      fprintf(stdout, 'freq span %f Hz - %f Hz\n', min(f), max(f));
      fprintf(stdout, 'binres = %f Hz\n', f(2)-f(1));
    end
  
    % if showplot, only show a plot and then break from loop
    if showplot
      figure
      plot(f/1e3, 10*log10(psd));
      grid on;
      xlabel('Frequency (in kHz)');
      ylabel('Power/Frequency (dB/Hz)');
      title(['Noise PSD Comparison - ' num2str(cnt)], 'Interpreter', 'none'); 
      break 
    else

      % write psd data to file
      fname = [outputname_prefix '-' num2str(cnt) '.psd'];
      fprintf (stdout, 'writing %s\n', fname); 
      fod = fopen(fname, 'w');
      if fod ~= -1
        fwrite(fod, psd); 
        fclose(fod);
      end
    end
    
    % inc the count for output file numbering
    cnt = cnt + 1;
  end


catch err
  fclose(fid);
  
  rethrow (err);
end
