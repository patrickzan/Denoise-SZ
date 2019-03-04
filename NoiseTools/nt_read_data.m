function [p,data]=nt_read_data(fname)
%[p,data]=nt_read_data(fname) - read data from file
%
% 
%  fname: file to read
%  
% 
VERBOSE=1;

if nargin < 1 ; error('!'); end
if nargin < 2; filetype=[]; end
if nargin < 3; variable=[]; end

if ~isa(fname, 'char'); 
    error('filename is not char string');
end
if ~exist(fname,'file')==2;
    error(['file >', data, '< not found']);
end

% standard fields
p.fname=fname;
p.read_with=[];
p.sr=[];

% intercept .mat files
if numel(fname)>4 & fname(end-3:end)=='.mat'
    p.read_with='.mat file';
    if VERBOSE; disp('mat file'); end
    % list variables in file, ask user to choose one
    S=whos('-file',fname);
    var_names=char(S.name);
    var_sizes=round([S.bytes]/1024)';
    a=repmat(' (', size(var_names,1),1);
    b=cellstr(num2str(var_sizes, '%9d'));
    b=char(b);
    c=[repmat(' Mb)', size(var_names,1),1)];
    i=listdlg('liststring',cellstr([var_names,a,b,c]),...
        'name', ['Select variable in file ',fname], ...
        'listsize', [600 300], ...
        'OKstring','Select',...
        'PromptString','select:');
    if isempty(i); data=[]; return; end
    if nargout>1;
        data=load(fname,deblank(var_names(i,:)));
        % if it's a structure, list fields, ask user to choose one
        while isstruct(data);
            if numel(data)>1
                i=listdlg('liststring',cellstr(S),...
                    'name', ['Select element of strucure array ',var_names(i,:)], ...
                    'listsize', [600 300], ...
                    'OKstring','Select',...
                    'PromptString','select:');
                if i ; data=data(i); end
            end
            S=fieldnames(data);
            i=listdlg('liststring',cellstr(S),...
                'name', ['Select field in struct ',var_names(i,:)], ...
                'listsize', [600 300], ...
                'OKstring','Select',...
                'PromptString','select:');
            if i ; data=getfield(data,S{i}); end
        end
    end
    return
end

% intercept Yokogawa files
if numel(fname)>4 & (fname(end-3:end)=='.con' | fname(end-3:end)=='.sqd')
    p.read_with='yokogawa 2013';
    p.acq_cond = getYkgwHdrAcqCond(fname);
    p.channel_info=getYkgwHdrChannel(fname);
    p.system_info=getYkgwHdrSystem(fname);
    p.event=getYkgwHdrEvent(fname);
    % read other info?
    p.sr=p.acq_cond.sample_rate;
    if nargout>1;
        data=getYkgwData(fname)';
    end
    return
end
   
         
% select file reader among those available
has_ft_reader=0; 
has_sopen=0;
if 2==exist('ft_read_header');
    has_ft_read_header=1;
else
    warning('function ft_read_header() not found: download FieldTrip and/or adjust path');
end
if 2==exist('sopen');
    has_sopen=1;
else
    warning('function sopen() not found: download BIOSIG and/or adjust path');
end
    
    
if has_ft_read_header
    isftReadable=0;
    try
        % readable by FieldTrip?
        h=ft_read_header(fname);
        isftReadable=1;
    catch
        ; % can't read
    end
end
if ~isftReadable & has_sopen
    isBiosigReadable=0;
    try
        % readable by FieldTrip?
        h=sopen(fname);
        isBiosigReadable=1;
        sclose(h);
    catch
        ; % can't read
    end
end
    
if isftReadable
    if VERBOSE; disp('read with FieldTrip'); end
    h=ft_read_header(fname);    
    p.header=h;
    p.read_with='FieldTrip';
    p.sr=h.Fs;
    if nargout>1;
        data=ft_read_data(fname)';
    end
elseif isBiosigReadable
    if VERBOSE; disp('read with Biosig'); end
    h=sopen(fname);
    p.header=h;
    p.read_with('BIOSIG');
    p.sr=h.SampleRate;
    if nargout>1;
        data=sread(fname)';
    end
    sclose(h);
else
    ismatfile=0;
    try
        % .mat file?
        S=whos('-file',data);
        if numel(S)>1
            if nargout==2
                [p,data]=nt_read_data([fname,'.mat']);
            else
                [p,data]=nt_read_data([fname,'.mat']);
            end
        end
    catch
        disp(['File >',fname,'< is not a matlab file, and FieldTrip and BIOSIG can''t read it']);
    end
end
    

