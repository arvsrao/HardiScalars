%
%  function ret=write_mrtrix_tracks(tracks,filename)
%

function ret=write_mrtrix_tracks(tracks,filename)

datatype = lower(tracks.datatype);
byteorder = datatype(end-1:end);

if strcmp(byteorder, 'le')
  f = fopen (filename, 'w', 'l');
  datatype = datatype(1:end-2);
elseif strcmp(byteorder, 'be')
  f = fopen (filename, 'w', 'b');
  datatype = datatype(1:end-2);
else
  disp ('unexpected data type - aborting')
  return;
end


if (f<1) 
  disp (['error opening ' filename ]);
  return
end

fprintf(f,'mrtrix tracks\n');
count=13;
keys=fieldnames(tracks);
for i=1:length(keys)
value=getfield(tracks,keys{i});
if( ~strcmp(keys{i},'data') & ~strcmp(keys{i},'seed') & ~strcmp(keys{i},'include') & ~strcmp(keys{i},'offset'))
header=[keys{i},': ',value,'\n'];
count=count+length(header)-1;
fprintf(f,header);
end
if strcmp(keys{i},'seed')
    header=['roi: seed',value,'\n'];
count=count+length(header)-1;
fprintf(f,header);
end
if strcmp(keys{i},'include')
    header=['roi: include',value,'\n'];
count=count+length(header)-1;
fprintf(f,header);
end
if strcmp(keys{i},'offset')
    offset=str2num(value);
end
keys{i}
end


%offset=count;
%offset=offset+length(int2str(offset))+9+4;
header_offset=['file: . ',int2str(offset),'\n'];
fprintf(f,header_offset);
fprintf(f,'END\n');
fclose(f);

if strcmp(byteorder, 'le')
  f = fopen (filename, 'a', 'l');
  datatype = datatype(1:end-2);
elseif strcmp(byteorder, 'be')
  f = fopen (filename, 'a', 'b');
  datatype = datatype(1:end-2);
else
  disp ('unexpected data type - aborting')
  return;
end

if (f<1) 
  disp (['error opening ' filename ]);
  return
end

position=ftell(f);
buffer=offset-position;
fseek (f, offset, -1);
while(position-offset<0)
fwrite(f,0,'char');
position=ftell(f);
end

data=[];
for i=1:length(tracks.data)
    data=[data;tracks.data{i};NaN,NaN,NaN];
end
disp('writing data');
data=reshape(data',1,prod(size(data)));
fwrite(f, data, datatype);
fclose(f);
ret=1;
end
