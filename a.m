path ='C:\Users\anna9\OneDrive\Υπολογιστής\speech_to_text\isolation\';
iso = {};
for i=1:9
    files = dir(fullfile(path,'*.waV'));
    %fn = files(i).name;
    [filepath,name,ext] = fileparts(files(i).name);
    %[x,fs] = audioread(fn);
    iso{i} = {files(i),name}; 
end
