delete(instrfindall);
s1=serial('COM5','BaudRate',115200);
s1.BytesAvailableFcnMode = 'byte';
s1.InputBufferSize = 200;
s1.BytesAvailableFcnCount=200;

try
 fopen(s1);
catch err
 fclose(instrfind);
 error('err')
end

while 1
 y=fread(s1);
 stem(y,'.');
 %yy=str2double(y);
 y=char(y)
 plot(y);
 drawnow
end

fclose(s1)


hhhh