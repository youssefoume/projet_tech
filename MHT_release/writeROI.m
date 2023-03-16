function writeROI(data,filename)
%ROI cell:data
fid = fopen(filename,'wt');
for i=1:length(data)
    
    fprintf(fid,'%g\t',data(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
end