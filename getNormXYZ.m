function dataNorm=getNormXYZ(data)
%gets the norm (or sum vector) of a matrix of three columns (x,y,z)
if size(data,2)~=3
    error('x y z');
end
for i=1:size(data,1)
    dataNorm(i,1)=norm(data(i,:));
end