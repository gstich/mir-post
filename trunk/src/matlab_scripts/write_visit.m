function write_visit(file,X,Y,V)

nx = size(X,1);
ny = size(X,2);

dlmwrite(file,'VARIABLES = "X","Y","var"','delimiter','');
dlmwrite(file,['ZONE I=',num2str(nx),', J=', num2str(ny), ', F=POINT'],'delimiter','','-append');

% Modify the format
DATA = zeros(ny,3,nx);
DATA(:,1,:) = X';
DATA(:,2,:) = Y';
DATA(:,3,:) = V';

%dlmwrite(file,DATA,'delimiter',' ','-append','precision','%12.12f');


for j=1:ny
    for i=1:nx
        dlmwrite(file,DATA(j,:,i),'delimiter',' ','-append','precision','%12.12f');
    end
    disp(j);
end