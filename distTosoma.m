function value = distTosoma(data_xyz,connct,type,label)
% calculate path distance from every skeletal points to the estimated centre of soma from SWC file
distlist = 0;

if sum(type==1)== 3
    %find rows of peripheral soma points 
    row_soma = find(type==1);
    row_somap = row_soma(connct(row_soma)==1);
    data_xyz(row_somap,:) = [];
    connct(row_somap) = [];
    label(row_somap) = [];   
    type(row_somap) = [];
%     connct(connct~=-1)=connct(connct~=-1)-2;
    connct(connct==-1) = 1;

elseif sum(type==1) == 1
    connct(connct==-1) = 1;
end

for i = 2:length(connct)
    row_parent  = find(label==connct(i));
    disttoparent = norm(data_xyz(i,:)-data_xyz(row_parent,:));
    distitosoma = distlist(row_parent)+ disttoparent;
    distlist = [distlist; distitosoma];
end

%axon postive dendrite negative sign
value = distlist.*(-sign(type-2.5));

end

