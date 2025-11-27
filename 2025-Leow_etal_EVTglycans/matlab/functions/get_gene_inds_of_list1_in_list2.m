function [inds] = get_gene_inds_of_list1_in_list2(list1,list2)
N_list1 = numel(list1);
inds   = zeros(N_list1,1);
for i=1:N_list1
%     i
    temp = find(strcmpi(list1{i},list2));
    if ~isempty(temp) & numel(temp)==1
%         list1{i}
%         temp
        inds(i)=temp;
%         disp([list1{i} ' found once'])
    elseif  numel(temp)>1
        disp([list1{i} ' duplicated'])
%         inds(i)=-2;
    else
%         disp([list1{i} ' not found'])
    end
end
end

