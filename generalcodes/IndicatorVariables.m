%function IndicatorVariables
function GroupSplSaiSir=IndicatorVariables(saisir,positions1,positions2)
% input: saisir - data set
%        index - positions in the samples name used as an indicator
% output: GroupSplSaiSir - group matrix in a saisir structure
% Genus group
if nargin == 1
    positions1 = 1:size(saisir.i,2);
    Group=create_group1(saisir,positions1);
elseif nargin == 2
    Group=create_group1(saisir,positions1);
else % (nargin=3)
Group=create_group1(saisir,positions1,positions2);
end
Max=max(Group.d);
[n m]=size(Group.d);
GroupSplitted=zeros(n,Max);
for i=1:Max
    test=find(Group.d==i);
    GroupSplitted(test,i)=1;
end

GroupSplSaiSir.d = GroupSplitted;
GroupSplSaiSir.v = Group.g.i;
if nargin == 3
    GroupSplSaiSir.i = Group.i(:,positions1:positions2);
else    
    GroupSplSaiSir.i = Group.i(:,positions1);
end
GroupSplSaiSir.g = Group.d;
end
