close all;
clear;

S = fileread('Claranet.gml');
nodes = regexp(S, 'node.*?id (?<id>\d+).*?label\s*"(?<label>[^"]*)"', 'names');
edges = regexp(S, 'edge.*?source\s*(?<source>\d+).*?target\s*(?<target>\d+)', 'names');
all_ids = {nodes.id};
all_names = {nodes.label};
all_sources = {edges.source};
all_targets = {edges.target};
[source_found, s] = ismember(all_sources, all_ids);
nfidx = find(~source_found);
if ~isempty(nfidx)
   error('Source ids not found in node list, starting from "%s"', edges(nfidx(1).source));
end
[target_found, t] = ismember(all_targets, all_ids);
nfidx = find(~target_found);
if ~isempty(nfidx)
   error('Target ids not found in node list, starting from "%s"', edges(nfidx(1).target));
end
EdgeTable = table([s.', t.'], ones(length(s),1), 'VariableNames', {'EndNodes' 'Weight'});
NodeTable = table(all_ids.', 'VariableNames',{'Name'});
elabels =1:size(EdgeTable,1)
% for i=1:size(EdgeTable,1)
%     elabels =[elabels int2str(i)];
% end
G = graph(EdgeTable,NodeTable);
h=plot(G);
labeledge(h,1:numedges(G),elabels);
EdgeMatrix=EdgeTable.EndNodes;
EdgeMatrix = [EdgeMatrix, rand(size(EdgeMatrix,1),1)*1000];
AdjacencyMatrix = zeros(size(nodes, 2));
for i=1:1:size(EdgeMatrix,1)
    AdjacencyMatrix(EdgeMatrix(i,1),EdgeMatrix(i,2)) = 1;
    AdjacencyMatrix(EdgeMatrix(i,2),EdgeMatrix(i,1)) = 1;
end
save('adjacencymatrix1','AdjacencyMatrix');
save('links1','EdgeMatrix');

tic
% Operational research topic 3£ºNetwork Tomography
%------Construct netCostMatrix------:
%links = csvread('links.csv', 1, 0);
%adjacencymatrix = csvread('adjacencymatrix.csv');
links = load('links1.mat','EdgeMatrix').EdgeMatrix;
adjacencymatrix = load('adjacencymatrix1.mat', 'AdjacencyMatrix').AdjacencyMatrix;
nodes = uniquetol(links(:,1:2));
true_links = links(:,3);
netCostMatrix = zeros(length(nodes),length(nodes));
for i = 1:length(links)
    netCostMatrix(links(i,1),links(i,2)) = links(i,3);
    netCostMatrix(links(i,2),links(i,1)) = links(i,3);
end
for i = 1:length(nodes)
    for j = 1:length(nodes)
        if netCostMatrix(i,j) == 0
            netCostMatrix(i,j) = inf;
        end
    end
end
%------------Monitor placement-------------:
monitors = [];
for i = 1:length(nodes)
    degree = 0;
    for j = 1:length(nodes)
        degree = degree + adjacencymatrix(i,j);
    end
    if degree == 2
        monitors(end+1) = i;
    end
end
monitors = [1 2 3 5 6 7 9 10 12 14];
k = 1;
 highlight(h,monitors, 'NodeColor','red');
 
%------------Call kShortestPath------------:
[shortestPaths, totalCosts] = kShortestPath(netCostMatrix, monitors, links,k,length(links));

%------Construct measurement matrix A------:
A = zeros(length(shortestPaths),length(links));
for i = 1:length(shortestPaths)
    for j = 1:length(shortestPaths{i})-1
        for m = 1:length(links)
            if min(shortestPaths{i}(j),shortestPaths{i}(j+1)) == links(m,1) && max(shortestPaths{i}(j),shortestPaths{i}(j+1)) == links(m,2)
                A(i,m)=1;
            end
        end
    end
end

[R,jb] = rref(A');
A = A(jb,:);
colorset = rand(size(A,1),3);
for i=1:1:size(A,1)
    edgepath=find(A(i,:)==1);
    figure;
    h=plot(G);
    highlight(h,'Edges',edgepath, 'EdgeColor',colorset(i,:),'linewidth',2);
end

totalCosts = totalCosts(jb);

%------Add measurement noise------:
totalCosts = totalCosts + randn(1,length(totalCosts));

%----------Display results--------:
pred_links = pinv(A'*A)*A'*totalCosts';
mae = sum(abs(pred_links-true_links))/length(pred_links);
independent_path_num = length(jb);
fprintf('monitors =');
disp(monitors);
fprintf('mean_absolute_error =');
disp(mae);
fprintf('total_select_path_num =');
disp(length(shortestPaths));
fprintf('independent_path_num =');
disp(independent_path_num);
toc


%%%% Find perfect cut set for specified scapegoat links. 
R=A;
Ls=[9];
[lm]=Perfectcut(R,Ls)
figure;
h=plot(G);
labeledge(h,1:numedges(G),elabels);
highlight(h,'Edges',Ls, 'EdgeColor','blue','linewidth',3);
highlight(h,'Edges',lm, 'EdgeColor','red','linewidth',3);