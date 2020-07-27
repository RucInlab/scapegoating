function [] =  NetworkTomography()
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
monitors = [0 2 5 6 8 13]+1;
k = 5;

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
